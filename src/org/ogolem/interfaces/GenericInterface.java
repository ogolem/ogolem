/*
Copyright (c) 2013-2014, J. M. Dieterich
                2015, J. M. Dieterich and B. Hartke
                2017-2020, J. M. Dieterich and B. Hartke
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * All advertising materials mentioning features or use of this software
      must display the following acknowledgement:

      This product includes software of the ogolem.org project developed by
      J. M. Dieterich and B. Hartke (Christian-Albrechts-University Kiel, Germany)
      and contributors.

    * Neither the name of the ogolem.org project, the University of Kiel
      nor the names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
package org.ogolem.interfaces;

import java.io.File;
import org.ogolem.core.AbstractLocOpt;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CartesianFullBackend;
import org.ogolem.core.Constants;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Gradient;
import org.ogolem.core.StreamGobbler;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * A generic interface for ogolem. Quick'n'dirty way of enabling the use of other codes then the
 * ones we have explicit (and optimized) interfaces for.
 *
 * <p>Specifications: * we need the environment variable OGO_GENERALCMD to be set * OGO_GENERALOPTS
 * are optional but will be appended to CMD if existing. * we create a directory for each task and
 * write in it a file input.xyz * we call the command in that directory * we expect as output a file
 * output.xyz in standard xyz file format and in the second row (comment row) we expect as the
 * second token the energy in kJ/mol.
 *
 * <p>Happy scripting! :-)
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class GenericInterface extends AbstractLocOpt implements CartesianFullBackend {

  private static final long serialVersionUID = (long) 20130927;
  private final String cmd;
  private final String opts;

  public GenericInterface(final GlobalConfig globconf) throws Exception {
    super(globconf);
    final String tmp = System.getenv("OGO_GENERALCMD");
    if (tmp == null) {
      throw new Exception("OGO_GENERALCMD not specified in environment.");
    } else {
      cmd = tmp;
    }
    final String tmp2 = System.getenv("OGO_GENERALOPTS");
    if (tmp2 == null) {
      System.err.println("WARNING: OGO_GENERALOPTS not specified in environment.");
      opts = "";
    } else {
      opts = tmp2;
    }
  }

  private GenericInterface(final GenericInterface orig) {
    super(orig);
    this.cmd = orig.cmd;
    this.opts = orig.opts;
  }

  @Override
  public GenericInterface copy() {
    return new GenericInterface(this);
  }

  @Override
  public String myID() {
    return "Generic interface calling: " + cmd;
  }

  @Override
  public String myIDandMethod() {
    return "Generic interface calling: " + cmd;
  }

  @Override
  public String getMethodID() {
    return "Generic interface calling: " + cmd;
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      final long lID,
      final CartesianCoordinates cartes,
      final boolean[][] constraints,
      final boolean isConstricted,
      final BondInfo bonds)
      throws Exception {

    // create a directory
    final String dirName = "generic-" + lID + "-opt";
    final String xyzFile = dirName + File.separator + "input.xyz";
    final String xyzRes = dirName + File.separator + "output.xyz";
    OutputPrimitives.createAFolder(dirName);

    // put the xyz info there
    final String[] printXYZ = cartes.createPrintableCartesians();
    OutputPrimitives.writeOut(xyzFile, printXYZ, false);

    // execute the command
    final Runtime rt = Runtime.getRuntime();
    final Process proc = rt.exec(cmd + " " + opts, null, new File(dirName));

    // any error message?
    final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

    // any output?
    final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

    // kick them off
    errorGobbler.start();
    outputGobbler.start();

    // any error???
    if (proc.waitFor() != 0) {
      throw new Exception(cmd + " returns non-zero return value (local optimization).");
    } // cmd should(!) have completed normally...

    // read the output back in
    final CartesianCoordinates res = cartes.copy();
    try {
      final String[] resDat = InputPrimitives.readFileIn(xyzRes);
      final String[] atoms = res.getAllAtomTypes();
      final double[][] xyz = res.getAllXYZCoord();
      final int noAtoms = Integer.parseInt(resDat[0].trim());
      if (noAtoms != res.getNoOfAtoms()) {
        throw new Exception("Wrong number of atoms in result xyz.");
      }
      final double e =
          Double.parseDouble(resDat[1].trim().split("\\s+")[1]) * Constants.KJTOHARTREE;
      res.setEnergy(e);
      for (int i = 2; i < noAtoms + 2; i++) {
        final String[] line = resDat[i].trim().split("\\s+");
        // check the atom type
        if (!line[0].equalsIgnoreCase(atoms[i - 2])) {
          throw new Exception("Wrong atom in output. " + line[0] + " should be " + atoms[i - 2]);
        }
        // parse the coordinates
        xyz[0][i - 2] = Double.parseDouble(line[1]) * Constants.ANGTOBOHR;
        xyz[1][i - 2] = Double.parseDouble(line[2]) * Constants.ANGTOBOHR;
        xyz[2][i - 2] = Double.parseDouble(line[3]) * Constants.ANGTOBOHR;
      }
    } catch (Exception e) {
      throw new Exception("Failure to read output of generic code " + cmd, e);
    }

    // cleanup
    ManipulationPrimitives.remove(dirName);

    return res;
  }

  @Override
  public void gradientCalculation(
      long lID,
      int iIteration,
      double[] xyz1D,
      String[] saAtomTypes,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      int iNoOfAtoms,
      float[] faCharges,
      short[] iaSpins,
      BondInfo bonds,
      final Gradient gradient,
      final boolean hasRigidEnv) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public double energyCalculation(
      long lID,
      int iIteration,
      double[] xyz1D,
      String[] saAtomTypes,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      int iNoOfAtoms,
      float[] faCharges,
      short[] iaSpins,
      BondInfo bonds,
      final boolean hasRigidEnv) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public CartesianFullBackend getBackend() {
    return new GenericInterface(this);
  }
}
