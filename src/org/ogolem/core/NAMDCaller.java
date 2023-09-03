/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.core;

import java.io.File;
import java.io.IOException;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Calling NAMD to do the local optimization.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
final class NAMDCaller extends AbstractLocOpt {

  // the ID
  private static final long serialVersionUID = (long) 20101108;

  private final int noIterations;
  private final String outName;

  NAMDCaller(final GlobalConfig globconf) {
    super(globconf);
    this.noIterations = globconf.maxIterLocOpt;
    this.outName = globconf.outputFolder;
  }

  private NAMDCaller(final NAMDCaller orig) {
    super(orig);
    this.noIterations = orig.noIterations;
    this.outName = orig.outName;
  }

  @Override
  public NAMDCaller copy() {
    return new NAMDCaller(this);
  }

  @Override
  public String myIDandMethod() {
    return "NAMD";
  }

  @Override
  public String myID() {
    return "NAMD";
  }

  @Override
  public Molecule localOptimization(Molecule mStartMolecule) {
    // this caller is not capable of molecular optimization due to the amber binary ff format
    System.out.println(
        "WARNING: NAMD caller can't locally optimize a molecule. Returning non-optimized molecule.");
    mStartMolecule.setEnergy(FixedValues.NONCONVERGEDENERGY);

    return mStartMolecule;
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      long lID,
      CartesianCoordinates cartes,
      final boolean[][] constraints,
      final boolean isConstricted,
      final BondInfo bonds)
      throws ConvergenceException {

    final String sNAMDFolder = "namd" + lID;

    // first create the NAMD folder where everything else will take place
    try {
      OutputPrimitives.createAFolder(sNAMDFolder);
    } catch (IOException e) {
      throw new ConvergenceException("Error creating NAMD folder!", e);
    }

    final float[] faCharges = cartes.getAllCharges();
    final short[] iaSpins = cartes.getAllSpins();

    /*
     * write the output aka input for the to be called program
     */
    try {
      // yes, indeed, it IS the amber coordinates format
      AmberCaller.writeAmberCoordinates(cartes.getAllXYZCoord(), sNAMDFolder);
    } catch (IOException e) {
      throw new ConvergenceException(
          "Problem in writing geometry for namd input (local optimization).", e);
    }

    final String sParameterFile =
        System.getProperty("user.dir") + System.getProperty("file.separator") + outName + ".prmtop";
    final String sCopiedPrmFile =
        System.getProperty("user.dir")
            + System.getProperty("file.separator")
            + sNAMDFolder
            + System.getProperty("file.separator")
            + outName
            + ".prmtop";
    try {
      // yup, we use the amber prmtop
      OutputPrimitives.createLink(sParameterFile, sCopiedPrmFile);
    } catch (IOException e) {
      throw new ConvergenceException(
          "Problem in copying prmtop for namd input (local optimization).", e);
    }

    // now we need to indeed write a configuration file
    try {
      writeNAMDConfiguration(sNAMDFolder, outName, noIterations);
    } catch (Exception e) {
      throw new ConvergenceException("Problem in writing namd input (local optimization).", e);
    }

    // we can just call namd2 with the config file, that's it
    Process proc = null;
    String[] saOutput;
    try {
      Runtime rt = Runtime.getRuntime();
      String[] saCommand = new String[] {"namd2", "configuration.conf"};
      File dir = new File(sNAMDFolder);
      proc = rt.exec(saCommand, null, dir);

      // any error message?
      StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      int iExitValue = proc.waitFor();
      if (iExitValue != 0) {
        if (DEBUG) {
          final String[] out = outputGobbler.getData();
          final String[] err = errorGobbler.getData();
          System.out.println("DEBUG: Out coming");
          for (final String s : out) {
            System.out.println(s);
          }
          System.out.println("DEBUG: Err coming");
          for (final String s : err) {
            System.out.println(s);
          }
          System.out.println("DEBUG: Done");
        }
        throw new ConvergenceException(
            "NAMD returns non-zero return value (local optimization): " + iExitValue + ".");
      } else {
        // namd should(!) have completed normally...
        saOutput = outputGobbler.getData();
      }

    } catch (Exception e) {
      System.err.println(e.toString());
      throw new ConvergenceException("NAMD has a problem (local optimization).", e);
    } finally {
      if (proc != null) {
        proc.destroy();
      }
    }

    if (DEBUG) {
      for (final String s : saOutput) {
        System.out.println("DEBUG: " + s);
      }
    }

    /*
     * read the output of NAMD
     */
    CartesianCoordinates optCartes;
    try {
      optCartes =
          Input.ReadXYZNAMDOutput(
              sNAMDFolder,
              cartes.getAllAtomTypes(),
              cartes.getNoOfAtoms(),
              cartes.getNoOfMolecules(),
              cartes.getAllAtomsPerMol());
    } catch (Exception e) {
      throw new ConvergenceException("Problem in reading the output of NAMD.", e);
    }

    /*
     * read the energy in
     */
    if (!saOutput[saOutput.length - 1].equalsIgnoreCase("Program finished.")) {
      throw new ConvergenceException("NAMD output didn't contain a program finished line.");
    } else {
      // now we search for the energy
      final String[] energyLine = saOutput[saOutput.length - 11].trim().split("\\s+");

      if (DEBUG) {
        for (final String s : energyLine) {
          System.out.println("DEBUG: energyLine " + s);
        }
        System.out.println("DEBUG: Trying to parse " + energyLine[11]);
      }

      // parse it
      double energy = FixedValues.NONCONVERGEDENERGY;
      try {
        energy = Double.parseDouble(energyLine[11]) * Constants.KCALTOHARTREE;
        if (DEBUG) {
          System.out.println("DEBUG: Energy is " + energy);
        }
      } catch (Exception e) {
        throw new ConvergenceException("NAMDCaller: Couldn't parse the energy.", e);
      }

      // set it in
      optCartes.setEnergy(energy);
    }

    /*
     * clean up
     */
    try {
      ManipulationPrimitives.remove(sNAMDFolder);
    } catch (Exception e) {
      System.err.println("Problem cleaning NAMD files up." + e.toString());
    }
    optCartes.setAllCharges(faCharges);
    optCartes.setAllSpins(iaSpins);

    return optCartes;
  }

  private static void writeNAMDConfiguration(
      final String folder, final String runName, final int iterations) throws IOException {

    final String[] input = new String[14];

    input[0] = "minimization on";
    input[1] = "temperature 0";
    input[2] = "numsteps " + iterations;
    input[3] = "exclude scaled1-4";
    input[4] = "1-4scaling 1.0";
    input[5] = "scnb 1";
    input[6] = "switching off";
    input[7] = "cutoff 998";
    input[8] = "pairlistdist 999";
    input[9] = "amber on";
    input[10] = "parmfile " + runName + ".prmtop";
    input[11] = "ambercoor coordinates.inpcrd";
    input[12] = "outputname opt";
    input[13] = "binaryoutput no";

    final String configPath = folder + System.getProperty("file.separator") + "configuration.conf";

    OutputPrimitives.writeOut(configPath, input, false);
  }
}
