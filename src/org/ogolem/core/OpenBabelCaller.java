/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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

import static org.ogolem.core.Constants.*;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * This calls openbabel for local geometry optimization. Openbabel DOES NOT (at least at the moment)
 * support charged atoms. The routine will conserve the charges though. Anyhow, it is highly
 * discouraged at the moment to use openbabel as an backend for systems containing charges!
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class OpenBabelCaller extends AbstractLocOpt {

  // the ID
  private static final long serialVersionUID = (long) 20150421;

  static enum FORCEFIELD {
    GHEMICAL,
    MMFF94,
    MMFF94s,
    UFF
  };

  private final FORCEFIELD forceField;
  private String[] saConnects;
  private final int iNoOfIterations;
  private final boolean noCaching;

  /**
   * The constructor for usage as a local optimizing engine.
   *
   * @param globconf A complete configuration set.
   */
  OpenBabelCaller(GlobalConfig globconf, FORCEFIELD whichFF, final boolean noCache) {
    super(globconf);

    if (DEBUG) {
      System.out.println("DEBUG: Configuring openBabel with " + noCache + " " + whichFF.name());
    }

    this.noCaching = noCache;
    this.forceField = whichFF;
    this.saConnects = null;
    this.iNoOfIterations = globconf.maxIterLocOpt;
  }

  private OpenBabelCaller(OpenBabelCaller orig) {
    super(orig);
    this.noCaching = orig.noCaching;
    this.forceField = orig.forceField;
    this.saConnects = null;
    this.iNoOfIterations = orig.iNoOfIterations;
  }

  @Override
  public OpenBabelCaller copy() {
    return new OpenBabelCaller(this);
  }

  @Override
  public String myID() {
    return "OpenBabel";
  }

  @Override
  public String myIDandMethod() {
    return "OpenBabel: " + forceField.name();
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      long lID,
      CartesianCoordinates cartes,
      boolean[][] constraints,
      boolean isConstricted,
      final BondInfo bonds)
      throws ConvergenceException {
    String sGeometryFile = "openbabel" + lID + ".mdl";

    float[] faCharges = cartes.getAllCharges();
    short[] iaSpins = cartes.getAllSpins();

    /*
     * write the output aka the cartesian coordinates as an mdl-format file
     */
    if (saConnects == null || noCaching) {
      // first time usage of this object, create the connectivity
      doConnects(bonds, cartes.getNoOfAtoms());
      if (DEBUG) {
        System.out.println("DEBUG: Finished bonding in openBabel.");
      }
    }

    // put the things together to an mdl file
    String[] saMDLContent = new String[saConnects.length + cartes.getNoOfAtoms() + 12];
    saMDLContent[0] = "null";
    saMDLContent[1] = " OGOLEM";
    saMDLContent[2] = "";

    // now a line with the number of atoms and bonds, rest like openbabel would do it
    saMDLContent[3] = " 0 0 0 0 0 999 V3000";
    saMDLContent[4] = "M  V30 BEGIN CTAB";
    saMDLContent[5] = "M  V30 COUNTS " + cartes.getNoOfAtoms() + " " + saConnects.length + " 0 0 0";
    saMDLContent[6] = "M  V30 BEGIN ATOM";

    // now the cartesian coordinates
    DecimalFormat format = new DecimalFormat("0.00000");
    String[] saAtoms = cartes.getAllAtomTypes();
    double[][] daCoords = cartes.getAllXYZCoord();
    for (int i = 7; i < cartes.getNoOfAtoms() + 7; i++) {
      saMDLContent[i] =
          "M  V30 "
              + (i - 6)
              + " "
              + saAtoms[i - 7]
              + " "
              + format.format(daCoords[0][i - 7] * BOHRTOANG)
              + " "
              + format.format(daCoords[1][i - 7] * BOHRTOANG)
              + " "
              + format.format(daCoords[2][i - 7] * BOHRTOANG)
              + " 0";
    }

    saMDLContent[cartes.getNoOfAtoms() + 7] = "M  V30 END ATOM";
    saMDLContent[cartes.getNoOfAtoms() + 8] = "M  V30 BEGIN BOND";

    // now append the bond information
    for (int i = 0; i < saConnects.length; i++) {
      saMDLContent[i + cartes.getNoOfAtoms() + 9] = saConnects[i];
    }

    // end the mdl format
    saMDLContent[saMDLContent.length - 3] = "M  V30 END BOND";
    saMDLContent[saMDLContent.length - 2] = "M  V30 END CTAB";
    saMDLContent[saMDLContent.length - 1] = "M  END";

    try {
      OutputPrimitives.writeOut(sGeometryFile, saMDLContent, false);
    } catch (Exception e) {
      throw new ConvergenceException(
          "Problem in writing geometry for openbabel input (local optimization).", e);
    }

    /*
     * call babel
     */
    String[] saOutput;
    String[] saErrOut;
    try {
      Runtime rt = Runtime.getRuntime();
      Process proc =
          rt.exec(
              new String[] {
                "obminimize", "-ff", forceField.name(), "-n", "" + iNoOfIterations, sGeometryFile
              });

      // any error message?
      StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error? THIS NOT DOABLE with openbabel unfortunately... :-(
      proc.waitFor();
      saOutput = outputGobbler.getData();
      saErrOut = errorGobbler.getData();

    } catch (Exception e) {
      e.printStackTrace(System.err);
      throw new ConvergenceException("openbabel has a problem (local optimization).", e);
    }

    /*
     * Read the output
     */
    try {
      cartes =
          TranslateOutToCartes(
              cartes.getNoOfAtoms(),
              cartes.getNoOfMolecules(),
              cartes.getAllAtomsPerMol(),
              saOutput,
              saErrOut,
              cartes.getAllAtomTypes());
    } catch (Exception e) {
      throw new ConvergenceException("Failure in reading output of openbabel.", e);
    }

    /*
     * clean up
     */
    try {
      ManipulationPrimitives.remove(sGeometryFile);
    } catch (Exception e) {
      System.err.println("Problem cleaning openbabel file up." + e.toString());
    }

    /*
     * return it
     */
    cartes.setAllCharges(faCharges);
    cartes.setAllSpins(iaSpins);

    return cartes;
  }

  private CartesianCoordinates TranslateOutToCartes(
      int iNoOfAtoms,
      int iNoOfMolecules,
      int[] iaAtsPerMol,
      String[] saOutput,
      String[] saErrOut,
      String[] saAtoms)
      throws CastException {

    // generate the new cartesian
    final CartesianCoordinates cartes =
        new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaAtsPerMol);
    cartes.setAllAtomTypes(saAtoms);

    // beginning of geometry is typically at pos 2
    int geomCount = 2;
    int atomCount = 0;
    int geomEnd = iNoOfAtoms + 2;

    final double[] coords = new double[3];
    while (geomCount < geomEnd) {
      // loop over lines
      final String s = saOutput[geomCount].trim();
      if (s.equalsIgnoreCase("no hcounts")) {
        geomEnd++;
        geomCount++;
        continue;
      }

      // actual line
      final String[] sa = s.split("\\s+");
      int coordOffset = 5;
      // try to figure the coordOffset out
      try {
        final int i = Integer.parseInt(sa[4]);
      } catch (Exception e) {
        // Exception perfectly normal, ignore. May happen in case of e.g.
        // amino acids which somehow get an extra column
        coordOffset = 6;
      }
      try {
        final double p = Double.parseDouble(sa[1]);
        coordOffset = 1;
      } catch (Exception e) {
        // may happen
      }
      try {
        coords[0] = Double.parseDouble(sa[coordOffset]) * ANGTOBOHR;
        coords[1] = Double.parseDouble(sa[coordOffset + 1]) * ANGTOBOHR;
        coords[2] = Double.parseDouble(sa[coordOffset + 2]) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException(e);
      }

      cartes.setXYZCoordinatesOfAtom(coords, atomCount);

      geomCount++;
      atomCount++;
    }

    // the energy is one row above the "convergence reached" info
    final String sTemp = saErrOut[saErrOut.length - 3];
    // cut the first columns off
    final String[] sa = sTemp.split("\\s+");

    double dEnergy = FixedValues.NONCONVERGEDENERGY;
    try {
      dEnergy = Double.parseDouble(sa[1]) * KJTOHARTREE;
    } catch (Exception e) {
      throw new CastException("Couldn't cast the energy of babel.", e);
    }
    cartes.setEnergy(dEnergy);

    return cartes;
  }

  private void doConnects(final BondInfo bonds, final int noAtoms) {

    final List<String> conns = new ArrayList<>(noAtoms);
    int connCount = 0;
    for (int i = 0; i < noAtoms; i++) {
      for (int j = i + 1; j < noAtoms; j++) {
        if (DEBUG) {
          System.out.println("DEBUG: bonds array: " + i + "  " + j + "  " + bonds.hasBond(i, j));
        }
        if (bonds.hasBond(i, j)) {
          // there is a bond aka connect
          connCount++;
          // assume single bond as we for the time being have not more info
          final String sC = "M  V30 " + connCount + " 1 " + (i + 1) + " " + (j + 1);
          conns.add(sC);
        }
      }
    }

    // transfer the connect arraylist to the string array
    saConnects = new String[conns.size()];

    for (int i = 0; i < saConnects.length; i++) {
      saConnects[i] = conns.get(i);
      if (DEBUG) {
        System.out.println("DEBUG: Connections: " + saConnects[i]);
      }
    }
  }
}
