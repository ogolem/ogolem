/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2012, J. M. Dieterich
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
package org.ogolem.switches;

import static org.ogolem.core.Constants.*;

import java.util.Iterator;
import java.util.LinkedList;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.StreamGobbler;
// physical constants

/**
 * Uses tinker for the local optimization of a set of cartesian coordinates.
 *
 * @author Johannes Dieterich
 * @version 2012-06-24
 */
final class TinkerLocOpt implements LocalOptimization {

  // for the time being we are restricted to mm3, but that might of course change in the future
  private static final String sWhichParameters = "mm3.prm";

  private final double dBlowBondFac;

  private final String sLocProgram;

  private final String sTinkerOptions;

  private boolean bParamsExist = false;

  private final boolean bDebug;

  TinkerLocOpt(final int iWhichMethod, final double dBlowFactorBonds, final boolean bDebugThis) {
    this.dBlowBondFac = dBlowFactorBonds;
    if (iWhichMethod == 0) {
      this.sLocProgram = "minimize";
      this.sTinkerOptions = " 0.1";
    } else if (iWhichMethod == 1) {
      this.sLocProgram = "newton";
      this.sTinkerOptions = " a a 0.01";
    } else if (iWhichMethod == 2) {
      this.sLocProgram = "optimize";
      this.sTinkerOptions = " 0.01";
    } else {
      this.sLocProgram = "minimize";
      this.sTinkerOptions = " 0.1";
    }
    bDebug = bDebugThis;
  }

  @Override
  public CartesianCoordinates locOptThis(final CartesianCoordinates startCartes, final int iID) {

    // TODO since this is not working
    if (true) {
      System.err.println("THIS IS NOT WORKING (TINKER LCOOPT)");
      return null;
    }

    /*
     * first we need to figure out what atomic types we have, for that we
     * need to create the bonding information
     */
    final boolean[][] baBonds =
        LittleHelpers.bondingInfo(startCartes, dBlowBondFac).getFullBondMatrix();

    if (bDebug) {
      System.out.println("DEBUG: Operating with bonding blow factor: " + dBlowBondFac);
      // we print the bonding info out
      System.out.println("DEBUG: Bonding info coming from TinkerLocOpt:");
      for (int i = 0; i < baBonds.length; i++) {
        String sBonds = " ";
        for (int j = 0; j < baBonds[0].length; j++) {
          sBonds += baBonds[i][j] + "  ";
        }
        System.out.println("DEBUG: " + sBonds);
      }
    }

    final int[][] iaConnectivities = assignConnects(startCartes.getAllAtomTypes(), baBonds);

    final int[] iaAtomTypes = assignAtomTypes(startCartes.getAllAtomTypes(), iaConnectivities);

    /*
     * create the input and put parameter files in the correct spot
     */
    final String sTinkerInput = "tinker" + iID + ".xyz";
    final String sTinkerBasis = "tinker" + iID;
    final String sTinkerKeyword = "tinker" + iID + ".key";

    final String[] saTinkerInp = new String[startCartes.getNoOfAtoms() + 1];

    saTinkerInp[0] = startCartes.getNoOfAtoms() + " created by OGOLEM";

    // create the first half
    final String[] saAtoms = startCartes.getAllAtomTypes();
    final double[][] daXYZ = startCartes.getAllXYZCoordsCopy();
    for (int i = 1; i < startCartes.getNoOfAtoms() + 1; i++) {
      saTinkerInp[i] =
          "\t"
              + i
              + "\t"
              + saAtoms[i - 1]
              + "\t"
              + daXYZ[0][i - 1]
              + "\t"
              + daXYZ[1][i - 1]
              + "\t"
              + daXYZ[2][i - 1]
              + "\t";
    }

    // create the second half, atom type and connectivity
    for (int i = 1; i < startCartes.getNoOfAtoms() + 1; i++) {
      saTinkerInp[i] += "\t" + iaAtomTypes[i - 1];
      for (final int iConn : iaConnectivities[i - 1]) {
        saTinkerInp[i] += "\t" + (iConn + 1);
      }
    }

    // write it out
    try {
      Output.printMiscToFile(sTinkerInput, saTinkerInp);
    } catch (Exception e) {
      System.err.println("WARNING: Couldn't write tinker inputfile.");
      e.printStackTrace(System.err);
      return null;
    }

    /*
     * create the keyword file for tinker
     */
    String[] saKeywordContent = new String[1];
    saKeywordContent[0] = "PARAMETERS " + sWhichParameters;
    try {
      Output.printMiscToFile(sTinkerKeyword, saKeywordContent);
    } catch (Exception e) {
      System.err.println("WARNING: Couldn't write tinker keyfile.");
      e.printStackTrace(System.err);
      return null;
    }

    /*
     * check whether the parameter file is present
     */
    if (!bParamsExist) {
      final boolean bParams = SwitchesInput.isFilePresents(sWhichParameters);
      if (bParams) {
        bParamsExist = true;
      } else {
        System.err.println("ERROR: No parameter file present. Aborting the locopt now.");
        return null;
      }
    }

    /*
     * debug, anyone?
     */
    if (bDebug) {
      final String sCartesianFile = "debug_geometry" + iID + ".xyz";
      final String[] saCartes = startCartes.createPrintableCartesians();

      final String sParamsFile = "debug_atomtypes" + iID + ".prm";
      final String[] saParams = new String[iaAtomTypes.length];
      for (int i = 0; i < saParams.length; i++) {
        saParams[i] = " " + iaAtomTypes[i];
      }

      try {
        Output.printMiscToFile(sCartesianFile, saCartes);
        Output.printMiscToFile(sParamsFile, saParams);
      } catch (Exception e) {
        System.err.println(
            "WARNING: Couldn't write out either cartesians or atom types for "
                + "debugging."
                + e.toString());
      }
    }

    /*
     * call tinker
     */
    String[] saTinkerStream;
    try {
      Runtime rt = Runtime.getRuntime();
      Process proc = rt.exec(sLocProgram + " " + sTinkerInput + sTinkerOptions);

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
        System.err.println("WARNING: Tinker returns non-zero return value.");
        return null;
      } else {
        // tinker should(!) have completed normally...
        saTinkerStream = outputGobbler.getData();
      }

    } catch (Exception e) {
      System.err.println("WARNING: Tinker has a problem (local optimization). " + e.toString());
      return null;
    }

    if (saTinkerStream.length == 0) {
      System.err.println("Tinker stream does not contain any info.");
      return null;
    }

    /*
     * read tinkers output
     */
    CartesianCoordinates endCartes = new CartesianCoordinates(startCartes);
    String sTinkerOutput = sTinkerInput;
    try {
      endCartes =
          SwitchesInput.readXYZTinker(
              sTinkerOutput + "_2",
              endCartes.getNoOfAtoms(),
              endCartes.getNoOfMolecules(),
              endCartes.getAllAtomsPerMol());
    } catch (Exception e) {
      System.err.println("Problem in reading the output of tinker." + e.toString());
      return null;
    }

    /*
     * get the energy
     */
    double dEnergy = 1E7;
    for (int i = 0; i < saTinkerStream.length; i++) {
      if (saTinkerStream[i].contains("Final Function Value :")) {
        String sTemp;
        sTemp = saTinkerStream[i];
        sTemp = sTemp.substring(sTemp.indexOf(":") + 1);
        sTemp = sTemp.trim();
        try {
          dEnergy = Double.parseDouble(sTemp) * KCALTOHARTREE;
        } catch (Exception e) {
          System.err.println("Failure to cast the energy of tinker." + e.toString());
          return null;
        }
        break;
      } else if (i == saTinkerStream.length - 1) {
        // there is no line saying "Final Function value"
        System.err.println("WARNING: Tinker didn't converge!");
        return null;
      }
    }
    endCartes.setEnergy(dEnergy);

    /*
     * clean up
     */
    if (!bDebug) {
      try {
        SwitchesInput.removeTinkerFiles(sTinkerBasis);
      } catch (Exception e) {
        System.err.println("WARNING: Problem cleaning tinker files up." + e.toString());
      }
    }

    /*
     * return it
     */
    return endCartes;
  }

  private static int[][] assignConnects(final String[] saAtoms, final boolean[][] baBonds) {

    final int[][] iaConnects = new int[saAtoms.length][];

    for (int i = 0; i < saAtoms.length; i++) {

      // collect connectivity
      final LinkedList<Integer> llTempConn = new LinkedList<>();

      final boolean[] baPartBonds = baBonds[i];

      for (int j = 0; j < baPartBonds.length; j++) {
        final boolean bCurrBond = baPartBonds[j];
        if (bCurrBond && j != i) {
          llTempConn.add(j);
        }
      }

      // now put connectivity in a usable form
      final Iterator<Integer> itConns = llTempConn.iterator();
      final int[] iaConns = new int[llTempConn.size()];

      int it = 0;
      while (itConns.hasNext()) {
        iaConns[it] = itConns.next();
        it++;
      }

      iaConnects[i] = iaConns;
    }

    return iaConnects;
  }

  private static int[] assignAtomTypes(final String[] saAtoms, final int[][] iaConnects) {

    final int[] iaAtomTypes = new int[saAtoms.length];

    // now assign the correct atom type (in mm3 context)
    AtomTypeFinder:
    for (int i = 0; i < saAtoms.length; i++) {
      final String sCurrAtom = saAtoms[i];
      final int[] iaCurrConns = iaConnects[i];

      if (sCurrAtom.equalsIgnoreCase("H")) {
        /*
         * possibilities: 5, 21, 23, 24, 28, 44, 48, 73, 124
         */
        if (iaCurrConns.length == 1) {
          // this SHOULD be the connectivity of the H atom ;-)
          if (!saAtoms[iaCurrConns[0]].equalsIgnoreCase("N")
              && !saAtoms[iaCurrConns[0]].equalsIgnoreCase("O")
              && !saAtoms[iaCurrConns[0]].equalsIgnoreCase("S")) {
            // easy atom type
            // standard hydrogen
            iaAtomTypes[i] = 5;
            // ignores others (enol: 73, acetylene: 124)
            continue AtomTypeFinder;
          } else if (saAtoms[iaCurrConns[0]].equalsIgnoreCase("O")) {
            // OH alcohol hydrogen
            iaAtomTypes[i] = 21;
            continue AtomTypeFinder;
          } else if (saAtoms[iaCurrConns[0]].equalsIgnoreCase("N")) {
            // NH amine/imine hydrogen
            iaAtomTypes[i] = 23;
            // ignores amide (28) and ammonium (48)!
            continue AtomTypeFinder;
          } else if (saAtoms[iaCurrConns[0]].equalsIgnoreCase("S")) {
            // SH thiol hydrogen
            iaAtomTypes[i] = 44;
            continue AtomTypeFinder;
          }
        } else {
          System.err.println(
              "ERROR: We are detection a "
                  + iaCurrConns.length
                  + " bonded hydrogen. Using dummy. Handle with care.");
          iaAtomTypes[i] = 112;
          continue AtomTypeFinder;
        }
      } else if (sCurrAtom.equalsIgnoreCase("C")) {
        /*
         * possibilities: 1, 2, 3, 4, 22, 29(not needed), 30(not needed),
         * 38, 50, 56, 57, 58, 67, 68, 71, 106, 113(not needed),
         * 114(not needed), 160, 161, 162
         */
        if (iaCurrConns.length == 4) {
          /*
           * remaining possibilities: 1, 22 (cyclopronane, not needed), 56 (cyclobutane,
           * not needed)
           */
          iaAtomTypes[i] = 1;
          continue AtomTypeFinder;
        } else if (iaCurrConns.length == 3) {
          /*
           * remaining possibilities: 2, 3, 29 (C radical, not needed), 38 (cyclopropene,
           * not needed), 50, 57 (cyclobutene, not needed), 58 (cyclobutanone, not needed),
           * 67 (cyclopropanone, not needed), 71 (ketonium carbon, not needed),
           * 113 (CH ferrocene, not needed), 114 (CC ferrocene, not needed),
           * 160/161/162 (nucleic acid, not needed)
           */
          if (saAtoms[iaCurrConns[0]].equalsIgnoreCase("O")) {
            if (saAtoms[iaCurrConns[iaCurrConns[0]]].length() == 1) {
              // the oxygen is really just connected to a C
              iaAtomTypes[i] = 3;
              continue AtomTypeFinder;
            }
          } else if (saAtoms[iaCurrConns[1]].equalsIgnoreCase("O")) {
            if (saAtoms[iaCurrConns[iaCurrConns[1]]].length() == 1) {
              // the oxygen is really just connected to a C
              iaAtomTypes[i] = 3;
              continue AtomTypeFinder;
            }
          } else if (saAtoms[iaCurrConns[2]].equalsIgnoreCase("O")) {
            if (saAtoms[iaCurrConns[iaCurrConns[2]]].length() == 1) {
              // the oxygen is really just connected to a C
              iaAtomTypes[i] = 3;
              continue AtomTypeFinder;
            }
          }

          // REMARK: 50 doesn't work since it is not flexible enough (not all torsions available)
          if (saAtoms[iaCurrConns[0]].equalsIgnoreCase("C")
                  && saAtoms[iaCurrConns[1]].equalsIgnoreCase("C")
              || saAtoms[iaCurrConns[0]].equalsIgnoreCase("C")
                  && saAtoms[iaCurrConns[2]].equalsIgnoreCase("C")
              || saAtoms[iaCurrConns[1]].equalsIgnoreCase("C")
                  && saAtoms[iaCurrConns[2]].equalsIgnoreCase("C")) {
            iaAtomTypes[i] = 2;
            continue AtomTypeFinder;
          }
          // TODO 2,3, 50 distinguish
        } else if (iaCurrConns.length == 2) {
          /*
           * remaining possibilities: 4, 30 (carbonium ion, not needed), 68 (allene,
           * not needed), 106 (ketene, not needed),
           */
          iaAtomTypes[i] = 4;
          continue AtomTypeFinder;
        } else if (iaCurrConns.length == 1) {
          /*
           * remaining possibilities: nothing -> dummy
           */
          iaAtomTypes[i] = 112;
          continue AtomTypeFinder;
        } else {
          System.err.println(
              "ERROR: We are detecting carbon with "
                  + iaCurrConns.length
                  + " bonds. Using dummy. Handle with care.");
          iaAtomTypes[i] = 112;
          continue AtomTypeFinder;
        }
      } else if (sCurrAtom.equalsIgnoreCase("O")) {
        /*
         * possibilities: 6, 7, 41, 47, 49, 69, 70, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
         * 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
         * 115, 116, 117, 118, 119, 120, 121, 145, 148, 149, 159
         */
        iaAtomTypes[i] = 6;
        continue AtomTypeFinder;
        // TODO we just are "sure" that ATM there is just that kind
      } else if (sCurrAtom.equalsIgnoreCase("N")) {
        /*
         * possibilities: 8, 9, 10, 39, 40, 43, 45, 46, 72, 107, 108, 109, 110, 111, 143,
         * 144, 146, 150, 155, 164
         */
        iaAtomTypes[i] = 9;
        continue AtomTypeFinder;
        // REMARK: 107 doesn't work for azo since there are not enough torsions
        // TODO we just assume Nsp2 ATM
      } else if (sCurrAtom.equalsIgnoreCase("S")) {
        /*
         * possibilities: 15, 16, 17, 18, 42, 74, 104, 105, 154
         */
        iaAtomTypes[i] = 15;
        continue AtomTypeFinder;
        // TODO there is nothing like this ATM
      } else {
        /*
         * since for the time being our limiting factor is the mndo
         * OMx parametrization, the above atom types should do ATM
         */
        System.err.println(
            "ERROR: We do not know the atom type " + sCurrAtom + " using dummy. Handle with care.");
        System.err.println(
            "Please contact the author, he might be willing to fix this for you. :-)");
        continue AtomTypeFinder;
      }
    }

    return iaAtomTypes;
  }
}
