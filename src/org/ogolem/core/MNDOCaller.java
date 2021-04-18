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

import static org.ogolem.core.Constants.BOHRTOANG;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Calls the MNDO package for local optimizations using semiempirical methods.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class MNDOCaller extends AbstractLocOpt {

  // the ID
  private static final long serialVersionUID = (long) 20150421;

  static enum METHOD {
    MNDOd,
    OM3,
    PM3,
    OM2,
    OM1,
    AM1,
    MNDOC,
    MNDO,
    MINDO3,
    CNDO2,
    SCCDFTB,
    SCCDFTBJorgensen
  };

  private final METHOD whichMethod;

  private final boolean bCOSMOWater;

  private int[] iaAtomicNosGeom;

  /**
   * The constructor for usage as a local optimizing engine.
   *
   * @param globconf A complete configuration set.
   */
  MNDOCaller(final GlobalConfig globconf, METHOD whichMethod, boolean bWaterSolv) {
    super(globconf);
    this.whichMethod = whichMethod;
    this.bCOSMOWater = bWaterSolv;
  }

  private MNDOCaller(final MNDOCaller orig) {
    super(orig);
    this.whichMethod = orig.whichMethod;
    this.bCOSMOWater = orig.bCOSMOWater;
  }

  @Override
  public MNDOCaller copy() {
    return new MNDOCaller(this);
  }

  @Override
  public String myID() {
    return "MNDO";
  }

  @Override
  public String myIDandMethod() {

    String method = "";
    switch (whichMethod) {
      case MNDOd:
        method = "MNDO/d";
        break;
      case OM3:
        method = "OM3";
        break;
      case PM3:
        method = "PM3";
        break;
      case OM2:
        method = "OM2";
        break;
      case OM1:
        method = "OM1";
        break;
      case AM1:
        method = "AM1";
        break;
      case MNDOC:
        method = "MNDOC";
        break;
      case MNDO:
        method = "MNDO";
        break;
      case MINDO3:
        method = "MINDO/3";
        break;
      case CNDO2:
        method = "CNDO/2";
        break;
      case SCCDFTB:
        method = "SCC-DFTB";
        break;
      case SCCDFTBJorgensen:
        method = "SCC-DFTB w/ Jorgensen-corr.";
        break;
      default:
        break;
    }

    if (bCOSMOWater) {
      method += " with COSMO water";
    }

    return "MNDO: " + method;
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      long id,
      CartesianCoordinates cartes,
      boolean[][] constraints,
      boolean isConstricted,
      final BondInfo bonds)
      throws Exception {

    final boolean bMolecule = (cartes.getNoOfMolecules() == 1);

    final String sMNDOInput = "mndo" + id + ".inp";
    final String sMNDOOutput = "mndo" + id + ".out";
    final String sMNDOBasis = "mndo" + id;

    final float[] faCharges = cartes.getAllCharges();
    final short[] iaSpins = cartes.getAllSpins();

    // get the atomic numbers together
    int[] iaCurrAtNos = new int[faCharges.length];
    if (!bMolecule && iaAtomicNosGeom == null) {
      // get them fresh
      final String[] saAtoms = cartes.getAllAtomTypes();

      for (int i = 0; i < saAtoms.length; i++) {
        iaCurrAtNos[i] = AtomicProperties.giveAtomicNumber(saAtoms[i]);
      }

      iaAtomicNosGeom = iaCurrAtNos;
    } else if (!bMolecule && iaAtomicNosGeom != null) {
      iaCurrAtNos = iaAtomicNosGeom;
    } else if (bMolecule) {
      // get them fresh
      final String[] saAtoms = cartes.getAllAtomTypes();
      iaCurrAtNos = new int[saAtoms.length];
      for (int i = 0; i < saAtoms.length; i++) {
        iaCurrAtNos[i] = AtomicProperties.giveAtomicNumber(saAtoms[i]);
      }
    } else {
      System.err.println("ERROR: In MNDO atomic numbers. Please notify the author.");
    }

    /*
     * create a folder for this run
     */
    try {
      OutputPrimitives.createAFolder(sMNDOBasis);
    } catch (Exception e) {
      throw new ConvergenceException(
          "WARNING: Problem in creating folder for mndo input (local optimization).", e);
    }

    /*
     * write the output aka input for the to be called program
     */
    try {
      writeMNDOInput(
          sMNDOBasis,
          sMNDOInput,
          whichMethod,
          cartes.getAllXYZCoord(),
          iaCurrAtNos,
          cartes.getTotalCharge(),
          cartes.getTotalSpin(),
          bCOSMOWater);
    } catch (Exception e) {
      throw new ConvergenceException(
          "WARNING: Problem in writing geometry for MNDO input (local optimization).", e);
    }

    /*
     * call mndo
     */
    Process proc;
    try {
      Runtime rt = Runtime.getRuntime();
      String sCommand = "mndo.sh " + sMNDOInput + " " + sMNDOOutput;
      String[] saEnvp = new String[0];
      File dir = new File(sMNDOBasis);
      proc = rt.exec(sCommand, saEnvp, dir);

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
        throw new ConvergenceException(
            "WARNING: MNDO returns non-zero return value (local optimization).");
      } else {
        // mndo should(!) have completed normally...
      }

    } catch (Exception e) {
      throw new ConvergenceException("WARNING: MNDO has a problem (local optimization).", e);
    }

    /*
     * read mndo's output in
     */
    String[] saOutput;
    try {
      saOutput = Input.ReadFile(sMNDOBasis + System.getProperty("file.separator") + sMNDOOutput);
    } catch (Exception e) {
      throw new ConvergenceException("WARNING: Couldn't read in MNDO output.", e);
    }

    /*
     * translate mndo's output
     */
    try {
      cartes =
          readCartes(
              saOutput,
              cartes.getNoOfAtoms(),
              cartes.getNoOfMolecules(),
              cartes.getAllAtomsPerMol(),
              cartes.getAllAtomTypes());
    } catch (Exception e) {
      throw new ConvergenceException("WARNING: Problem in translating the output of MNDO.", e);
    }

    /*
     * clean up
     */
    try {
      ManipulationPrimitives.remove(sMNDOBasis);
    } catch (Exception e) {
      System.err.println("WARNING: Problem cleaning MNDO files up." + e.toString());
    }

    /*
     * return it
     */
    cartes.setAllCharges(faCharges);
    cartes.setAllSpins(iaSpins);

    return cartes;
  }

  private CartesianCoordinates readCartes(
      final String[] saOutput,
      final int iNoOfAtoms,
      final int iNoOfMolecules,
      final int[] iaAtsPerMol,
      final String[] saAtoms)
      throws CastException {

    final CartesianCoordinates cartes =
        new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaAtsPerMol);
    cartes.setAllAtomTypes(saAtoms);

    double dEnergy = 0.0;
    boolean bEnergyFound = false;
    int iGeomStart = 0;
    boolean bGeomFound = false;

    for (int i = saOutput.length - 1; i >= 0; i--) {
      final String sTemp = saOutput[i].trim();
      if (sTemp.startsWith("FINAL HEAT OF FORMATION")) {
        String sTemp2 = sTemp.substring(24).trim();
        // now the end off
        sTemp2 = sTemp2.substring(0, sTemp2.indexOf(" ")).trim();
        try {
          dEnergy = Double.parseDouble(sTemp2) * Constants.KCALTOHARTREE;
        } catch (Exception e) {
          throw new CastException("WARNING: Couldn't cast energy", e);
        }
        bEnergyFound = true;
      } else if (sTemp.startsWith("FINAL CARTESIAN GRADIENT NORM")) {
        iGeomStart = i + 8;
        bGeomFound = true;
      }

      if (bEnergyFound && bGeomFound) {
        break;
      }
    }

    if (!bEnergyFound) {
      throw new CastException("WARNING: Couldn't find energy.");
    }

    if (iGeomStart == 0) {
      throw new CastException("WARNING: Couldn't find optimized geometry.");
    }

    // put the energy in
    cartes.setEnergy(dEnergy);

    // then the coordinates (translate angstrom to bohr)
    for (int i = 0; i < iNoOfAtoms; i++) {
      String sTemp = saOutput[i + iGeomStart].trim();
      final double[] daXYZ = new double[3];

      // first get rid of the two columns in the beginning
      sTemp = sTemp.substring(sTemp.indexOf(" ")).trim();
      sTemp = sTemp.substring(sTemp.indexOf(" ")).trim();

      // then get the coordinates
      for (int j = 0; j < 3; j++) {
        final int iIndex = sTemp.indexOf(" ");
        final String sCurrCoord = sTemp.substring(0, iIndex).trim();
        sTemp = sTemp.substring(iIndex).trim();

        // try casting it
        try {
          daXYZ[j] = Double.parseDouble(sCurrCoord) * Constants.ANGTOBOHR;
        } catch (Exception e) {
          throw new CastException("WARNING: Couldn't cast coordinate.", e);
        }

        // now get rid of the star (just for the first two needed)
        if (j < 2) {
          sTemp = sTemp.substring(sTemp.indexOf(" ")).trim();
        }
      }

      cartes.setXYZCoordinatesOfAtom(daXYZ, i);
    }
    return cartes;
  }

  private static void writeMNDOInput(
      final String sFolder,
      final String sInputFile,
      final METHOD iWhichMethod,
      final double[][] daXYZ,
      final int[] iaAtomicNos,
      final int iCharge,
      final int iSpin,
      final boolean bCOSMOWater)
      throws IOException {

    // method ID to MNDO method ID
    int iMethodID;
    switch (iWhichMethod) {
      case MNDOd:
        iMethodID = -10;
        break;
      case OM3:
        iMethodID = -8;
        break;
      case PM3:
        iMethodID = -7;
        break;
      case OM2:
        iMethodID = -6;
        break;
      case OM1:
        iMethodID = -5;
        break;
      case AM1:
        iMethodID = -2;
        break;
      case MNDOC:
        iMethodID = -1;
        break;
      case MNDO:
        iMethodID = 0;
        break;
      case MINDO3:
        iMethodID = 1;
        break;
      case CNDO2:
        iMethodID = 2;
        break;
      case SCCDFTB:
        iMethodID = 5;
        break;
      case SCCDFTBJorgensen:
        iMethodID = 6;
        break;
      default:
        iMethodID = -8;
        System.err.println(
            "ERROR: MNDO method to method translation. " + "Contact the author. Using OM3 now.");
    }

    final String[] saInput = new String[iaAtomicNos.length + 5];
    saInput[0] = "iop=" + iMethodID + " igeom=1 jop=0 iform=1 mplib=0 ";

    // water solvation using the COSMO modell
    if (bCOSMOWater) {
      saInput[0] += "icosmo=1 +";
    } else {
      // no solvation
      saInput[0] += "+";
    }

    /*
     * MNDO's definition of spin is (unfortunately...) slightly different
     * from normal.
     * 1 means two single occupied orbitals forming in total a singlet
     * 0 means a normal singlet state
     * the rest is as usual
     * WE DO NOT AT THE MOMENT USE THE SPECIAL SINGLET FEATURE!
     */
    int iMult;
    if (iSpin == 0) {
      iMult = 0;
    } else {
      iMult = iSpin + 1;
    }
    saInput[1] = "kharge=" + iCharge + " imult=" + iMult;
    saInput[2] = "Automatically created by OGOLEM";
    saInput[3] = " ";

    // coordinates
    final DecimalFormat form = new DecimalFormat("0.00000000");
    for (int i = 0; i < iaAtomicNos.length; i++) {
      saInput[i + 4] =
          " "
              + iaAtomicNos[i]
              + "   "
              + form.format(daXYZ[0][i] * BOHRTOANG)
              + " 1 "
              + form.format(daXYZ[1][i] * BOHRTOANG)
              + " 1 "
              + form.format(daXYZ[2][i] * BOHRTOANG)
              + " 1 ";
    }

    // and an empty line in the end is important...
    saInput[saInput.length - 1] = "";

    // write it out
    final String sMNDOInp = sFolder + System.getProperty("file.separator") + sInputFile;
    OutputPrimitives.writeOut(sMNDOInp, saInput, false);
  }
}
