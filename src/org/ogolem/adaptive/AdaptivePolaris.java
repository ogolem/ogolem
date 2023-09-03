/*
Copyright (c) 2010-2014, J. M. Dieterich
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
package org.ogolem.adaptive;

import static org.ogolem.core.Constants.KCALTOHARTREE;

import java.io.File;
import java.util.ArrayList;
import java.util.Locale;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.StreamGobbler;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Makes OGOLEM go nuclear. Highly specific implementation.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
final class AdaptivePolaris extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20150727;
  private static final String[] DONNEE_HEADER = {
    "*******************************************",
    "*                                         *",
    "*    Fichier de donnees pour Polaris      *",
    "*                                         *",
    "*     Types d'analyse disponibles         *",
    "*                                         *",
    "*    sgp : Single point                   *",
    "*    tra :  analyse d une trajectoire     *",
    "*    rel : convergence sur dynamique      *",
    "*    nco : interaction a Ncorps           *",
    "*    opt : optimisation                   *",
    "*    vib : frequence de vibration         *",
    "*    dyE : dynamique a energie constante  *",
    "*    dyT : dynamique a T constante        *",
    "*    QUD : analyse de quenching           *",
    "*                                         *",
    "*******************************************",
    "",
    "Type de calcul : nco       ",
    "",
    "---------------------------------",
    "|   repertoire champ de force   |",
    "---------------------------------",
    ""
  };

  private static final String[] DONNEE_MIDDLE = {
    "",
    "--------------------------------------------",
    "|   Liste des sous-elements a considerer   |",
    "--------------------------------------------",
    "",
    "MRES",
    "WION",
    "",
    "-----------------------------",
    "|   Repertoire de travail   |",
    "-----------------------------",
    ""
  };

  private final String[] filesToCopy;

  private final String[] donneeData;

  private final String[] metData;

  private final String[] repData;

  private final boolean bUseTotalEnergy;

  private long callCounter = 0;

  AdaptivePolaris(final boolean fitTotalEnergy) {

    bUseTotalEnergy = fitTotalEnergy;

    final String sRefParams1 = "MET.FFM";
    final String sRefParams2 = "REP.FFM";
    final String sInpReference = "donnee-ref.aux";

    // which files are to be copied?
    final String sAuxFile = "adaptive-polaris.aux";
    filesToCopy = readFilesToCopy(sAuxFile);

    // read the file-stubs in
    donneeData = readDoineeRef(sInpReference);
    metData = readRefIn("OGOPARAM", sRefParams1);
    repData = readRefIn("OGOPARAM", sRefParams2);
  }

  AdaptivePolaris(final AdaptivePolaris orig) {
    filesToCopy = orig.filesToCopy.clone();
    donneeData = orig.donneeData.clone();
    metData = orig.metData.clone();
    repData = orig.repData.clone();
    bUseTotalEnergy = orig.bUseTotalEnergy;
  }

  @Override
  public AdaptivePolaris copy() {
    return new AdaptivePolaris(this);
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {

    callCounter++;

    final String sFolder = "polaris" + params.getID() + "_" + callCounter;
    final String sTotalFolder =
        System.getProperty("user.dir") + System.getProperty("file.separator") + sFolder;
    final String sMetFile = sFolder + System.getProperty("file.separator") + "MET.FFM";
    final String sRepFile = sFolder + System.getProperty("file.separator") + "REP.FFM";
    final String sDoneeFile = sFolder + System.getProperty("file.separator") + "donnee";

    // create the folder
    try {
      OutputPrimitives.createAFolder(sFolder);
    } catch (Exception e) {
      System.err.println("ERROR: Couldn't create folder. " + e.toString());
      return FixedValues.NONCONVERGEDENERGY;
    }

    // copy the files
    for (final String file : filesToCopy) {
      try {
        final String newFile = sFolder + System.getProperty("file.separator") + file;
        OutputPrimitives.copyFile(file, newFile);
      } catch (Exception e) {
        System.err.println("ERROR: Couldn't copy file. " + e.toString());
      }
    }

    // copy the specific geometry files
    try {
      OutputPrimitives.copyFile(
          "GEOM.xyz_" + geomID, sFolder + System.getProperty("file.separator") + "GEOM.xyz");
      OutputPrimitives.copyFile(
          "GEOM.sle_" + geomID, sFolder + System.getProperty("file.separator") + "GEOM.sle");
    } catch (Exception e) {
      System.err.println("ERROR: Couldn't copy geometry file. " + e.toString());
    }

    // glue the donee and the other files
    final String[] saMetRepData = glueMetRepTogether(metData, repData, params.getAllParamters());
    final String[] saDoneeData = glueDonneeTogether(donneeData, sTotalFolder, "GEOM");

    // write the files
    try {
      OutputPrimitives.writeOut(sMetFile, saMetRepData[0], false);
      OutputPrimitives.writeOut(sRepFile, saMetRepData[1], false);
      OutputPrimitives.writeOut(sDoneeFile, saDoneeData, false);
    } catch (Exception e) {
      System.err.println("ERROR: Couldn't write data to file." + e.toString());
      return FixedValues.NONCONVERGEDENERGY;
    }

    // call the polaris command
    try {
      Runtime rt = Runtime.getRuntime();
      String sPolarisCmd = System.getenv("OGO_POLARISCMD");
      if (sPolarisCmd == null) {
        // option not exported, falling back to default
        sPolarisCmd = "Polaris";
      }
      File dir = new File(sFolder);
      Process proc = rt.exec(new String[] {sPolarisCmd}, null, dir);

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
        System.err.println("Polaris returns non-zero return value (local optimization).");
        return FixedValues.NONCONVERGEDENERGY;
      } else {
        // polaris should(!) have completed normally...
      }
    } catch (Exception e) {
      e.printStackTrace(System.err);
      System.err.println("Polaris has a problem (local optimization).");
      return FixedValues.NONCONVERGEDENERGY;
    }

    // read the energy
    final double dEnergy = grepEnergy(sFolder + System.getProperty("file.separator") + "GEOM.out");

    // delete the folder
    try {
      ManipulationPrimitives.remove(sFolder);
    } catch (Exception e) {
      System.err.println("ERROR: Couldn't remove folder. " + e.toString());
    }

    return dEnergy;
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds,
      final double[] grad) {

    System.out.println(
        "INFO: You are using a numerical gradient. This is "
            + "very slow und potentially has a bad convergence pattern. "
            + "Consider using no local optimization or a gradient free one.");

    // numerical gradient, but nothing else possible
    return NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, grad);
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      final ArrayList<CartesianCoordinates> refCartes, final String sMethod) {

    // the number of parameters can be directly accessed from the
    final int iNoOfParams = metData.length - 1 + repData.length - 1;
    final int[] iaParamsPerAt = {iNoOfParams};
    final String[] keys = {"unknown"};

    final AdaptiveParameters paramStub =
        new AdaptiveParameters(iNoOfParams, -1, keys, iaParamsPerAt, sMethod);

    return paramStub;
  }

  @Override
  public double[][] minMaxBordersForParams(final AdaptiveParameters params) {

    System.out.println(
        "WARNING: This interface is not to be used without "
            + "custom parameter borders. Using extremely big borders now. "
            + "This will be highly inefficient.");

    final int iNoOfParams = params.getNumberOfParamters();

    final double[][] daBorders = new double[2][iNoOfParams];

    for (int i = 0; i < iNoOfParams; i++) {
      daBorders[0][i] = -1000.0;
      daBorders[1][i] = 1000.0;
    }

    return daBorders;
  }

  private String[] readFilesToCopy(final String auxFile) {

    // read file in
    String[] sa;
    try {
      sa = Input.ReadFile(auxFile);
    } catch (Exception e) {
      System.err.println("ERROR: Cannot read aux file! " + e.toString());
      return new String[0];
    }

    if (!sa[0].trim().equalsIgnoreCase("###OGOLEMAUX###")) {
      System.err.println("ERROR: No valid aux file!");
      return new String[0];
    }

    final String[] toCopy = new String[sa.length - 1];
    for (int i = 0; i < toCopy.length; i++) {
      toCopy[i] = sa[i + 1];
    }

    return toCopy;
  }

  private String[] readDoineeRef(final String file) {

    // read file in
    String[] sa;
    try {
      sa = Input.ReadFile(file);
    } catch (Exception e) {
      System.err.println("ERROR: Cannot read file " + file + ". " + e.toString());
      return new String[0];
    }

    return sa;
  }

  private String[] readRefIn(final String paramKey, final String file) {

    // read file in
    String[] sa;
    try {
      sa = Input.ReadFile(file);
    } catch (Exception e) {
      System.err.println("ERROR: Cannot read file " + file + ". " + e.toString());
      return new String[0];
    }

    String s = "";
    for (String t : sa) {
      s += t + System.getProperty("line.separator");
    }

    // now chop
    final String[] sa2 = s.split(paramKey);

    return sa2;
  }

  private String[] glueDonneeTogether(
      final String[] donneeReference, final String sTotalFolder, final String sXYZ) {

    // a lot of copying required ;-)
    final String[] saCompl =
        new String[DONNEE_HEADER.length + DONNEE_MIDDLE.length + donneeReference.length + 9];

    int iOffset = 0;
    System.arraycopy(DONNEE_HEADER, 0, saCompl, iOffset, DONNEE_HEADER.length);
    iOffset += DONNEE_HEADER.length;
    saCompl[iOffset] = sTotalFolder;
    iOffset++;
    System.arraycopy(DONNEE_MIDDLE, 0, saCompl, iOffset, DONNEE_MIDDLE.length);
    iOffset += DONNEE_MIDDLE.length;
    saCompl[iOffset] = sTotalFolder;
    iOffset++;
    saCompl[iOffset] = "";
    iOffset++;
    saCompl[iOffset] = "--------------------------";
    iOffset++;
    saCompl[iOffset] = "|   Fichiers a analyser  |";
    iOffset++;
    saCompl[iOffset] = "--------------------------";
    iOffset++;
    saCompl[iOffset] = "";
    iOffset++;
    saCompl[iOffset] = sXYZ;
    iOffset++;
    saCompl[iOffset] = "";
    iOffset++;
    System.arraycopy(donneeReference, 0, saCompl, iOffset, donneeReference.length);

    return saCompl;
  }

  private String[] glueMetRepTogether(
      final String[] metReference, final String[] repReference, final double[] parameters) {

    final String[] data = new String[2];

    int iParamCounter = 0;
    // first MET.FFM
    data[0] = metReference[0];
    for (int i = 1; i < metReference.length; i++) {
      data[0] += String.format(Locale.US, "%10.4f", parameters[iParamCounter]) + metReference[i];
      iParamCounter++;
    }

    // then REP.FFM
    data[1] = repReference[0];
    for (int i = 1; i < repReference.length; i++) {
      data[1] += String.format(Locale.US, "%10.4f", parameters[iParamCounter]) + repReference[i];
      iParamCounter++;
    }

    return data;
  }

  private double grepEnergy(final String file) {

    // read file in
    String[] data;
    try {
      data = Input.ReadFile(file);
    } catch (Exception e) {
      System.err.println("ERROR: Couldn't read file in. " + e.toString());
      return FixedValues.NONCONVERGEDENERGY;
    }

    // loop through
    double energy = FixedValues.NONCONVERGEDENERGY;
    for (int i = data.length - 1; i >= 0; i--) {
      final String s = data[i].trim();
      if (bUseTotalEnergy) {
        if (s.contains("Energie totale")) {
          final String[] sa = s.split("\\s+");
          try {
            energy = Double.parseDouble(sa[4]) * KCALTOHARTREE;
            return energy;
          } catch (Exception e) {
            System.err.println("ERROR: Couldn't parse energy. " + e.toString());
            return FixedValues.NONCONVERGEDENERGY;
          }
        }
      } else {
        if (s.contains("E(interaction)")) {
          final String[] sa = s.split("\\s+");
          try {
            energy = Double.parseDouble(sa[2]) * KCALTOHARTREE;
            return energy;
          } catch (Exception e) {
            System.err.println("ERROR: Couldn't parse energy. " + e.toString());
            return FixedValues.NONCONVERGEDENERGY;
          }
        }
      }
    }

    System.err.println("ERROR: No energy found in Polaris output.");
    return energy;
  }
}
