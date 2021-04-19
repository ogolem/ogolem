/*
Copyright (c) 2010-2014, J. M. Dieterich
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
import java.util.ArrayList;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Calls the MD parts of Tinker for a pseudo local optimization.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
final class TinkerMDCaller extends AbstractLocOpt {

  private static final long serialVersionUID = (long) 20101215;

  static enum METHOD {
    pss,
    anneal,
    monte,
    sniffer
  };

  private METHOD whichProg;
  private final String mdProg;
  private final boolean useSolvation;

  private final String[] saTinkerSecondHalf;
  private final String sWhichParameters;

  TinkerMDCaller(GlobalConfig conf, METHOD whichMethod, boolean useSolv) {
    super(conf);
    this.whichProg = whichMethod;

    switch (whichMethod) {
      case pss:
        mdProg = "pss";
        break;
      case anneal:
        mdProg = "anneal";
        break;
      case monte:
        mdProg = "monte";
        break;
      case sniffer:
        mdProg = "sniffer";
        break;
      default:
        mdProg = "pss";
        break;
    }

    this.useSolvation = useSolv;

    final boolean[][] bonds = conf.geoConf.bonds.getFullBondMatrix();
    saTinkerSecondHalf = new String[bonds.length + 1];
    saTinkerSecondHalf[0] = " ";

    String[] saAuxInput;
    final String sWhichAuxFile = conf.outputFolder + "-tinkermd.aux";
    try {
      saAuxInput = Input.ReadFile(sWhichAuxFile);
    } catch (Exception e) {
      System.err.println("Error in reading the mandatory aux-file!" + e.toString());
      throw new RuntimeException("Error in reading mandatory aux-file in TinkerMDCaller.", e);
    }
    if (!saAuxInput[0].trim().equalsIgnoreCase("###OGOLEMAUX###")) {
      System.err.println("ERROR: The read in auxiliary file is not valid!");
      throw new RuntimeException(
          "Auxiliary input file for Tinker not valid. Not starting with OGOLEMAUX ID.");
    }

    // read in, which parameters are to be taken
    sWhichParameters = saAuxInput[1].trim();

    for (int i = 2; i < bonds.length + 2; i++) {
      // get the force field parameter code in every single line
      saAuxInput[i] = saAuxInput[i].trim();
      saAuxInput[i] = saAuxInput[i].substring(saAuxInput[i].indexOf(" "));
      saAuxInput[i] = saAuxInput[i].trim();
      saTinkerSecondHalf[i - 1] = saAuxInput[i] + "\t";
    }

    // now we act on the boolean[][] bond information
    for (int i = 0; i < bonds.length; i++) {
      for (int j = 0; j < bonds.length; j++) {
        if (i == j) {
          // not specified in tinker
        } else {
          if (bonds[i][j]) {
            // bond, add that to the connectivity info
            int k = j + 1;
            saTinkerSecondHalf[i + 1] = saTinkerSecondHalf[i + 1] + k + "\t";
          } else {
            // no bond, go on
          }
        }
      }
    }
  }

  private TinkerMDCaller(TinkerMDCaller orig) {
    super(orig);
    this.whichProg = orig.whichProg;
    this.mdProg = orig.mdProg;
    this.useSolvation = orig.useSolvation;
    this.saTinkerSecondHalf = orig.saTinkerSecondHalf.clone();
    this.sWhichParameters = orig.sWhichParameters;
  }

  @Override
  public TinkerMDCaller copy() {
    return new TinkerMDCaller(this);
  }

  @Override
  public String myIDandMethod() {
    return "Tinker MD: " + mdProg;
  }

  @Override
  public String myID() {
    return "Tinker MD";
  }

  @Override
  public Molecule localOptimization(Molecule mStartMolecule) {

    System.out.println(
        "WARNING: Tinker MD caller does not support molecular "
            + "optimization. Returning non-optimized molecule.");

    mStartMolecule.setEnergy(FixedValues.NONCONVERGEDENERGY);

    return mStartMolecule;
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      long lID,
      CartesianCoordinates cartes,
      final boolean[][] baConstraints,
      final boolean isConstricted,
      final BondInfo bonds)
      throws ConvergenceException {

    final String folder = "tinkermd" + lID;

    final float[] faCharges = cartes.getAllCharges();
    final short[] iaSpins = cartes.getAllSpins();

    // create the folder
    try {
      OutputPrimitives.createAFolder(folder);
    } catch (Exception e) {
      throw new ConvergenceException(e);
    }

    // copy the parameter file
    try {
      OutputPrimitives.createLink(
          sWhichParameters, folder + System.getProperty("file.separator") + sWhichParameters);
    } catch (Exception e) {
      throw new ConvergenceException(e);
    }

    // create the xyz file
    final String sTinkerInput =
        folder + System.getProperty("file.separator") + "tinker" + lID + ".xyz";
    final String[] saTinkerInp = new String[cartes.getNoOfAtoms() + 1];

    saTinkerInp[0] = cartes.getNoOfAtoms() + " created by OGOLEM";

    // create the first half
    String[] saAtoms = cartes.getAllAtomTypes();
    double[][] daXYZ = cartes.getAllXYZCoord();
    for (int i = 1; i < cartes.getNoOfAtoms() + 1; i++) {
      saTinkerInp[i] =
          "\t"
              + i
              + "\t"
              + saAtoms[i - 1]
              + "\t"
              + daXYZ[0][i - 1] * Constants.BOHRTOANG
              + "\t"
              + daXYZ[1][i - 1] * Constants.BOHRTOANG
              + "\t"
              + daXYZ[2][i - 1] * Constants.BOHRTOANG
              + "\t";
    }

    for (int i = 1; i < cartes.getNoOfAtoms() + 1; i++) {
      saTinkerInp[i] += saTinkerSecondHalf[i];
    }

    // write it out
    try {
      OutputPrimitives.writeOut(sTinkerInput, saTinkerInp, false);
    } catch (Exception e) {
      throw new ConvergenceException("Error whilst writing the manipulated tinker input out", e);
    }

    // create the tinker key file
    final String sTinkerKeyword =
        folder + System.getProperty("file.separator") + "tinker" + lID + ".key";
    final ArrayList<String> alKeywords = new ArrayList<>();
    alKeywords.add("PARAMETERS " + sWhichParameters);

    if (useSolvation) {
      // XXX which solvation to use?
      alKeywords.add("SOLVATE ASP");
    }

    for (int i = 0; i < baConstraints[0].length; i++) {
      if (baConstraints[0][i] && baConstraints[1][i] && baConstraints[2][i]) {
        alKeywords.add(
            "RESTRAIN-POSITION "
                + (i + 1)
                + " "
                + daXYZ[0][i] * Constants.BOHRTOANG
                + "  "
                + daXYZ[1][i] * Constants.BOHRTOANG
                + "  "
                + daXYZ[2][i] * Constants.BOHRTOANG
                + "  10000");
        // we use a higher (aka steeper) potential here to turn the restraint more into a constraint
      } else if (!baConstraints[0][i] && !baConstraints[1][i] && !baConstraints[2][i]) {
        continue;
      } else {
        System.err.println(
            "WARNING: Tinker caller only allows complete constraints/restraints on an atom.");
      }
    }

    final String[] saKeywordContent = new String[alKeywords.size()];
    for (int i = 0; i < alKeywords.size(); i++) {
      saKeywordContent[i] = alKeywords.get(i);
    }

    try {
      OutputPrimitives.writeOut(sTinkerKeyword, saKeywordContent, false);
    } catch (Exception e) {
      throw new ConvergenceException("Failure in writing the tinker key file.", e);
    }

    /*
     * call the specified tinker program
     */
    String[] saTinkerStream;
    try {
      Runtime rt = Runtime.getRuntime();

      Process proc = null;
      String s = System.getenv("OGO_TINKERMDPATH");
      String s2 = System.getenv("OGO_TINKERMDOPT");
      // fallbacks
      if (s == null || s2 == null) {
        s = "monte";
        s2 = " 1000 C 3.0 500 0.01";
        whichProg = METHOD.monte;
      } else {
        s += mdProg;
        // pad with whitespace
        s2 = " " + s2;
      }

      proc = rt.exec(s + " " + sTinkerInput + s2);

      // any error message?
      final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      final int iExitValue = proc.waitFor();
      if (iExitValue != 0) {
        throw new ConvergenceException("Tinker MD returns non-zero return value.");
      } else {
        // tinker should(!) have completed normally...
        saTinkerStream = outputGobbler.getData();
      }

    } catch (Exception e) {
      throw new ConvergenceException("Tinker MD has a problem (local optimization).", e);
    }

    if (saTinkerStream.length == 0) {
      System.err.println("Tinker MD stream does not contain any info.");
      throw new ConvergenceException();
    }

    // grep the energy first, is easiest...
    final double energy = readEnergy(saTinkerStream);

    // now read the new cartesian coordinates, slightly more difficult...
    String outGeoFile = null;
    if (whichProg == METHOD.pss || whichProg == METHOD.anneal) {
      final File fold = new File(folder);
      final String[] fileList = fold.list();
      int lastNo = -1;
      for (String file : fileList) {
        // get the number, take care of the case that it is no number
        final String[] part = file.split("\\.");
        try {
          int i = Integer.parseInt(part[1]);
          lastNo = Math.max(lastNo, i);
        } catch (Exception e) {
          // happens for e.g. .xyz .key and so on, just fine
          continue;
        }
      }

      if (lastNo < 0) throw new ConvergenceException("No snapshot found.");

      // now we know the filename
      outGeoFile = folder + System.getProperty("file.separator") + "tinker" + lID + "." + lastNo;
    } else if (whichProg == METHOD.monte || whichProg == METHOD.sniffer) {
      outGeoFile = folder + System.getProperty("file.separator") + "tinker" + lID + ".xyz_2";
    }

    if (outGeoFile == null) throw new ConvergenceException("No file known.");

    try {
      cartes =
          Input.ReadXYZTinker(
              outGeoFile,
              cartes.getNoOfAtoms(),
              cartes.getNoOfMolecules(),
              cartes.getAllAtomsPerMol());
    } catch (Exception e) {
      throw new ConvergenceException("Problem in reading the output of tinker.", e);
    }
    // due to the wonderful tinker output, this is absolutely required
    cartes.setAllAtomTypes(saAtoms);

    cartes.setEnergy(energy);

    /*
     * clean up
     */
    try {
      ManipulationPrimitives.remove(folder);
    } catch (Exception e) {
      System.err.println("Problem cleaning tinker files up." + e.toString());
    }

    /*
     * return it
     */
    cartes.setAllCharges(faCharges);
    cartes.setAllSpins(iaSpins);

    return cartes;
  }

  private double readEnergy(final String[] stream) throws ConvergenceException {
    double energy = FixedValues.NONCONVERGEDENERGY;

    if (whichProg == METHOD.pss) {
      // we try reading pss output
      for (int i = stream.length - 1; i >= 0; i--) {
        final String s = stream[i].trim();
        if (s.startsWith("Final Function Value")) {
          final String[] sa = s.split("\\s+");
          try {
            energy = Double.parseDouble(sa[6]) * Constants.KCALTOHARTREE;
          } catch (Exception e) {
            throw new ConvergenceException("Couldn't parse energy.", e);
          }
          break;
        }
      }
    } else if (whichProg == METHOD.anneal) {
      // we try reading anneal output
      for (int i = stream.length - 1; i >= 0; i--) {
        final String s = stream[i].trim();
        if (s.startsWith("Total Energy")) {
          final String[] sa = s.split("\\s+");
          try {
            energy = Double.parseDouble(sa[2]) * Constants.KCALTOHARTREE;
          } catch (Exception e) {
            throw new ConvergenceException("Couldn't parse energy.", e);
          }
          break;
        }
      }

    } else if (whichProg == METHOD.monte) {
      // we try reading monte output
      for (int i = stream.length - 1; i >= 0; i--) {
        final String s = stream[i].trim();
        if (s.startsWith("Global Minimum Energy Value :")) {
          final String[] sa = s.split("\\s+");
          try {
            energy = Double.parseDouble(sa[5]) * Constants.KCALTOHARTREE;
          } catch (Exception e) {
            throw new ConvergenceException("Couldn't parse energy.", e);
          }
          break;
        }
      }
    } else if (whichProg == METHOD.sniffer) {
      // we try reading sniffer output
      for (int i = stream.length - 1; i >= 0; i--) {
        final String s = stream[i].trim();
        if (s.startsWith("Final Function Value")) {
          final String[] sa = s.split("\\s+");
          try {
            energy = Double.parseDouble(sa[4]) * Constants.KCALTOHARTREE;
          } catch (Exception e) {
            throw new ConvergenceException("Couldn't parse energy.", e);
          }
          break;
        }
      }
    } else {
      System.err.println(
          "ERROR: Unknown method choice for Tinker MD. I do "
              + "not know how to handly the streaming output. :-(");
    }

    return energy;
  }
}
