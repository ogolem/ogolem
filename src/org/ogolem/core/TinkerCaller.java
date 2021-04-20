/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
              2017-2021, J. M. Dieterich and B. Hartke
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

import java.util.ArrayList;
import org.ogolem.helpers.Tuple;
import org.ogolem.io.OutputPrimitives;

/**
 * This calls the tinker program suite for local optimization. Since force fields in general suck
 * badly when it comes to automatic generation of input files, the user is requested to specify a
 * file named @code{YOURJOBNAME.aux} in which the choice of force field atom types is specified for
 * all atoms in the system. This might sound (and is) a little inconvenient but it is so far the
 * best choice. Besides, there needs to be the actual parameter file in the work directory. Also,
 * Tinker can be used as a backend. But this is highly discouraged for normal application as the
 * calling overhead is significant in comparison to the force field execution time!
 *
 * @author Johannes Dieterich, Bernd Hartke
 * @version 2021-04-19
 */
class TinkerCaller extends AbstractLocOpt implements CartesianFullBackend {
  // TODO delete always all files...
  // the ID
  private static final long serialVersionUID = (long) 20130712;
  private static final boolean DEBUG = false;

  static enum METHOD {
    minimize,
    optimize,
    newton,
    custom
  };

  private final String sLocProgram;

  private final String sTinkerOptions;

  private final boolean useSolvation;

  private String[] saTinkerSecondHalf;

  private String sWhichParameters;

  private final boolean useFileOut;

  private final boolean uberCustom;

  private final int noAtomsCluster;
  private final int noAtomsEnvironment;
  private final int noAtomsTotal;
  private final boolean envIsRigid;

  /**
   * The constructor for usage as a local optimizing engine.
   *
   * @param globconf A complete configuration set.
   */
  TinkerCaller(
      GlobalConfig globconf,
      METHOD whichMethod,
      boolean solvation,
      boolean fileout,
      boolean customKey,
      boolean customBonds) {
    super(globconf);
    switch (whichMethod) {
      case minimize:
        sLocProgram = "minimize";
        sTinkerOptions = " 0.1";
        break;
      case newton:
        sLocProgram = "newton";
        sTinkerOptions = " a a 0.01";
        break;
      case optimize:
        sLocProgram = "optimize";
        sTinkerOptions = " 0.01";
        break;
      case custom:
        sLocProgram = null;
        sTinkerOptions = null;
        break;
      default:
        // should never end up here
        sLocProgram = null;
        sTinkerOptions = null;
        break;
    }
    this.useSolvation = solvation;
    this.useFileOut = fileout;
    this.uberCustom = customKey;
    this.sWhichParameters = null;

    final Geometry gTmp = new Geometry(globconf.geoConf);

    this.noAtomsCluster = gTmp.getNumberOfAtoms();
    this.noAtomsEnvironment = (gTmp.containsEnvironment()) ? gTmp.getEnvironment().atomsInEnv() : 0;
    this.noAtomsTotal = noAtomsCluster + noAtomsEnvironment;
    this.envIsRigid =
        (gTmp.containsEnvironment()) ? gTmp.getEnvironment().isEnvironmentRigid() : false;

    saTinkerSecondHalf = new String[noAtomsTotal + 1];
    saTinkerSecondHalf[0] = " ";

    String[] saAuxInput;
    final String sWhichAuxFile = globconf.outputFolder + "-tinker.aux";

    try {
      saAuxInput = Input.ReadFile(sWhichAuxFile);
    } catch (Exception e) {
      System.err.println("Error in reading the mandatory aux-file!" + e.toString());
      saAuxInput = null;
    }

    if (saAuxInput == null) {
      throw new Error("Mandatory aux file could not be parsed");
    }

    if (!saAuxInput[0].trim().equalsIgnoreCase("###OGOLEMAUX###")) {
      System.err.println("ERROR: The read in auxiliary file is not valid!");
      throw new RuntimeException(
          "Auxiliary input file for Tinker not valid. Not starting with OGOLEMAUX ID.");
    }

    // read in, which parameters are to be taken
    sWhichParameters = saAuxInput[1].trim();

    if (customBonds && customKey) {
      for (int i = 2; i < noAtomsTotal + 2; i++) {
        final String[] sa = saAuxInput[i].trim().split("\\s", 2);
        saTinkerSecondHalf[i - 1] = sa[1];
      }
    } else {
      for (int i = 2; i < noAtomsTotal + 2; i++) {
        // get the force field parameter code in every single line
        final String[] sa = saAuxInput[i].trim().split("\\s", 2);
        saTinkerSecondHalf[i - 1] = sa[1] + "\t";
      }

      // once more ugly java code! ;-)
      CartesianCoordinates cartes;
      BondInfo bonds;
      if (gTmp.containsEnvironment()) {
        final Tuple<CartesianCoordinates, BondInfo> tup =
            gTmp.getCartesiansAndBondsWithEnvironment();
        cartes = tup.getObject1();
        bonds = tup.getObject2();
      } else {
        cartes = gTmp.getCartesians();
        bonds = gTmp.getBondInfo();
      }
      final boolean[][] bondMat = bonds.getFullBondMatrix();
      final String[] atoms = cartes.getAllAtomTypes();

      // now we act on the boolean[][] bond information
      for (int i = 0; i < bondMat.length; i++) {

        if (atoms[i].equalsIgnoreCase("DM")) {
          // real dummy for TIP5P
          continue;
        }

        for (int j = 0; j < bondMat.length; j++) {

          if (atoms[j].equalsIgnoreCase("DM")) {
            // real dummy for TIP5P
            continue;
          }

          if ((atoms[j].equalsIgnoreCase("LP") && !atoms[i].equalsIgnoreCase("O"))
              || (atoms[j].equalsIgnoreCase("M") && !atoms[i].equalsIgnoreCase("O"))
              || (atoms[i].equalsIgnoreCase("LP") && !atoms[j].equalsIgnoreCase("O"))
              || (atoms[i].equalsIgnoreCase("M") && !atoms[j].equalsIgnoreCase("O"))) {
            // dummy atom for TIP4P/5P and no connection to the O
            continue;
          }

          // not specified in tinker
          if (i == j) {
            continue;
          }
          if (bondMat[i][j]) {
            // bond, add that to the connectivity info
            int k = j + 1;
            saTinkerSecondHalf[i + 1] += k + "\t";
          }
        }
      }
    }
  }

  private TinkerCaller(final TinkerCaller orig) {
    super(orig);
    // shallow copy should be enough
    this.sLocProgram = orig.sLocProgram;
    this.sTinkerOptions = orig.sTinkerOptions;
    this.sWhichParameters = orig.sWhichParameters;
    this.useSolvation = orig.useSolvation;
    this.useFileOut = orig.useFileOut;
    this.uberCustom = orig.uberCustom;
    if (orig.saTinkerSecondHalf != null) {
      this.saTinkerSecondHalf = orig.saTinkerSecondHalf.clone();
    }
    this.noAtomsCluster = orig.noAtomsCluster;
    this.noAtomsEnvironment = orig.noAtomsEnvironment;
    this.noAtomsTotal = orig.noAtomsTotal;
    this.envIsRigid = orig.envIsRigid;
  }

  @Override
  public TinkerCaller copy() {
    return new TinkerCaller(this);
  }

  @Override
  public String myIDandMethod() {
    if (sLocProgram == null) return "Tinker: custom";
    else return "Tinker: " + sLocProgram;
  }

  @Override
  public String myID() {
    if (sLocProgram == null) return "Tinker: custom";
    else return "Tinker: " + sLocProgram;
  }

  @Override
  public Molecule localOptimization(Molecule mStartMolecule) {
    // this caller is not capable of molecular optimization due to the tinker aux file format
    System.out.println(
        "WARNING: Tinker caller can't locally optimize a molecule. Returning non-optimized molecule.");
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

    String sTinkerInput = "tinker" + lID + ".xyz";
    String sTinkerBasis = "tinker" + lID;
    String sTinkerKeyword = "tinker" + lID + ".key";

    float[] faCharges = cartes.getAllCharges();
    short[] iaSpins = cartes.getAllSpins();

    String[] saTinkerInp = new String[cartes.getNoOfAtoms() + 1];

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

    if (!uberCustom) {
      /*
       * create the keyword file for tinker, restraints used
       */
      ArrayList<String> alKeywords = new ArrayList<>();
      alKeywords.add("PARAMETERS " + sWhichParameters);

      if (useSolvation) {
        // XXX which solvation to use?
        alKeywords.add("SOLVATE ASP");
      }

      if (envIsRigid) {
        // mark the environment atoms as inactive
        alKeywords.add("INACTIVE -" + (this.noAtomsCluster + 1) + " " + this.noAtomsTotal);
      }

      for (int i = 0; i < this.noAtomsCluster; i++) {
        if (baConstraints[0][i] && baConstraints[1][i] && baConstraints[2][i]) {
          alKeywords.add("INACTIVE " + (i + 1));
        } else if (!baConstraints[0][i] && !baConstraints[1][i] && !baConstraints[2][i]) {
          continue;
        } else {
          System.err.println(
              "WARNING: Tinker caller only allows complete constraints/restraints on an atom.");
        }
      }

      String[] saKeywordContent = new String[alKeywords.size()];
      for (int i = 0; i < alKeywords.size(); i++) {
        saKeywordContent[i] = alKeywords.get(i);
      }

      try {
        OutputPrimitives.writeOut(sTinkerKeyword, saKeywordContent, false);
      } catch (Exception e) {
        throw new ConvergenceException("Failure in writing the tinker key file.", e);
      }
    }

    /*
     * call the specified tinker program
     */
    String[] saTinkerStream;
    try {
      Runtime rt = Runtime.getRuntime();

      Process proc;
      if (sLocProgram == null) {
        // then we call a custom command with custom options
        String s = System.getenv("OGO_TINKERCMD");
        String s2 = System.getenv("OGO_TINKEROPT");
        // fallbacks
        if (s == null) s = "minimize";
        if (s2 == null) {
          s2 = " 0.1";
        } else {
          // pad with whitespace
          s2 = " " + s2;
        }
        proc = rt.exec(s + " " + sTinkerInput + s2);
      } else {
        proc = rt.exec(sLocProgram + " " + sTinkerInput + sTinkerOptions);
      }

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
        throw new ConvergenceException("Tinker returns non-zero return value.");
      } else {
        // tinker should(!) have completed normally...
        saTinkerStream = outputGobbler.getData();
      }

    } catch (Exception e) {
      throw new ConvergenceException("Tinker has a problem (local optimization).", e);
    }

    if (useFileOut) {
      try {
        saTinkerStream = Input.ReadFile(sTinkerBasis + ".out");
      } catch (Exception e) {
        throw new ConvergenceException("Couldn't read out file set through fileout.", e);
      }
    }

    if (saTinkerStream.length == 0) {
      System.err.println("Tinker stream does not contain any info.");
      throw new ConvergenceException();
    }

    /*
     * read tinkers output
     */
    String sTinkerOutput = sTinkerInput;
    try {
      cartes =
          Input.ReadXYZTinker(
              sTinkerOutput + "_2",
              cartes.getNoOfAtoms(),
              cartes.getNoOfMolecules(),
              cartes.getAllAtomsPerMol(),
              cartes.getReferenceEnvironmentCopy());
    } catch (Exception e) {
      throw new ConvergenceException("Problem in reading the output of tinker.", e);
    }
    // due to the wonderful tinker output, this is absolutely required
    cartes.setAllAtomTypes(saAtoms);

    /*
     * get the energy
     */
    double dEnergy = FixedValues.NONCONVERGEDENERGY;
    for (int i = saTinkerStream.length - 1; i >= 0; i--) {
      if (saTinkerStream[i].contains("Final Function Value :")) {
        final String[] sa = saTinkerStream[i].trim().split("\\s+");
        try {
          dEnergy = Double.parseDouble(sa[4]) * Constants.KCALTOHARTREE;
        } catch (Exception e) {
          if (DEBUG) {
            System.out.println("DEBUG: output");
            for (final String s : saTinkerStream) {
              System.out.println(s);
            }
          }
          throw new ConvergenceException("Failure to cast the energy of tinker.", e);
        }
        break;
      }
    }

    if (dEnergy == FixedValues.NONCONVERGEDENERGY) {
      if (DEBUG) {
        if (saTinkerStream.length > 6) {
          for (int i = saTinkerStream.length - 6; i < saTinkerStream.length; i++) {
            System.err.println(saTinkerStream[i]);
          }
        }
      }
      // there was no line saying "Final Function value"
      throw new ConvergenceException("Tinker didn't converge!");
    }

    cartes.setEnergy(dEnergy);

    /*
     * clean up
     */
    try {
      Input.RemoveTinkerFiles(sTinkerBasis);
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

  @Override
  public String getMethodID() {
    return "tinker analyze with: " + sWhichParameters;
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
      final BondInfo bonds,
      final Gradient gradient,
      final boolean hasRigidEnv) {
    final Gradient numGrad =
        NumericalGradients.numericalGradient(
            lID,
            iIteration,
            xyz1D,
            saAtomTypes,
            atomNos,
            atsPerMol,
            energyparts,
            iNoOfAtoms,
            faCharges,
            iaSpins,
            bonds,
            this,
            hasRigidEnv);
    gradient.copyDataIn(numGrad);
  }

  @Override
  public double energyCalculation(
      long lID,
      int iIteration,
      double[] xyz1D,
      String[] saAtoms,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      int iNoOfAtoms,
      float[] faCharges,
      short[] iaSpins,
      final BondInfo bonds,
      final boolean hasRigidEnv) {

    final String sTinkerBasis = "tinker" + lID + "_" + iIteration;
    final String sTinkerInput = sTinkerBasis + ".xyz";
    final String sTinkerKeyword = sTinkerBasis + ".key";

    final String[] saTinkerInp = new String[iNoOfAtoms + 1];

    saTinkerInp[0] = iNoOfAtoms + " created by OGOLEM";

    // create the first half
    for (int i = 1; i < iNoOfAtoms + 1; i++) {
      saTinkerInp[i] =
          "\t"
              + i
              + "\t"
              + saAtoms[i - 1]
              + "\t"
              + xyz1D[(i - 1)] * Constants.BOHRTOANG
              + "\t"
              + xyz1D[iNoOfAtoms + (i - 1)] * Constants.BOHRTOANG
              + "\t"
              + xyz1D[2 * iNoOfAtoms + (i - 1)] * Constants.BOHRTOANG
              + "\t";
    }

    for (int i = 1; i < iNoOfAtoms + 1; i++) {
      saTinkerInp[i] += saTinkerSecondHalf[i];
    }

    // write it out
    try {
      OutputPrimitives.writeOut(sTinkerInput, saTinkerInp, false);
    } catch (Exception e) {
      e.printStackTrace(System.err);
      cleanUp(sTinkerBasis);
      return FixedValues.NONCONVERGEDENERGY;
    }

    if (!uberCustom) {
      /*
       * create the keyword file for tinker, restraints used
       */
      ArrayList<String> alKeywords = new ArrayList<>();
      alKeywords.add("PARAMETERS " + sWhichParameters);

      if (useSolvation) {
        // XXX which solvation to use?
        alKeywords.add("SOLVATE ASP");
      }

      String[] saKeywordContent = new String[alKeywords.size()];
      for (int i = 0; i < alKeywords.size(); i++) {
        saKeywordContent[i] = alKeywords.get(i);
      }

      try {
        OutputPrimitives.writeOut(sTinkerKeyword, saKeywordContent, false);
      } catch (Exception e) {
        e.printStackTrace(System.err);
        cleanUp(sTinkerBasis);
        return FixedValues.NONCONVERGEDENERGY;
      }
    }

    /*
     * call the specified tinker program
     */
    String[] saTinkerStream;
    try {
      final Runtime rt = Runtime.getRuntime();

      final Process proc = rt.exec("analyze " + sTinkerInput + " E");

      // any error message or output?
      final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");
      final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      final int exitVal = proc.waitFor();
      if (exitVal != 0) {
        throw new ConvergenceException("Tinker returns non-zero return value.");
      } else {
        // tinker should(!) have completed normally...
        saTinkerStream = outputGobbler.getData();
      }

    } catch (Exception e) {
      e.printStackTrace(System.err);
      cleanUp(sTinkerBasis);
      return FixedValues.NONCONVERGEDENERGY;
    }

    if (useFileOut) {
      try {
        saTinkerStream = Input.ReadFile(sTinkerBasis + ".out");
      } catch (Exception e) {
        e.printStackTrace(System.err);
        cleanUp(sTinkerBasis);
        return FixedValues.NONCONVERGEDENERGY;
      }
    }

    if (saTinkerStream.length == 0) {
      System.err.println("Tinker stream does not contain any info.");
      cleanUp(sTinkerBasis);
      return FixedValues.NONCONVERGEDENERGY;
    }

    /*
     * get the energy
     */
    double dEnergy = FixedValues.NONCONVERGEDENERGY;
    for (int i = saTinkerStream.length - 1; i >= 0; i--) {
      if (saTinkerStream[i].trim().startsWith("Total Potential Energy :")) {
        final String[] sa = saTinkerStream[i].trim().split("\\s+");
        try {
          dEnergy = Double.parseDouble(sa[4]) * Constants.KCALTOHARTREE;
        } catch (Exception e) {
          dEnergy = FixedValues.NONCONVERGEDENERGY;
        }
        break;
      }
    }

    if (dEnergy == FixedValues.NONCONVERGEDENERGY) {
      if (DEBUG) {
        if (saTinkerStream.length > 6) {
          for (int i = saTinkerStream.length - 6; i < saTinkerStream.length; i++) {
            System.err.println(saTinkerStream[i]);
          }
        }
      }
      // there was no line saying "Total Potential Energy"
    }

    /*
     * clean up
     */
    cleanUp(sTinkerBasis);

    return dEnergy;
  }

  @Override
  public CartesianFullBackend getBackend() {
    return new TinkerCaller(this);
  }

  private static void cleanUp(final String tinkerBasis) {
    try {
      Input.RemoveTinkerFiles(tinkerBasis);
    } catch (Exception e) {
      System.err.println("Problem cleaning tinker files up." + e.toString());
    }
  }
}
