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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * This calls the DFTBplus program for local optimizations.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class DFTBplusCaller extends AbstractLocOpt {

  // the ID
  private static final long serialVersionUID = (long) 20150421;

  static enum WHICHDRIVER {
    steepestdescent,
    conjugategradient
  };

  private final WHICHDRIVER whichDriver;
  private final boolean withSCC;
  private int[] atomNumbers = null;
  private String[] dftbAtomTypes = null;
  private final int maxIter;
  private final double maxForce;
  private final String cmd;
  private final boolean custom;
  private final String[] auxInp;

  /**
   * The constructor. Needs the global config to work.
   *
   * @param globConf the global configuration
   * @param whichMethod the method as an integer choice (valid:0/1/2/3/99, catch all to 0)
   */
  public DFTBplusCaller(
      final GlobalConfig globConf,
      final WHICHDRIVER whichMethod,
      final boolean doSCC,
      final boolean custom) {
    super(globConf);
    this.whichDriver = whichMethod;
    this.withSCC = doSCC;
    this.custom = custom;
    this.maxIter = globConf.maxIterLocOpt;
    this.maxForce = globConf.threshLocOptGradient;
    final String tmp = System.getenv("OGO_DFTBPCMD");
    if (tmp == null) {
      cmd = "dftb+";
    } else {
      cmd = tmp;
    }

    // read aux file if custom
    String[] tmpInp = null;
    if (custom) {
      try {
        final String auxFile = globConf.outputFolder + "-dftbp.aux";
        if (DEBUG) {
          System.out.println("DEBUG: Reading in dftbp aux file " + auxFile);
        }
        final String[] tmpIn = InputPrimitives.readFileIn(auxFile);
        if (DEBUG) {
          System.out.println("DEBUG: We have read: ");
          for (final String s : tmpIn) {
            System.out.println(s);
          }
        }
        String s = "";
        final String lineSep = System.getProperty("line.separator");
        for (final String sx : tmpIn) {
          s += sx + lineSep;
        }
        if (DEBUG) {
          System.out.println("DEBUG: Again, we read: ");
          System.out.println(s);
        }
        // lets split it
        final String[] tmpInpX = s.split("OGOCOORDS");
        if (DEBUG) {
          System.out.println("DEBUG: After the split: ");
          System.out.println(tmpInpX[0]);
          System.out.println("DEBUG: part 2");
          System.out.println(tmpInpX[1]);
        }
        tmpInp = new String[2];
        if (tmpInpX.length == 1) {
          // apparently the OGOCOORDS was in the first line
          tmpInp[0] = "";
          tmpInp[1] = tmpInpX[0];
        } else if (tmpInpX.length == 2) {
          tmpInp[0] = tmpInpX[0];
          tmpInp[1] = tmpInpX[1];
        } else {
          throw new Exception(
              "Your stub contains more than one OGOCOORDS identifier. We will let this break!");
        }
      } catch (Exception e) {
        System.err.println(
            "ERROR: Couldn't read in dftb+ aux file "
                + globConf.outputFolder
                + "-dftbp.aux"
                + ". This will break horribly!");
        e.printStackTrace(System.err);
        tmpInp = new String[2];
      }
    }
    this.auxInp = tmpInp;
  }

  private DFTBplusCaller(final DFTBplusCaller orig) {
    super(orig);
    this.whichDriver = orig.whichDriver;
    this.withSCC = orig.withSCC;
    this.maxIter = orig.maxIter;
    this.maxForce = orig.maxForce;
    this.cmd = orig.cmd;
    this.custom = orig.custom;
    this.auxInp = (orig.auxInp == null) ? null : orig.auxInp.clone();
  }

  @Override
  public DFTBplusCaller copy() {
    return new DFTBplusCaller(this);
  }

  @Override
  public String myID() {
    return "dftb+";
  }

  @Override
  public String myIDandMethod() {
    String driver = "";
    if (whichDriver == WHICHDRIVER.steepestdescent) {
      driver = "steepest descent";
    } else if (whichDriver == WHICHDRIVER.conjugategradient) {
      driver = "conjugate gradient";
    }

    if (withSCC) {
      driver += "w/ SCC";
    } else {
      driver += "w/o SCC";
    }

    return "dftb+ " + driver;
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      long lID,
      CartesianCoordinates cartes,
      boolean[][] baConstraints,
      boolean isConstricted,
      final BondInfo bonds)
      throws ConvergenceException {

    final String sFolder = "dftbp_" + lID;

    final float[] faCharges = cartes.getAllCharges();
    final short[] iaSpins = cartes.getAllSpins();

    /*
     * write the output aka input for the to be called program
     */
    try {
      WriteOutput(
          sFolder,
          cartes.getAllXYZCoord(),
          cartes.getAllAtomTypes(),
          cartes.getTotalCharge(),
          baConstraints);
    } catch (Exception e) {
      throw new ConvergenceException(
          "Problem in writing geometry for dftbp input (local optimization).", e);
    }

    /*
     * call dftb+
     */
    try {
      final Runtime rt = Runtime.getRuntime();
      final File dir = new File(sFolder);
      final Process proc = rt.exec(cmd, null, dir);

      // any error message?
      final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      final int exitValue = proc.waitFor();
      if (exitValue != 0) {
        throw new ConvergenceException("DFTB+ returns non-zero return value (local optimization).");
      } // else: dftb+ should(!) have completed normally...

    } catch (Exception e) {
      throw new ConvergenceException("DFTB+ has a problem (local optimization).", e);
    }

    /*
     * read dftb+'s output
     */
    try {
      cartes =
          Input.ReadXYZDFTBPlusOutput(
              sFolder,
              cartes.getNoOfAtoms(),
              cartes.getNoOfMolecules(),
              cartes.getAllAtomsPerMol(),
              iaSpins,
              faCharges);
    } catch (Exception e) {
      throw new ConvergenceException("Problem in reading the output of dftb+.", e);
    }

    /*
     * clean up
     */
    try {
      ManipulationPrimitives.remove(sFolder);
    } catch (Exception e) {
      System.err.println("Problem cleaning dftb+ files up." + e.toString());
    }

    /*
     * return it
     */
    return cartes;
  }

  private void WriteOutput(
      String sFolder,
      double[][] daXYZ,
      String[] saAtoms,
      int iTotalCharge,
      boolean[][] baConstraints)
      throws Exception {

    if (dftbAtomTypes == null) {
      // we create it fresh
      dftbAtomTypes = getAtomTypes(saAtoms);
      atomNumbers = getAtomNumbers(saAtoms, dftbAtomTypes);
    } // else: we assume it was already done once

    if (custom) {
      writeCustomInput(sFolder, daXYZ, atomNumbers, dftbAtomTypes);
      return;
    }

    try {
      WriteDFTBPInput(
          sFolder,
          daXYZ,
          atomNumbers,
          dftbAtomTypes,
          iTotalCharge,
          whichDriver,
          withSCC,
          baConstraints,
          maxIter,
          maxForce);
    } catch (Exception e) {
      throw e;
    }
  }

  private String[] getAtomTypes(final String[] atoms) {

    final List<String> types = new LinkedList<>();
    // always add the 0'th atom
    types.add(atoms[0]);
    for (int i = 1; i < atoms.length; i++) {
      if (!types.contains(atoms[i])) {
        types.add(atoms[i]);
      }
    }

    final String[] saTypes = new String[types.size()];
    for (int i = 0; i < types.size(); i++) {
      saTypes[i] = types.get(i);
    }

    return saTypes;
  }

  private int[] getAtomNumbers(final String[] atoms, final String[] atomTypes) {

    final int[] numbers = new int[atoms.length];

    for (int i = 0; i < atoms.length; i++) {
      for (int j = 0; j < atomTypes.length; j++) {
        if (atoms[i].equalsIgnoreCase(atomTypes[j])) {
          numbers[i] = j + 1;
          break;
        }
      }
    }

    return numbers;
  }

  private void writeCustomInput(
      final String folder, final double[][] xyz, final int[] whichAtoms, final String[] typeNames)
      throws Exception {

    final int noOfAtoms = whichAtoms.length;

    if (DEBUG) {
      System.out.println("DEBUG: I have the following two Strings for aux input.");
      System.out.println("DEBUG: STRING ONE:");
      System.out.println(auxInp[0]);
      System.out.println("DEBUG: STRING TWO:");
      System.out.println(auxInp[1]);
      System.out.println("DEBUG: STRINGS OVER: NUMBER OF ATOMS " + noOfAtoms);
    }

    final List<String> output = new ArrayList<>(noOfAtoms + 8);
    output.add(auxInp[0]);
    output.add("Geometry = {");
    String typeNameString = "TypeNames = {";
    for (final String typeName : typeNames) {
      typeNameString += "\"" + typeName + "\" ";
    }
    typeNameString += "}";
    output.add(typeNameString);

    // the coordinates
    output.add("TypesAndCoordinates [Bohr] = {");
    for (int i = 0; i < noOfAtoms; i++) {
      output.add(whichAtoms[i] + "\t" + xyz[0][i] + "\t" + xyz[1][i] + "\t" + +xyz[2][i]);
    }
    output.add("}");

    output.add("}");

    output.add(auxInp[1]);

    OutputPrimitives.createAFolder(folder, false);
    // put the input file into the folder
    final String inpFile = folder + System.getProperty("file.separator") + "dftb_in.hsd";
    OutputPrimitives.writeOut(inpFile, output, false);
  }

  private static void WriteDFTBPInput(
      final String folder,
      final double[][] xyz,
      final int[] whichAtoms,
      final String[] typeNames,
      final int totalCharge,
      final WHICHDRIVER whichDriver,
      final boolean withSCC,
      final boolean[][] constraints,
      final int maxIter,
      final double dMaxForce)
      throws IOException {

    final int noOfAtoms = whichAtoms.length;

    final List<String> constrInp = new LinkedList<>();
    for (int j = 0; j < noOfAtoms; j++) {
      if (constraints[0][j]) {
        constrInp.add("" + (j + 1) + " 1.0  0.0 0.0");
      }
      if (constraints[1][j]) {
        constrInp.add("" + (j + 1) + " 0.0  1.0 0.0");
      }
      if (constraints[2][j]) {
        constrInp.add("" + (j + 1) + " 0.0  0.0 1.0");
      }
    }

    int noOfConstr = constrInp.size();
    if (noOfConstr > 0) {
      // two additional lines for opening and closing tags
      noOfConstr += 2;
    }

    final String[] output = new String[noOfAtoms + 22 + typeNames.length + noOfConstr];
    output[0] = "Geometry = {";
    output[1] = "TypeNames = {";
    for (final String typeName : typeNames) {
      output[1] += "\"" + typeName + "\" ";
    }
    output[1] += "}";

    // the coordinates
    output[2] = "TypesAndCoordinates [Bohr] = {";
    for (int i = 3; i < noOfAtoms + 3; i++) {
      output[i] =
          whichAtoms[i - 3] + "\t" + xyz[0][i - 3] + "\t" + xyz[1][i - 3] + "\t" + +xyz[2][i - 3];
    }
    output[noOfAtoms + 3] = "}";

    output[noOfAtoms + 4] = "}";

    // driver
    if (whichDriver == WHICHDRIVER.steepestdescent) {
      output[noOfAtoms + 5] = "Driver = SteepestDescent {";
    } else if (whichDriver == WHICHDRIVER.conjugategradient) {
      output[noOfAtoms + 5] = "Driver = ConjugateGradient {";
    } else {
      output[noOfAtoms + 5] = "Driver = SteepestDescent {";
    }

    output[noOfAtoms + 6] = "MaxSteps = " + maxIter;

    output[noOfAtoms + 7] = "MaxForceComponent = " + dMaxForce;

    output[noOfAtoms + 8] = "MovedAtoms = 1:-1";

    int constrOffset = noOfAtoms + 9;
    if (noOfConstr > 0) {
      // opening tag
      output[constrOffset] = "Constraints = {";
      constrOffset++;

      for (final String s : constrInp) {
        output[constrOffset] = s;
        constrOffset++;
      }

      output[constrOffset] = "}";
    }

    output[noOfAtoms + 9 + noOfConstr] = "}";

    // hamiltonian
    output[noOfAtoms + 10 + noOfConstr] = "Hamiltonian = DFTB {";
    if (withSCC) {
      output[noOfAtoms + 11 + noOfConstr] = "SCC = yes";
    } else {
      output[noOfAtoms + 11 + noOfConstr] = "SCC = no";
    }

    // slater koster files
    output[noOfAtoms + 12 + noOfConstr] = "SlaterKosterFiles = Type2FileNames {";
    output[noOfAtoms + 13 + noOfConstr] =
        "Prefix = \""
            + System.getProperty("user.dir")
            + System.getProperty("file.separator")
            + "\" ";
    output[noOfAtoms + 14 + noOfConstr] = "Separator = \"-\" ";
    output[noOfAtoms + 15 + noOfConstr] = "Suffix = \".skf\"";
    output[noOfAtoms + 16 + noOfConstr] = "LowerCaseTypeName = No";
    output[noOfAtoms + 17 + noOfConstr] = "}";

    // charge
    output[noOfAtoms + 18 + noOfConstr] = "Charge = " + totalCharge + ".0";

    // maximum angular momentum
    output[noOfAtoms + 19 + noOfConstr] = "MaxAngularMomentum = {";
    for (int i = 0; i < typeNames.length; i++) {
      // XXX this might need a more advanced logic!
      String sRest = "\"p\"";
      if (typeNames[i].equalsIgnoreCase("H")) {
        sRest = "\"s\"";
      }
      output[i + 20 + noOfAtoms + noOfConstr] = typeNames[i] + " = " + sRest;
    }
    output[noOfAtoms + 20 + typeNames.length + noOfConstr] = "}";

    // don't forget closing }
    output[noOfAtoms + 21 + typeNames.length + noOfConstr] = "}";

    OutputPrimitives.createAFolder(folder, false);

    // put the input file into the folder
    final String sInputFile = folder + System.getProperty("file.separator") + "dftb_in.hsd";
    OutputPrimitives.writeOut(sInputFile, output, false);
  }
}
