/*
Copyright (c) 2012-2014, J. M. Dieterich
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
import java.text.DecimalFormat;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Interface to the VASP program package. Limited to selective optimizations and cubic cells.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class VASPCaller extends AbstractLocOpt {

  // the ID
  private static final long serialVersionUID = (long) 20121024;
  private final String whichFunctional;
  private final String[] copyableFiles;
  private final String[] poscarStub;
  private final double[][] cell;
  private final int whichReference;
  private final String refAtom;
  private final double[] refCoords;
  private final int[] order;
  private final int[] offPerType;
  private final boolean symmetrize;
  private final double mirrorPlanePos;
  private boolean printedWarning = false;
  private final DecimalFormat f = new DecimalFormat("#.################");

  VASPCaller(
      final GlobalConfig globconf,
      final String cellInfo,
      final String orderInfo,
      final boolean symmetrize,
      final double mirrorPlanePos,
      final int whichReference,
      final String refAtom,
      final double[] refCoords) {

    super(globconf);
    this.whichFunctional = "CUSTOM";
    this.symmetrize = symmetrize;
    this.mirrorPlanePos = mirrorPlanePos;
    this.whichReference = whichReference;
    this.refAtom = refAtom;
    this.refCoords = refCoords;

    // read in which files should be copied over
    final String auxFile = globconf.outputFolder + "-vasp.aux";
    String[] sa = null;
    try {
      sa = Input.ReadFile(auxFile);
    } catch (Exception e) {
      System.err.println(
          "ERROR: Couldn't read in VASP file list. This will fail! " + e.getMessage());
    }
    this.copyableFiles = sa;

    // read in the poscar stub
    String[] poscar;
    try {
      poscar = Input.ReadFile("POSCAR");
    } catch (Exception e) {
      System.err.println("Error in reading the custom input file!" + e.toString());
      throw new RuntimeException("Error reading in custom input file!", e);
    }

    // custom method: split the input
    String s = "";
    for (final String t : poscar) {
      s += t + System.getProperty("line.separator");
    }
    this.poscarStub = s.split("OGOCOORDS");

    // read the cell size from the input
    this.cell = new double[3][3];
    cell[0][0] = 1.0;
    cell[1][1] = 1.0;
    cell[2][2] = 1.0;
    try {
      final String[] sa2 = cellInfo.trim().split("\\,");
      int c = 0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          // parse and convert A to a.u.
          cell[i][j] = Double.parseDouble(sa2[c].trim()) * Constants.ANGTOBOHR;
          c++;
        }
      }
    } catch (Exception e) {
      System.err.println(
          "ERROR: Your cell input does not have the correct formatting. This will fail! "
              + e.getMessage());
    }

    // now the information on order and number of types
    final String[] sa3 = orderInfo.split(",");
    System.out.println(orderInfo);
    this.order = new int[sa3.length];
    this.offPerType = new int[sa3.length];
    int c = 0;
    for (int i = 0; i < sa3.length; i++) {
      final String[] tmp = sa3[i].split("\\.");
      order[i] = AtomicProperties.giveAtomicNumber(tmp[0].trim());
      offPerType[i] = c;
      c += Integer.parseInt(tmp[1].trim());
    }
  }

  private VASPCaller(final VASPCaller orig) {
    super(orig);
    this.whichFunctional = orig.whichFunctional;
    this.symmetrize = orig.symmetrize;
    this.mirrorPlanePos = orig.mirrorPlanePos;
    this.copyableFiles = orig.copyableFiles.clone();
    this.poscarStub = orig.poscarStub.clone();
    this.cell = orig.cell.clone();
    this.whichReference = orig.whichReference;
    this.refAtom = orig.refAtom;
    this.refCoords = orig.refCoords.clone();
    this.order = orig.order;
    this.offPerType = orig.offPerType;
    this.printedWarning = orig.printedWarning;
  }

  @Override
  public VASPCaller copy() {
    return new VASPCaller(this);
  }

  @Override
  public String myIDandMethod() {
    return "VASP: " + whichFunctional;
  }

  @Override
  public String myID() {
    return "VASP: " + whichFunctional;
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      final long lID,
      final CartesianCoordinates cartes,
      final boolean[][] constraints,
      final boolean isConstricted,
      final BondInfo bonds)
      throws InitIOException, IOException, ConvergenceException, CastException {

    // create our new directory
    final String folder = "vasp" + lID;
    OutputPrimitives.createAFolder(folder);

    // copy all files there that we know of
    final String sep = System.getProperty("file.separator");
    for (final String file : copyableFiles) {

      if (file.trim().equals("POSCAR")) {
        System.err.println("WARNING: Please do not specify POSCAR to be copied. Thank you!");
        continue;
      }

      final String newFile = folder + sep + file.trim();
      OutputPrimitives.createLink(file.trim(), newFile);
    }

    // setup the POSCAR file
    final short[] atomNos = cartes.getAllAtomNumbers();
    final double[][] xyz = cartes.getAllXYZCoord();
    final int noAtoms = cartes.getNoOfAtoms();
    final int noEnv = cartes.getAtomsPerMol(cartes.getNoOfMolecules() - 1);

    /*final String[] types = cartes.getAllAtomTypes();
    System.out.println("" + noAtoms);
    System.out.println("");
    for(int i = 0; i < noAtoms; i++){
        System.out.println(types[i] + "\t" + xyz[0][i]*Constants.BOHRTOANG + "\t" + xyz[1][i]*Constants.BOHRTOANG + "\t" + xyz[2][i]*Constants.BOHRTOANG);
    }*/

    final String[] poscar =
        (symmetrize) ? new String[(noAtoms - noEnv) * 2 + noEnv + 2] : new String[noAtoms + 2];
    poscar[0] = poscarStub[0];
    poscar[poscar.length - 1] = poscarStub[1];

    if (this.whichReference >= 0) {
      // check where our atom is now
      final String[] atoms = cartes.getAllAtomTypes();
      if (!atoms[whichReference].equalsIgnoreCase(this.refAtom)) {
        System.err.println(
            "ERROR: Vasp Caller has a wrong atom for move reference "
                + refAtom
                + "\t"
                + atoms[whichReference]);
      }
      final double[] move = new double[3];
      move[0] = xyz[0][whichReference] - refCoords[0];
      move[1] = xyz[1][whichReference] - refCoords[1];
      move[2] = xyz[2][whichReference] - refCoords[2];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < xyz[0].length; j++) {
          xyz[i][j] -= move[i];
        }
      }

      System.out.println("DEBUG: After operation...");
      System.out.println(
          atoms[whichReference]
              + "\t"
              + xyz[0][whichReference] * Constants.BOHRTOANG
              + "\t"
              + xyz[1][whichReference] * Constants.BOHRTOANG
              + "\t"
              + xyz[2][whichReference] * Constants.BOHRTOANG);
    }

    final int[] countTypes = new int[order.length];
    final int[] wherePoints = new int[atomNos.length];
    for (int i = 0; i < atomNos.length; i++) {
      final int where = pointMe(atomNos[i], countTypes);
      wherePoints[i] = where;
      poscar[where] =
          formatOutput(
              xyz[0][i],
              xyz[1][i],
              xyz[2][i],
              constraints[0][i],
              constraints[1][i],
              constraints[2][i]);
    }

    // add (if symmetrize is needed) another set of the same, just mirrored
    if (symmetrize) {
      for (int i = 0; i < (noAtoms - noEnv); i++) {
        final int where = pointMe(atomNos[i], countTypes);
        final double newZ = mirrorPlanePos + (mirrorPlanePos - xyz[2][i]);
        poscar[where] =
            formatOutput(
                xyz[0][i],
                xyz[1][i],
                newZ,
                constraints[0][i],
                constraints[1][i],
                constraints[2][i]);
      }
    }

    OutputPrimitives.writeOut(folder + sep + "POSCAR", poscar, false);

    // run VASP
    Process proc = null;
    try {

      String vaspCmd = System.getenv("OGO_VASPCMD");
      if (vaspCmd == null) {
        vaspCmd = "vasp";
      }

      Runtime rt = Runtime.getRuntime();
      File dir = new File(folder);
      proc = rt.exec(vaspCmd, null, dir);

      // any error message?
      StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      final int exitValue = proc.waitFor();
      if (exitValue != 0) {
        throw new ConvergenceException(
            "WARNING: VASP returns non-zero return value (local optimization): "
                + exitValue
                + ". For ID "
                + lID);
      } // else: VASP should(!) have completed normally...

    } catch (Exception e) {
      System.err.println(e.toString());
      throw new ConvergenceException("VASP has a problem (local optimization). For ID " + lID, e);
    } finally {
      if (proc != null) {
        proc.destroy();
      }
    }

    // read output in
    final String[] output = Input.ReadFile(folder + sep + "OUTCAR");
    final CartesianCoordinates optCart = new CartesianCoordinates(cartes);

    // where do the coordinates start?
    int start = 0;
    for (int i = output.length - 1; i >= 0; i--) {
      if (output[i].trim().startsWith("POSITION")) {
        start = i + 2;
        break;
      }
    }

    // nasty...
    start--; // to get rid of the plus one from pointMe()
    final double[][] xyzOpt = optCart.getAllXYZCoord();
    for (int i = 0; i < atomNos.length; i++) {
      final String[] dat = output[start + wherePoints[i]].trim().split("\\s+");
      xyzOpt[0][i] = Double.parseDouble(dat[0]) * Constants.ANGTOBOHR;
      xyzOpt[1][i] = Double.parseDouble(dat[1]) * Constants.ANGTOBOHR;
      xyzOpt[2][i] = Double.parseDouble(dat[2]) * Constants.ANGTOBOHR;
      // no explicit switch for symmetrize needed
    }

    // energy
    double e = FixedValues.NONCONVERGEDENERGY;
    for (int i = output.length - 1; i >= 0; i--) {
      if (output[i]
          .trim()
          .startsWith("FREE ENERGIE OF THE")) { // yup, there is a spelling error in the VASP output
        e = Double.parseDouble(output[i + 2].trim().split("\\s+")[4]) * Constants.EVTOHARTREE;
        break;
      }
    }
    optCart.setEnergy(e);

    final String[] types = cartes.getAllAtomTypes();
    System.out.println("" + noAtoms);
    System.out.println("" + optCart.getEnergy());
    for (int i = 0; i < noAtoms; i++) {
      System.out.println(
          types[i]
              + "\t"
              + xyzOpt[0][i] * Constants.BOHRTOANG
              + "\t"
              + xyzOpt[1][i] * Constants.BOHRTOANG
              + "\t"
              + xyzOpt[2][i] * Constants.BOHRTOANG);
    }

    // clean directory out
    ManipulationPrimitives.remove(folder);

    return optCart;
  }

  private int pointMe(final int atomNo, final int[] counter) {

    for (int i = 0; i < order.length; i++) {
      if (atomNo == order[i]) {
        counter[i]++;
        return offPerType[i]
            + counter[
                i]; // actually, the 1 too much is good as it cures poscar[0] being our first stub
        // half :-)
      }
    }

    System.err.println(
        "ERROR: The atom number you specified "
            + atomNo
            + " is not in the order you specified. This is wrong!");
    return 0;
  }

  private static double[] findOrigin(final double[][] xyz) {

    final double[] origin = new double[3];
    for (int i = 0; i < 3; i++) {
      double d = Double.MAX_VALUE;
      for (final double p : xyz[i]) {
        d = Math.min(d, p);
      }
      origin[i] = d;
    }

    return origin;
  }

  private String formatOutput(
      final double x,
      final double y,
      final double z,
      final boolean xc,
      final boolean yc,
      final boolean zc) { // , final double[] origin){

    if (!printedWarning) {
      System.err.println("WARNING: Currently, we assume a cubic cell.");
      System.err.println("WARNING: Currently, we asume a selective optimization.");
      System.err.println("WARNING: Contact author(s) if you want to complain about this.");
      printedWarning = true;
      // XXX better both once up for it
    }

    // transfer from absolute to relative coordinates
    // final double rx = (x-origin[0]) / cell[0][0];
    // final double ry = (y-origin[1]) / cell[1][1];
    // final double rz = (z-origin[2]) / cell[2][2];
    final double rx = x / cell[0][0];
    final double ry = y / cell[1][1];
    final double rz = z / cell[2][2];

    // format
    final StringBuffer buff = new StringBuffer(80);
    buff.append("  ")
        .append(f.format(rx))
        .append("  ")
        .append(f.format(ry))
        .append("  ")
        .append(f.format(rz));
    if (xc) {
      buff.append(" F");
    } else {
      buff.append(" T");
    }

    if (yc) {
      buff.append("  F");
    } else {
      buff.append("  T");
    }

    if (zc) {
      buff.append("  F");
    } else {
      buff.append("  T");
    }

    return buff.toString();
  }
}
