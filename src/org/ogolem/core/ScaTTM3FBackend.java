/*
Copyright (c) 2012-2014, J. M. Dieterich
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

import java.io.File;
import java.util.Locale;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.properties.DipoleMoment;
import org.ogolem.scaTTM3F.TTM3F;

/**
 * Calls the Scala-TTM3F force field for highly exact water clusters.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class ScaTTM3FBackend implements CartesianFullBackend {

  private static final long serialVersionUID = (long) 20150215;
  private static final boolean ANALYSISOUT = false;
  private static final boolean TESTGRADIENT = false;
  private static final boolean TRAPNANS = true;
  // typical length of a hydrogen bond in water is 197 pm (according to wiki), ergo this is would
  // really be way too short already
  private static final double CUTDIST = 0.9;
  private static final double CUTDISTSQ = CUTDIST * CUTDIST; // this is in angstrom (sq)!
  private static final double CUTENERGY = 0.1; // add some hartree
  private static final double CUTGRADIENT = 0.1; // add some hartree/bohr

  private final boolean simpleExtraTerm;
  private final boolean DEBUG;
  private final double CUTOFF;
  private final int noWat;
  private final int noAts;
  private final double[][] coordCache;
  private final double[][] gradCache;
  private final double[][] gradCache2;
  private final TTM3F ttm3f;
  private final boolean printDipols;
  private final double[] dipolsCache;
  private final String outFolder;
  private final String dipolLog;
  private final String totDipLog;
  private final String[] dipolsOut;
  private final String[] totDipolOut;
  private long count;

  public ScaTTM3FBackend(
      final int noWaters,
      final double cut,
      final boolean printDipols,
      final String outFolder,
      final boolean simpleCutoffTerm) {
    this.DEBUG = (GlobalConfig.DEBUGLEVEL > 2);
    this.CUTOFF = cut;
    this.noWat = noWaters;
    this.noAts = 3 * noWat;
    this.coordCache = new double[3 * noWat][3];
    this.gradCache = new double[3 * noWat][3];
    this.gradCache2 = new double[3 * noWat][3];
    this.ttm3f = new TTM3F(noWat);
    this.printDipols = printDipols;
    this.dipolsCache = new double[3 * noWat];
    this.dipolsOut = (printDipols) ? new String[noWat + 2] : null;
    this.outFolder = outFolder;
    this.dipolLog = outFolder + File.separator + "dipols-scaTTM3F.log";
    this.totDipLog = outFolder + File.separator + "totdipol-scaTTM3F.trj";
    this.totDipolOut = (printDipols) ? new String[3] : null;
    this.count = 0l;
    this.simpleExtraTerm = simpleCutoffTerm;
  }

  @Override
  public ScaTTM3FBackend copy() {
    return new ScaTTM3FBackend(
        this.noWat, this.CUTOFF, this.printDipols, this.outFolder, this.simpleExtraTerm);
  }

  @Override
  public String getMethodID() {
    String s = "TTM3F (Scala version):\n";
    s += "\t cutoff energy:  " + this.CUTOFF + "\n";
    s += "\t getting dipols: " + this.printDipols + "\n";
    if (this.printDipols) {
      s += "\t dipol output:   " + this.outFolder + "\n";
    }
    s += "\t simple anti-cold-fusion term: " + this.simpleExtraTerm;

    return s;
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

    assert (noAts * 3 == xyz1D.length);
    count++;

    // zero the partial contributions
    for (int i = 0; i < energyparts.length; i++) {
      energyparts[i] = 0.0;
    }

    moveCoords(xyz1D);

    double e;
    try {
      e = ttm3f.ttm3f(coordCache, gradCache, energyparts, false, true, dipolsCache) * KCALTOHARTREE;
    } catch (Exception ex) {
      System.out.println("WARNING: Exception in TTM3F FF: " + ex.toString());
      if (DEBUG) {
        ex.printStackTrace(System.err);
        System.err.println("The geometry causing this was: ");
        printCoords();
      }
      e = FixedValues.NONCONVERGEDENERGY;
    }

    final double[][] gradMat = gradient.getTotalGradient();

    if (Double.isInfinite(e)
        || Double.isNaN(e)
        || e < -CUTOFF
        || Math.abs(e) > FixedValues.NONCONVERGEDENERGY) {
      if (ANALYSISOUT && Double.isNaN(e)) {
        System.out.println("Found energy to be " + e + " for " + lID);
        System.out.println("Coord cache content: ");
        printCoords();
        System.out.println("Actual coord content:");
        for (final double d : xyz1D) {
          System.out.println(" " + d);
        }
      }
      e = FixedValues.NONCONVERGEDENERGY;
    }

    if (printDipols) {
      dipolLogging();
    }

    /*for(int i = 0; i < noAts; i++){
              for(int j = 0; j < 3; j++){
                  final double d = gradCache[i][j];
                  //if(Double.isInfinite(d) || Double.isNaN(d) || d < -CUTOFF) gradCache[i][j] = FixedValues.NONCONVERGEDGRADIENT;
    if(Double.isInfinite(d) || Double.isNaN(d) || Math.abs(d) > FixedValues.NONCONVERGEDGRADIENT){gradCache[i][j] = FixedValues.NONCONVERGEDGRADIENT;}
              }
          }*/

    moveGrad(gradMat);

    if (simpleExtraTerm) {

      final boolean[][] bondMat = bonds.getFullBondMatrix();
      boolean anythingFound = false;
      for (int i = 0; i < noAts - 1; i++) {
        for (int j = i + 1; j < noAts; j++) {

          if (bondMat[i][j]) {
            continue;
          } // this would not help if there would be a water collapsing in itself w/ an aritifically
          // low energy: never seen that though.

          final double dX = coordCache[i][0] - coordCache[j][0];
          final double dY = coordCache[i][1] - coordCache[j][1];
          final double dZ = coordCache[i][2] - coordCache[j][2];
          final double distSq = dX * dX + dY * dY + dZ * dZ;

          if (distSq <= CUTDISTSQ) {

            anythingFound = true;

            // increase energy
            final double distInv = Math.sqrt(distSq) - CUTDIST;
            final double divProdX = dX * distInv;
            final double divProdY = dY * distInv;
            final double divProdZ = dZ * distInv;
            e += CUTENERGY * distInv * distInv;

            if (DEBUG) {
              System.out.println(
                  "DEBUG: distance (sq) too small. Distance: "
                      + distSq
                      + " vs "
                      + CUTDISTSQ
                      + " for "
                      + i
                      + " "
                      + j
                      + " for id "
                      + lID);
            }

            // set gradient
            gradMat[0][i] = CUTGRADIENT * divProdX;
            gradMat[0][j] = -(CUTGRADIENT * divProdX);
            gradMat[1][i] = CUTGRADIENT * divProdY;
            gradMat[1][j] = -(CUTGRADIENT * divProdY);
            gradMat[2][i] = CUTGRADIENT * divProdZ;
            gradMat[2][j] = -(CUTGRADIENT * divProdZ);

            // increase gradient
            /*gradMat[0][i] += CUTGRADIENT * divProdX;
            gradMat[0][j] -= CUTGRADIENT * divProdX;
            gradMat[1][i] += CUTGRADIENT * divProdY;
            gradMat[1][j] -= CUTGRADIENT * divProdY;
            gradMat[2][i] += CUTGRADIENT * divProdZ;
            gradMat[2][j] -= CUTGRADIENT * divProdZ*/
          }
        }
      }

      if (DEBUG && anythingFound) {
        System.out.println("DEBUG: Energy is " + e + " for a problematic case...");
      }
    }

    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < noAts; i++) {
        final double d = gradMat[j][i];
        if (Double.isInfinite(d)
            || Double.isNaN(d)
            || Math.abs(d) > FixedValues.NONCONVERGEDGRADIENT) {
          if (DEBUG) {
            System.out.println(
                "DEBUG: Trying to sanitize "
                    + j
                    + " "
                    + i
                    + " is "
                    + d
                    + " not "
                    + FixedValues.NONCONVERGEDGRADIENT);
          }
          gradMat[j][i] = FixedValues.NONCONVERGEDGRADIENT;
        }
      }
    }

    gradient.setTotalEnergy(e);

    if (TESTGRADIENT) {
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
      final double[][] num = numGrad.getTotalGradient();
      for (int i = 0; i < noAts; i++) {
        for (int j = 0; j < 3; j++) {
          if (Math.abs(num[j][i] - gradMat[j][i]) > 1E-7)
            System.out.println(
                "DEBUG GRADIENTDIFF Diff "
                    + num[j][i]
                    + "\t"
                    + gradMat[j][i]
                    + "\t"
                    + (gradMat[j][i] / num[j][i])
                    + "\tat\t"
                    + i
                    + "\t"
                    + j);
          else {
            System.out.println("GRADIENTDIFF FINE");
          }
        }
      }
      System.out.println("GEOMETRY");
      for (int i = 0; i < noAts; i++) {
        System.out.println(
            " " + coordCache[i][0] + "\t" + coordCache[i][1] + "\t" + coordCache[i][2]);
      }

      System.out.println("NUMERICAL GRADIENT");
      for (int i = 0; i < noAts; i++) {
        System.out.println(" " + num[0][i] + "\t" + num[1][i] + "\t" + num[2][i]);
      }
    }
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
      final BondInfo bonds,
      final boolean hasRigidEnv) {

    assert (xyz1D.length == 3 * noAts);
    count++;

    // zero the partial contributions
    for (int i = 0; i < energyparts.length; i++) {
      energyparts[i] = 0.0;
    }

    moveCoords(xyz1D);

    double e;
    try {
      e = ttm3f.ttm3f(coordCache, gradCache2, energyparts, true, true, dipolsCache) * KCALTOHARTREE;
    } catch (Exception ex) {
      System.out.println("WARNING: Exception in TTM3F FF: " + ex.toString());
      if (DEBUG) {
        ex.printStackTrace(System.err);
        System.err.println("The geometry causing this was: ");
        printCoords();
      }
      return FixedValues.NONCONVERGEDENERGY;
    }

    if (printDipols) {
      dipolLogging();
    }

    if (simpleExtraTerm) {

      final boolean[][] bondMat = bonds.getFullBondMatrix();
      for (int i = 0; i < noAts - 1; i++) {
        for (int j = i + 1; j < noAts; j++) {

          if (bondMat[i][j]) {
            continue;
          } // this would not help if there would be a water collapsing in itself w/ an aritifically
          // low energy: never seen that though.

          final double dX = coordCache[i][0] - coordCache[j][0];
          final double dY = coordCache[i][1] - coordCache[j][1];
          final double dZ = coordCache[i][2] - coordCache[j][2];
          final double distSq = dX * dX + dY * dY + dZ * dZ;

          if (distSq <= CUTDISTSQ) {
            // increase energy
            final double distInv = Math.sqrt(distSq) - CUTDIST;
            e += CUTENERGY * distInv * distInv;
          }
        }
      }
    }

    if (Double.isInfinite(e)
        || Double.isNaN(e)
        || e < -CUTOFF
        || Math.abs(e) > FixedValues.NONCONVERGEDENERGY) {
      if (ANALYSISOUT && Double.isNaN(e)) {
        System.out.println("Found energy to be " + e + " for " + lID);
        System.out.println("Coord cache content: ");
        printCoords();
        System.out.println("Actual coord content:");
        for (final double d : xyz1D) {
          System.out.println(" " + d);
        }
      }
      return FixedValues.NONCONVERGEDENERGY;
    }

    return e;
  }

  private void moveCoords(final double[] coords) {

    if (TRAPNANS) {
      for (int i = 0; i < coords.length; i++) {
        if (Double.isNaN(coords[i])) {
          System.out.println("ERROR: NaN in coords[] " + i + ". Obviously fatal.");
          throw new RuntimeException("NaN in coords[] " + i + ". Obviously fatal.");
        }
      }
    }

    for (int i = 0; i < noWat; i++) {
      coordCache[i][0] = coords[3 * i] * BOHRTOANG; // oxygen
      coordCache[i][1] = coords[noAts + 3 * i] * BOHRTOANG; // oxygen
      coordCache[i][2] = coords[2 * noAts + 3 * i] * BOHRTOANG; // oxygen
      coordCache[noWat + 2 * i][0] = coords[3 * i + 1] * BOHRTOANG; // hydrogen 1
      coordCache[noWat + 2 * i][1] = coords[noAts + 3 * i + 1] * BOHRTOANG; // hydrogen 1
      coordCache[noWat + 2 * i][2] = coords[2 * noAts + 3 * i + 1] * BOHRTOANG; // hydrogen 1
      coordCache[noWat + 2 * i + 1][0] = coords[3 * i + 2] * BOHRTOANG; // hydrogen 2
      coordCache[noWat + 2 * i + 1][1] = coords[noAts + 3 * i + 2] * BOHRTOANG; // hydrogen 2
      coordCache[noWat + 2 * i + 1][2] = coords[2 * noAts + 3 * i + 2] * BOHRTOANG; // hydrogen 2
    }

    if (DEBUG) {
      System.out.println(noAts + "\n");
      int c = 0;
      for (int i = 0; i < noAts; i++) {
        c++;
        String s = (c <= noWat) ? "O" : "H";
        for (int j = 0; j < 3; j++) {
          s += "\t" + coordCache[i][j];
        }
        System.out.println(s);
      }
    }
  }

  private void printCoords() {

    System.err.println(noAts);
    System.err.println();
    for (int i = 0; i < noWat; i++) {
      final double oX = coordCache[i][0]; // oxygen
      final double oY = coordCache[i][1]; // oxygen
      final double oZ = coordCache[i][2]; // oxygen
      final double h1X = coordCache[noWat + 2 * i][0]; // hydrogen 1
      final double h1Y = coordCache[noWat + 2 * i][1]; // hydrogen 1
      final double h1Z = coordCache[noWat + 2 * i][2]; // hydrogen 1
      final double h2X = coordCache[noWat + 2 * i + 1][0]; // hydrogen 2
      final double h2Y = coordCache[noWat + 2 * i + 1][1]; // hydrogen 2
      final double h2Z = coordCache[noWat + 2 * i + 1][2]; // hydrogen 2
      System.err.println("O\t" + oX + "\t" + oY + "\t" + oZ);
      System.err.println("H\t" + h1X + "\t" + h1Y + "\t" + h1Z);
      System.err.println("H\t" + h2X + "\t" + h2Y + "\t" + h2Z);
    }
  }

  private void moveGrad(final double[][] grad) {

    for (int i = 0; i < noWat; i++) {
      for (int j = 0; j < 3; j++) {
        grad[j][3 * i] = gradCache[i][j] * KCALTOHARTREE * BOHRTOANG;
        grad[j][3 * i + 1] = gradCache[noWat + 2 * i][j] * KCALTOHARTREE * BOHRTOANG;
        grad[j][3 * i + 2] = gradCache[noWat + 2 * i + 1][j] * KCALTOHARTREE * BOHRTOANG;
      }
    }

    if (TRAPNANS) {
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < noWat; i++) {
          if (Double.isNaN(grad[j][i])) {
            System.out.println("ERROR: NaN in grad[][] " + j + "  " + i + ". Obviously fatal.");
            throw new RuntimeException("NaN in grad[][] " + j + "  " + i + ". Obviously fatal.");
          }
        }
      }
    }
  }

  public DipoleMoment getLastDipol() {

    // compute total dipole moment
    final double[] totDipol = new double[3];
    for (int i = 0; i < noWat; i++) {
      totDipol[0] += dipolsCache[3 * i];
      totDipol[1] += dipolsCache[3 * i + 1];
      totDipol[2] += dipolsCache[3 * i + 2];
    }

    return new DipoleMoment(totDipol);
  }

  private void dipolLogging() {

    if (DEBUG) {
      // get individual dipols only if debug level is high
      dipolsOut[0] = "" + count;
      dipolsOut[1] = "";
      for (int i = 0; i < noWat; i++) {
        dipolsOut[i + 2] =
            "\t"
                + (i + 1)
                + "\t"
                + dipolsCache[3 * i]
                + "\t"
                + dipolsCache[3 * i + 1]
                + "\t"
                + dipolsCache[3 * i + 2];
      }
      try {
        OutputPrimitives.writeOut(dipolLog, dipolsOut, true);
      } catch (Exception ex) {
        System.err.println("WARNING: Failure to write log entry " + count + " to file " + dipolLog);
        ex.printStackTrace(System.err);
      }
    }

    // always compute and log the total dipole moment for IR prediction purposes
    final double[] totDipol = new double[3];
    for (int i = 0; i < noWat; i++) {
      totDipol[0] += dipolsCache[3 * i];
      totDipol[1] += dipolsCache[3 * i + 1];
      totDipol[2] += dipolsCache[3 * i + 2];
    }

    // bring the output into the same format as the trj output in MDSystem in order to use the
    // analyzation tools for that.
    totDipolOut[0] = "" + 1; // one "atom" (whatever that is...)
    totDipolOut[1] = "" + count;
    // first fields zero, total dipole moment in place of the velocity
    totDipolOut[2] =
        "XX" // dummy :-)
            + " "
            + String.format(Locale.US, "%14.7f", 0.0)
            + " "
            + String.format(Locale.US, "%14.7f", 0.0)
            + " "
            + String.format(Locale.US, "%14.7f", 0.0)
            + " "
            + String.format(Locale.US, "%14.7f", totDipol[0])
            + " "
            + String.format(Locale.US, "%14.7f", totDipol[1])
            + " "
            + String.format(Locale.US, "%14.7f", totDipol[2]);

    // attach it to our logfile
    try {
      OutputPrimitives.writeOut(totDipLog, totDipolOut, true);
    } catch (Exception ex) {
      System.err.println("WARNING: Failure to write log entry " + count + " to file " + totDipLog);
      ex.printStackTrace(System.err);
    }
  }
}
