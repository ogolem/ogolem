/*
Copyright (c) 2012-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
              2017-2020, J. M. Dieterich and B. Hartke
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

import static org.ogolem.core.Constants.BOHRTOANG;
import static org.ogolem.core.Constants.EVTOHARTREE;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.skalevala.Parameters;
import org.ogolem.skalevala.Runot;
import org.ogolem.skalevala.Vainamoinen;

/**
 * Calls the sKalevala in an adaptive manner.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class AdaptiveSkalevalaCaller extends AbstractAdaptiveBackend {

  private static final long serialVersionUID = (long) 20201228;

  private final int runotID;
  private final Runot runot;
  private final boolean useCaching;
  private Parameters p;
  private org.ogolem.skalevala.CartesianCoordinates c;

  public AdaptiveSkalevalaCaller(
      final int runotID,
      final boolean isInAdaptive,
      final AdaptiveParameters params,
      final boolean useCaches) {
    this.runotID = runotID;
    this.runot = new Vainamoinen(runotID);
    this.useCaching = useCaches;
    // STFU
    org.ogolem.skalevala.Configuration.printTimings_$eq(false);
    if (!isInAdaptive) {
      // translate different parameter definitions
      final String[] atoms = params.getAllKeysCopy();
      int max = 0;
      for (final String atom : atoms) {
        final int no = AtomicProperties.giveAtomicNumber(atom);
        max = Math.max(no, max);
      }
      p = runot.createParameterStub(max);
      // copy data in
      copyParamData(p, params);
    }
  }

  public AdaptiveSkalevalaCaller(final AdaptiveSkalevalaCaller orig) {
    this.runotID = orig.runotID;
    this.runot = orig.runot.clone();
    this.useCaching = orig.useCaching;
    if (p != null) this.p = orig.p.copy();
    if (c != null) this.c = orig.c.clone();
  }

  @Override
  public AdaptiveSkalevalaCaller copy() {
    return new AdaptiveSkalevalaCaller(this);
  }

  @Override
  public String getMethodID() {
    return "sKalevala powered energy with " + Vainamoinen.nimi(runotID);
  }

  @Override
  public double energyCalculation(
      final long lID,
      final int iteration,
      final double[] xyz1d,
      final String[] atomTypes,
      final short[] atomNos,
      final int[] atsPerMol,
      final double[] energyparts,
      final int noOfAtoms,
      final float[] faCharges,
      final short[] iaSpins,
      final BondInfo bonds,
      final boolean hasRigidEnv) {

    float fCharge = 0;
    int spin = 0;
    for (int i = 0; i < noOfAtoms; i++) {
      fCharge += faCharges[i];
      spin += iaSpins[i];
    }
    final int charge = (int) fCharge;

    final var symBonds = bonds.getFullBondMatrix();
    final var symBuffer = symBonds.underlyingStorageBuffer();
    int symBondsIdx = 0;
    final var fullBondMatrix = new boolean[noOfAtoms][noOfAtoms];
    for (int i = 0; i < noOfAtoms; i++) {
      fullBondMatrix[i][i] = true; // mark self-bonds
      for (int j = i + i; j < noOfAtoms; j++) {
        fullBondMatrix[i][j] = symBuffer[symBondsIdx];
        fullBondMatrix[j][i] = symBuffer[symBondsIdx];
        symBondsIdx++;
      }
    }

    final org.ogolem.skalevala.CartesianCoordinates ca =
        new org.ogolem.skalevala.CartesianCoordinates(
            cartes2Cartes(xyz1d, noOfAtoms), atomTypes, atomNos, spin, charge, fullBondMatrix);

    final double e = runot.energy(ca, p);

    return e * EVTOHARTREE;
  }

  @Override
  public void gradientCalculation(
      final long lID,
      final int iteration,
      final double[] xyz1d,
      final String[] atomTypes,
      final short[] atomNos,
      final int[] atsPerMol,
      final double[] energyparts,
      final int noOfAtoms,
      final float[] faCharges,
      final short[] iaSpins,
      final BondInfo bonds,
      final Gradient gradient,
      final boolean hasRigidEnv) {

    float fCharge = 0;
    int spin = 0;
    for (int i = 0; i < noOfAtoms; i++) {
      fCharge += faCharges[i];
      spin += iaSpins[i];
    }
    final int charge = (int) fCharge;

    final var symBonds = bonds.getFullBondMatrix();
    final var symBuffer = symBonds.underlyingStorageBuffer();
    int symBondsIdx = 0;
    final var fullBondMatrix = new boolean[noOfAtoms][noOfAtoms];
    for (int i = 0; i < noOfAtoms; i++) {
      fullBondMatrix[i][i] = true; // mark self-bonds
      for (int j = i + i; j < noOfAtoms; j++) {
        fullBondMatrix[i][j] = symBuffer[symBondsIdx];
        fullBondMatrix[j][i] = symBuffer[symBondsIdx];
        symBondsIdx++;
      }
    }

    final org.ogolem.skalevala.CartesianCoordinates ca =
        new org.ogolem.skalevala.CartesianCoordinates(
            cartes2Cartes(xyz1d, noOfAtoms), atomTypes, atomNos, spin, charge, fullBondMatrix);

    final org.ogolem.skalevala.Gradient g = runot.gradient(ca, p);

    gradient2Gradient(g, noOfAtoms, gradient);
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      int geomID,
      final BondInfo bonds) {

    if (!useCaching || c == null) {

      final var noOfAtoms = cartes.getNoOfAtoms();
      final var symBonds = bonds.getFullBondMatrix();
      final var symBuffer = symBonds.underlyingStorageBuffer();
      int symBondsIdx = 0;
      final var fullBondMatrix = new boolean[noOfAtoms][noOfAtoms];
      for (int i = 0; i < noOfAtoms; i++) {
        fullBondMatrix[i][i] = true; // mark self-bonds
        for (int j = i + i; j < noOfAtoms; j++) {
          fullBondMatrix[i][j] = symBuffer[symBondsIdx];
          fullBondMatrix[j][i] = symBuffer[symBondsIdx];
          symBondsIdx++;
        }
      }

      // initialize c fresh
      c = createCartesStub(cartes, fullBondMatrix);
    }

    if (p == null) {
      // initialize the local sKalevala parameter object
      final short[] ia = cartes.getAllAtomNumbers();
      int maxNum = 0;
      for (final int i : ia) {
        maxNum = Math.max(i, maxNum);
      }
      p = runot.createParameterStub(maxNum);
    }

    copyParamData(p, params);
    cartes2Cartes(c, cartes);

    final double e = runot.energy(c, p);

    return e * EVTOHARTREE;
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      int geomID,
      final BondInfo bonds,
      final double[] grad) {

    if (!useCaching || c == null) {

      final var noOfAtoms = cartes.getNoOfAtoms();
      final var symBonds = bonds.getFullBondMatrix();
      final var symBuffer = symBonds.underlyingStorageBuffer();
      int symBondsIdx = 0;
      final var fullBondMatrix = new boolean[noOfAtoms][noOfAtoms];
      for (int i = 0; i < noOfAtoms; i++) {
        fullBondMatrix[i][i] = true; // mark self-bonds
        for (int j = i + i; j < noOfAtoms; j++) {
          fullBondMatrix[i][j] = symBuffer[symBondsIdx];
          fullBondMatrix[j][i] = symBuffer[symBondsIdx];
          symBondsIdx++;
        }
      }

      // initialize c fresh
      c = createCartesStub(cartes, fullBondMatrix);
    }

    if (p == null) {
      // initialize the local sKalevala parameter object
      final short[] ia = cartes.getAllAtomNumbers();
      int maxNum = 0;
      for (final int i : ia) {
        maxNum = Math.max(i, maxNum);
      }
      p = runot.createParameterStub(maxNum);
    }

    copyParamData(p, params);
    cartes2Cartes(c, cartes);

    System.err.println("ERROR: Routine for parameter gradient with sKalevala not implemented yet.");
    return 1000.0;
    // TODO
  }

  @Override
  public double[][] minMaxBordersForParams(final AdaptiveParameters params) {

    final String[] keys = params.getAllKeysCopy();
    final int noOfParams = params.getNumberOfParamters();
    final double[][] borders = new double[2][noOfParams];

    int pc = 0;
    for (final String key : keys) {
      final int noForKey = params.getAmountOfParametersForKey(key);

      // first are the s/p/d ionic potentials
      borders[0][pc] = 0.0;
      borders[1][pc] = 100.0;
      pc++;
      borders[0][pc] = 0.0;
      borders[1][pc] = 100.0;
      pc++;
      borders[0][pc] = 0.0;
      borders[1][pc] = 100.0;
      pc++;

      // then the exponents for s/p/d
      borders[0][pc] = 1.0;
      borders[1][pc] = 5.0;
      pc++;
      borders[0][pc] = 1.0;
      borders[1][pc] = 5.0;
      pc++;
      borders[0][pc] = 1.0;
      borders[1][pc] = 5.0;
      pc++;

      // then k_s k_p and k_d
      borders[0][pc] = 0.3;
      borders[1][pc] = 2.0;
      pc++;
      borders[0][pc] = 0.3;
      borders[1][pc] = 2.0;
      pc++;
      borders[0][pc] = 0.3;
      borders[1][pc] = 2.0;
      pc++;

      // a/b/c for electronic interaction
      borders[0][pc] = 0.1;
      borders[1][pc] = 4.0;
      pc++;
      borders[0][pc] = 0.1;
      borders[1][pc] = 4.0;
      pc++;
      borders[0][pc] = 0.1;
      borders[1][pc] = 4.0;
      pc++;

      // delta and epsilon for nuclear interaction
      borders[0][pc] = 0.1;
      borders[1][pc] = 5.0;
      pc++;
      borders[0][pc] = 0.1;
      borders[1][pc] = 5.0;
      pc++;

      // the additional rest (e.g. LJ parameters)
      for (int i = pc; i < noForKey; i++) {
        // XXX this is very generic and therefore inefficient!
        borders[0][pc] = -100.0;
        borders[1][pc] = 100.0;
        pc++;
      }
    }

    return borders;
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      final ArrayList<CartesianCoordinates> refCartes, final String sMethod) {

    // figure out which atoms we are dealing with
    final List<String> allAtoms = new ArrayList<>();
    for (final CartesianCoordinates cartes : refCartes) {
      final String[] atoms = cartes.getAllAtomTypes();
      for (final String atom : atoms) {
        if (!allAtoms.contains(atom)) {
          allAtoms.add(atom);
        }
      }
    }

    // check maximum number
    int maxNum = 0;
    for (final String atom : allAtoms) {
      maxNum = Math.max(maxNum, AtomicProperties.giveAtomicNumber(atom));
    }

    // initialize the local sKalevala parameter object
    p = runot.createParameterStub(maxNum);

    // now we need to know how many parameters we need per atom
    final int[] noParams = p.getParameterDimsPerAtom();
    int tot = 0;
    for (final int i : noParams) {
      tot += i;
    }

    final int noAtoms = allAtoms.size();
    final int[] paramsPerAt = new int[noAtoms];
    final String[] atoms = new String[noAtoms];
    for (int i = 0; i < noAtoms; i++) {
      atoms[i] = allAtoms.get(i);
      paramsPerAt[i] = tot;
    }

    final AdaptiveParameters paramStub =
        new AdaptiveParameters(noAtoms * tot, -1, atoms, paramsPerAt, sMethod);

    return paramStub;
  }

  private static void copyParamData(final Parameters p, final AdaptiveParameters adapt) {

    final String[] keys = adapt.getAllKeysCopy();

    for (final String key : keys) {
      final int no = AtomicProperties.giveAtomicNumber(key);
      final double[] vals = adapt.getParametersForKey(key);
      p.putValuesIn(vals, no);
    }
  }

  private static org.ogolem.skalevala.CartesianCoordinates createCartesStub(
      final CartesianCoordinates cartes, final boolean[][] bonding) {

    final double[][] coords = new double[cartes.getNoOfAtoms()][3];
    final org.ogolem.skalevala.CartesianCoordinates c =
        new org.ogolem.skalevala.CartesianCoordinates(
            coords,
            cartes.getAllAtomTypes(),
            cartes.getAllAtomNumbers(),
            cartes.getTotalSpin(),
            cartes.getTotalCharge(),
            bonding);

    return c;
  }

  private static void cartes2Cartes(
      final org.ogolem.skalevala.CartesianCoordinates c, final CartesianCoordinates cartes) {

    final int noAtoms = cartes.getNoOfAtoms();
    final double[][] xyz = cartes.getAllXYZCoord();
    final double[][] cxyz = c.coordinates();
    for (int i = 0; i < noAtoms; i++) {
      cxyz[i][0] = xyz[0][i];
      cxyz[i][1] = xyz[1][i];
      cxyz[i][2] = xyz[2][i];
    }
  }

  private static double[][] cartes2Cartes(final double[] xyz1d, final int noAtoms) {

    final double[][] cxyz = new double[noAtoms][3];

    for (int i = 0; i < noAtoms; i++) {
      cxyz[i][0] = xyz1d[i];
      cxyz[i][1] = xyz1d[noAtoms + i];
      cxyz[i][2] = xyz1d[2 * noAtoms + i];
    }

    return cxyz;
  }

  private static void gradient2Gradient(
      final org.ogolem.skalevala.Gradient grad, final int noOfAtoms, final Gradient gradient) {

    final double[][] gradMat = gradient.getTotalGradient();

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < noOfAtoms; j++) {
        gradMat[i][j] = grad.gradient().r(i, j) * EVTOHARTREE * BOHRTOANG;
      }
    }

    gradient.setTotalEnergy(grad.energy() * EVTOHARTREE);
  }
}
