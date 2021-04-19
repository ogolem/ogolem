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
package org.ogolem.adaptive;

import static org.ogolem.math.FastFunctions.pow;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.helpers.Machine;

/**
 * Adaptive Lennard Jones force field. Will not only take into account classical 6/12 Terms but also
 * more, depending on the users wish. Initial code to use a Stillinger-Weber-Gong potential for
 * three body terms (untested and buggy).
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public final class AdaptiveLJFF extends AbstractAdaptiveBackend {

  // XXX implement Horner scheme for the pow() calls

  // the ID
  private static final long serialVersionUID = (long) 20200321;
  private static final boolean DEBUG = false;

  private final boolean easyMix;
  private final AdaptiveParameters params;

  private final double machPrec;

  private final int startPot;
  private final int endPot;
  private final int potIncr;

  private final double closeCutBlow;
  private final double farCutBlow;

  private final boolean discard2Body;
  private final boolean use3Body;

  private final boolean useCaching;
  private int[] paramOffsetCache;
  private int[] paramOffsetCache2b1;
  private int[][] paramOffsetCache3b;
  // TODO debug 3body
  // TODO allow coord gradient to be w/o 2body contrib

  /**
   * Constructor
   *
   * @param isInAdaptive
   * @param easyLBMix
   * @param startPot
   * @param endPot
   * @param potIncr
   * @param parameters
   * @param use3Body
   * @param closeCutBlow
   * @param distCutBlow
   * @param useCache
   * @param discard2Body
   */
  public AdaptiveLJFF(
      final boolean isInAdaptive,
      final boolean easyLBMix,
      final int startPot,
      final int endPot,
      final int potIncr,
      final AdaptiveParameters parameters,
      final boolean use3Body,
      final double closeCutBlow,
      final double distCutBlow,
      final boolean useCache,
      final boolean discard2Body) {

    assert (startPot > 0);
    assert (endPot >= startPot);
    assert (potIncr >= 0);

    this.machPrec = Machine.calcMachinePrecision();
    this.easyMix = easyLBMix;
    this.startPot = startPot;
    this.endPot = endPot;
    this.potIncr = potIncr;
    this.discard2Body = discard2Body;
    this.use3Body = use3Body;
    this.closeCutBlow = closeCutBlow;
    this.farCutBlow = distCutBlow;
    this.useCaching = useCache;

    if (!isInAdaptive) {
      params = parameters;
    } else {
      // we simply do not need them
      params = null;
    }
  }

  private AdaptiveLJFF(final AdaptiveLJFF orig) {
    this.easyMix = orig.easyMix;
    this.discard2Body = orig.discard2Body;
    this.use3Body = orig.use3Body;
    this.closeCutBlow = orig.closeCutBlow;
    this.farCutBlow = orig.farCutBlow;
    this.machPrec = orig.machPrec;
    this.endPot = orig.endPot;
    this.potIncr = orig.potIncr;
    this.startPot = orig.startPot;
    this.useCaching = orig.useCaching;

    if (orig.params != null) {
      this.params = new AdaptiveParameters(orig.params);
    } else {
      this.params = null;
    }

    if (orig.paramOffsetCache != null) {
      this.paramOffsetCache = orig.paramOffsetCache.clone();
    } else {
      this.paramOffsetCache = null;
    }

    if (orig.paramOffsetCache2b1 != null) {
      this.paramOffsetCache2b1 = orig.paramOffsetCache2b1.clone();
    } else {
      this.paramOffsetCache2b1 = null;
    }

    // shallow copy should be sufficient, I hope.
    if (orig.paramOffsetCache3b != null) {
      this.paramOffsetCache3b = orig.paramOffsetCache3b.clone();
    } else {
      this.paramOffsetCache3b = null;
    }
  }

  @Override
  public AdaptiveLJFF copy() {
    return new AdaptiveLJFF(this);
  }

  @Override
  public String getMethodID() {
    return "Adaptive LJFF";
  }

  @Override
  public void gradientCalculation(
      final long lID,
      final int iIteration,
      final double[] daXYZ1D,
      final String[] saAtomTypes,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      final int iNoOfAtoms,
      final float[] faCharges,
      final short[] iaSpins,
      final BondInfo bonds,
      final Gradient gradient,
      final boolean hasRigidEnv) {

    if (useCaching && paramOffsetCache == null) {
      // initialize caches
      initializeCaches(params, saAtomTypes, use3Body);
    }

    // the matrix for the gradient
    gradient.zeroGradient();
    final double[][] gradientMat = gradient.getTotalGradient();

    // get in O(N)
    final double[] radii = new double[iNoOfAtoms];
    for (int i = 0; i < iNoOfAtoms; i++) {
      radii[i] = AtomicProperties.giveRadius(atomNos[i]);
    }

    double[] daFirstParams = null;
    double[] daSecParams;
    double[] daParams = params.getAllParamters();
    int counter = -1;

    double twoBodyEnergy = 0.0;

    int lastOffset = 0;
    for (int i = 0; i < atsPerMol.length - 1; i++) {
      lastOffset += atsPerMol[i];
    }

    final int firstLoopAtomNo = (hasRigidEnv) ? lastOffset - 1 : iNoOfAtoms - 1;
    for (int i = 0; i < firstLoopAtomNo; i++) {
      if (easyMix) {
        // first part of parameters
        daFirstParams = params.getParametersForKey(saAtomTypes[i]);
        if (daFirstParams == null) {
          continue;
        }
      }

      final double rad1 = radii[i];

      final double xi = daXYZ1D[i];
      final double yi = daXYZ1D[iNoOfAtoms + i];
      final double zi = daXYZ1D[2 * iNoOfAtoms + i];

      double gradXI = gradientMat[0][i];
      double gradYI = gradientMat[1][i];
      double gradZI = gradientMat[2][i];

      for (int j = i + 1; j < iNoOfAtoms; j++) {

        counter++;

        final double rad2 = radii[j];

        // check whether one should use a cutoff
        final double dDistX = xi - daXYZ1D[j];
        final double dDistY = yi - daXYZ1D[iNoOfAtoms + j];
        final double dDistZ = zi - daXYZ1D[2 * iNoOfAtoms + j];
        final double dDist = Math.sqrt(dDistX * dDistX + dDistY * dDistY + dDistZ * dDistZ);
        final double dAddedRadii = rad1 + rad2;
        if (dDist < (closeCutBlow * dAddedRadii)) {
          // use cutoff. speculative... ;-)
          gradientMat[0][i] += FixedValues.NONCONVERGEDGRADIENT;
          gradientMat[0][j] -= FixedValues.NONCONVERGEDGRADIENT;
          gradientMat[1][i] += FixedValues.NONCONVERGEDGRADIENT;
          gradientMat[1][j] -= FixedValues.NONCONVERGEDGRADIENT;
          gradientMat[2][i] += FixedValues.NONCONVERGEDGRADIENT;
          gradientMat[2][j] -= FixedValues.NONCONVERGEDGRADIENT;
          continue;
        } else if (dDist > (farCutBlow * dAddedRadii)) {
          continue;
        }

        int offset = 0;
        if (easyMix) {
          // second part of parameters
          daSecParams = params.getParametersForKey(saAtomTypes[j]);
          if (daSecParams == null) {
            continue;
          }

          // mix them
          daParams = mixParams(daFirstParams, daSecParams);
        } else {
          if (!useCaching) {
            String sPair = saAtomTypes[i] + saAtomTypes[j];
            daParams = params.getParametersForKey(sPair);
            if (daParams == null) {
              sPair = saAtomTypes[j] + saAtomTypes[i];
              daParams = params.getParametersForKey(sPair);
              if (daParams == null) {
                // ok this IS a problem
                System.err.println(
                    "WARNING: We can't find parameters for the pair "
                        + saAtomTypes[i]
                        + " and "
                        + saAtomTypes[j]
                        + ".");
                continue;
              }
            }
          } else {
            offset = paramOffsetCache[counter];
          }
        }

        final double dDistInv = 1.0 / dDist;
        // divide the distances in all three dimensions by the total distance to cover for the
        // coordinate system
        final double dDivProdX = dDistX * dDistInv;
        final double dDivProdY = dDistY * dDistInv;
        final double dDivProdZ = dDistZ * dDistInv;

        final double epsFac = 4 * daParams[offset];

        // calculate the contributions
        int iCounter = offset + 1;
        double dTempDeriv = 0.0;
        for (int iPot = startPot; iPot <= endPot; iPot = iPot + potIncr) {
          final double dSign = Math.signum(daParams[iCounter]);
          final double sigTerm = pow((daParams[iCounter] * dDistInv), iPot);
          final double twoBodyContr = epsFac * dSign * sigTerm;
          twoBodyEnergy += twoBodyContr;
          dTempDeriv += twoBodyContr * dDistInv * (-iPot);
          iCounter++;
        }

        gradXI += dTempDeriv * dDivProdX;
        gradientMat[0][j] -= dTempDeriv * dDivProdX;
        gradYI += dTempDeriv * dDivProdY;
        gradientMat[1][j] -= dTempDeriv * dDivProdY;
        gradZI += dTempDeriv * dDivProdZ;
        gradientMat[2][j] -= dTempDeriv * dDivProdZ;
      }

      gradientMat[0][i] = gradXI;
      gradientMat[1][i] = gradYI;
      gradientMat[2][i] = gradZI;
    }

    /*
     * Now the Stillinger-Weber-Gong term
     */

    if (use3Body && saAtomTypes.length > 2) {

      // XXXX we do not use hasRigidEnv for 3body!

      // put the 1D coordinates into 3D ones - XXX remove
      final double[][] daXYZ = new double[3][iNoOfAtoms];
      System.arraycopy(daXYZ1D, 0, daXYZ[0], 0, iNoOfAtoms);
      System.arraycopy(daXYZ1D, iNoOfAtoms, daXYZ[1], 0, iNoOfAtoms);
      System.arraycopy(daXYZ1D, iNoOfAtoms * 2, daXYZ[2], 0, iNoOfAtoms);

      int off2b1;
      int off2b2;
      int off2b3;
      int off3b;
      counter = -1;
      int c = -1;
      daParams = params.getAllParamters();

      for (int k = 0; k < iNoOfAtoms; k++) {
        for (int j = 0; j < k; j++) {

          c++;
          if (!useCaching) {
            off2b1 = posOfKey2sin3(saAtomTypes, k, j, params);
          } else {
            off2b1 = paramOffsetCache2b1[c];
          }

          final double dDistJKX = daXYZ[0][k] - daXYZ[0][j];
          final double dDistJKY = daXYZ[1][k] - daXYZ[1][j];
          final double dDistJKZ = daXYZ[2][k] - daXYZ[2][j];
          final double dDistJKSq = dDistJKX * dDistJKX + dDistJKY * dDistJKY + dDistJKZ * dDistJKZ;

          final double dDistJK = Math.sqrt(dDistJKSq);

          // ok, check whether the distance is above the equilibrium distance for this pair
          if (dDistJK > daParams[off2b1]) {
            continue;
          }

          final double dRJKRCC = 1.0 / (dDistJK - daParams[off2b1]);
          final double dGRJKR = daParams[off2b1 + 1] * dRJKRCC;
          final double dFJK = Math.exp(dGRJKR);

          for (int i = 0; i < j; i++) {

            counter++;

            if (!useCaching) {
              off3b = posOfKey3in3(saAtomTypes, i, j, k, params);
              off2b2 = posOfKey2sin3(saAtomTypes, i, j, params);
              off2b3 = posOfKey2sin3(saAtomTypes, i, k, params);
            } else {
              off3b = paramOffsetCache3b[2][counter];
              off2b2 = paramOffsetCache3b[0][counter];
              off2b3 = paramOffsetCache3b[1][counter];
            }

            final double dIJX = daXYZ[0][i] - daXYZ[0][j];
            final double dIJY = daXYZ[1][i] - daXYZ[1][j];
            final double dIJZ = daXYZ[2][i] - daXYZ[2][j];

            final double dDistIJSq = dIJX * dIJX + dIJY * dIJY + dIJZ * dIJZ;

            final double dIKX = daXYZ[0][i] - daXYZ[0][k];
            final double dIKY = daXYZ[1][i] - daXYZ[1][k];
            final double dIKZ = daXYZ[2][i] - daXYZ[2][k];

            final double dDistIKSq = dIKX * dIKX + dIKY * dIKY + dIKZ * dIKZ;

            // ok, check whether the distance is above the equilibrium distance for this pair
            if (dDistIJSq > daParams[off2b2] * daParams[off2b2]) {
              continue;
            }
            if (dDistIKSq > daParams[off2b3] * daParams[off2b3]) {
              continue;
            }

            final double dDistIJ = Math.sqrt(dDistIJSq);
            final double dDistIK = Math.sqrt(dDistIKSq);

            // finally the real derivative
            final double dRIJRCC = 1.0 / (dDistIJ - daParams[off2b2]);
            final double dGRIJR = daParams[off2b2 + 1] * dRIJRCC;
            final double dFIJ = Math.exp(dGRIJR);

            final double dRIKRCC = 1.0 / (dDistIK - daParams[off2b3]);
            final double dGRIKR = daParams[off2b3 + 1] * dRIKRCC;
            final double dFIK = Math.exp(dGRIKR);

            final double dCTIJK = (dDistIJSq + dDistJKSq - dDistIKSq) / (2.0 * dDistIJ * dDistJK);
            final double dCTIKJ = (dDistIKSq + dDistJKSq - dDistIJSq) / (2.0 * dDistIK * dDistJK);
            final double dCTJIK = (dDistIJSq + dDistIKSq - dDistJKSq) / (2.0 * dDistIJ * dDistIK);

            final double dXJIK = dCTJIK + 1.0 / 3.0;
            final double dXJIK2 = dXJIK * dXJIK;

            final double dYJIK = dCTJIK + daParams[off3b];
            final double dYJIK2 = dYJIK * dYJIK + daParams[off3b + 1];

            final double dXIJK = dCTIJK + 1.0 / 3.0;
            final double dXIJK2 = dXIJK * dXIJK;

            final double dYIJK = dCTIJK + daParams[off3b];
            final double dYIJK2 = dYIJK * dYIJK + daParams[off3b + 1];

            final double dXIKJ = dCTIKJ + 1.0 / 3.0;
            final double dXIKJ2 = dXIKJ * dXIKJ;

            final double dYIKJ = dCTIKJ + daParams[off3b];
            final double dYIKJ2 = dYIKJ * dYIKJ + daParams[off3b + 1];

            final double dV31 = dFIJ * dFIK * dXJIK2 * dYJIK2;
            final double dV32 = dFIJ * dFJK * dXIJK2 * dYIJK2;
            final double dV33 = dFIK * dFJK * dXIKJ2 * dYIKJ2;

            final double dPJIK =
                daParams[off3b + 2]
                    * dFIJ
                    * dFIK
                    * 2.0
                    * dXJIK
                    * (dYJIK2 + dYJIK)
                    * dDistJK
                    / dDistIJ
                    / dDistIK;

            final double dPIJK =
                daParams[off3b + 2]
                    * dFIJ
                    * dFJK
                    * 2.0
                    * dXIJK
                    * (dYIJK2 + dYIJK)
                    * dDistIK
                    / dDistIJ
                    / dDistJK;

            final double dPIKJ =
                daParams[off3b + 2]
                    * dFIK
                    * dFJK
                    * 2.0
                    * dXIKJ
                    * (dYIKJ2 + dYIKJ)
                    * dDistIJ
                    / dDistIK
                    / dDistJK;

            final double dDJIKIJ = -daParams[off3b + 2] * dV31 * dGRIJR * dRIJRCC + dPJIK * dCTIJK;
            final double dDJIKIK = -daParams[off3b + 2] * dV31 * dGRIKR * dRIKRCC + dPJIK * dCTIKJ;
            final double dDJIKJK = -dPJIK;

            final double dDIJKIJ = -daParams[off3b + 2] * dV32 * dGRIJR * dRIJRCC + dPIJK * dCTJIK;
            final double dDIJKIK = -dPIJK;
            final double dDIJKJK = -daParams[off3b + 2] * dV32 * dGRJKR * dRJKRCC + dPIJK * dCTIKJ;

            final double dDIKJIJ = -dPIKJ;
            final double dDIKJIK = -daParams[off3b + 2] * dV33 * dGRIKR * dRIKRCC + dPIKJ * dCTJIK;
            final double dDJKJJK = -daParams[off3b + 2] * dV33 * dGRJKR * dRJKRCC + dPIKJ * dCTIJK;

            final double dDSumIJ = dDJIKIJ + dDIJKIJ + dDIKJIJ;
            final double dDSumIK = dDJIKIK + dDIJKIK + dDIKJIK;
            final double dDSumJK = dDJIKJK + dDIJKJK + dDJKJJK;

            gradientMat[0][i] += dDSumIJ * dIJX / dDistIJ + dDSumIK * dIKX / dDistIK;
            gradientMat[0][j] -= dDSumIJ * dIJX / dDistIJ + dDSumJK * dDistJKX / dDistJK;
            gradientMat[0][k] -= dDSumIK * dIKX / dDistIK - dDSumJK * dDistJKX / dDistJK;

            gradientMat[1][i] += dDSumIJ * dIJY / dDistIJ + dDSumIK * dIKY / dDistIK;
            gradientMat[1][j] -= dDSumIJ * dIJY / dDistIJ + dDSumJK * dDistJKY / dDistJK;
            gradientMat[1][k] -= dDSumIK * dIKY / dDistIK - dDSumJK * dDistJKY / dDistJK;

            gradientMat[2][i] += dDSumIJ * dIJZ / dDistIJ + dDSumIK * dIKZ / dDistIK;
            gradientMat[2][j] -= dDSumIJ * dIJZ / dDistIJ + dDSumJK * dDistJKZ / dDistJK;
            gradientMat[2][k] -= dDSumIK * dIKZ / dDistIK - dDSumJK * dDistJKZ / dDistJK;

            // TODO derivative of SWG term... test it!
          }
        }
      }
    }

    gradient.setGradientTotal(gradientMat);

    if (use3Body && saAtomTypes.length > 2) {
      // one last energy evaluation
      // XXX remove this!
      final double dEnergy =
          energyCalculation(
              lID,
              iIteration,
              daXYZ1D,
              saAtomTypes,
              atomNos,
              atsPerMol,
              energyparts,
              iNoOfAtoms,
              faCharges,
              iaSpins,
              bonds,
              hasRigidEnv);
      gradient.setTotalEnergy(dEnergy);
    } else {
      gradient.setTotalEnergy(twoBodyEnergy);
    }

    if (DEBUG) {
      // calculate a numerical gradient
      final Gradient numericalGrad =
          NumericalGradients.numericalGradient(
              lID,
              iIteration,
              daXYZ1D,
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

      final double[][] daNumGrad = numericalGrad.getTotalGradient();

      for (int i = 0; i < daNumGrad.length; i++) {
        for (int j = 0; j < daNumGrad[0].length; j++) {
          if (Math.abs(daNumGrad[i][j] - gradientMat[i][j]) >= 1E-5) {
            // we have a noticable difference
            System.err.println(
                "DEBUG: Analytical vs numerical gradient is: "
                    + gradientMat[i][j]
                    + " vs "
                    + daNumGrad[i][j]);
          } else {
            System.out.append("DEBUG: Analytical vs numerical gradient was fine. :-)");
          }
        }
      }
    }
  }

  @Override
  public double energyCalculation(
      final long lID,
      final int iIteration,
      final double[] daXYZ1D,
      final String[] saAtomTypes,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      final int iNoOfAtoms,
      final float[] faCharges,
      final short[] iaSpins,
      final BondInfo bonds,
      final boolean hasRigidEnv) {

    double dEnergy = 0.0;
    double d2Body = 0.0;

    if (DEBUG) {
      for (final String s : saAtomTypes) {
        System.out.println("DEBUG: " + s);
      }
    }

    // get in O(N)
    final double[] radii = new double[iNoOfAtoms];
    for (int i = 0; i < iNoOfAtoms; i++) {
      radii[i] = AtomicProperties.giveRadius(atomNos[i]);
    }

    double[] daFirstParams = null;
    double[] daSecParams;
    double[] daParams = params.getAllParamters();
    int counter = -1;

    if (useCaching && paramOffsetCache == null) {
      // initialize caches
      initializeCaches(params, saAtomTypes, use3Body);
    }

    int lastOffset = 0;
    for (int i = 0; i < atsPerMol.length - 1; i++) {
      lastOffset += atsPerMol[i];
    }

    final int firstLoopAtomNo = (hasRigidEnv) ? lastOffset - 1 : iNoOfAtoms - 1;
    for (int i = 0; i < firstLoopAtomNo; i++) {
      if (easyMix) {
        // first part of parameters
        daFirstParams = params.getParametersForKey(saAtomTypes[i]);
        if (daFirstParams == null) {
          continue;
        }
      }

      final double rad1 = radii[i];

      final double xi = daXYZ1D[i];
      final double yi = daXYZ1D[iNoOfAtoms + i];
      final double zi = daXYZ1D[2 * iNoOfAtoms + i];

      for (int j = i + 1; j < iNoOfAtoms; j++) {

        counter++;

        final double rad2 = radii[j];

        final double dX = xi - daXYZ1D[j];
        final double dY = yi - daXYZ1D[iNoOfAtoms + j];
        final double dZ = zi - daXYZ1D[2 * iNoOfAtoms + j];

        final double dDist = Math.sqrt(dX * dX + dY * dY + dZ * dZ);

        // check whether one should use a cutoff
        final double dAddedRadii = rad1 + rad2;

        if (dDist < (closeCutBlow * dAddedRadii)) {
          // use cutoff.
          d2Body += FixedValues.NONCONVERGEDENERGY;
          continue;
        } else if (dDist > (farCutBlow * dAddedRadii)) {
          continue;
        }

        int offset = 0;
        if (easyMix) {
          // second part of parameters
          daSecParams = params.getParametersForKey(saAtomTypes[j]);
          if (daSecParams == null) {
            continue;
          }

          // mix them
          daParams = mixParams(daFirstParams, daSecParams);
        } else {
          if (!useCaching) {
            String sPair = saAtomTypes[i] + saAtomTypes[j];
            daParams = params.getParametersForKey(sPair);
            if (daParams == null) {
              sPair = saAtomTypes[j] + saAtomTypes[i];
              daParams = params.getParametersForKey(sPair);
              if (daParams == null) {
                // ok this IS a problem
                System.err.println(
                    "WARNING: We can't find parameters for the pair "
                        + saAtomTypes[i]
                        + " and "
                        + saAtomTypes[j]
                        + ".");
                continue;
              }
            }
          } else {
            offset = paramOffsetCache[counter];
          }
        }

        // calculate the contributions
        final double epsFac = 4 * daParams[offset];
        final double distInv = 1.0 / dDist;

        int iCounter = offset + 1;
        for (int iPot = startPot; iPot <= endPot; iPot = iPot + potIncr) {
          final double dSign = Math.signum(daParams[iCounter]);
          d2Body += epsFac * dSign * pow((daParams[iCounter] * distInv), iPot);
          iCounter++;
        }
      }
    }

    if (!discard2Body) {
      dEnergy += d2Body;
    }

    /*
     * Stillinger-Weber-Gong term
     */
    if (use3Body && iNoOfAtoms > 2) {

      // put the 1D coordinates into 3D ones - XXX remove
      final double[][] daXYZ = new double[3][iNoOfAtoms];
      System.arraycopy(daXYZ1D, 0, daXYZ[0], 0, iNoOfAtoms);
      System.arraycopy(daXYZ1D, iNoOfAtoms, daXYZ[1], 0, iNoOfAtoms);
      System.arraycopy(daXYZ1D, iNoOfAtoms * 2, daXYZ[2], 0, iNoOfAtoms);

      int off2b1;
      int off2b2;
      int off2b3;
      int off3b;
      counter = -1;
      int c = -1;
      daParams = params.getAllParamters();

      for (int k = 0; k < iNoOfAtoms; k++) {
        for (int j = 0; j < k; j++) {

          c++;
          if (!useCaching) {
            off2b1 = posOfKey2sin3(saAtomTypes, k, j, params);
          } else {
            off2b1 = paramOffsetCache2b1[c];
          }

          final double dJKX = daXYZ[0][k] - daXYZ[0][j];
          final double dJKY = daXYZ[1][k] - daXYZ[1][j];
          final double dJKZ = daXYZ[2][k] - daXYZ[2][j];

          final double dDistJKSq = dJKX * dJKX + dJKY * dJKY + dJKZ * dJKZ;
          final double dDistJK = Math.sqrt(dDistJKSq);

          // ok, check whether the distance is above the equilibrium distance for this pair
          if (dDistJK > daParams[off2b1]) {
            continue;
          }

          final double dRJKRCC = 1.0 / (Math.sqrt(dDistJKSq) - daParams[off2b1]);
          final double dGRJKR = daParams[off2b1 + 1] * dRJKRCC;
          final double dFJK = Math.exp(dGRJKR);

          for (int i = 0; i < j; i++) {

            counter++;

            if (!useCaching) {
              off3b = posOfKey3in3(saAtomTypes, i, j, k, params);
              off2b2 = posOfKey2sin3(saAtomTypes, i, j, params);
              off2b3 = posOfKey2sin3(saAtomTypes, i, k, params);
            } else {
              off3b = paramOffsetCache3b[2][counter];
              off2b2 = paramOffsetCache3b[0][counter];
              off2b3 = paramOffsetCache3b[1][counter];
            }

            final double dIJX = daXYZ[0][i] - daXYZ[0][j];
            final double dIJY = daXYZ[1][i] - daXYZ[1][j];
            final double dIJZ = daXYZ[2][i] - daXYZ[2][j];

            final double dDistIJSq = dIJX * dIJX + dIJY * dIJY + dIJZ * dIJZ;

            final double dIKX = daXYZ[0][i] - daXYZ[0][k];
            final double dIKY = daXYZ[1][i] - daXYZ[1][k];
            final double dIKZ = daXYZ[2][i] - daXYZ[2][k];

            final double dDistIKSq = dIKX * dIKX + dIKY * dIKY + dIKZ * dIKZ;

            // ok, check whether the distance is above the equilibrium distance for this pair
            if (dDistIJSq > daParams[off2b2] * daParams[off2b2]) {
              continue;
            }
            if (dDistIKSq > daParams[off2b3] * daParams[off2b3]) {
              continue;
            }

            final double dDistIJ = Math.sqrt(dDistIJSq);
            final double dDistIK = Math.sqrt(dDistIKSq);

            // now we can finally start to evaluate the contributions
            final double dFIJ = Math.exp(daParams[off2b2 + 1] / (dDistIJ - daParams[off2b2]));
            final double dFIK = Math.exp(daParams[off2b3 + 1] / (dDistIK - daParams[off2b3]));

            final double dCTIJK = (dDistIJSq + dDistJKSq - dDistIKSq) / (2.0 * dDistIJ * dDistJK);
            final double dCTIKJ = (dDistIKSq + dDistJKSq - dDistIJSq) / (2.0 * dDistIK * dDistJK);
            final double dCTJIK = (dDistIJSq + dDistIKSq - dDistJKSq) / (2.0 * dDistIJ * dDistIK);

            final double t1 = dCTJIK + 1.0 / 3.0;
            final double dXJIK = t1 * t1;
            final double t2 = dCTJIK + daParams[off3b];
            final double dYJIK = t2 * t2 + daParams[off3b + 1];

            final double t3 = dCTIJK + 1.0 / 3.0;
            final double dXIJK = t3 * t3;
            final double t4 = dCTIJK + daParams[off3b];
            final double dYIJK = t4 * t4 + daParams[off3b + 1];

            final double t5 = dCTIKJ + 1.0 / 3.0;
            final double dXIKJ = t5 * t5;
            final double t6 = dCTIKJ + daParams[off3b];
            final double dYIKJ = t6 * t6 + daParams[off3b + 1];

            dEnergy +=
                daParams[off3b + 2]
                    * (dFIJ * dFIK * dXJIK * dYJIK
                        + dFIJ * dFJK * dXIJK * dYIJK
                        + dFIK * dFJK * dXIKJ * dYIKJ);
          }
        }
      }
    }

    return dEnergy;
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {
    return energyOfStructWithParams(
        cartes,
        params,
        geomID,
        bonds,
        cartes.containedEnvType() == CartesianCoordinates.ENVTYPE.RIGID);
  }

  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds,
      final boolean hasRigidEnv) {

    final int iNoOfAtoms = cartes.getNoOfAtoms();
    final String[] saAtoms = cartes.getAllAtomTypes();
    final double[][] daXYZ = cartes.getAllXYZCoord();

    double dEnergy = 0.0;
    double d2Body = 0.0;

    if (DEBUG) {
      for (final String s : saAtoms) {
        System.out.println("DEBUG: " + s);
      }
    }

    final short[] atomNos = cartes.getAllAtomNumbers();
    // get in O(N)
    final double[] radii = new double[iNoOfAtoms];
    for (int i = 0; i < iNoOfAtoms; i++) {
      radii[i] = AtomicProperties.giveRadius(atomNos[i]);
    }

    double[] daFirstParams = null;
    double[] daSecParams;
    double[] daParams = params.getAllParamters();
    int counter = -1;

    if (useCaching && paramOffsetCache == null) {
      // initialize caches
      initializeCaches(params, saAtoms, use3Body);
    }

    final int[] atsPerMol = cartes.getAllAtomsPerMol();
    int lastOffset = 0;
    for (int i = 0; i < atsPerMol.length - 1; i++) {
      lastOffset += atsPerMol[i];
    }

    final int firstLoopAtomNo = (hasRigidEnv) ? lastOffset - 1 : iNoOfAtoms - 1;
    for (int i = 0; i < firstLoopAtomNo; i++) {
      if (easyMix) {
        // first part of parameters
        daFirstParams = params.getParametersForKey(saAtoms[i]);
        if (daFirstParams == null) {
          continue;
        }
      }

      final double rad1 = radii[i];

      final double xi = daXYZ[0][i];
      final double yi = daXYZ[1][i];
      final double zi = daXYZ[2][i];

      for (int j = i + 1; j < iNoOfAtoms; j++) {

        counter++;

        final double rad2 = radii[j];

        final double dX = xi - daXYZ[0][j];
        final double dY = yi - daXYZ[1][j];
        final double dZ = zi - daXYZ[2][j];

        final double dDist = Math.sqrt(dX * dX + dY * dY + dZ * dZ);

        // check whether one should use a cutoff
        final double dAddedRadii = rad1 + rad2;

        if (dDist < (closeCutBlow * dAddedRadii)) {
          // use cutoff.
          d2Body += FixedValues.NONCONVERGEDENERGY;
          continue;
        } else if (dDist > (farCutBlow * dAddedRadii)) {
          continue;
        }

        int offset = 0;
        if (easyMix) {
          // second part of parameters
          daSecParams = params.getParametersForKey(saAtoms[j]);
          if (daSecParams == null) {
            continue;
          }

          // mix them
          daParams = mixParams(daFirstParams, daSecParams);
        } else {
          if (!useCaching) {
            String sPair = saAtoms[i] + saAtoms[j];
            daParams = params.getParametersForKey(sPair);
            if (daParams == null) {
              sPair = saAtoms[j] + saAtoms[i];
              daParams = params.getParametersForKey(sPair);
              if (daParams == null) {
                // ok this IS a problem
                System.err.println(
                    "WARNING: We can't find parameters for the pair "
                        + saAtoms[i]
                        + " and "
                        + saAtoms[j]
                        + ".");
                continue;
              }
            }
          } else {
            offset = paramOffsetCache[counter];
          }
        }

        // calculate the contributions
        final double epsFac = 4 * daParams[offset];
        final double distInv = 1.0 / dDist;

        int iCounter = offset + 1;
        for (int iPot = startPot; iPot <= endPot; iPot = iPot + potIncr) {
          final double dSign = Math.signum(daParams[iCounter]);
          d2Body += epsFac * dSign * pow((daParams[iCounter] * distInv), iPot);
          iCounter++;
        }
      }
    }

    if (!discard2Body) {
      dEnergy += d2Body;
    }

    /*
     * Stillinger-Weber-Gong term
     */
    if (use3Body && cartes.getNoOfAtoms() > 2) {
      int off2b1;
      int off2b2;
      int off2b3;
      int off3b;
      counter = -1;
      int c = -1;
      daParams = params.getAllParamters();

      // XXX we do not use hasRigidEnv for 3 body

      for (int k = 0; k < iNoOfAtoms; k++) {
        for (int j = 0; j < k; j++) {

          c++;
          if (!useCaching) {
            off2b1 = posOfKey2sin3(saAtoms, k, j, params);
          } else {
            off2b1 = paramOffsetCache2b1[c];
          }

          final double dJKX = daXYZ[0][k] - daXYZ[0][j];
          final double dJKY = daXYZ[1][k] - daXYZ[1][j];
          final double dJKZ = daXYZ[2][k] - daXYZ[2][j];

          final double dDistJKSq = dJKX * dJKX + dJKY * dJKY + dJKZ * dJKZ;
          final double dDistJK = Math.sqrt(dDistJKSq);

          // ok, check whether the distance is above the equilibrium distance for this pair
          if (dDistJK > daParams[off2b1]) {
            continue;
          }

          final double dRJKRCC = 1.0 / (Math.sqrt(dDistJKSq) - daParams[off2b1]);
          final double dGRJKR = daParams[off2b1 + 1] * dRJKRCC;
          final double dFJK = Math.exp(dGRJKR);

          for (int i = 0; i < j; i++) {

            counter++;

            if (!useCaching) {
              off3b = posOfKey3in3(saAtoms, i, j, k, params);
              off2b2 = posOfKey2sin3(saAtoms, i, j, params);
              off2b3 = posOfKey2sin3(saAtoms, i, k, params);
            } else {
              off3b = paramOffsetCache3b[2][counter];
              off2b2 = paramOffsetCache3b[0][counter];
              off2b3 = paramOffsetCache3b[1][counter];
            }

            final double dIJX = daXYZ[0][i] - daXYZ[0][j];
            final double dIJY = daXYZ[1][i] - daXYZ[1][j];
            final double dIJZ = daXYZ[2][i] - daXYZ[2][j];

            final double dDistIJSq = dIJX * dIJX + dIJY * dIJY + dIJZ * dIJZ;

            final double dIKX = daXYZ[0][i] - daXYZ[0][k];
            final double dIKY = daXYZ[1][i] - daXYZ[1][k];
            final double dIKZ = daXYZ[2][i] - daXYZ[2][k];

            final double dDistIKSq = dIKX * dIKX + dIKY * dIKY + dIKZ * dIKZ;

            // ok, check whether the distance is above the equilibrium distance for this pair
            if (dDistIJSq > daParams[off2b2] * daParams[off2b2]) {
              continue;
            }
            if (dDistIKSq > daParams[off2b3] * daParams[off2b3]) {
              continue;
            }

            final double dDistIJ = Math.sqrt(dDistIJSq);
            final double dDistIK = Math.sqrt(dDistIKSq);

            // now we can finally start to evaluate the contributions
            final double dFIJ = Math.exp(daParams[off2b2 + 1] / (dDistIJ - daParams[off2b2]));
            final double dFIK = Math.exp(daParams[off2b3 + 1] / (dDistIK - daParams[off2b3]));

            final double dCTIJK = (dDistIJSq + dDistJKSq - dDistIKSq) / (2.0 * dDistIJ * dDistJK);
            final double dCTIKJ = (dDistIKSq + dDistJKSq - dDistIJSq) / (2.0 * dDistIK * dDistJK);
            final double dCTJIK = (dDistIJSq + dDistIKSq - dDistJKSq) / (2.0 * dDistIJ * dDistIK);

            final double t1 = dCTJIK + 1.0 / 3.0;
            final double dXJIK = t1 * t1;
            final double t2 = dCTJIK + daParams[off3b];
            final double dYJIK = t2 * t2 + daParams[off3b + 1];

            final double t3 = dCTIJK + 1.0 / 3.0;
            final double dXIJK = t3 * t3;
            final double t4 = dCTIJK + daParams[off3b];
            final double dYIJK = t4 * t4 + daParams[off3b + 1];

            final double t5 = dCTIKJ + 1.0 / 3.0;
            final double dXIKJ = t5 * t5;
            final double t6 = dCTIKJ + daParams[off3b];
            final double dYIKJ = t6 * t6 + daParams[off3b + 1];

            dEnergy +=
                daParams[off3b + 2]
                    * (dFIJ * dFIK * dXJIK * dYJIK
                        + dFIJ * dFJK * dXIJK * dYIJK
                        + dFIK * dFJK * dXIKJ * dYIKJ);
          }
        }
      }
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

    // XXX just a numerical gradient ATM
    return NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, grad);
  }

  @Override
  public double[][] minMaxBordersForParams(final AdaptiveParameters params) {

    final String[] saKeys = params.getAllKeysCopy();
    final int iNoOfParams = params.getNumberOfParamters();

    final double[][] daBorders = new double[2][iNoOfParams];
    int iCounter = 0;
    int iParamCounter = 0;

    for (final String key : saKeys) {
      final int iNoForKey = params.getAmountOfParametersForKey(key);

      if (key.length() <= 4 && iNoForKey != 3) {

        if (iNoForKey == 2) {
          // this is the second part of the 3 body term
          // first parameter is a cutoff distance (in bohr)
          daBorders[0][iParamCounter] = 20;
          daBorders[1][iParamCounter] = 30;
          iParamCounter++;
          // second parameter is the gamma
          daBorders[0][iParamCounter] = -10;
          daBorders[1][iParamCounter] = 10;
          iParamCounter++;
        } else {
          // is a two-body set
          for (int j = 0; j < iNoForKey; j++) {
            if (iCounter == 0) {
              // first is an epsilon: speculative
              daBorders[0][iParamCounter] = 0.00001;
              daBorders[1][iParamCounter] = 0.00090;
              iCounter++;
              iParamCounter++;
            } else {
              // the sigmas are of course VERY speculative
              daBorders[0][iParamCounter] = -15.0;
              daBorders[1][iParamCounter] = 15.0;
              iCounter++;
              iParamCounter++;
            }

            if (iCounter >= iNoForKey) {
              iCounter = 0;
            }
          }
        }
      } else if (key.length() <= 6 && key.length() > 4 && iNoForKey == 3) {
        // is a three-body set
        for (int j = 0; j < iNoForKey; j++) {
          /*
           * since we somewhat try to resemble angles, 2pi should be a
           * minimum, we add something on top and done
           */
          daBorders[0][iParamCounter] = -10.0;
          daBorders[1][iParamCounter] = +10.0;
          iParamCounter++;
        }
      } else {
        System.err.println(
            "ERROR: Unknown case detected for key " + key + " no of params " + iNoForKey);
        for (int i = 0; i < iNoForKey; i++) {
          /*
           * Use huge default values
           */
          daBorders[0][iParamCounter] = -1000.0;
          daBorders[1][iParamCounter] = +1000.0;
          iParamCounter++;
        }
      }
    }

    return daBorders;
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      final ArrayList<CartesianCoordinates> refCartes, final String sMethod) {

    // loop through all reference geometries and figure a list of non-redundant atom tyes out
    final LinkedList<String> llAtoms = new LinkedList<>();

    final Iterator<CartesianCoordinates> itRefGeoms = refCartes.iterator();

    while (itRefGeoms.hasNext()) {
      final CartesianCoordinates cartesTemp = itRefGeoms.next();
      final String[] saAtoms = cartesTemp.getAllAtomTypes();

      for (final String saAtom : saAtoms) {
        if (!llAtoms.contains(saAtom)) {
          // add it to the list
          llAtoms.add(saAtom);
        }
      }
    }

    /*
     * triples
     */
    int iNoOfTriples = 0;
    if (use3Body) {
      // remember the symmetry of the term: HeNeKr is same as NeKrHe and so on
      for (int i = 0; i < llAtoms.size(); i++) {
        for (int j = 0; j <= i; j++) {
          for (int k = 0; k <= j; k++) {
            iNoOfTriples++;
          }
        }
      }
    }

    /*
     * pairs
     */

    String[] saAtoms;
    int[] iaParamsPerAt;

    int iParamSum = 0;

    if (easyMix) {

      if (!use3Body) {
        saAtoms = new String[llAtoms.size()];
        iaParamsPerAt = new int[llAtoms.size()];
      } else {
        // because of the additional two body param of the 3body term
        saAtoms = new String[llAtoms.size() * 2 + iNoOfTriples];
        iaParamsPerAt = new int[llAtoms.size() * 2 + iNoOfTriples];
      }

      // for each atom we need a set of parameters
      for (int i = 0; i < saAtoms.length; i++) {
        String sTempAtom = llAtoms.get(i);
        // independent of which atom, it is always an epsilon value and a set of sigmas
        iaParamsPerAt[i] = (int) Math.floor((endPot - startPot) / (double) potIncr) + 2;
        iParamSum += iaParamsPerAt[i];
        saAtoms[i] = sTempAtom;
      }
    } else {
      int iNoOfPairs = 0;
      for (int i = 0; i < llAtoms.size(); i++) {
        for (int j = 0; j <= i; j++) {
          iNoOfPairs++;
        }
      }

      if (!use3Body) {
        saAtoms = new String[iNoOfPairs];
        iaParamsPerAt = new int[iNoOfPairs];
      } else {
        // since we have the additional 2-Body parameters of the 3-Body term, crude but true
        saAtoms = new String[iNoOfPairs * 2 + iNoOfTriples];
        iaParamsPerAt = new int[iNoOfPairs * 2 + iNoOfTriples];
      }

      // we need a set of parameters for each pair
      final int iNoOfAtoms = llAtoms.size();

      int iCounter = 0;
      for (int i = 0; i < iNoOfAtoms; i++) {
        for (int j = i; j < iNoOfAtoms; j++) {
          final String sPair = llAtoms.get(i) + llAtoms.get(j);
          saAtoms[iCounter] = sPair;
          // independent of which pair, it is always an epsilon value and a set of sigmas
          iaParamsPerAt[iCounter] =
              (int) Math.floor((double) (endPot - startPot) / (double) potIncr) + 2;
          iParamSum += iaParamsPerAt[iCounter];
          iCounter++;
          if (use3Body) {
            saAtoms[iCounter] = sPair.toLowerCase();
            // R_ij, gamma_ij
            iaParamsPerAt[iCounter] = 2;
            iParamSum += 2;
            iCounter++;
          }
        }
      }
    }

    /*
     * again triples
     */
    if (use3Body) {

      int iCounter = saAtoms.length - iNoOfTriples;

      // remember the symmetry of the term: HeNeKr is same as NeKrHe and so on
      for (int i = 0; i < llAtoms.size(); i++) {
        for (int j = 0; j <= i; j++) {
          for (int k = 0; k <= j; k++) {
            final String sTriple = llAtoms.get(i) + llAtoms.get(j) + llAtoms.get(k);
            saAtoms[iCounter] = sTriple;
            // c_0, c_1, lambda
            iaParamsPerAt[iCounter] = 3;
            iParamSum += 3;
            iCounter++;
          }
        }
      }
    }

    final AdaptiveParameters paramStub =
        new AdaptiveParameters(iParamSum, -1, saAtoms, iaParamsPerAt, sMethod);

    return paramStub;
  }

  private double[] mixParams(final double[] daFirstParams, final double[] daSecParams) {

    final double[] daMixed = new double[daFirstParams.length];

    // standard Lorentz-Berthelot mixing

    // first one is always a LJ epsilon
    daMixed[0] = Math.sqrt(daFirstParams[0] * daSecParams[0]);

    // the rest are sigmas
    for (int i = 1; i < daFirstParams.length; i++) {
      daMixed[i] = 0.5 * (daFirstParams[i] + daSecParams[i]);
    }

    return daMixed;
  }

  private void initializeCaches(
      final AdaptiveParameters params, final String[] atoms, final boolean use3b) {

    final int noOfAtoms = atoms.length;
    this.paramOffsetCache = new int[(noOfAtoms * noOfAtoms - noOfAtoms) / 2];

    int counter = 0;
    for (int i = 0; i < noOfAtoms - 1; i++) {
      for (int j = i + 1; j < noOfAtoms; j++) {
        paramOffsetCache[counter] = posOfKey(atoms, i, j, params);
        counter++;
      }
    }

    if (!use3b) {
      return;
    }

    // how many 3b? I am too bored to calculate this now...
    counter = 0;
    int c = 0;
    for (int i = 0; i < atoms.length; i++) {
      for (int j = 0; j < i; j++) {
        c++;
        for (int k = 0; k < j; k++) {
          counter++;
        }
      }
    }

    this.paramOffsetCache3b = new int[3][counter];
    this.paramOffsetCache2b1 = new int[c];

    counter = 0;
    c = 0;
    for (int i = 0; i < atoms.length; i++) {
      for (int j = 0; j < i; j++) {
        paramOffsetCache2b1[c] = posOfKey2sin3(atoms, i, j, params);
        c++;
        for (int k = 0; k < j; k++) {
          // 2b part
          paramOffsetCache3b[0][counter] = posOfKey2sin3(atoms, k, j, params);
          paramOffsetCache3b[1][counter] = posOfKey2sin3(atoms, k, i, params);
          // 3b part
          paramOffsetCache3b[2][counter] = posOfKey3in3(atoms, i, j, k, params);
          counter++;
        }
      }
    }
  }

  private int posOfKey(
      final String[] atoms, final int i, final int j, final AdaptiveParameters params) {

    // get parameter position for this pair
    String pair = atoms[i] + atoms[j];
    int pos = params.getStartPointForKey(pair);
    if (pos >= 0) {
      return pos;
    }

    pair = atoms[j] + atoms[i];
    pos = params.getStartPointForKey(pair);
    if (pos >= 0) {
      return pos;
    }

    System.err.println("WARNING: No parameters for key " + pair + ". Returning wrong ones!");
    return 0;
  }

  private int posOfKey2sin3(
      final String[] atoms, final int i, final int j, final AdaptiveParameters params) {

    String d = (atoms[i] + atoms[j]).toLowerCase();
    int pos = params.getStartPointForKey(d);
    if (pos >= 0) {
      return pos;
    }

    d = (atoms[j] + atoms[i]).toLowerCase();
    pos = params.getStartPointForKey(d);
    if (pos >= 0) {
      return pos;
    }

    System.err.println("WARNING: No parameters for key " + d + ". Returning wrong ones!");
    return 0;
  }

  private int posOfKey3in3(
      final String[] atoms,
      final int i,
      final int j,
      final int k,
      final AdaptiveParameters params) {

    String t = atoms[i] + atoms[j] + atoms[k];
    int pos = params.getStartPointForKey(t);
    if (pos >= 0) {
      return pos;
    }

    t = atoms[i] + atoms[k] + atoms[j];
    pos = params.getStartPointForKey(t);
    if (pos >= 0) {
      return pos;
    }

    t = atoms[j] + atoms[i] + atoms[k];
    pos = params.getStartPointForKey(t);
    if (pos >= 0) {
      return pos;
    }

    t = atoms[j] + atoms[k] + atoms[i];
    pos = params.getStartPointForKey(t);
    if (pos >= 0) {
      return pos;
    }

    t = atoms[k] + atoms[i] + atoms[j];
    pos = params.getStartPointForKey(t);
    if (pos >= 0) {
      return pos;
    }

    t = atoms[k] + atoms[j] + atoms[i];
    pos = params.getStartPointForKey(t);
    if (pos >= 0) {
      return pos;
    }

    // we are through with all possibilities :-(
    System.err.println("WARNING: No parameters for key " + t + ". Returning wrong ones!");
    return 0;
  }
}
