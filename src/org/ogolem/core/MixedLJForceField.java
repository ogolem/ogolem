/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015-2021, J. M. Dieterich and B. Hartke
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

import java.util.Arrays;
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * This provides a backend for mixed atom type Lennard-Jones calculations. Uses Lorentz-Berthelot
 * combination rules.
 *
 * @author Johannes Dieterich
 * @version 2021-04-22
 */
public class MixedLJForceField implements CartesianFullBackend {

  // the ID
  private static final long serialVersionUID = (long) 20210422;

  private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

  private final boolean cache;
  private double[] eps;
  private double[] sig;

  public MixedLJForceField(final boolean useCaching) {
    this.cache = useCaching;
  }

  @Override
  public MixedLJForceField copy() {
    return new MixedLJForceField(cache);
  }

  @Override
  public String getMethodID() {
    return "LJ Force Field for heterogenous cases";
  }

  @Override
  public double energyCalculation(
      final long lID,
      final int iIteration,
      final double[] xyz1D,
      final String[] saAtomTypes,
      final short[] atomNos,
      final int[] atsPerMol,
      final double[] energyparts,
      final int iNoOfAtoms,
      final float[] faCharges,
      final short[] iaSpins,
      final BondInfo bonds,
      final boolean hasRigidEnv) {

    // zero the partial contributions
    Arrays.fill(energyparts, 0.0);

    // some cutoff constants
    final double t1 = 1.0 / 0.64;
    final double t1Sq = t1 * t1;
    final double t1Hex = t1Sq * t1Sq * t1Sq;
    final double t112 = t1Hex * t1Hex;

    // get all LJ parameters in O(N)
    if (!cache || eps == null) {
      eps = new double[iNoOfAtoms];
      sig = new double[iNoOfAtoms];
      for (int i = 0; i < iNoOfAtoms; i++) {
        if (atomNos[i] == 0) {
          continue;
        } // dummy
        eps[i] = AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[i]);
        sig[i] = AtomicProperties.giveLennardJonesSigma(saAtomTypes[i]);
      }
    }

    final var vConst2 = DoubleVector.broadcast(SPECIES, 10000.0);

    // get the squared distances between all atoms
    double dPotEnergyAdded = 0.0;

    int lastOffset = 0;
    for (int i = 0; i < atsPerMol.length - 1; i++) {
      lastOffset += atsPerMol[i];
    }

    final int firstLoopAtomNo = (hasRigidEnv) ? lastOffset - 1 : iNoOfAtoms - 1;
    for (int i = 0; i < firstLoopAtomNo; i++) {

      if (atomNos[i] == 0) {
        continue;
      } // dummy

      final double dEpsilon1 = eps[i];
      final var vEps1 = DoubleVector.broadcast(SPECIES, dEpsilon1);
      final double dSigma1 = sig[i];
      final var vSig1 = DoubleVector.broadcast(SPECIES, dSigma1);
      final double x0 = xyz1D[i];
      final var vX0 = DoubleVector.broadcast(SPECIES, x0);
      final double y0 = xyz1D[i + iNoOfAtoms];
      final var vY0 = DoubleVector.broadcast(SPECIES, y0);
      final double z0 = xyz1D[i + 2 * iNoOfAtoms];
      final var vZ0 = DoubleVector.broadcast(SPECIES, z0);

      var vPotEnergy = DoubleVector.zero(SPECIES);

      final int loopBound = SPECIES.loopBound(iNoOfAtoms - i - 1);
      int j = i + 1;
      for (; j < loopBound + i + 1; j += SPECIES.length()) {
        final var vEps2 = DoubleVector.fromArray(SPECIES, eps, j);
        final var vSig2 = DoubleVector.fromArray(SPECIES, sig, j);

        // mixing
        final var vEps = vEps2.mul(vEps1).lanewise(VectorOperators.SQRT);
        final var vSig = vSig1.add(vSig2).mul(0.5);

        // the cutoff distance
        final var vSeam = vSig.mul(0.64);
        final var vSeamSq = vSeam.mul(vSeam);

        final var vDX = vX0.sub(DoubleVector.fromArray(SPECIES, xyz1D, j));
        final var vDXSq = vDX.mul(vDX);
        final var vDY = vY0.sub(DoubleVector.fromArray(SPECIES, xyz1D, j + iNoOfAtoms));
        final var vDXYSq = vDY.fma(vDY, vDXSq);
        final var vDZ = vZ0.sub(DoubleVector.fromArray(SPECIES, xyz1D, j + 2 * iNoOfAtoms));
        final var vDistSq = vDZ.fma(vDZ, vDXYSq);
        final var vComp = vDistSq.lt(vSeamSq);
        final boolean anyLower = vComp.anyTrue();
        final var vInvRPow2 = vSig.mul(vSig).div(vDistSq);
        final var vInvRPow6 = vInvRPow2.mul(vInvRPow2).mul(vInvRPow2);
        final var vInvRPow12 = vInvRPow6.mul(vInvRPow6);
        var vVecTmp =
            vInvRPow12
                .sub(vInvRPow6)
                .mul(vEps)
                .mul(4); // XXX I bet one could do this with bitshifts
        final var vEP = DoubleVector.fromArray(SPECIES, energyparts, j);
        if (anyLower) {
          // the unlikely case
          final var vDist = vDistSq.lanewise(VectorOperators.SQRT);
          final var vConst1 = vEps.mul(4 * (t112 - t1Hex)).sub(10000.0).div(vSeam);
          final var vCutTmp = vDist.fma(vConst1, vConst2);
          vVecTmp = vVecTmp.blend(vCutTmp, vComp);
        }
        final var vEPAdd = vEP.add(vVecTmp);
        vEPAdd.intoArray(energyparts, j);
        vPotEnergy = vPotEnergy.add(vVecTmp);
      }

      final double redEnergy = vPotEnergy.reduceLanes(VectorOperators.ADD);
      energyparts[i] += redEnergy;
      dPotEnergyAdded += redEnergy;

      for (; j < iNoOfAtoms; j++) {

        if (atomNos[j] == 0) {
          continue;
        } // dummy

        // the Lennard-Jones parameters, set two
        final double dEpsilon2 = eps[j];
        final double dSigma2 = sig[j];

        final double dEpsilon = Math.sqrt(dEpsilon1 * dEpsilon2);
        final double dSigma = 0.5 * (dSigma1 + dSigma2);

        // the cutoff distance
        final double dSeam = 0.64 * dSigma;
        final double dSeamSquared = dSeam * dSeam;

        final double dDistX = x0 - xyz1D[j];
        final double dDistY = y0 - xyz1D[j + iNoOfAtoms];
        final double dDistZ = z0 - xyz1D[j + 2 * iNoOfAtoms];
        final double dDistSquared = dDistX * dDistX + dDistY * dDistY + dDistZ * dDistZ;
        if (dDistSquared > dSeamSquared) {
          final double dInvRPow2 = dSigma * dSigma / dDistSquared;
          final double dInvRPow6 = dInvRPow2 * dInvRPow2 * dInvRPow2;
          final double dInvRPow12 = dInvRPow6 * dInvRPow6;

          /*
           * The used Lennard-Jones potential is of the form
           * V(r_{ij})= 4\cdot\epsilon\cdot \left[\left(\dfrac{\sigma}{r_{ij}}\right)^{12}-\left(\dfrac{\sigma}{r_{ij}}^{6}\right)\right].
           */
          final double tmp = 4.0 * dEpsilon * (dInvRPow12 - dInvRPow6);
          energyparts[i] += tmp;
          energyparts[j] += tmp;
          dPotEnergyAdded += tmp;
        } else {
          // System.out.println("Atoms are too close together");

          // more constants... needed for cutting of the potential
          final double invSeam = 1.0 / dSeam;
          final double dConst1 = (4.0 * dEpsilon * (t112 - t1Hex) - 10000.0) * invSeam;
          final double dConst2 = 10000.0;

          final double dDist = Math.sqrt(dDistSquared);
          final double tmp = dConst1 * dDist + dConst2;
          energyparts[i] += tmp;
          energyparts[j] += tmp;
          dPotEnergyAdded += tmp;
        }
      }
    }

    return dPotEnergyAdded;
  }

  /**
   * same as above in blue... This provides the gradient, derived with exactly the same methods as
   * above.
   */
  @Override
  public void gradientCalculation(
      final long lID,
      final int iIteration,
      final double[] xyz1D,
      final String[] saAtomTypes,
      final short[] atomNos,
      final int[] atsPerMol,
      final double[] energyparts,
      final int iNoOfAtoms,
      final float[] faCharges,
      final short[] iaSpins,
      final BondInfo bonds,
      final Gradient gradient,
      final boolean hasRigidEnv) {

    // zero the partial contributions
    Arrays.fill(energyparts, 0.0);

    // the matrix for the gradient
    gradient.zeroGradient();
    final double[][] daGradientMat = gradient.getTotalGradient();

    // some cutoff constants
    final double t1 = 1.0 / 0.64;
    final double t1Sq = t1 * t1;
    final double t1Hex = t1Sq * t1Sq * t1Sq;
    final double t112 = t1Hex * t1Hex;

    final double t112_Hex = t112 - t1Hex;
    final var vT112_Hex = DoubleVector.broadcast(SPECIES, t112_Hex);

    // final var vT1Hex = DoubleVector.broadcast(SPECIES, t1Hex);
    // final var vT112 = DoubleVector.broadcast(SPECIES, t112);

    // get all LJ parameters in O(N)
    if (!cache || eps == null) {
      eps = new double[iNoOfAtoms];
      sig = new double[iNoOfAtoms];
      for (int i = 0; i < iNoOfAtoms; i++) {
        if (atomNos[i] == 0) {
          continue;
        } // dummy
        eps[i] = AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[i]);
        sig[i] = AtomicProperties.giveLennardJonesSigma(saAtomTypes[i]);
      }
    }

    final var vConst2 = DoubleVector.broadcast(SPECIES, 10000.0);
    final var one = DoubleVector.broadcast(SPECIES, 1);

    // calculate all pair distances
    double dPotEnergyAdded = 0.0;

    int lastOffset = 0;
    for (int i = 0; i < atsPerMol.length - 1; i++) {
      lastOffset += atsPerMol[i];
    }

    final int firstLoopAtomNo = (hasRigidEnv) ? lastOffset - 1 : iNoOfAtoms - 1;
    for (int i = 0; i < firstLoopAtomNo; i++) {

      if (atomNos[i] == 0) {
        continue;
      } // dummy

      final double dEpsilon1 = eps[i];
      final var vEps1 = DoubleVector.broadcast(SPECIES, dEpsilon1);
      final double dSigma1 = sig[i];
      final var vSig1 = DoubleVector.broadcast(SPECIES, dSigma1);
      final double x0 = xyz1D[i];
      final var vX0 = DoubleVector.broadcast(SPECIES, x0);
      final double y0 = xyz1D[i + iNoOfAtoms];
      final var vY0 = DoubleVector.broadcast(SPECIES, y0);
      final double z0 = xyz1D[i + 2 * iNoOfAtoms];
      final var vZ0 = DoubleVector.broadcast(SPECIES, z0);

      var vGradXI = DoubleVector.zero(SPECIES);
      var vGradYI = DoubleVector.zero(SPECIES);
      var vGradZI = DoubleVector.zero(SPECIES);

      var vPotEnergy = DoubleVector.zero(SPECIES);

      final int loopBound = SPECIES.loopBound(iNoOfAtoms - i - 1);
      int j = i + 1;
      for (; j < loopBound + i + 1; j += SPECIES.length()) {
        final var vEps2 = DoubleVector.fromArray(SPECIES, eps, j);
        final var vSig2 = DoubleVector.fromArray(SPECIES, sig, j);

        // mixing
        final var vEps = vEps2.mul(vEps1).lanewise(VectorOperators.SQRT);
        final var vSig = vSig1.add(vSig2).mul(0.5);

        // the cutoff distance
        final var vSeam = vSig.mul(0.64);
        final var vSeamSq = vSeam.mul(vSeam);

        final var vDX = vX0.sub(DoubleVector.fromArray(SPECIES, xyz1D, j));
        final var vDXSq = vDX.mul(vDX);
        final var vDY = vY0.sub(DoubleVector.fromArray(SPECIES, xyz1D, j + iNoOfAtoms));
        final var vDXYSq = vDY.fma(vDY, vDXSq);
        final var vDZ = vZ0.sub(DoubleVector.fromArray(SPECIES, xyz1D, j + 2 * iNoOfAtoms));
        final var vDistSq = vDZ.fma(vDZ, vDXYSq);
        final var vDist = vDistSq.lanewise(VectorOperators.SQRT);
        final var vDistInv = one.div(vDist);
        final var vDivProdX = vDX.mul(vDistInv);
        final var vDivProdY = vDY.mul(vDistInv);
        final var vDivProdZ = vDZ.mul(vDistInv);
        final var vComp = vDistSq.lt(vSeamSq);
        final boolean anyLower = vComp.anyTrue();
        if (anyLower) {
          // any atoms too close - the unlikely case
          if (GlobalConfig.DEBUGLEVEL > 0)
            System.err.println(
                "WARNING: LJ gradient: Atoms too close together, we take the cutoff potential.");

          final var vCst1 = vEps.mul(4).mul(vT112_Hex).sub(10000);
          final var vConst1 = vCst1.div(vSeam);
          final var vTmp = vConst1.mul(vDist).add(vConst2);

          // regular energy contributions
          final var vInvRPow2 = vSig.mul(vSig).div(vDistSq);
          final var vInvRPow6 = vInvRPow2.mul(vInvRPow2).mul(vInvRPow2);
          final var vInvRPow12 = vInvRPow6.mul(vInvRPow6);
          final var vVecTmpUncorr =
              vInvRPow12
                  .sub(vInvRPow6)
                  .mul(vEps)
                  .mul(4); // XXX I bet one could do this with bitshifts
          // blend vVecTmp w/ the cutoff potential
          final var vVecTmp = vVecTmpUncorr.blend(vTmp, vComp);

          final var vEP = DoubleVector.fromArray(SPECIES, energyparts, j);
          final var vEPAdd = vEP.add(vVecTmp);
          vEPAdd.intoArray(energyparts, j);
          vPotEnergy = vPotEnergy.add(vVecTmp);
          // dTemp = -48.0 * dEpsilon * dInvRPow12 * dDistInv + 24.0 * dEpsilon * dInvRPow6 *
          // dDistInv;
          final var vDistGrad = vEps.mul(vDistInv).mul(24);
          final var vGrad2 = vInvRPow12.mul(vDistGrad).mul(-2);
          final var vGradTmp = vInvRPow6.fma(vDistGrad, vGrad2);
          final var vNegGradTmp = vGradTmp.neg();

          // don't bother blending the gradient - even on the scalar path we're leaving them be
          // so on the vectorpath, they'll be what they are - but the energy will be sane and cutoff

          vGradXI = vGradTmp.fma(vDivProdX, vGradXI);
          final var vGOJ = DoubleVector.fromArray(SPECIES, daGradientMat[0], j);
          vNegGradTmp.fma(vDivProdX, vGOJ).intoArray(daGradientMat[0], j);

          vGradYI = vGradTmp.fma(vDivProdY, vGradYI);
          final var vG1J = DoubleVector.fromArray(SPECIES, daGradientMat[1], j);
          vNegGradTmp.fma(vDivProdY, vG1J).intoArray(daGradientMat[1], j);

          vGradZI = vGradTmp.fma(vDivProdZ, vGradZI);
          final var vG2J = DoubleVector.fromArray(SPECIES, daGradientMat[2], j);
          vNegGradTmp.fma(vDivProdZ, vG2J).intoArray(daGradientMat[2], j);
        } else {
          final var vInvRPow2 = vSig.mul(vSig).div(vDistSq);
          final var vInvRPow6 = vInvRPow2.mul(vInvRPow2).mul(vInvRPow2);
          final var vInvRPow12 = vInvRPow6.mul(vInvRPow6);
          final var vVecTmp =
              vInvRPow12
                  .sub(vInvRPow6)
                  .mul(vEps)
                  .mul(4); // XXX I bet one could do this with bitshifts
          final var vEP = DoubleVector.fromArray(SPECIES, energyparts, j);
          final var vEPAdd = vEP.add(vVecTmp);
          vEPAdd.intoArray(energyparts, j);
          vPotEnergy = vPotEnergy.add(vVecTmp);
          // dTemp = -48.0 * dEpsilon * dInvRPow12 * dDistInv + 24.0 * dEpsilon * dInvRPow6 *
          // dDistInv;
          final var vDistGrad = vEps.mul(vDistInv).mul(24);
          final var vGrad2 = vInvRPow12.mul(vDistGrad).mul(-2);
          final var vGradTmp = vInvRPow6.fma(vDistGrad, vGrad2);
          final var vNegGradTmp = vGradTmp.neg();

          vGradXI = vGradTmp.fma(vDivProdX, vGradXI);
          final var vGOJ = DoubleVector.fromArray(SPECIES, daGradientMat[0], j);
          vNegGradTmp.fma(vDivProdX, vGOJ).intoArray(daGradientMat[0], j);

          vGradYI = vGradTmp.fma(vDivProdY, vGradYI);
          final var vG1J = DoubleVector.fromArray(SPECIES, daGradientMat[1], j);
          vNegGradTmp.fma(vDivProdY, vG1J).intoArray(daGradientMat[1], j);

          vGradZI = vGradTmp.fma(vDivProdZ, vGradZI);
          final var vG2J = DoubleVector.fromArray(SPECIES, daGradientMat[2], j);
          vNegGradTmp.fma(vDivProdZ, vG2J).intoArray(daGradientMat[2], j);
        }
      }

      final double redEnergy = vPotEnergy.reduceLanes(VectorOperators.ADD);
      energyparts[i] += redEnergy;
      dPotEnergyAdded += redEnergy;

      double gradXI = daGradientMat[0][i] + vGradXI.reduceLanes(VectorOperators.ADD);
      double gradYI = daGradientMat[1][i] + vGradYI.reduceLanes(VectorOperators.ADD);
      double gradZI = daGradientMat[2][i] + vGradZI.reduceLanes(VectorOperators.ADD);

      for (; j < iNoOfAtoms; j++) {

        if (atomNos[j] == 0) {
          continue;
        } // dummy

        // the Lennard-Jones parameters, set two
        final double dEpsilon2 = eps[j];
        final double dSigma2 = sig[j];

        final double dEpsilon = Math.sqrt(dEpsilon1 * dEpsilon2);
        final double dSigma = 0.5 * (dSigma1 + dSigma2);

        // the cutoff distance
        final double dSeam = 0.64 * dSigma;
        final double dSeamSquared = dSeam * dSeam;

        final double dDistX = x0 - xyz1D[j];
        final double dDistY = y0 - xyz1D[j + iNoOfAtoms];
        final double dDistZ = z0 - xyz1D[j + 2 * iNoOfAtoms];
        final double dDistSquared = dDistX * dDistX + dDistY * dDistY + dDistZ * dDistZ;

        final double dDist = Math.sqrt(dDistSquared);
        final double dDistInv = 1.0 / dDist;
        // divide the distances in all three dimensions by the total distance to cover for the
        // coordinate system
        final double dDivProdX = dDistX * dDistInv;
        final double dDivProdY = dDistY * dDistInv;
        final double dDivProdZ = dDistZ * dDistInv;

        double dTemp;
        if (dDistSquared > dSeamSquared) {
          final double dInvRPow2 = dSigma * dSigma / dDistSquared;
          final double dInvRPow6 = dInvRPow2 * dInvRPow2 * dInvRPow2;
          final double dInvRPow12 = dInvRPow6 * dInvRPow6;

          /*
           * The used Lennard-Jones potential is of the form
           * V(r_{ij})= 4\cdot\epsilon\cdot \left[\left(\dfrac{\sigma}{r_{ij}}\right)^{12}-\left(\dfrac{\sigma}{r_{ij}}^{6}\right)\right].
           */

          final double tmp = 4.0 * dEpsilon * (dInvRPow12 - dInvRPow6);
          energyparts[i] += tmp;
          energyparts[j] += tmp;
          dPotEnergyAdded += tmp;
          dTemp = -48.0 * dEpsilon * dInvRPow12 * dDistInv + 24.0 * dEpsilon * dInvRPow6 * dDistInv;
        } else {
          if (GlobalConfig.DEBUGLEVEL > 0)
            System.err.println(
                "WARNING: LJ gradient: Atoms too close together, we take the cutoff potential.");

          // more constants... needed for cutting of the potential
          final double dConst1 = (4.0 * dEpsilon * (t112 - t1Hex) - 10000.0) / dSeam;
          final double dConst2 = 10000.0;

          final double tmp = dConst1 * dDist + dConst2;
          energyparts[i] += tmp;
          energyparts[j] += tmp;
          dPotEnergyAdded += tmp;
          dTemp = dConst1;
        }
        gradXI += dTemp * dDivProdX;
        daGradientMat[0][j] -= dTemp * dDivProdX;
        gradYI += dTemp * dDivProdY;
        daGradientMat[1][j] -= dTemp * dDivProdY;
        gradZI += dTemp * dDivProdZ;
        daGradientMat[2][j] -= dTemp * dDivProdZ;
      }

      daGradientMat[0][i] = gradXI;
      daGradientMat[1][i] = gradYI;
      daGradientMat[2][i] = gradZI;
    }

    gradient.setTotalEnergy(dPotEnergyAdded);
  }
}
