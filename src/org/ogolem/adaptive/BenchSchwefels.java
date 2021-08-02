/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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
package org.ogolem.adaptive;

import java.util.ArrayList;
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;

/**
 * A benchmark function, Schwefels function in n-D. Optimal solution: 420.9687. Requires a single
 * reference geometry with a single dummy atom. Benchmarks are always ugly. ;-)
 *
 * @author Johannes Dieterich
 * @version 2021-07-21
 */
public final class BenchSchwefels extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20210721;
  private static final boolean DEBUG = false;
  private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;
  private final int iDims;

  public BenchSchwefels(final int iDimensionality) {
    assert (iDimensionality > 0);
    this.iDims = iDimensionality;
  }

  @Override
  public BenchSchwefels copy() {
    return new BenchSchwefels(iDims);
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {
    // f(x(0)...x(n)) = 418.9829*n + sum(-x(i)*sin(sqrt(abs(x(i))))

    // we of course do not actually use the cartesian coordinates, just the
    // parameters
    final double[] daParams = params.getAllParamters();

    final int upperBound = SPECIES.loopBound(iDims);
    final var vCutVal = DoubleVector.broadcast(SPECIES, 500);
    final var vCutCons = DoubleVector.broadcast(SPECIES, 0.02);

    var vRes = DoubleVector.broadcast(SPECIES, 0.0);
    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vPI = DoubleVector.fromArray(SPECIES, daParams, i);
      // first check if we need the cutoff
      final var vPIAbs = vPI.abs();
      final var vCutMask = vPIAbs.compare(VectorOperators.GT, vCutVal);
      final var anyCut = vCutMask.anyTrue();
      if (anyCut) {
        // cutoff
        final var vCut = vPI.mul(vPI).mul(vCutCons).neg(); // to avoid the sub() below

        // regular
        final var vPAbsSq = vPIAbs.sqrt();
        final var vResLocUn = vPI.mul(vPAbsSq.lanewise(VectorOperators.SIN));

        // blend
        final var vResLoc = vResLocUn.blend(vCut, vCutMask);
        vRes = vRes.sub(vResLoc);

      } else {
        final var vPAbsSq = vPIAbs.sqrt();
        final var vResLoc = vPI.mul(vPAbsSq.lanewise(VectorOperators.SIN));
        vRes = vRes.sub(vResLoc);
      }
    }

    double dEnergy = 418.9829 * iDims + vRes.reduceLanes(VectorOperators.ADD);
    for (; i < iDims; i++) {
      if (daParams[i] > 500.0 || daParams[i] < -500.0) {
        // take the cutoff potential
        dEnergy += 0.02 * daParams[i] * daParams[i];
      } else {
        dEnergy -= daParams[i] * Math.sin(Math.sqrt(Math.abs(daParams[i])));
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
      final double[] daGrad) {

    // we of course do not actually use the cartesian coordinates, just the
    // parameters

    final int iNoOfParams = params.getNumberOfParamters();

    final double[] daParams = params.getAllParamters();

    final int upperBound = SPECIES.loopBound(iDims);
    final var vCutVal = DoubleVector.broadcast(SPECIES, 500);
    final var vCutCons = DoubleVector.broadcast(SPECIES, 0.02);
    final var vCutCons2 = DoubleVector.broadcast(SPECIES, 0.04);
    final var vTwo = DoubleVector.broadcast(SPECIES, 2);

    var vRes = DoubleVector.broadcast(SPECIES, 0.0);
    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vPI = DoubleVector.fromArray(SPECIES, daParams, i);
      // first check if we need the cutoff
      final var vPIAbs = vPI.abs();
      final var vCutMask = vPIAbs.compare(VectorOperators.GT, vCutVal);
      final var anyCut = vCutMask.anyTrue();
      if (anyCut) {
        // cutoff
        final var vCut = vPI.mul(vPI).mul(vCutCons).neg(); // to avoid the sub() below
        final var vCutGrad = vPI.mul(vCutCons2);

        // regular
        final var vPAbsR = vPIAbs.sqrt();
        final var vT1 = vPAbsR.lanewise(VectorOperators.SIN);
        final var vResLocUn = vPI.mul(vT1);

        // blend
        final var vResLoc = vResLocUn.blend(vCut, vCutMask);
        vRes = vRes.sub(vResLoc);

        final var vNumSec = vPIAbs.neg().mul(vPAbsR.lanewise(VectorOperators.COS));
        final var vDenomSec = vPAbsR.mul(vTwo);

        final var vSec = vNumSec.div(vDenomSec);
        final var vAnalyGradUn = vT1.neg().add(vSec);

        // blend
        final var vAnalyGrad = vAnalyGradUn.blend(vCutGrad, vCutMask);
        vAnalyGrad.intoArray(daGrad, i);

      } else {
        final var vPAbsR = vPIAbs.sqrt();
        final var vT1 = vPAbsR.lanewise(VectorOperators.SIN);
        final var vResLoc = vPI.mul(vT1);
        vRes = vRes.sub(vResLoc);

        final var vNumSec = vPIAbs.neg().mul(vPAbsR.lanewise(VectorOperators.COS));
        final var vDenomSec = vPAbsR.mul(vTwo);

        final var vSec = vNumSec.div(vDenomSec);
        final var vAnalyGrad = vT1.neg().add(vSec);
        vAnalyGrad.intoArray(daGrad, i);
      }
    }

    double dEnergy = 418.9829 * iDims + vRes.reduceLanes(VectorOperators.ADD);
    for (; i < iNoOfParams; i++) {

      if (daParams[i] > 500.0 || daParams[i] < -500.0) {
        // take the deviation of the cutoff potential
        daGrad[i] = 0.04 * daParams[i];
        dEnergy += 0.02 * daParams[i] * daParams[i];
      } else {
        /*
         * f'(x(i)) = -(sqrt(x^2) cos(sqrt(|x|)))/(2 sqrt(|x|))-sin(sqrt(|x|))
         * done using wolframalpha
         */

        final double pi = daParams[i];
        final double piAbs = Math.abs(pi);
        final double piAbsR = Math.sqrt(piAbs);
        final double t1 = Math.sin(piAbsR);

        final double dNumSecTerm = -piAbs * Math.cos(piAbsR);
        final double dDenomSecTerm = 2.0 * piAbsR;

        final double dAnalyticalVal = -t1 + dNumSecTerm / dDenomSecTerm;

        daGrad[i] = dAnalyticalVal;
        dEnergy -= daParams[i] * t1;
      }
    }

    if (DEBUG) {
      final double[] numGrad = new double[daGrad.length];

      final double numFunc =
          NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, numGrad);
      System.out.println("DEBUG OUTPUT FOR NUMGRAD");
      System.out.println("Func " + numFunc + "\t" + dEnergy + "\t" + (numFunc - dEnergy));
      for (int x = 0; x < daGrad.length; x++) {
        System.out.println(numGrad[x] + "\t" + daGrad[x] + "\t" + (numGrad[x] - daGrad[x]));
        daGrad[x] = numGrad[x];
      }

      return numFunc;
    }

    return dEnergy;
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      final ArrayList<CartesianCoordinates> refCartes, final String sMethod) {

    final String[] saAtoms = {"XX"};
    final int[] iaParamsPerAt = {iDims};
    final AdaptiveParameters paramStub =
        new AdaptiveParameters(iDims, -1, saAtoms, iaParamsPerAt, sMethod);

    return paramStub;
  }

  @Override
  public double[][] minMaxBordersForParams(final AdaptiveParameters params) {

    final double[][] daMinMax = new double[2][iDims];

    for (int i = 0; i < iDims; i++) {
      daMinMax[0][i] = -500;
      daMinMax[1][i] = 500;
    }

    return daMinMax;
  }
}
