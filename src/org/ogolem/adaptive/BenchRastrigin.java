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
 * Rastrigins benchmark function in nD. Global minimum: f(x)=0; xi=0.
 *
 * @author Johannes Dieterich
 * @version 2021-07-20
 */
public final class BenchRastrigin extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20210720;
  public static final boolean DEBUG = false;

  private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

  private final int iDims;

  public BenchRastrigin(final int dims) {
    this.iDims = dims;
  }

  @Override
  public BenchRastrigin copy() {
    return new BenchRastrigin(iDims);
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {

    final double[] daParams = params.getAllParamters();

    final int upperBound = SPECIES.loopBound(iDims);
    final double twoPI = 2 * Math.PI;
    final var vTwoPI = DoubleVector.broadcast(SPECIES, twoPI);
    final var vCutVal = DoubleVector.broadcast(SPECIES, 5.12);
    final var vTen = DoubleVector.broadcast(SPECIES, 10);

    var vRes = DoubleVector.broadcast(SPECIES, 0.0);
    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vPI = DoubleVector.fromArray(SPECIES, daParams, i);
      // first check if we need the cutoff
      final var vPIAbs = vPI.abs();
      final var vCutMask = vPIAbs.compare(VectorOperators.GT, vCutVal);
      final var anyCut = vCutMask.anyTrue();
      if (anyCut) {
        final var vPSq = vPI.mul(vPI);

        // cutoff part
        final var vCut = vPSq.mul(vTen);

        // regular parts
        final var vPC = vPI.mul(vTwoPI).lanewise(VectorOperators.COS).mul(vTen);
        final var vResLocUn = vPSq.sub(vPC);

        // blend
        final var vResLoc = vResLocUn.blend(vCut, vCutMask);

        vRes = vRes.add(vResLoc);
      } else {
        final var vPSq = vPI.mul(vPI);
        final var vPC = vPI.mul(vTwoPI).lanewise(VectorOperators.COS).mul(vTen);
        final var vResLoc = vPSq.sub(vPC);
        vRes = vRes.add(vResLoc);
      }
    }

    double dEnergy = 10 * iDims + vRes.reduceLanes(VectorOperators.ADD);
    for (; i < iDims; i++) {
      if (daParams[i] > 5.12 || daParams[i] < -5.12) {
        // take the cutoff potential
        dEnergy += 10.0 * daParams[i] * daParams[i];
      } else {
        dEnergy += daParams[i] * daParams[i] - 10.0 * Math.cos(twoPI * daParams[i]);
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

    final double[] daParams = params.getAllParamters();
    final double twoPI = 2 * Math.PI;
    final double twentyPI = 20 * Math.PI;

    final int upperBound = SPECIES.loopBound(iDims);
    final var vTwoPI = DoubleVector.broadcast(SPECIES, twoPI);
    final var vTwentyPI = DoubleVector.broadcast(SPECIES, twentyPI);
    final var vCutVal = DoubleVector.broadcast(SPECIES, 5.12);
    final var vTwo = DoubleVector.broadcast(SPECIES, 2);
    final var vTen = DoubleVector.broadcast(SPECIES, 10);

    var vRes = DoubleVector.broadcast(SPECIES, 0.0);
    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vPI = DoubleVector.fromArray(SPECIES, daParams, i);
      // first check if we need the cutoff
      final var vPIAbs = vPI.abs();
      final var vCutMask = vPIAbs.compare(VectorOperators.GT, vCutVal);
      final var anyCut = vCutMask.anyTrue();
      if (anyCut) {
        final var vPSq = vPI.mul(vPI);

        // cutoff part
        final var vCut = vPSq.mul(vTen);
        final var vCutGrad = vPI.mul(vTen);

        // regular parts
        final var vT1 = vPI.mul(vTwoPI);
        final var vPC = vT1.lanewise(VectorOperators.COS).mul(vTen);
        final var vResLocUn = vPSq.sub(vPC);

        final var vT2 = vT1.lanewise(VectorOperators.SIN).mul(vTwentyPI);
        final var vGradUn = vPI.mul(vTwo).add(vT2);

        // blend
        final var vResLoc = vResLocUn.blend(vCut, vCutMask);
        vRes = vRes.add(vResLoc);
        final var vGrad = vGradUn.blend(vCutGrad, vCutMask);
        vGrad.intoArray(daGrad, i);
      } else {
        final var vT1 = vPI.mul(vTwoPI);
        final var vPSq = vPI.mul(vPI);
        final var vPC = vT1.lanewise(VectorOperators.COS).mul(vTen);
        final var vResLoc = vPSq.sub(vPC);
        vRes = vRes.add(vResLoc);
        final var vT2 = vT1.lanewise(VectorOperators.SIN).mul(vTwentyPI);
        final var vGrad = vPI.mul(vTwo).add(vT2);
        vGrad.intoArray(daGrad, i);
      }
    }

    double dEnergy = 10 * iDims + vRes.reduceLanes(VectorOperators.ADD);
    for (; i < iDims; i++) {
      if (daParams[i] > 5.12 || daParams[i] < -5.12) {
        // take the deviation of the cutoff potential
        daGrad[i] = 10.0 * daParams[i];
        dEnergy += 10.0 * daParams[i] * daParams[i];
      } else {
        /*
         * f'(x(i)) = 2x_i + 20pi*sin(2pi x_i)
         */
        final double t1 = twoPI * daParams[i];
        daGrad[i] = 2.0 * daParams[i] + twentyPI * Math.sin(t1);
        dEnergy += daParams[i] * daParams[i] - 10.0 * Math.cos(t1);
      }
    }

    if (DEBUG) {
      // calculate a numerical gradient and compare it to the analytical one
      final double[] daNumGrad = new double[iDims];
      final double numE =
          NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, daNumGrad);

      for (int x = 0; x < daNumGrad.length; x++) {
        if (Math.abs(daNumGrad[i] - daGrad[x]) >= 1E-7) {
          // we have a noticable difference
          System.err.println(
              "DEBUG: Analytical vs numerical gradient is: " + daGrad[x] + " vs " + daNumGrad[x]);
        } else {
          System.out.println("DEBUG: Analytical vs numerical gradient was fine. :-)");
        }
      }
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
      daMinMax[0][i] = -5.12;
      daMinMax[1][i] = 5.12;
    }

    return daMinMax;
  }
}
