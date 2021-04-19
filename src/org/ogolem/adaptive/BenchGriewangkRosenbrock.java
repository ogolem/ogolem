/*
Copyright (c) 2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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

import java.util.ArrayList;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;

/**
 * Composite Griewangk-Rosenbrock function.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29 TODO somewhat this function seems to be potentially broken (not providing the
 *     correct minima after a first test, it seems)
 */
final class BenchGriewangkRosenbrock extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20150727;
  public static final boolean DEBUG = false;

  private final int iDims;
  private final double mult;

  BenchGriewangkRosenbrock(final int dims) {
    this.iDims = dims;
    this.mult = 10.0 / (dims - 1);
  }

  @Override
  public BenchGriewangkRosenbrock copy() {
    return new BenchGriewangkRosenbrock(iDims);
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {

    final double[] daParams = params.getAllParamters();
    final int iNoOfParams = daParams.length;
    assert (iNoOfParams == iDims);

    double dFuncVal = 0.0;
    for (int i = 0; i < iNoOfParams - 1; i++) {

      final double pi = daParams[i];
      final double piSq = pi * pi;
      final double t1 = piSq - daParams[i + 1];
      final double t2 = pi - 1;
      final double t3 = t2 * t2;

      final double dSi = 100.0 * t1 * t1 + t3;
      dFuncVal += dSi * 0.00025 - Math.cos(dSi);
    }

    dFuncVal *= mult;
    dFuncVal += 10.0;

    return dFuncVal;
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds,
      final double[] grad) {

    final double[] daParams = params.getAllParamters();
    final int iNoOfParams = daParams.length;
    assert (iNoOfParams == iDims);

    for (int i = 0; i < grad.length; i++) grad[i] = 0.0;

    double dFuncVal = 0.0;
    for (int i = 0; i < iNoOfParams - 1; i++) {
      final double pi = daParams[i];
      final double pj = daParams[i + 1];
      final double piSq = pi * pi;
      final double t1 = piSq - pj;
      final double t1Sq = t1 * t1;
      final double t2 = pi - 1;
      final double t3 = t2 * t2;

      final double dSi = 100.0 * t1Sq + t3;

      final double sinSi = Math.sin(dSi);
      final double sinSiT = 4000 * sinSi + 1;

      grad[i] += mult * 0.0005 * (200 * pi * (piSq - pj) + pi - 1) * sinSiT;
      // grad[i  ] += mult*(0.00025*(400*pi*t1+2*t2)+(400*pi*t1+2*t2)*sinSi);
      grad[i + 1] -= mult * 0.05 * t1 * sinSiT;
      // grad[i+1] += mult*(0.05*(pj-piSq)-200*t1*sinSi);
      dFuncVal += dSi * 0.00025 - Math.cos(dSi);
    }

    dFuncVal *= mult;
    dFuncVal += 10.0;

    if (DEBUG) {
      final double[] numGrad = new double[grad.length];

      final double numFunc =
          NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, numGrad);
      System.out.println("DEBUG OUTPUT FOR NUMGRAD");
      System.out.println("Func " + numFunc + "\t" + dFuncVal + "\t" + (numFunc - dFuncVal));
      for (int i = 0; i < grad.length; i++) {
        System.out.println(numGrad[i] + "\t" + grad[i] + "\t" + (numGrad[i] - grad[i]));
      }
    }

    return dFuncVal;
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
      daMinMax[0][i] = -5.00;
      daMinMax[1][i] = 5.00;
    }

    return daMinMax;
  }
}
