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

import static java.lang.Math.*;

import java.util.ArrayList;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;

/**
 * Ackleys benchmark function in nD. Global minimum: f(x)=0, x(i)=0.0.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
final class BenchAckley extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20150727;
  private static final boolean DEBUG = false;
  private static final double TWOPI = PI * 2;
  private final boolean correctGrad;
  private final int iDims;
  private final double dimsInv;

  BenchAckley(final int dims, final boolean corrGrad) {
    this.iDims = dims;
    this.dimsInv = 1.0 / dims;
    this.correctGrad = corrGrad;
  }

  @Override
  public BenchAckley copy() {
    return new BenchAckley(iDims, correctGrad);
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {

    final double[] p = params.getAllParamters();
    final int no = p.length;

    double a = 0.0;
    double b = 0.0;
    for (int i = 0; i < no; i++) {
      final double pi = p[i];
      if (pi > 32.768 || pi < -32.768) {
        a += FixedValues.NONCONVERGEDENERGY;
        b += FixedValues.NONCONVERGEDENERGY;
        continue;
      }
      a += pi * pi;
      b += cos(TWOPI * pi);
    }

    final double e = -20 * exp(-0.2 * sqrt(a * dimsInv)) - exp(b * dimsInv) + 20 + E;

    return e;
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds,
      final double[] daGrad) {

    final double[] p = params.getAllParamters();
    assert (p.length == params.getNumberOfParamters());
    final int no = params.getNumberOfParamters();

    double a = 0.0;
    double b = 0.0;
    for (int i = 0; i < no; i++) {
      /*
       * done using wolframalpha :-)
       */
      final double pi = p[i];

      if (pi > 32.768) {
        a += FixedValues.NONCONVERGEDENERGY;
        b += FixedValues.NONCONVERGEDENERGY;
        // daGrad[i] -= FixedValues.NONCONVERGEDGRADIENT;
        continue;
      } else if (pi < -32.768) {
        a += FixedValues.NONCONVERGEDENERGY;
        b += FixedValues.NONCONVERGEDENERGY;
        // daGrad[i] += FixedValues.NONCONVERGEDGRADIENT;
        continue;
      }

      final double pSq = pi * pi;
      a += pSq;
      final double trigi = TWOPI * pi;
      final double bi = cos(trigi);
      b += bi;

      if (correctGrad) {
        daGrad[i] = TWOPI * sin(trigi); // we use this as a cache
      } else {
        final double pSqDims = sqrt(pSq * dimsInv);

        final double firstNum = 4 * exp(-0.2 * pSqDims) * pi;
        final double firstDenom = pSqDims * iDims;
        final double firstTerm = firstNum / firstDenom;
        final double secNum = TWOPI * sin(trigi) * exp(dimsInv * bi);
        final double secTerm = secNum * dimsInv;

        daGrad[i] = firstTerm + secTerm;
      }
    }

    final double t1 = exp(b * dimsInv);
    final double t2 = sqrt(a * dimsInv);
    final double t2Inv = 1 / (iDims * t2);
    final double t3 = exp(-0.2 * t2);

    if (correctGrad) {
      for (int i = 0; i < no; i++) {
        final double secNum = daGrad[i] * t1;
        final double secTerm = secNum * dimsInv;
        daGrad[i] = 4 * p[i] * t3 * t2Inv + secTerm;
      }
    }

    final double e = -20 * t3 - t1 + 20 + E;

    if (DEBUG) {
      final double[] numGrad = new double[daGrad.length];

      final double numFunc =
          NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, numGrad);
      System.out.println("DEBUG OUTPUT FOR NUMGRAD");
      System.out.println("Func " + numFunc + "\t" + e + "\t" + (numFunc - e));
      for (int i = 0; i < daGrad.length; i++) {
        System.out.println(numGrad[i] + "\t" + daGrad[i] + "\t" + (numGrad[i] - daGrad[i]));
      }
    }

    return e;
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
      daMinMax[0][i] = -32.768;
      daMinMax[1][i] = 32.768;
    }

    return daMinMax;
  }
}
