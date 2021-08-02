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
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;

/**
 * Ackleys benchmark function in nD. Global minimum: f(x)=0, x(i)=0.0.
 *
 * @author Johannes Dieterich
 * @version 2021-07-21
 */
public final class BenchAckley extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20150727;
  private static final boolean DEBUG = false;
  private static final double TWOPI = PI * 2;
  private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

  private final int iDims;
  private final double dimsInv;

  public BenchAckley(final int dims) {
    this.iDims = dims;
    this.dimsInv = 1.0 / dims;
  }

  @Override
  public BenchAckley copy() {
    return new BenchAckley(iDims);
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {

    final double[] p = params.getAllParamters();

    final int upperBound = SPECIES.loopBound(iDims);
    final var vCutVal = DoubleVector.broadcast(SPECIES, 32.768);
    final var vCutCons = DoubleVector.broadcast(SPECIES, FixedValues.NONCONVERGEDENERGY);
    final var vTwoPI = DoubleVector.broadcast(SPECIES, 2 * Math.PI);

    var vA = DoubleVector.broadcast(SPECIES, 0.0);
    var vB = DoubleVector.broadcast(SPECIES, 0.0);
    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vPI = DoubleVector.fromArray(SPECIES, p, i);
      // first check if we need the cutoff
      final var vPIAbs = vPI.abs();
      final var vCutMask = vPIAbs.compare(VectorOperators.GT, vCutVal);
      final var anyCut = vCutMask.anyTrue();
      if (anyCut) {
        final var vPISqUn = vPI.mul(vPI);
        final var vPISq = vPISqUn.blend(vCutCons, vCutMask);
        vA = vPISq.add(vA);

        final var vBTemUn = vPI.mul(vTwoPI).lanewise(VectorOperators.COS);
        final var vBTem = vBTemUn.blend(vCutCons, vCutMask);
        vB = vBTem.add(vB);
      } else {
        vA = vPI.fma(vPI, vA);
        final var vBTem = vPI.mul(vTwoPI).lanewise(VectorOperators.COS);
        vB = vBTem.add(vB);
      }
    }

    double a = vA.reduceLanes(VectorOperators.ADD);
    double b = vB.reduceLanes(VectorOperators.ADD);
    for (; i < iDims; i++) {
      final double pi = p[i];
      if (pi > 32.768 || pi < -32.768) {
        a += FixedValues.NONCONVERGEDENERGY;
        b += FixedValues.NONCONVERGEDENERGY;
        continue;
      }
      a = Math.fma(pi, pi, a);
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

    final int upperBound = SPECIES.loopBound(iDims);
    final var vCutVal = DoubleVector.broadcast(SPECIES, 32.768);
    final var vCutCons = DoubleVector.broadcast(SPECIES, FixedValues.NONCONVERGEDENERGY);
    final var vTwoPI = DoubleVector.broadcast(SPECIES, 2 * Math.PI);

    var vA = DoubleVector.broadcast(SPECIES, 0.0);
    var vB = DoubleVector.broadcast(SPECIES, 0.0);
    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vPI = DoubleVector.fromArray(SPECIES, p, i);
      // first check if we need the cutoff
      final var vPIAbs = vPI.abs();
      final var vCutMask = vPIAbs.compare(VectorOperators.GT, vCutVal);
      final var anyCut = vCutMask.anyTrue();
      if (anyCut) {
        final var vPISqUn = vPI.mul(vPI);
        final var vPISq = vPISqUn.blend(vCutCons, vCutMask);
        vA = vPISq.add(vA);

        final var vTrigi = vPI.mul(vTwoPI);
        final var vBTemUn = vTrigi.lanewise(VectorOperators.COS);
        final var vBTem = vBTemUn.blend(vCutCons, vCutMask);
        vB = vBTem.add(vB);

        final var vGradTmp = vTrigi.lanewise(VectorOperators.SIN).mul(vTwoPI);
        vGradTmp.intoArray(daGrad, i);

      } else {
        vA = vPI.fma(vPI, vA);
        final var vTrigi = vPI.mul(vTwoPI);
        final var vBTem = vTrigi.lanewise(VectorOperators.COS);
        vB = vBTem.add(vB);
        final var vGradTmp = vTrigi.lanewise(VectorOperators.SIN).mul(vTwoPI);
        vGradTmp.intoArray(daGrad, i);
      }
    }

    double a = vA.reduceLanes(VectorOperators.ADD);
    double b = vB.reduceLanes(VectorOperators.ADD);
    for (; i < iDims; i++) {
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

      a = Math.fma(pi, pi, a);
      final double trigi = TWOPI * pi;
      final double bi = cos(trigi);
      b += bi;

      daGrad[i] = TWOPI * sin(trigi); // we use this as a cache
    }

    final double t1 = exp(b * dimsInv);
    final double t2 = sqrt(a * dimsInv);
    final double t2Inv = (Math.abs(t2) < 1e-10) ? 0.0 : 1 / (iDims * t2);
    final double t3 = exp(-0.2 * t2);

    final var vT1 = DoubleVector.broadcast(SPECIES, t1);
    final var vT2Inv = DoubleVector.broadcast(SPECIES, t2Inv);
    final var vT3 = DoubleVector.broadcast(SPECIES, t3);
    final var vFour = DoubleVector.broadcast(SPECIES, 4);
    final var vDimsInv = DoubleVector.broadcast(SPECIES, dimsInv);

    int x = 0;
    for (; x < upperBound; x += SPECIES.length()) {
      final var vGradLoad = DoubleVector.fromArray(SPECIES, daGrad, x);
      final var vSecNum = vGradLoad.mul(vT1);
      final var vSecTerm = vSecNum.mul(vDimsInv);
      final var vP = DoubleVector.fromArray(SPECIES, p, x);
      final var vGrad = vP.mul(vFour).mul(vT3).mul(vT2Inv).add(vSecTerm);

      vGrad.intoArray(daGrad, x);
    }

    for (; x < iDims; x++) {
      final double secNum = daGrad[x] * t1;
      final double secTerm = secNum * dimsInv;

      daGrad[x] = 4 * p[x] * t3 * t2Inv + secTerm;
    }

    final double e = -20 * t3 - t1 + 20 + E;

    if (DEBUG) {
      final double[] numGrad = new double[daGrad.length];

      final double numFunc =
          NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, numGrad);
      System.out.println("DEBUG OUTPUT FOR NUMGRAD");
      System.out.println("Func " + numFunc + "\t" + e + "\t" + (numFunc - e));
      for (int y = 0; y < daGrad.length; y++) {
        System.out.println(numGrad[y] + "\t" + daGrad[y] + "\t" + (numGrad[y] - daGrad[y]));
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
