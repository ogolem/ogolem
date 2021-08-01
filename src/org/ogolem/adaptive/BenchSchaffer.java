/*
Copyright (c) 2010, J. M. Dieterich and B. Hartke
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
 * Schaffers F7 benchmark function in N-D. Global minimum: f(x,y)=0: x=0,y=0.
 *
 * @author Johannes Dieterich
 * @version 2021-07-18
 */
public final class BenchSchaffer extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20210718;
  private static final boolean DEBUG = false;

  private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

  private final int iDims;

  public BenchSchaffer(final int dims) {
    assert (dims > 1);
    this.iDims = dims;
  }

  @Override
  public BenchSchaffer copy() {
    return new BenchSchaffer(iDims);
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {

    final double[] daParams = params.getAllParamters();

    final int upperBound = SPECIES.loopBound(iDims - 1);
    final var oneFifth = DoubleVector.broadcast(SPECIES, 1.0 / 5.0);
    final var one = DoubleVector.broadcast(SPECIES, 1.0);
    final var fifty = DoubleVector.broadcast(SPECIES, 50.0);

    var vRes = DoubleVector.broadcast(SPECIES, 0.0);
    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vPI = DoubleVector.fromArray(SPECIES, daParams, i);
      final var vPN = DoubleVector.fromArray(SPECIES, daParams, i + 1);
      final var vPISq = vPI.mul(vPI);
      final var vPNSq = vPN.mul(vPN);
      final var vSI = vPISq.add(vPNSq).sqrt();
      final var vSIR = vSI.sqrt();
      final var vT1Arg = vSI.pow(oneFifth).mul(fifty);
      final var vT1 = vT1Arg.lanewise(VectorOperators.SIN);
      final var vT11 = vT1.mul(vT1).add(one);
      final var vSIR_T11 = vSIR.mul(vT11);
      vRes = vRes.add(vSIR_T11);
    }

    double dEnergy = vRes.reduceLanes(VectorOperators.ADD);
    for (; i < iDims - 1; i++) {
      final double pi = daParams[i];
      final double piSq = pi * pi;
      final double piN = daParams[i + 1];
      final double piNSq = piN * piN;
      final double si = Math.sqrt(piSq + piNSq);
      final double siR = Math.sqrt(si);
      final double t1 = Math.sin(50 * Math.pow(si, 1.0 / 5.0));
      dEnergy += siR * (t1 * t1 + 1);
    }

    dEnergy /= (iDims - 1.0);
    dEnergy *= dEnergy;

    return dEnergy;
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds,
      final double[] grad) {

    final double[] daParams = params.getAllParamters();

    final double norm = 1 / (iDims - 1.0);

    final var vOne = DoubleVector.broadcast(SPECIES, 1.0);
    final var oneFifth = DoubleVector.broadcast(SPECIES, 1.0 / 5.0);
    final var v13_20th = DoubleVector.broadcast(SPECIES, 13.0 / 20.0);
    final var one = DoubleVector.broadcast(SPECIES, 1.0);
    final var fifty = DoubleVector.broadcast(SPECIES, 50.0);
    final var two = DoubleVector.broadcast(SPECIES, 2.0);

    final int upperBound = SPECIES.loopBound(iDims - 1);
    var vRes = DoubleVector.broadcast(SPECIES, 0.0);
    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vPI = DoubleVector.fromArray(SPECIES, daParams, i);
      final var vPN = DoubleVector.fromArray(SPECIES, daParams, i + 1);
      final var vPISq = vPI.mul(vPI);
      final var vPNSq = vPN.mul(vPN);
      final var vPIS = vPISq.add(vPNSq);
      final var vSI = vPIS.sqrt();
      final var vSIR = vSI.sqrt();
      final var vT1 = vSI.pow(oneFifth).mul(fifty);
      final var vT2s = vT1.lanewise(VectorOperators.SIN);
      final var vT2sSq = vT2s.mul(vT2s);

      final var vT2sSq1 = vT2sSq.add(one);
      final var vT3 = vSIR.mul(vT2sSq1);
      vRes = vRes.add(vT3);

      final var vSI34 = vSIR.mul(vSIR).mul(vSIR);
      final var vSecDivPrior = vOne.div(vSI34.mul(two));
      final var vSecDivCorr =
          vSI34.abs().lt(1e-15).and(vPI.abs().lt(1e-10)).and(vPN.abs().lt(1e-10));
      final var vSecDiv = vSecDivPrior.blend(0.0, vSecDivCorr);
      final var vSecX = vPI.mul(vSecDiv);
      final var vSecY = vPN.mul(vSecDiv);

      final var vThird = vT2sSq.mul(vSecDiv);
      final var vThirdX = vPI.mul(vThird);
      final var vThirdY = vPN.mul(vThird);

      final var vSI1320 = vPIS.pow(v13_20th);
      final var vT2C = vT1.lanewise(VectorOperators.COS);
      final var vFourthPrior = vT2s.mul(20).mul(vT2C).div(vSI1320);
      // final var vFourthCorr = vPIS.abs().lt(1e-20);
      // we "should" be blending w/ vFourthCorr, however, the reason
      // these are close to zero are the parameters - hence, save
      // that eval
      final var vFourth = vFourthPrior.blend(0.0, vSecDivCorr);
      // vFourthCorr == vSecDivCorr?
      final var vFourthX = vPI.mul(vFourth);
      final var vFourthY = vPN.mul(vFourth);

      final var vGradI = vSecX.add(vThirdX).add(vFourthX);
      final var vGradN = vSecY.add(vThirdY).add(vFourthY); // XXX why not just blend here?

      // putting this into the array is a bit more difficult due to overlap
      // namely: assuming 4-wide vectors, we'd look at
      // grad [i]   [i+1]   [i+2]   [i+3]   [i+4]
      //   +gI[0]  +gI[1]  +gI[2]  +gI[3]
      //           +gN[0]  +gN[1]  +gN[2]  +gN[3]

      final var vGradNShift = vGradN.broadcast(0).slice(vGradN.length() - 1, vGradN);

      final var vGradIP = DoubleVector.fromArray(SPECIES, grad, i);
      final var vGradIRes = vGradI.add(vGradIP).add(vGradNShift);
      vGradIRes.intoArray(grad, i);

      // if we didn't do the little shift trick above, this would be the correct sequence
      // final var vGradNP = DoubleVector.fromArray(SPECIES, grad, i + 1);
      // final var vGradNRes = vGradN.add(vGradNP);
      // vGradNRes.intoArray(grad, i + 1); // XXX MAYBE WE SHOULD BLEND HERE
      grad[i + SPECIES.length()] += vGradN.lane(vGradN.length() - 1);
    }

    double dEnergy = vRes.reduceLanes(VectorOperators.ADD);
    for (; i < iDims - 1; i++) {

      final double pi = daParams[i];
      final double piSq = pi * pi;
      final double piN = daParams[i + 1];
      final double piNSq = piN * piN;
      final double si = Math.sqrt(piSq + piNSq);
      final double siR = Math.sqrt(si);
      final double t1 = 50 * Math.pow(si, 1.0 / 5.0);
      final double t2s = Math.sin(t1);
      final double t2sSq = t2s * t2s;

      final double t3 = siR * (t2sSq + 1);

      final double si34 = siR * siR * siR;
      final double secDiv =
          (Math.abs(si34) < 1e-20 && Math.abs(pi) < 1e-10 && Math.abs(piN) < 1e-10)
              ? 0.0
              : 1 / (2 * si34);
      final double secX = pi * secDiv;
      final double secY = piN * secDiv;

      final double third = t2sSq * secDiv;
      final double thirdX = pi * third;
      final double thirdY = piN * third;

      final double pis = piSq + piNSq;
      final double si1320 = Math.pow((pis), 13.0 / 20.0);
      final double t2c = Math.cos(t1);
      // the last parameter may give troubles if we're not checking it individually
      // effectively a very low pi may compensate for a slightly higher piN
      final double fourth =
          (Math.abs(pis) < 1e-20 && Math.abs(pi) < 1e-10 && Math.abs(piN) < 1e-10)
              ? 0.0
              : 20 * t2s * t2c / si1320;
      final double fourthX = pi * fourth;
      final double fourthY = piN * fourth;

      grad[i] += (secX + thirdX + fourthX);
      grad[i + 1] += (secY + thirdY + fourthY);
      dEnergy += t3;
    }

    dEnergy *= norm;
    final double norm2 = dEnergy * 2 * norm;
    final var vNorm2 = DoubleVector.broadcast(SPECIES, norm2);
    final int upperBound2 = SPECIES.loopBound(iDims);
    int x = 0;
    for (; x < upperBound2; x += SPECIES.length()) {
      final var vGradLoad = DoubleVector.fromArray(SPECIES, grad, x);
      final var vGrad = vGradLoad.mul(vNorm2);
      vGrad.intoArray(grad, x);
    }

    for (; x < iDims; x++) grad[x] *= norm2;

    dEnergy *= dEnergy;

    if (DEBUG) {
      // calculate a numerical gradient and compare it to the analytical one
      final double[] daNumGrad = new double[iDims];
      final double numE =
          NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, daNumGrad);

      for (int y = 0; y < daNumGrad.length; y++) {
        if (Math.abs(daNumGrad[y] - grad[y]) >= 1E-7) {
          // we have a noticable difference
          System.err.println(
              "DEBUG: Analytical vs numerical gradient is: " + grad[y] + " vs " + daNumGrad[y]);
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
      daMinMax[0][i] = -100.0;
      daMinMax[1][i] = 100.0;
    }

    return daMinMax;
  }
}
