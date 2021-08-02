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
 * Lunacek's bi-Rastrigin function. Global minimum at xi = 2.5. Fake minimum at xi=-2.5.
 *
 * @author Johannes Dieterich
 * @version 2021-07-22
 */
public final class BenchLunacek extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20210722;
  public static final boolean bDebug = false;
  private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

  private final int iDims;

  // some things that are costly to evaluate
  private static final double dD = 1.0;
  private static final double dMu0 = 2.5;
  private final double dS;
  private final double dMu1;

  public BenchLunacek(final int dims) {
    this.iDims = dims;
    // commented, this is the bbob setting (be aware that this *might* change the gradient
    // expression)
    // this.dS = 1.0 - 1.0/(2*Math.sqrt((iDims+20))-8.2);
    this.dS = 0.7;
    this.dMu1 = -Math.sqrt((dMu0 * dMu0 - dD) / dS);
  }

  @Override
  public BenchLunacek copy() {
    return new BenchLunacek(iDims);
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {

    final double[] daParams = params.getAllParamters();

    final int upperBound = SPECIES.loopBound(iDims);
    final var vCutVal = DoubleVector.broadcast(SPECIES, 5.0);
    final var vMu0 = DoubleVector.broadcast(SPECIES, dMu0);
    final var vMu1 = DoubleVector.broadcast(SPECIES, dMu1);
    final var vTwoPI = DoubleVector.broadcast(SPECIES, 2.0 * Math.PI);

    var vSphere1 = DoubleVector.broadcast(SPECIES, 0.0);
    var vSphere2 = DoubleVector.broadcast(SPECIES, 0.0);
    var vRastrigin = DoubleVector.broadcast(SPECIES, 0.0);

    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vPI = DoubleVector.fromArray(SPECIES, daParams, i);
      // first check if we need the cutoff
      final var vPIAbs = vPI.abs();
      final var vCutMask = vPIAbs.compare(VectorOperators.GT, vCutVal);
      final var anyCut = vCutMask.anyTrue();
      if (anyCut) {
        // cutoff potential
        final var vCutPot = vPI.mul(vPI).neg();

        // regular & blend
        final var vMu0P = vPI.sub(vMu0);
        final var vMu1P = vPI.sub(vMu1);

        vSphere1 = vMu0P.fma(vMu0P, vSphere1);
        vSphere2 = vMu1P.fma(vMu1P, vSphere2);

        final var vTRas = vMu0P.mul(vTwoPI).lanewise(VectorOperators.COS).blend(vCutPot, vCutMask);
        vRastrigin = vTRas.add(vRastrigin);
      } else {
        final var vMu0P = vPI.sub(vMu0);
        final var vMu1P = vPI.sub(vMu1);

        vSphere1 = vMu0P.fma(vMu0P, vSphere1);
        vSphere2 = vMu1P.fma(vMu1P, vSphere2);

        final var vTRas = vMu0P.mul(vTwoPI).lanewise(VectorOperators.COS);
        vRastrigin = vTRas.add(vRastrigin);
      }
    }

    double dSpherePart1 = vSphere1.reduceLanes(VectorOperators.ADD);
    double dSpherePart2 = vSphere2.reduceLanes(VectorOperators.ADD);
    double dRastriginPart = vRastrigin.reduceLanes(VectorOperators.ADD);
    for (; i < iDims; i++) {

      // boundaries
      if (daParams[i] > 5.0 || daParams[i] < -5.0) {
        // needs to be negative so that it gets positive in the end...
        dRastriginPart -= daParams[i] * daParams[i];
        continue;
      }

      final double mu0P = daParams[i] - dMu0;
      final double mu1P = daParams[i] - dMu1;
      dSpherePart1 = Math.fma(mu0P, mu0P, dSpherePart1);
      dSpherePart2 = Math.fma(mu1P, mu1P, dSpherePart2);

      dRastriginPart += Math.cos(2.0 * Math.PI * (daParams[i] - dMu0));
    }

    final double dSphereContrib = Math.min(dSpherePart1, (dD * iDims + dS * dSpherePart2));

    final double dTotal = dSphereContrib + 10.0 * iDims - 10.0 * dRastriginPart;

    return dTotal;
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds,
      final double[] grad) {

    final double[] daParams = params.getAllParamters();

    final int upperBound = SPECIES.loopBound(iDims);
    final var vCutVal = DoubleVector.broadcast(SPECIES, 5.0);
    final var vMu0 = DoubleVector.broadcast(SPECIES, dMu0);
    final var vMu1 = DoubleVector.broadcast(SPECIES, dMu1);
    final var vTwoPI = DoubleVector.broadcast(SPECIES, 2.0 * Math.PI);
    final var vTwo = DoubleVector.broadcast(SPECIES, 2.0);
    final var vTwentyPI = DoubleVector.broadcast(SPECIES, 20.0 * Math.PI);
    final var vCutGrad = DoubleVector.broadcast(SPECIES, FixedValues.NONCONVERGEDGRADIENT);

    var vSphere1 = DoubleVector.broadcast(SPECIES, 0.0);
    var vSphere2 = DoubleVector.broadcast(SPECIES, 0.0);
    var vRastrigin = DoubleVector.broadcast(SPECIES, 0.0);

    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vPI = DoubleVector.fromArray(SPECIES, daParams, i);
      // first check if we need the cutoff
      final var vPIAbs = vPI.abs();
      final var vCutMask = vPIAbs.compare(VectorOperators.GT, vCutVal);
      final var anyCut = vCutMask.anyTrue();
      if (anyCut) {
        // cutoff potential
        final var vCutPot = vPI.mul(vPI).neg();

        // regular & blend
        final var vMu0P = vPI.sub(vMu0);
        final var vMu1P = vPI.sub(vMu1);

        vSphere1 = vMu0P.fma(vMu0P, vSphere1);
        vSphere2 = vMu1P.fma(vMu1P, vSphere2);

        final var vT1 = vMu0P.mul(vTwoPI);
        final var vTCos = vT1.lanewise(VectorOperators.COS).blend(vCutPot, vCutMask);
        vRastrigin = vTCos.add(vRastrigin);
        final var vGrad =
            vT1.lanewise(VectorOperators.SIN).mul(vTwentyPI).blend(vCutGrad, vCutMask);
        vGrad.intoArray(grad, i);
      } else {
        final var vMu0P = vPI.sub(vMu0);
        final var vMu1P = vPI.sub(vMu1);

        vSphere1 = vMu0P.fma(vMu0P, vSphere1);
        vSphere2 = vMu1P.fma(vMu1P, vSphere2);

        final var vT1 = vMu0P.mul(vTwoPI);
        final var vTCos = vT1.lanewise(VectorOperators.COS);
        vRastrigin = vTCos.add(vRastrigin);
        final var vGrad = vT1.lanewise(VectorOperators.SIN).mul(vTwentyPI);
        vGrad.intoArray(grad, i);
      }
    }

    double dSpherePart1 = vSphere1.reduceLanes(VectorOperators.ADD);
    double dSpherePart2 = vSphere2.reduceLanes(VectorOperators.ADD);
    double dRastriginPart = vRastrigin.reduceLanes(VectorOperators.ADD);
    for (; i < iDims; i++) {

      // boundaries
      if (daParams[i] > 5.0 || daParams[i] < -5.0) {
        // needs to be negative so that it gets positive in the end...
        dRastriginPart -= daParams[i] * daParams[i];
        grad[i] += FixedValues.NONCONVERGEDGRADIENT;
        continue;
      }

      final double mu0P = daParams[i] - dMu0;
      final double mu1P = daParams[i] - dMu1;
      dSpherePart1 = Math.fma(mu0P, mu0P, dSpherePart1);
      dSpherePart2 = Math.fma(mu1P, mu1P, dSpherePart2);

      final double t1 = 2 * Math.PI * (daParams[i] - dMu0);
      dRastriginPart += Math.cos(t1);
      grad[i] = +20 * Math.PI * Math.sin(t1); // *2.0*10.0 because of the *-10.0 in dTotal
    }

    // grad[] contains the Rastrigin gradients now, add the sphere part
    double dSphereContrib;
    final double t1 = dD * iDims + dS * dSpherePart2;
    if (dSpherePart1 < t1) {
      dSphereContrib = dSpherePart1;

      int x = 0;
      for (; x < upperBound; x += SPECIES.length()) {
        final var vPI = DoubleVector.fromArray(SPECIES, daParams, x);
        final var vGradLoad = DoubleVector.fromArray(SPECIES, grad, x);
        final var vTmp = vPI.sub(vMu0);
        final var vGrad = vTmp.fma(vTwo, vGradLoad);
        vGrad.intoArray(grad, x);
      }

      for (; x < iDims; x++) grad[x] += 2 * (daParams[x] - dMu0);
    } else {
      dSphereContrib = t1;

      final var vS = DoubleVector.broadcast(SPECIES, dS);

      int x = 0;
      for (; x < upperBound; x += SPECIES.length()) {
        final var vPI = DoubleVector.fromArray(SPECIES, daParams, x);
        final var vGradLoad = DoubleVector.fromArray(SPECIES, grad, x);
        final var vTmp = vPI.sub(vMu1).mul(vS);
        final var vGrad = vTmp.fma(vTwo, vGradLoad);
        vGrad.intoArray(grad, x);
      }

      for (; x < iDims; x++) grad[x] += 2 * dS * (daParams[x] - dMu1);
    }

    final double dTotal = dSphereContrib + 10.0 * iDims - 10.0 * dRastriginPart;

    return dTotal;
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
