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
 * Lunacek's bi-Rastrigin function. Global minimum at xi = 2.5. Fake minimum at xi=-2.5.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
final class BenchLunacek extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20150727;
  public static final boolean bDebug = false;

  private final int iDims;

  // some things that are costly to evaluate
  private static final double dD = 1.0;
  private static final double dMu0 = 2.5;
  private final double dS;
  private final double dMu1;

  BenchLunacek(final int dims) {
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
    final int iNoOfParams = daParams.length;

    double dSpherePart1 = 0.0;
    double dSpherePart2 = 0.0;
    double dRastriginPart = 0.0;
    for (int i = 0; i < iNoOfParams; i++) {

      // boundaries
      if (daParams[i] > 5.0 || daParams[i] < -5.0) {
        // needs to be negative so that it gets positive in the end...
        dRastriginPart -= daParams[i] * daParams[i];
        continue;
      }

      final double mu0P = daParams[i] - dMu0;
      final double mu1P = daParams[i] - dMu1;
      dSpherePart1 += mu0P * mu0P;
      dSpherePart2 += mu1P * mu1P;

      dRastriginPart += Math.cos(2.0 * Math.PI * (daParams[i] - dMu0));
    }

    final double dSphereContrib = Math.min(dSpherePart1, (dD * iNoOfParams + dS * dSpherePart2));

    final double dTotal = dSphereContrib + 10.0 * iNoOfParams - 10.0 * dRastriginPart;

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
    final int iNoOfParams = daParams.length;

    double dSpherePart1 = 0.0;
    double dSpherePart2 = 0.0;
    double dRastriginPart = 0.0;
    for (int i = 0; i < iNoOfParams; i++) {

      // boundaries
      if (daParams[i] > 5.0 || daParams[i] < -5.0) {
        // needs to be negative so that it gets positive in the end...
        dRastriginPart -= daParams[i] * daParams[i];
        grad[i] += FixedValues.NONCONVERGEDGRADIENT;
        continue;
      }

      final double mu0P = daParams[i] - dMu0;
      final double mu1P = daParams[i] - dMu1;
      dSpherePart1 += mu0P * mu0P;
      dSpherePart2 += mu1P * mu1P;

      final double t1 = 2 * Math.PI * (daParams[i] - dMu0);
      dRastriginPart += Math.cos(t1);
      grad[i] = +20 * Math.PI * Math.sin(t1); // *2.0*10.0 because of the *-10.0 in dTotal
    }

    // grad[] contains the Rastrigin gradients now, add the sphere part
    double dSphereContrib;
    final double t1 = dD * iNoOfParams + dS * dSpherePart2;
    if (dSpherePart1 < t1) {
      dSphereContrib = dSpherePart1;
      for (int i = 0; i < iNoOfParams; i++) grad[i] += 2 * (daParams[i] - dMu0);
    } else {
      dSphereContrib = t1;
      for (int i = 0; i < iNoOfParams; i++) grad[i] += 2 * dS * (daParams[i] - dMu1);
    }

    final double dTotal = dSphereContrib + 10.0 * iNoOfParams - 10.0 * dRastriginPart;

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
