/*
Copyright (c) 2020, J. M. Dieterich and B. Hartke
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
 * Calls an external code to be an "adaptivable".
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class AdaptiveExternalCaller extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20200725;

  private final int noParams;

  AdaptiveExternalCaller(final int noParameters) {
    this.noParams = noParameters;
  }

  AdaptiveExternalCaller(final AdaptiveExternalCaller orig) {
    this.noParams = orig.noParams;
  }

  @Override
  public AdaptiveExternalCaller copy() {
    return new AdaptiveExternalCaller(this);
  }

  @Override
  public double energyOfStructWithParams(
      CartesianCoordinates cartes, AdaptiveParameters params, int geomID, BondInfo bonds) {
    throw new RuntimeException("not supported by external adaptivable.");
  }

  @Override
  public double gradientOfStructWithParams(
      CartesianCoordinates cartes,
      AdaptiveParameters params,
      int geomID,
      BondInfo bonds,
      double[] grad) {
    throw new RuntimeException("not supported by external adaptivable.");
  }

  @Override
  public double[][] minMaxBordersForParams(AdaptiveParameters params) {

    // really should never be used, instead specify custom parameter bounds in the input.

    final double[][] minMax = new double[2][noParams];

    for (int i = 0; i < noParams; i++) {
      minMax[0][i] = -1000.0;
      minMax[1][i] = 1000.0;
    }

    return minMax;
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      ArrayList<CartesianCoordinates> refCartes, String sMethod) {

    final String[] saAtoms = {"XX"};
    final int[] iaParamsPerAt = {noParams};
    final AdaptiveParameters paramStub =
        new AdaptiveParameters(noParams, -1, saAtoms, iaParamsPerAt, sMethod);

    return paramStub;
  }
}
