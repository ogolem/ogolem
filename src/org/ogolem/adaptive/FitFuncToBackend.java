/*
Copyright (c) 2015-2020, J. M. Dieterich and B. Hartke
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

import org.ogolem.adaptive.genericfitness.GenericFitnessFunction;
import org.ogolem.generic.GenericBackend;

/**
 * An adapter from adaptivable to backend.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
final class FitFuncToBackend implements GenericBackend<Double, AdaptiveParameters> {

  private static final long serialVersionUID = (long) 2020429;

  private final GenericFitnessFunction fitfunc;
  private final boolean normParams;
  private final double[] minBorders;
  private final double[] maxBorders;

  private AdaptiveParameters params;

  FitFuncToBackend(
      final GenericFitnessFunction fitfunc,
      final boolean normParams,
      final double[] minBorders,
      final double[] maxBorders) {
    this.fitfunc = fitfunc;
    this.normParams = normParams;
    this.minBorders = minBorders;
    this.maxBorders = maxBorders;
  }

  FitFuncToBackend(final FitFuncToBackend orig) {
    this.fitfunc = orig.fitfunc.copy();
    this.normParams = orig.normParams;
    this.maxBorders = orig.maxBorders.clone(); // shallow enough
    this.minBorders = orig.minBorders.clone();
  }

  @Override
  public GenericBackend<Double, AdaptiveParameters> copy() {
    return new FitFuncToBackend(this);
  }

  @Override
  public String getMyID() {
    return fitfunc.getClass().getName();
  }

  @Override
  public int numberOfActiveCoordinates(final AdaptiveParameters individual) {
    return individual.getNumberOfParamters();
  }

  @Override
  public double[] getActiveCoordinates(final AdaptiveParameters individual) {
    this.params = individual;
    return params.getAllParamters();
  }

  @Override
  public double gradient(final double[] currCoords, final double[] gradient, final int iteration) {

    paramsBackToBorders(currCoords);
    params.copyParameters(currCoords);
    final ParameterGradient grad = fitfunc.evaluateGradient(params);
    final double[] gg = grad.getGradient();

    System.arraycopy(gg, 0, gradient, 0, gg.length);

    return grad.getFunctionValue();
  }

  @Override
  public void resetToStable(final double[] coordinates) {
    paramsBackToBorders(coordinates);
  }

  @Override
  public void updateActiveCoordinates(
      final AdaptiveParameters individual, final double[] coordinates) {
    paramsBackToBorders(coordinates);
    individual.copyParameters(coordinates);
  }

  @Override
  public double fitness(final double[] currCoords, final int iteration) {
    paramsBackToBorders(currCoords);
    params.copyParameters(currCoords);
    return fitfunc.evaluateFitness(params);
  }

  @Override
  public AdaptiveParameters fitness(
      final AdaptiveParameters individual, final boolean forceOneEval) {

    final double funcVal = fitfunc.evaluateFitness(individual);
    individual.setFitness(funcVal);

    return individual;
  }

  private void paramsBackToBorders(final double[] params) {

    if (!normParams) return;

    for (int i = 0; i < params.length; i++) {
      if (params[i] < minBorders[i]) params[i] = minBorders[i];
      else if (params[i] > maxBorders[i]) params[i] = maxBorders[i];
    }
  }

  @Override
  public BOUNDSTYPE boundariesInRepresentation(final AdaptiveParameters individual) {
    return BOUNDSTYPE.ALL;
  }

  @Override
  public void bestEstimateBoundaries(
      final double[] currCoords, final double[] low, final double[] high) {

    System.arraycopy(minBorders, 0, low, 0, minBorders.length);
    System.arraycopy(maxBorders, 0, high, 0, minBorders.length);
  }
}
