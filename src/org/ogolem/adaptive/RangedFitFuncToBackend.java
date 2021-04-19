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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.adaptive.genericfitness.GenericFitnessFunction;
import org.ogolem.generic.Copyable;
import org.ogolem.generic.GenericBackend;

/**
 * A ranged fitness function. E.g., it can be used to optimize only a subset of parameters.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
class RangedFitFuncToBackend implements GenericBackend<Double, AdaptiveParameters> {

  private static final long serialVersionUID = (long) 2020429;

  private final GenericFitnessFunction fitfunc;
  private final boolean normParams;
  private final List<Range> ranges;
  private final double[] fullLowerParamBound;
  private final double[] fullUpperParamBound;
  private AdaptiveParameters paramCache;
  // private double[] fullGradCache;

  public RangedFitFuncToBackend(
      final GenericFitnessFunction fitFunc,
      final boolean normParams,
      final List<Range> ranges,
      final double[] fullMinBorders,
      final double[] fullMaxBorders) {

    System.err.println(
        "WARNING: THE FEATURE OF RANGED PARAMETER FITNESS FUNCTION IS NOT WELL TESTED!");
    System.err.println(
        "WARNING: THE FEATURE OF RANGED PARAMETER FITNESS FUNCTION IS NOT WELL TESTED!");
    System.err.println(
        "WARNING: THE FEATURE OF RANGED PARAMETER FITNESS FUNCTION IS NOT WELL TESTED!");
    System.err.println(
        "WARNING: THE FEATURE OF RANGED PARAMETER FITNESS FUNCTION IS NOT WELL TESTED!");
    System.err.println(
        "WARNING: THE FEATURE OF RANGED PARAMETER FITNESS FUNCTION IS NOT WELL TESTED!");

    this.fitfunc = fitFunc;
    this.ranges = ranges;
    this.fullLowerParamBound = fullMinBorders;
    this.fullUpperParamBound = fullMaxBorders;
    this.normParams = normParams;
  }

  RangedFitFuncToBackend(final RangedFitFuncToBackend orig) {
    this.fitfunc = orig.fitfunc.copy();
    this.ranges = new ArrayList<>();
    orig.ranges.forEach(
        (r) -> {
          ranges.add(r.copy());
        });
    this.fullLowerParamBound = orig.fullLowerParamBound.clone();
    this.fullUpperParamBound = orig.fullUpperParamBound.clone();
    this.normParams = orig.normParams;
  }

  @Override
  public GenericBackend<Double, AdaptiveParameters> copy() {
    return new RangedFitFuncToBackend(this);
  }

  @Override
  public String getMyID() {
    return "ranged: " + fitfunc.getClass().getName();
  }

  @Override
  public int numberOfActiveCoordinates(final AdaptiveParameters individual) {

    // loop over ranges and add up
    int noParams = 0;
    noParams = ranges.stream().map((r) -> (r.up - r.low + 1)).reduce(noParams, Integer::sum);

    return noParams;
  }

  @Override
  public double[] getActiveCoordinates(final AdaptiveParameters individual) {

    // shouldn't be too expensive, so do it again instead of caching
    final int noParams = numberOfActiveCoordinates(individual);

    final double[] redParams = new double[noParams];
    final double[] allParams = individual.getAllParamters();
    int off = 0;
    for (final Range r : ranges) {
      final int no = (r.up - r.low + 1);
      System.arraycopy(allParams, r.low, redParams, off, no);
      off += no;
    }

    this.paramCache = individual;
    // this.fullGradCache = new double[individual.getNumberOfParamters()];

    return redParams;
  }

  @Override
  public double fitness(final double[] currCoords, final int iteration) {

    updateActiveCoordinates(paramCache, currCoords);

    final double fit = fitfunc.evaluateFitness(paramCache);

    return fit;
  }

  @Override
  public double gradient(final double[] currCoords, final double[] gradient, final int iteration) {

    updateActiveCoordinates(paramCache, currCoords);
    final ParameterGradient grad = fitfunc.evaluateGradient(paramCache);
    final double[] fullGrad = grad.getTotalGradient();
    int off = 0;
    for (final Range r : ranges) {
      final int no = (r.up - r.low + 1);
      System.arraycopy(fullGrad, r.low, gradient, off, no);
      off += no;
    }

    return grad.getFunctionValue();
  }

  @Override
  public AdaptiveParameters fitness(
      final AdaptiveParameters individual, final boolean forceOneEval) {

    final double fit = fitfunc.evaluateFitness(individual);

    final AdaptiveParameters copy = individual.copy();
    copy.setFitness(fit);

    return copy;
  }

  @Override
  public void updateActiveCoordinates(
      final AdaptiveParameters individual, final double[] coordinates) {

    paramsBackToBorders(coordinates);
    final double[] allCoords = individual.getAllParamters();
    int off = 0;
    for (final Range r : ranges) {
      final int no = (r.up - r.low + 1);
      System.arraycopy(coordinates, off, allCoords, r.low, no);
      off += no;
    }
  }

  @Override
  public BOUNDSTYPE boundariesInRepresentation(final AdaptiveParameters individual) {
    return BOUNDSTYPE.ALL;
  }

  @Override
  public void bestEstimateBoundaries(
      final double[] currCoords, final double[] low, final double[] high) {

    int off = 0;
    for (final Range r : ranges) {
      final int no = (r.up - r.low + 1);
      System.arraycopy(fullLowerParamBound, r.low, low, off, no);
      System.arraycopy(fullUpperParamBound, r.low, high, off, no);
      off += no;
    }
  }

  @Override
  public void resetToStable(final double[] coordinates) {
    // no-op
  }

  private void paramsBackToBorders(final double[] params) {

    if (!normParams) return;

    int off = 0;
    for (final Range r : ranges) {
      final int no = (r.up - r.low + 1);
      for (int fullInd = r.low; fullInd < r.up + 1; fullInd++) {
        if (params[off] < fullLowerParamBound[fullInd]) params[off] = fullLowerParamBound[fullInd];
        if (params[off] > fullUpperParamBound[fullInd]) params[off] = fullUpperParamBound[fullInd];
        off++;
      }
    }
  }

  static class Range implements Serializable, Copyable {

    private static final long serialVersionUID = (long) 20150404;

    private final int low;
    private final int up;

    /**
     * Constructs a range.
     *
     * @param low the lower end of this range (inclusive)
     * @param up the upper end of this range (inclusive)
     */
    Range(final int low, final int up) {
      this.up = up;
      this.low = low;
    }

    Range(final Range orig) {
      this.up = orig.up;
      this.low = orig.low;
    }

    @Override
    public Range copy() {
      return new Range(this);
    }

    public int getLow() {
      return low;
    }

    public int getUp() {
      return up;
    }
  }
}
