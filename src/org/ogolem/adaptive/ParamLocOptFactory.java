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

import java.util.HashMap;
import org.ogolem.adaptive.genericfitness.GenericFitnessFunction;
import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.GenericLocOpt;
import org.ogolem.locopt.AbstractLocOptFactory;

/**
 * The adaptive parameter version of the local optimization factory.
 *
 * @author Johannes Dieterich
 * @version 2020-07-25
 */
class ParamLocOptFactory extends AbstractLocOptFactory<Double, AdaptiveParameters> {

  private static final long serialVersionUID = (long) 20150428;

  private final AdaptiveConf adapConf;
  private final double[] minBorders;
  private final double[] maxBorders;
  private final boolean normParams;

  ParamLocOptFactory(
      final double absConvThreshold,
      final int maxIter,
      final boolean normParams,
      final double[] minBorders,
      final double[] maxBorders,
      final HashMap<String, GenericBackend<Double, AdaptiveParameters>> backendDict,
      final AdaptiveConf adaptConf) {
    super(absConvThreshold, maxIter, backendDict);
    this.adapConf = adaptConf;
    this.normParams = normParams;
    this.minBorders = minBorders;
    this.maxBorders = maxBorders;
  }

  @Override
  protected GenericLocOpt<Double, AdaptiveParameters> buildSpecializedLocalOptimization(
      final String input) throws Exception {

    if (input.startsWith("external:")) {

      final String[] opts = tokenizeSecondLevel(input.substring("external:".length()).trim());

      GenericBackend<Double, AdaptiveParameters> backend = null;
      for (final String opt : opts) {
        if (opt.startsWith("backend=")) {
          backend = getBackend(opt.substring(8));
        } else {
          throw new RuntimeException("Unknown option " + opt + " in single point only!");
        }
      }

      checkSanityBackend(backend);

      return new ExternalLocOpt(backend);
    }

    // no specialized local optimization left

    return null;
  }

  @Override
  public GenericBackend<Double, AdaptiveParameters> parseBackend(final String backend)
      throws Exception {

    if (!backend.equalsIgnoreCase("alldefault")) {
      return null;
    }

    System.out.println(
        "INFO: we fall back to all default backend definition. This is an intermediate solution that will be removed soon!");
    final GenericFitnessFunction func = adapConf.getCartesianFitnessFunc();
    final GenericBackend<Double, AdaptiveParameters> back =
        new FitFuncToBackend(func, normParams, minBorders, maxBorders);

    return back;
  }
}
