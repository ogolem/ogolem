/*
Copyright (c) 2021, J. M. Dieterich and B. Hartke
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
package org.ogolem.microbenchmarks;

import java.util.Random;
import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.adaptive.BenchSchwefels;

/**
 * Benchmark the speed of computing Schwefel's benchmark function in n-D
 *
 * @author Johannes Dieterich
 * @version 2021-07-21
 */
final class SchwefelBench implements SingleMicroBenchmark {

  private final int noDims;
  private final BenchSchwefels schwefel;
  private final AdaptiveParameters p;

  SchwefelBench(final int noDims) {
    assert (noDims > 0);
    this.noDims = noDims;
    this.schwefel = new BenchSchwefels(noDims);
    this.p =
        new AdaptiveParameters(
            noDims, -1, new String[] {"XXX"}, new int[] {noDims}, "benchschwefel");

    final Random r = new Random(71); // seed to make sure it's always the same
    double[] params = p.getAllParamters();
    for (int i = 0; i < params.length; i++) {
      params[i] = -500 + r.nextDouble() * 1000;
    }
  }

  @Override
  public double runSingle() throws Exception {
    return schwefel.energyOfStructWithParams(null, p, 0, null);
  }

  @Override
  public String name() {
    return "Schwefel function benchmark in " + noDims + "D";
  }
}
