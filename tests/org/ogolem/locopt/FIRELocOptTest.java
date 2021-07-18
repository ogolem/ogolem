/*
Copyright (c) 2019, M. Dittner, B. Hartke
              2020, D. Behrens, M. Dittner, B. Hartke
              2021, D. Behrens, M. Dittner, J. M. Dieterich, B. Hartke
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

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
package org.ogolem.locopt;

import static org.junit.jupiter.api.Assertions.*;

import java.util.Random;
import org.junit.jupiter.api.Test;

/**
 * A test class for the FIRE.
 *
 * @author Dominik Behrens
 * @version 2021-07-17
 */
public class FIRELocOptTest {

  /** Testing the FIRE's optimize method. */
  @Test
  public void testOptimize() {
    System.out.println("optimize");

    // Generate harmonic Function
    final Random r = new Random();
    final double a = 0.0001;
    final double shift = 4.2;
    final double b = 0.0;
    final HarmonicFunction func = new HarmonicFunction(a, shift, b);

    // Some basic FIRE setup
    FIRELocOpt.FireConfig fc = new FIRELocOpt.FireConfig();
    fc.iMaxIterations = 1000;
    fc.dFMax = 1e-12;
    fc.dMaxMove = 0.2;
    final FIRELocOpt<Double, BasicOptimizableType> instance1 = new FIRELocOpt<>(fc, func);

    int noTrials = 100;
    int dim = 42;
    double numAcc = 1e-7;

    for (int trial = 0; trial < noTrials; trial++) {
      final double[] gen1 = new double[dim];

      // random init
      for (int x = 0; x < dim; x++) {
        gen1[x] = 42.0 * r.nextDouble() + shift;
      }

      final BasicOptimizableType ind1 = new BasicOptimizableType(gen1);

      final BasicOptimizableType res1 = instance1.optimize(ind1);

      // check
      assertTrue(Math.abs(res1.getFitness()) <= numAcc, "Fitness (I) wrong: " + res1.getFitness());

      final double[] res1Dat = res1.getGenomeAsDouble();
      for (int x = 0; x < dim; x++) {
        assertTrue(
            Math.abs(res1Dat[x] - shift) <= 10 * numAcc,
            "Error (I) from correct solution: "
                + Math.abs(res1Dat[x] - shift)
                + " in "
                + x); // make it a bit bigger
      }
    }
  }
}
