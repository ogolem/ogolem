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
package org.ogolem.adaptive;

import static org.junit.jupiter.api.Assertions.*;

import java.util.Arrays;
import java.util.Random;
import org.junit.jupiter.api.Test;

/**
 * A test class for Schaffer's F7 function
 *
 * @author Johannes Dieterich
 * @version 2021-07-17
 */
public class SchafferF7Test {

  @Test
  public void testMinimum() {

    final BenchSchaffer schaff = new BenchSchaffer(42);

    final AdaptiveParameters p = createStub(42);
    double[] params = p.getAllParamters();
    for (int i = 0; i < params.length; i++) {
      params[i] = 0.0; // known global minimum
    }

    final double fit = schaff.energyOfStructWithParams(null, p, 0, null);

    assertEquals(
        0.0, fit, 1e-10, "Function value of Schaffer F7 does not match global minimum. Is " + fit);
  }

  @Test
  public void testMinimum2() {

    final BenchSchaffer schaff = new BenchSchaffer(71);

    final AdaptiveParameters p = createStub(71);
    double[] params = p.getAllParamters();
    for (int i = 0; i < params.length; i++) {
      params[i] = 0.0; // known global minimum
    }

    final double[] grad = new double[71];
    final double fit = schaff.gradientOfStructWithParams(null, p, 0, null, grad);

    assertEquals(
        0.0, fit, 1e-10, "Function value of Schaffer F7 does not match global minimum. Is " + fit);

    for (int i = 0; i < grad.length; i++) {
      assertEquals(
          0.0,
          grad[i],
          1e-10,
          "Gradient value " + i + " of Schaffer F7 does not match minimum. Is " + grad[i]);
    }
  }

  @Test
  public void testRandomPoints() {

    final BenchSchaffer schaff = new BenchSchaffer(57);

    final AdaptiveParameters p = createStub(57);
    final Random r = new Random();

    final double[] grad = new double[57];
    final double[] numGrad = new double[57];

    for (int iter = 0; iter < 100; iter++) {

      final double[] params = p.getAllParamters();
      for (int i = 0; i < params.length; i++) {
        params[i] = -100.0 + r.nextDouble() * 200;
      }

      Arrays.fill(grad, 0.0);
      Arrays.fill(numGrad, 0.0);
      final double fit = schaff.energyOfStructWithParams(null, p, 0, null);
      final double fit2 = schaff.gradientOfStructWithParams(null, p, 0, null, grad);

      assertEquals(
          fit, fit2, 1e-10, "Mismatch between function value from function-only and gradient");

      // calculate a numerical gradient and compare it to the analytical one
      final double numE = NumericalGradients.calculateParamGrad(null, p, schaff, -1, null, numGrad);

      for (int i = 0; i < grad.length; i++) {
        final double gradAcc = (Math.abs(grad[i]) < 1.0) ? 5e-4 : 5e-4 * Math.abs(grad[i]);
        assertEquals(
            numGrad[i], grad[i], gradAcc, "Mismatch between numerical and analytical gradient.");
      }
    }
  }

  private AdaptiveParameters createStub(final int noDims) {

    final AdaptiveParameters p =
        new AdaptiveParameters(
            noDims, -1, new String[] {"XXX"}, new int[] {noDims}, "benchschafferf7");
    return p;
  }
}
