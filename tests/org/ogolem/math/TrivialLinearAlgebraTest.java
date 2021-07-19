/*
Copyright (c) 2020-2021, J. M. Dieterich and B. Hartke
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
package org.ogolem.math;

import static org.junit.jupiter.api.Assertions.*;

import contrib.jama.Matrix;
import java.util.Random;
import org.junit.jupiter.api.Test;

/**
 * Tests for trivial linear algebra
 *
 * @author Johannes Dieterich
 * @version 2021-07-17
 */
public class TrivialLinearAlgebraTest {

  @Test
  public void testMatMult() {
    System.out.println("matMult");

    final Matrix matA = new Matrix(3, 3);
    final Matrix matB = new Matrix(3, 1000);

    final Random r = new Random(42);

    final double[][] matADat = matA.getArray();
    final double[][] matBDat = matB.getArray();

    for (int i = 0; i < matADat.length; i++) {
      for (int j = 0; j < matADat[i].length; j++) {
        matADat[i][j] = r.nextDouble();
      }
    }

    for (int i = 0; i < matBDat.length; i++) {
      for (int j = 0; j < matBDat[i].length; j++) {
        matBDat[i][j] = r.nextDouble();
      }
    }

    // Jama reference path
    final Matrix matC = matA.times(matB);
    final double[][] matCDat = matC.getArray();

    // our implementation
    final double[][] matCOurDat = TrivialLinearAlgebra.matMult(matADat, matBDat);

    assertEquals(matCDat.length, matCOurDat.length);

    for (int i = 0; i < matCOurDat.length; i++) {
      assertArrayEquals(matCDat[i], matCOurDat[i], 1e-20);
    }
  }
}
