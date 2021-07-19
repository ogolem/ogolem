/*
Copyright (c) 2015-2021, J. M. Dieterich and B. Hartke
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
import org.junit.jupiter.api.Test;

/**
 * A test class for the LAPACK-style interface
 *
 * @author Johannes Dieterich
 * @version 2021-07-17
 */
public class LAPACKInterfaceTest {

  public LAPACKInterfaceTest() {}

  /** Test of invertMatrix method, of class LAPACKInterface. */
  @Test
  public void testInvertMatrix() {
    System.out.println("invertMatrix");

    final double[] matrix =
        new double[] {
          1, 0, 5,
          2, 1, 6,
          3, 4, 0
        };

    final double[][] matJama =
        new double[][] {
          new double[] {1, 2, 3},
          new double[] {0, 1, 4},
          new double[] {5, 6, 0}
        };
    final Matrix mat = new Matrix(matJama);
    final Matrix inv = mat.inverse();
    final double[][] invMat = inv.getArray();

    final int dim = 3;
    final double[] work = new double[dim * dim];
    final int[] piv = new int[dim];

    LAPACKInterface.invertMatrix(dim, matrix, work, piv);

    final double thresh = 1e-7;
    assertEquals(invMat[0][0], matrix[0], thresh);
    assertEquals(invMat[1][0], matrix[1], thresh);
    assertEquals(invMat[2][0], matrix[2], thresh);
    assertEquals(invMat[0][1], matrix[3], thresh);
    assertEquals(invMat[1][1], matrix[4], thresh);
    assertEquals(invMat[2][1], matrix[5], thresh);
    assertEquals(invMat[0][2], matrix[6], thresh);
    assertEquals(invMat[1][2], matrix[7], thresh);
    assertEquals(invMat[2][2], matrix[8], thresh);
  }
}
