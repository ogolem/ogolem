/*
Copyright (c) 2020, 2023 J. M. Dieterich and B. Hartke
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

/**
 * 3x3 matrix.
 *
 * @author Johannes Dieterich
 * @version 2023-06-28
 */
public final record Matrix3x3(
    double a00,
    double a01,
    double a02,
    double a10,
    double a11,
    double a12,
    double a20,
    double a21,
    double a22) {

  public final double determinant() {

    final double sub1 = sub2x2Det(a11, a12, a21, a22);
    final double sub2 = sub2x2Det(a10, a12, a20, a22);
    final double sub3 = sub2x2Det(a10, a11, a20, a21);
    final double tmp2 = -a01 * sub2;
    final double tmp12 = Math.fma(a00, sub1, tmp2);
    final double det = Math.fma(a02, sub3, tmp12);

    return det;
  }

  private static double sub2x2Det(
      final double a00, final double a01, final double a10, final double a11) {

    final double a101 = -a10 * a01;
    final double det = Math.fma(a00, a11, a101);

    return det;
  }

  public Matrix3x3 inverse() {

    final double det = determinant();
    if (Math.abs(det) < 1e-10) {
      System.out.println("Cannot invert 3x3 matrix - matrix is singular");
      return null;
    }

    // use method of cofactors:
    // https://metric.ma.ic.ac.uk/metric_public/matrices/inverses/inverses2.html
    final double a00minor = sub2x2Det(a11, a12, a21, a22);
    final double a01minor = -sub2x2Det(a10, a12, a20, a22);
    final double a02minor = sub2x2Det(a10, a11, a20, a21);
    final double a10minor = -sub2x2Det(a01, a02, a21, a22);
    final double a11minor = sub2x2Det(a00, a02, a20, a22);
    final double a12minor = -sub2x2Det(a00, a01, a20, a21);
    final double a20minor = sub2x2Det(a01, a02, a11, a12);
    final double a21minor = -sub2x2Det(a00, a02, a10, a12);
    final double a22minor = sub2x2Det(a00, a01, a10, a11);

    final double a00Inv = a00minor / det;
    final double a01Inv = a10minor / det;
    final double a02Inv = a20minor / det;
    final double a10Inv = a01minor / det;
    final double a11Inv = a11minor / det;
    final double a12Inv = a21minor / det;
    final double a20Inv = a02minor / det;
    final double a21Inv = a12minor / det;
    final double a22Inv = a22minor / det;

    Matrix3x3 inverse =
        new Matrix3x3(a00Inv, a01Inv, a02Inv, a10Inv, a11Inv, a12Inv, a20Inv, a21Inv, a22Inv);

    return inverse;
  }

  public Matrix3x3 multiply(final Matrix3x3 matB) {

    final double a00s = a00 * matB.a00;
    final double a001 = Math.fma(a01, matB.a10, a00s);
    final double a00New = Math.fma(a02, matB.a20, a001);
    final double a01s = a00 * matB.a01;
    final double a011 = Math.fma(a01, matB.a11, a01s);
    final double a01New = Math.fma(a02, matB.a21, a011);
    final double a02s = a00 * matB.a02;
    final double a021 = Math.fma(a01, matB.a12, a02s);
    final double a02New = Math.fma(a02, matB.a21, a021);
    final double a10s = a10 * matB.a00;
    final double a101 = Math.fma(a11, matB.a10, a10s);
    final double a10New = Math.fma(a12, matB.a20, a101);
    final double a11s = a10 * matB.a01;
    final double a111 = Math.fma(a11, matB.a11, a11s);
    final double a11New = Math.fma(a12, matB.a21, a111);
    final double a12s = a10 * matB.a02;
    final double a121 = Math.fma(a11, matB.a12, a12s);
    final double a12New = Math.fma(a12, matB.a21, a121);
    final double a20s = a20 * matB.a00;
    final double a201 = Math.fma(a21, matB.a10, a20s);
    final double a20New = Math.fma(a22, matB.a20, a201);
    final double a21s = a20 * matB.a01;
    final double a211 = Math.fma(a21, matB.a11, a21s);
    final double a21New = Math.fma(a22, matB.a21, a211);
    final double a22s = a20 * matB.a02;
    final double a221 = Math.fma(a21, matB.a12, a22s);
    final double a22New = Math.fma(a22, matB.a21, a221);

    return new Matrix3x3(a00New, a01New, a02New, a10New, a11New, a12New, a20New, a21New, a22New);
  }
}
