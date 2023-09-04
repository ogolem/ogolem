/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015-2023, J. M. Dieterich and B. Hartke
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

import java.util.Arrays;
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Performs some really trivial linear algebra.
 *
 * @author Johannes Dieterich
 * @version 2023-06-28
 */
public class TrivialLinearAlgebra {

  private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

  /**
   * Multiplications of two matrices where the resulting matrices is handed over to the method.
   *
   * @param matrixA
   * @param matrixB
   * @param res overridden on exit
   */
  public static void matMult(
      final double[][] matrixA, final double[][] matrixB, final double[][] res) {

    final int m = matrixA.length;
    final int n = matrixA[0].length;
    final int s = matrixB[0].length;

    matMult(matrixA, matrixB, res, m, n, s);
  }

  /**
   * Multiplications of two matrices where the resulting matrices is handed over to the method.
   *
   * @param matrixA
   * @param matrixB
   * @param res overridden on exit
   * @param m
   * @param n
   * @param s
   */
  public static void matMult(
      final double[][] matrixA,
      final double[][] matrixB,
      final double[][] res,
      final int m,
      final int n,
      final int s) {

    assert (res.length >= m);
    assert (res[0].length >= s);
    assert (matrixA.length >= m);
    assert (matrixA[0].length >= n);
    assert (matrixB.length >= n);
    assert (matrixB[0].length >= s);

    // in BLAS notation: M = m, N = s, K = n
    // we hence want to loop in order M=m (index i), K=n (index l), N=s (index j)

    final int upperBound = SPECIES.loopBound(s);

    for (int i = 0; i < m; i++) {

      final double[] rowRes = res[i];
      // zero the result
      Arrays.fill(rowRes, 0.0);

      for (int l = 0; l < n; l++) {
        final double scalA = matrixA[i][l];
        final double[] rowB = matrixB[l];

        int j = 0;
        if (upperBound > 0) {
          final var vScalA = DoubleVector.broadcast(SPECIES, scalA);
          for (; j < upperBound; j += SPECIES.length()) {
            final var vRowB = DoubleVector.fromArray(SPECIES, rowB, j);
            final var vRowRes = DoubleVector.fromArray(SPECIES, rowRes, j);
            final var vUpRes = vRowB.fma(vScalA, vRowRes);
            vUpRes.intoArray(rowRes, j);
          }
        }
        for (; j < s; j++) {
          rowRes[j] += scalA * rowB[j];
        }
      }
    }
  }

  public static void matMult(
      final double[][] matrixA,
      final double[][] matrixB,
      final double[][] res,
      final int m,
      final int n,
      final int s,
      final int offAX,
      final int offAY,
      final int offBX,
      final int offBY) {

    assert (res.length >= m);
    assert (res[0].length >= s);
    assert (matrixA.length >= m);
    assert (matrixA[0].length >= n);
    assert (matrixB.length >= n);
    assert (matrixB[0].length >= s);

    final int upperBound = SPECIES.loopBound(s);

    for (int i = 0; i < m; i++) {

      final double[] rowRes = res[i];
      // zero the result
      Arrays.fill(rowRes, 0.0);

      for (int l = 0; l < n; l++) {
        final double scalA = matrixA[i + offAX][l + offAY];
        final double[] rowB = matrixB[l + offBX];
        int j = 0;
        if (upperBound > 0) {
          final var vScalA = DoubleVector.broadcast(SPECIES, scalA);
          for (; j < upperBound; j += SPECIES.length()) {
            final var vRowB = DoubleVector.fromArray(SPECIES, rowB, j + offBY);
            final var vRowRes = DoubleVector.fromArray(SPECIES, rowRes, j);
            final var vUpRes = vRowB.fma(vScalA, vRowRes);
            vUpRes.intoArray(rowRes, j);
          }
        }
        for (; j < s; j++) {
          rowRes[j] = Math.fma(scalA, rowB[j + offBY], rowRes[j]);
        }
      }
    }
  }

  /**
   * This multiplies a matrix with another matrix.
   *
   * @param matrixA the first matrix.
   * @param matrixB The second Matrix.
   * @return The result.
   */
  public static double[][] matMult(final double[][] matrixA, final double[][] matrixB) {

    final int m = matrixA.length;
    final int s = matrixB[0].length;

    final double[][] matrixC = new double[m][s];
    matMult(matrixA, matrixB, matrixC, m, matrixA[0].length, s);

    return matrixC;
  }

  public static void matMult(
      final Matrix3x3 matrixA, final double[][] matrixB, final int s, final double[][] res) {

    assert (matrixA != null);
    assert (matrixB != null);
    assert (res != null);
    assert (s >= 0);
    assert (matrixB.length == 3);
    assert (res.length == 3);
    assert (matrixB[0].length >= s);
    assert (matrixB[1].length >= s);
    assert (matrixB[2].length >= s);
    assert (res[0].length >= s);
    assert (res[1].length >= s);
    assert (res[2].length >= s);

    final double[] res0 = res[0];
    final double[] res1 = res[1];
    final double[] res2 = res[2];

    final int upperBound = SPECIES.loopBound(s);
    if (upperBound == 0) {
      // special case this - think rotating H2O
      for (int j = 0; j < s; j++) {

        final double b0j = matrixB[0][j];

        final double res00j = matrixA.a00() * b0j;
        final double res10j = matrixA.a10() * b0j;
        final double res20j = matrixA.a20() * b0j;

        final double b1j = matrixB[1][j];

        final double res01j = Math.fma(matrixA.a01(), b1j, res00j);
        final double res11j = Math.fma(matrixA.a11(), b1j, res10j);
        final double res21j = Math.fma(matrixA.a21(), b1j, res20j);

        final double b2j = matrixB[2][j];

        res0[j] = Math.fma(matrixA.a02(), b2j, res01j);
        res1[j] = Math.fma(matrixA.a12(), b2j, res11j);
        res2[j] = Math.fma(matrixA.a22(), b2j, res21j);
      }

      return;
    }

    int j = 0;
    final var vA00 = DoubleVector.broadcast(SPECIES, matrixA.a00());
    final var vA10 = DoubleVector.broadcast(SPECIES, matrixA.a10());
    final var vA20 = DoubleVector.broadcast(SPECIES, matrixA.a20());
    final var vA01 = DoubleVector.broadcast(SPECIES, matrixA.a01());
    final var vA11 = DoubleVector.broadcast(SPECIES, matrixA.a11());
    final var vA21 = DoubleVector.broadcast(SPECIES, matrixA.a21());
    final var vA02 = DoubleVector.broadcast(SPECIES, matrixA.a02());
    final var vA12 = DoubleVector.broadcast(SPECIES, matrixA.a12());
    final var vA22 = DoubleVector.broadcast(SPECIES, matrixA.a22());

    for (; j < upperBound; j += SPECIES.length()) {
      final var vB00 = DoubleVector.fromArray(SPECIES, matrixB[0], j);
      final var vRes0 = vB00.mul(vA00);
      final var vRes1 = vB00.mul(vA10);
      final var vRes2 = vB00.mul(vA20);

      final var vB01 = DoubleVector.fromArray(SPECIES, matrixB[1], j);
      final var vRes0T = vB01.fma(vA01, vRes0);
      final var vRes1T = vB01.fma(vA11, vRes1);
      final var vRes2T = vB01.fma(vA21, vRes2);

      final var vB02 = DoubleVector.fromArray(SPECIES, matrixB[2], j);
      final var vRes0TT = vB02.fma(vA02, vRes0T);
      vRes0TT.intoArray(res0, j);
      final var vRes1TT = vB02.fma(vA12, vRes1T);
      vRes1TT.intoArray(res1, j);
      final var vRes2TT = vB02.fma(vA22, vRes2T);
      vRes2TT.intoArray(res2, j);
    }

    for (; j < s; j++) {

      final double b0j = matrixB[0][j];

      final double res00j = matrixA.a00() * b0j;
      final double res10j = matrixA.a10() * b0j;
      final double res20j = matrixA.a20() * b0j;

      final double b1j = matrixB[1][j];

      final double res01j = Math.fma(matrixA.a01(), b1j, res00j);
      final double res11j = Math.fma(matrixA.a11(), b1j, res10j);
      final double res21j = Math.fma(matrixA.a21(), b1j, res20j);

      final double b2j = matrixB[2][j];

      res0[j] = Math.fma(matrixA.a02(), b2j, res01j);
      res1[j] = Math.fma(matrixA.a12(), b2j, res11j);
      res2[j] = Math.fma(matrixA.a22(), b2j, res21j);
    }

    return;
  }

  /**
   * Calculates the dot product of two vectors.
   *
   * @param v1 array
   * @param v2 array
   * @return the dot product
   */
  public static double dotProduct(final double[] v1, final double[] v2) {

    assert (v1.length == v2.length);
    final int upperBound = SPECIES.loopBound(v1.length);

    var dVec = DoubleVector.zero(SPECIES);
    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var v1Vec = DoubleVector.fromArray(SPECIES, v1, i);
      final var v2Vec = DoubleVector.fromArray(SPECIES, v2, i);
      dVec = v1Vec.fma(v2Vec, dVec);
    }
    double d = dVec.reduceLanes(VectorOperators.ADD);

    for (; i < v1.length; i++) {
      d = Math.fma(v1[i], v2[i], d);
    }

    return d;
  }

  /**
   * Calculates the dot product of a vector with itself.
   *
   * @param n the number of elements to be used
   * @param v the vector to be used
   * @return the self dot product
   */
  public static double selfDotProduct(final int n, final double[] v) {

    assert (n <= v.length);
    if (n <= 0) return 0.0;
    final int upperBound = SPECIES.loopBound(n);

    var dVec = DoubleVector.zero(SPECIES);
    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vVec = DoubleVector.fromArray(SPECIES, v, i);
      dVec = vVec.fma(vVec, dVec);
    }
    double d = dVec.reduceLanes(VectorOperators.ADD);

    for (; i < n; i++) {
      d = Math.fma(v[i], v[i], d);
    }

    return d;
  }

  /**
   * Calculates the dot product of a vector with a shifted version of itself: dot(v[off1], v[off2])
   *
   * @param n the number of elements to be used
   * @param v the vector to be used
   * @param off1 first offset to use
   * @param off2 second offset to use
   * @return the shifted self dot product
   */
  public static final double shiftedSelfDotProduct(
      final int n, final double[] v, final int off1, final int off2) {

    assert (n + off1 <= v.length);
    assert (n + off2 <= v.length);

    if (n <= 0) return 0.0;
    final int upperBound = SPECIES.loopBound(n);

    var dVec = DoubleVector.zero(SPECIES);
    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vVec1 = DoubleVector.fromArray(SPECIES, v, off1 + i);
      final var vVec2 = DoubleVector.fromArray(SPECIES, v, off2 + i);
      dVec = vVec1.fma(vVec2, dVec);
    }
    double d = dVec.reduceLanes(VectorOperators.ADD);

    for (; i < n; i++) {
      d = Math.fma(v[off1 + i], v[off2 + i], d);
    }

    return d;
  }

  /**
   * Calculates the axpy inside a vector: v = alpha*x[offset] + v
   *
   * @param n the number of elements to scale
   * @param alpha the scaling factor
   * @param v the vector to scale (changed on exit)
   * @param offset the offset for the
   */
  public static void selfAxpy(final int n, final double alpha, final double[] v, final int offset) {

    assert (n + offset <= v.length);
    if (n <= 0 || alpha == 0) return;

    final int upperBound = SPECIES.loopBound(n);

    var alphaVec = DoubleVector.broadcast(SPECIES, alpha);

    int i = 0;
    for (; i < upperBound; i += SPECIES.length()) {
      final var vVecOff = DoubleVector.fromArray(SPECIES, v, offset + i);
      final var vVec = DoubleVector.fromArray(SPECIES, v, i);
      final var vRes = vVecOff.fma(alphaVec, vVec);
      vRes.intoArray(v, i);
    }

    for (; i < n; i++) {
      v[i] = Math.fma(alpha, v[offset + i], v[i]);
    }
  }

  /**
   * Calculates the cross product of two arrays with length 3.
   *
   * @param vectorOne array
   * @param vectorTwo array
   * @param result array
   */
  public static void crossProduct(
      final double[] vectorOne, final double[] vectorTwo, final double[] result) {

    assert (vectorOne.length == 3);
    assert (vectorTwo.length == 3);
    assert (result.length == 3);
    result[0] = vectorOne[1] * vectorTwo[2] - vectorOne[2] * vectorTwo[1];
    result[1] = vectorOne[2] * vectorTwo[0] - vectorOne[0] * vectorTwo[2];
    result[2] = vectorOne[0] * vectorTwo[1] - vectorOne[1] * vectorTwo[0];
  }

  /*
   * Computes the eigenvalues and vectors of a symmetric matrix using netlib.
   * @param matrix the symmetric matrix. only the UPPER part is required.
   * @param eVals the eigenvalues, changed on return.
   * @param eVecs the eigenvectors, changes on return.
   * @return an integer error code. if return code = 0: OK, otherwise there was a problem
   * TOTALLY UNTESTED AND PROBABLY NOT WORKING DUE TO INCOMPATIBILITIES!

  public static int symEigenValue(final double[][] matrix, final double[] eVals, final double[][] eVecs){

      assert(matrix != null);
      assert(matrix.length > 0);
      assert(matrix.length == matrix[0].length);

      final double[] mat = new double[matrix.length*matrix.length];
      for(int i = 0; i < matrix.length; i++){
          System.arraycopy(matrix[i], 0, mat, i*matrix.length, matrix[i].length);
      }

      final double[] work = new double[3*matrix.length-1];
      org.netlib.util.intW out = null;
      org.netlib.lapack.Dsyev.dsyev("V", "U", matrix.length, mat, -1, matrix.length,
              eVals, -1, work, -1, 3*matrix.length-1, out);

      if(out.val != 0){
          // some problem in the process
          System.out.println("WARNING: Symmetric eigenvalue decomposition failed. " + out.val);
          return out.val;
      }

      // copy over the eigenvalues
      for(int i = 0; i < matrix.length; i++){
          System.arraycopy(mat, i*matrix.length, eVecs[i], 0, matrix[i].length);
      }

      return 0;
  }*/

  /*
  public static int eigenValueDecomp(final double[][] matrix){

      assert(matrix != null);
      assert(matrix.length > 0);
      assert(matrix.length == matrix[0].length);

      final double[] mat = new double[matrix.length*matrix.length];
      for(int i = 0; i < matrix.length; i++){
          System.arraycopy(matrix[i], 0, mat, i*matrix.length, matrix[i].length);
      }

      final double[] evReal = new double[matrix.length];
      final double[] evImag = new double[matrix.length];
      Integer out;
      org.netlib.lapack.Dgeev.dgeev("N", "V", matrix.length, mat, -1, matrix.length,
              evReal, -1, evImag, -1, null, -1, dd, doubles4, i7, i8, doubles5, i9, i10, out);
  }*/
}
