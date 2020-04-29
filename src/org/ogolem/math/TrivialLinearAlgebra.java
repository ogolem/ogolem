/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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
 * Performs some really trivial linear algebra.
 * @author Johannes Dieterich
 * @version 2020-04-27
 */
public class TrivialLinearAlgebra {
    
    /**
     * Multiplications of two matrices where the resulting matrices is handed
     * over to the method.
     * @param matrixA
     * @param matrixB
     * @param res overridden on exit
     */
    public static void matMult(final double[][] matrixA, final double[][] matrixB,
            final double[][] res){
        
        final int m = matrixA.length;
        final int n = matrixA[0].length;
        final int s = matrixB[0].length;

        matMult(matrixA, matrixB, res, m, n, s);
    }
    
    /**
     * Multiplications of two matrices where the resulting matrices is handed
     * over to the method.
     * @param matrixA
     * @param matrixB
     * @param res overridden on exit
     * @param m 
     * @param n 
     * @param s 
     */
    public static void matMult(final double[][] matrixA, final double[][] matrixB,
            final double[][] res, final int m, final int n, final int s){
    
        assert(res.length >= m);
        assert(res[0].length >= s);
        assert(matrixA.length >= m);
        assert(matrixA[0].length >= n);
        assert(matrixB.length >= n);
        assert(matrixB[0].length >= s);

        // in BLAS notation: M = m, N = s, K = n
        // we hence want to loop in order M=m (index i), K=n (index l), N=s (index j)
        
        for (int i = 0; i < m; i++) {
            
            final double[] rowRes = res[i];
            // zero the result
            for (int j = 0; j < s; j++) {
                rowRes[j] = 0.0;
            }
            
            for (int l = 0; l < n; l++) {
                final double scalA = matrixA[i][l];
                final double[] rowB = matrixB[l];
                for (int j = 0; j < s; j++) {
                    rowRes[j] += scalA * rowB[j];
                }
            }
        }
    }
        
    public static void matMult(final double[][] matrixA, final double[][] matrixB,
            final double[][] res, final int m, final int n, final int s,
            final int offAX, final int offAY, final int offBX, final int offBY){
        
        assert(res.length >= m);
        assert(res[0].length >= s);
        assert(matrixA.length >= m);
        assert(matrixA[0].length >= n);
        assert(matrixB.length >= n);
        assert(matrixB[0].length >= s);

        for (int i = 0; i < m; i++) {
            
            final double[] rowRes = res[i];
            // zero the result
            for (int j = 0; j < s; j++) {
                rowRes[j] = 0.0;
            }
            
            for (int l = 0; l < n; l++) {
                final double scalA = matrixA[i+offAX][l+offAY];
                final double[] rowB = matrixB[l+offBX];
                for (int j = 0; j < s; j++) {
                    rowRes[j] += scalA * rowB[j+offBY];
                }
            }
        }
    }
    
    /**
     * This multiplies a matrix with another matrix.
     * @param matrixA the first matrix.
     * @param matrixB The second Matrix.
     * @return The result.
     */
    public static double[][] matMult(final double[][] matrixA, final double[][] matrixB){
        
        final int m = matrixA.length;
        final int s = matrixB[0].length;
        
        final double[][] matrixC = new double[m][s];
        matMult(matrixA, matrixB, matrixC, m, matrixA[0].length, s);
        
        return matrixC;
    }
    
    public static void matMult(final Matrix3x3 matrixA, final double[][] matrixB, final int s, final double[][] res){
        
        assert(matrixA != null);
        assert(matrixB != null);
        assert(res != null);
        assert(s >= 0);
        assert(matrixB.length == 3);
        assert(res.length == 3);
        assert(matrixB[0].length >= s);
        assert(matrixB[1].length >= s);
        assert(matrixB[2].length >= s);
        assert(res[0].length >= s);
        assert(res[1].length >= s);
        assert(res[2].length >= s);
        
        final double[] res0 = res[0];
        final double[] res1 = res[1];
        final double[] res2 = res[2];
        
        for(int j = 0; j < s; j++){
            
            final double b0j = matrixB[0][j];
            final double res00j = matrixA.a00 * b0j;
            final double res10j = matrixA.a10 * b0j;
            final double res20j = matrixA.a20 * b0j;
            
            res0[j] = res00j;
            res1[j] = res10j;
            res2[j] = res20j;
        }
        
        for(int j = 0; j < s; j++){
            
            final double b1j = matrixB[1][j];
            final double res01j = matrixA.a01 * b1j;
            final double res11j = matrixA.a11 * b1j;
            final double res21j = matrixA.a21 * b1j;
            
            res0[j] += res01j;
            res1[j] += res11j;
            res2[j] += res21j;
        }
        
        for(int j = 0; j < s; j++){
            
            final double b2j = matrixB[2][j];
            final double res02j = matrixA.a02 * b2j;
            final double res12j = matrixA.a12 * b2j;
            final double res22j = matrixA.a22 * b2j;
            
            res0[j] += res02j;
            res1[j] += res12j;
            res2[j] += res22j;
        }
    }
    
    /**
     * Calculates the dot product of two vectors.
     * @param v1 array
     * @param v2 array
     * @return the dot product
     */
    public static double dotProduct(final double[] v1, final double[] v2){
        
        assert(v1.length == v2.length);
        double d = 0.0;
        for(int i = 0; i < v1.length; i++){
        	d = Math.fma(v1[i], v2[i], d);
        }

        return d;
    }
    
    /**
     * Calculates the dot product of a vector with itself.
     * @param n the number of elements to be used
     * @param v the vector to be used
     * @return the self dot product
     */
    public static double selfDotProduct(final int n, final double[] v) {
    	
    	assert(n <= v.length);
    	double d = 0.0;
        for(int i = 0; i < n; i++){
        	d = Math.fma(v[i], v[i], d);
        }

        return d;
    }    
    
    /**
     * Calculates the cross product of two arrays with length 3.
     * @param vectorOne array
     * @param vectorTwo array
     * @param result array
     */
    public static void crossProduct(final double[] vectorOne, final double[] vectorTwo, final double[] result) {
        
        assert(vectorOne.length == 3);
        assert(vectorTwo.length == 3);
        assert(result.length == 3);
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
