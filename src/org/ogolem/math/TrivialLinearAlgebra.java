/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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
 * @version 2014-09-06
 */
public class TrivialLinearAlgebra {
    
    /**
     * Multiplications of two matrices where the resulting matrices is handed
     * over to the method.
     * @param matrixA
     * @param matrixB
     * @param res 
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
     * @param res 
     * @param m 
     * @param n 
     * @param s 
     */
    public static void matMult(final double[][] matrixA, final double[][] matrixB,
            final double[][] res, final int m, final int n, final int s){
    
        final double[] cache = new double[n]; // this is obviously responsible for some mem footprint
        matMult(matrixA, matrixB, res, m, n, s, cache);
    }
        
    /**
     * Multiplications of two matrices where the resulting matrices is handed
     * over to the method.
     * @param matrixA
     * @param matrixB
     * @param res 
     * @param m 
     * @param n 
     * @param s 
     * @param cache
     */
    public static void matMult(final double[][] matrixA, final double[][] matrixB,
            final double[][] res, final int m, final int n, final int s,
            final double[] cache){
        
        assert(res.length >= m);
        assert(res[0].length >= s);
        assert(matrixA.length >= m);
        assert(matrixA[0].length >= n);
        assert(matrixB.length >= n);
        assert(matrixB[0].length >= s);
        assert(cache.length >= n);

        for (int j = 0; j < s; j++) {
            for (int k = 0; k < n; k++) {cache[k] = matrixB[k][j];}
            for (int i = 0; i < m; i++) {

		final double[] rowA = matrixA[i];
		
                double d = 0.0;
                for (int k = 0; k < n; k++) {d += rowA[k] * cache[k];}
                
		res[i][j] = d;
            }
        }
    }
    
    public static void matMult(final double[][] matrixA, final double[][] matrixB,
            final double[][] res, final int m, final int n, final int s,
            final int offAX, final int offAY, final int offBX, final int offBY,
            final double[] cache){
        
        assert(res.length >= m);
        assert(res[0].length >= s);
        assert(matrixA.length >= m);
        assert(matrixA[0].length >= n);
        assert(matrixB.length >= n);
        assert(matrixB[0].length >= s);
        assert(cache.length >= m);

        for (int j = 0; j < s; j++) {
            for (int k = 0; k < n; k++) cache[k] = matrixB[k+offBX][j+offBY];
            for (int i = 0; i < m; i++) {

		final double[] rowA = matrixA[i+offAX]; 
		
                double d = 0.0;
                for (int k = 0; k < n; k++) d += rowA[k+offAY] * cache[k];
                
		res[i][j] = d;
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
            d += v1[i]*v2[i];
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
