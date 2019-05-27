/**
Copyright (c) 2014, J. M. Dieterich
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
 * A trivial interface for SOME BLAS-like functions (to be extended). May underneath call some
 * other (faster/more stable) implementation (like netlib/jlapack). Also may not be entirely
 * compatible with the original Fortran implementation (especially in edge cases!)
 * @author Johannes Dieterich
 * @version 2014-06-22
 */
public class BLASInterface {
    
    public static enum TRANSPOSE {YES, NO};
    
    private BLASInterface(){}; // no instantiation wanted or needed
    
    public static double ddot(final int n, final double[] dx, final int incx,
            final double[] dy, final int incy){
        
        if(incx == 1 && incy == 1 && n == dx.length){return TrivialLinearAlgebra.dotProduct(dx, dy);}
        
        double res = 0.0;
        int ix = 0;
        int iy = 0;
        for(int i = 0; i < n; i++){
            res += dx[ix]*dy[iy];
            ix += incx;
            iy += incy;
        }
        
        return res;
    }
    
    public static void daxpy(final int n, final double da, final double[] dx, final int incx,
            final double[] dy, final int incy){
        
        if(da == 0.0){return;}
        if(da == 1.0){
            dxpy(n,dx,incx,dy,incy);
            return;
        }
        
        if(incx == 1 && incy == 1){
            for(int i = 0; i < n; i++){
                dy[i] += da*dx[i];
            }
        } else {
            int ix = 0;
            int iy = 0;
            for(int i = 0; i < n; i++){
                dy[iy] += da*dx[ix];
                ix += incx;
                iy += incy;
            }
        }
    }
    
    /**
     * special non-BLAS interface for a daxpy with a = 1.0
     * @param n the number of elements
     * @param dx vector 1
     * @param incx increment for vector 1
     * @param dy vector 2
     * @param incy increment for vector 2
     */
    public static void dxpy(final int n, final double[] dx, final int incx,
            final double[] dy, final int incy){
        
        if(incx == 1 && incy == 1){
            for(int i = 0; i < n; i++){
                dy[i] += dx[i];
            }
        } else {
            int ix = 0;
            int iy = 0;
            for(int i = 0; i < n; i++){
                dy[iy] += dx[ix];
                ix += incx;
                iy += incy;
            }
        }
    }
    
    public static double dnrm2(final int n, final double[] dx, final int incx){
        
        if(n <= 0){return 0.0;}
        
        return org.netlib.blas.Dnrm2.dnrm2(n, dx, 0, incx);
    }
    
    public static void dcopy(final int n, final double[] dx, final int incx,
            final double[] dy, final int incy){
        
        if(incx == 1 && incy == 1){
            System.arraycopy(dx, 0, dy, 0, n);
            return;
        }
        
        org.netlib.blas.Dcopy.dcopy(n, dx, 0, incx, dy, 0, incy);
    }
    
    /*public static void dgemm(final TRANSPOSE fir, final TRANSPOSE sec, final int m, final int n,
            final int k, final double alpha, final double[][] a, final int lda, final double[][] b,
            final int ldb, final double beta, final double[][] c, final int ldc){
        
        final String t1 = (fir == TRANSPOSE.NO) ? "N" : "T";
        final String t2 = (sec == TRANSPOSE.NO) ? "N" : "T";
        final double[] aO = new double[m*n]; // actually not should eb a.length*a[0].length
        final double[] b0 = new double[];
        
        org.netlib.blas.Dgemm.dgemm(t1, t2, m, n, k, alpha, aO, 0, lda, bO, 0, ldb, beta, cO, 0, ldc);
    }*/
}
