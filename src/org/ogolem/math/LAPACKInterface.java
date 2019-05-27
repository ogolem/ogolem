/**
Copyright (c) 2015, J. M. Dieterich and B. Hartke
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

import org.netlib.util.intW;

/**
 * A trivial interface for SOME LAPACK-like functions (to be extended). May underneath call some
 * other (faster/more stable) implementation (like netlib/jlapack). Also may not be entirely
 * compatible with the original Fortran implementation (especially in edge cases!)
 * @author Johannes Dieterich
 * @version 2015-03-14
 */
public class LAPACKInterface {
    
    private LAPACKInterface(){}; // no instantiation wanted or needed
    
    
    public static void invertMatrix (final int dim, final double[] matrix, final double[] scr1,
            final int[] scr2){
        
        assert(scr1.length >= dim*dim);
        assert(matrix.length >= dim*dim);
        assert(scr2.length >= dim);
        
        intW info = new intW(0);
        // dgetrf_(&N,&N,A,&N,IPIV,&INFO);
        org.netlib.lapack.Dgetrf.dgetrf(dim,dim,matrix,0,dim,scr2,0,info);
        // dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
        org.netlib.lapack.Dgetri.dgetri(dim,matrix,0,dim,scr2,0,scr1,0,scr1.length,info);
    }
}
