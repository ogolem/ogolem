/**
Copyright (c) 2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.microbenchmarks;

import java.util.Random;

/**
 * Benchmark our internal, trivial matrix multiplication routine.
 * @author Johannes Dieterich
 * @version 2020-03-01
 */
class MatMultBenchmark implements SingleMicroBenchmark {

    private final int m;
    private final int n;
    private final int k;
    
    private final double[][] matA;
    private final double[][] matB;
    private final double[][] matC;
    
    MatMultBenchmark(final int m, final int n, final int k){
        assert(m > 0);
        assert(n > 0);
        assert(k > 0);
        this.m = m;
        this.n = n;
        this.k = k;
        
        final Random r = new Random(42);
        this.matA = new double[m][k];
        this.matB = new double[k][n];
        this.matC = new double[m][n];
        
        for(int i = 0; i < m; i++){
            for(int j = 0; j < k; j++){
                matA[i][j] = r.nextDouble();
            }
        }
        
        for(int i = 0; i < k; i++){
            for(int j = 0; j < n; j++){
                matB[i][j] = r.nextDouble();
            }
        }
    }
    
    @Override
    public double runSingle() throws Exception {
        org.ogolem.math.TrivialLinearAlgebra.matMult(matA, matB, matC);
        return matC[0][0];
    }

    @Override
    public String name() {
        return "matrix multiplication benchmark for sizes " + m + " / " + n + " / " + k;
    }
}
