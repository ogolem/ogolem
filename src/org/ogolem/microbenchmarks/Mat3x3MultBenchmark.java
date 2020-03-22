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
import org.ogolem.math.Matrix3x3;

/**
 * Benchmark our internal matrix multiplication routine specialized for 3x3
 * matrices.
 * @author Johannes Dieterich
 * @version 2020-03-01
 */
class Mat3x3MultBenchmark implements SingleMicroBenchmark {

    private final int n;
    
    private final Matrix3x3 matA;
    private final double[][] matB;
    private final double[][] matC;
    
    Mat3x3MultBenchmark(final int n){
        assert(n > 0);
        this.n = n;
        
        final Random r = new Random(42);
        this.matA = new Matrix3x3(r.nextDouble(), r.nextDouble(), r.nextDouble(),
        r.nextDouble(), r.nextDouble(), r.nextDouble(), r.nextDouble(), r.nextDouble(),
        r.nextDouble());
        this.matB = new double[3][n];
        this.matC = new double[3][n];
        
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < n; j++){
                matB[i][j] = r.nextDouble();
            }
        }
    }
    
    @Override
    public double runSingle() throws Exception {
        org.ogolem.math.TrivialLinearAlgebra.matMult(matA, matB, n, matC);
        return matC[0][0];
    }

    @Override
    public String name() {
        return "matrix 3x3 multiplication benchmark for sizes 3 / " + n + " / 3";
    }
}
