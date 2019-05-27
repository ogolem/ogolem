/**
Copyright (c) 2014, J. M. Dieterich
                2015, J. M. Dieterich and B. Hartke
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
package org.ogolem.locopt;

import java.util.Random;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * A test class for BOBYQA, the bound-constraint gradient-free optimizer.
 * @author Johannes Dieterich
 * @version 2015-03-28
 */
public class BOBYQALocOptTest {
    
    /**
     * Test of optimize method, of class PraxisLocOpt.
     */
    @Test
    public void testOptimize() {
        System.out.println("optimize");
        final Random r = new Random();
        final double a = 0.0001;
        final double shift = 4.2;
        final double b = 0.0;
        final double convergedTrust = 1e-6;
        final double initialTrust = 1e-1;
        final int maxIter = 9999;
        final int noTrials = 100;
        final int dim = 42;
        final HarmonicFunction func = new HarmonicFunction(a,shift,b, 0.0, 42*shift);
        
        final double[][] bounds = new double[2][dim];
        for(int x = 0; x < dim; x++){
            bounds[0][x] = 0.0;
            bounds[1][x] = 42.0+shift;
        }
        
        final BOBYQALocOpt<Double,BasicOptimizableType> instance1 = new BOBYQALocOpt<>(func,maxIter,initialTrust,convergedTrust);
        
        for(int trial = 0; trial < noTrials; trial++){
            final double[] gen1 = new double[dim];
        
            // random init
            for(int x = 0; x < dim; x++){
                final double d = 42.0*r.nextDouble()+shift;
                // normalize [0.0; 42+shift]
                gen1[x] = (d)/(42+shift);
            }
                    
            final BasicOptimizableType ind1 = new BasicOptimizableType(gen1);
        
            final BasicOptimizableType res1 = instance1.optimize(ind1);
            
            // check
            assertTrue("Fitness (I) wrong: " + res1.getFitness(), Math.abs(res1.getFitness()) <= convergedTrust);
            //SHIT! Does not work. So does the res1Dat stuff... CRAP!;
            
            final double[] res1Dat = res1.getGenomeAsDouble();
            for(int x = 0; x < dim; x++){
                assertTrue("Error (I) from correct solution: " + Math.abs(res1Dat[x]-shift) + " in " + x, Math.abs(res1Dat[x]-shift) <= 10*convergedTrust); // make it a bit bigger
            }
        }
    }
    
}
