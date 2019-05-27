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
package org.ogolem.locopt;

import org.ogolem.generic.ContinuousProblem;

/**
 * A basic optimizable type. Basically a wrapper for a double[].
 * @author Johannes Dieterich
 * @version 2014-12-15
 */
public class BasicOptimizableType extends ContinuousProblem<Double>{
    
    private static final long serialVersionUID = (long) 20141215;
    private double[] genome;

    BasicOptimizableType(final double[] genome){
        this.genome = genome;
    }
    
    @Override
    public ContinuousProblem<Double> clone() {
        return new BasicOptimizableType(genome.clone());
    }

    @Override
    public double[] getGenomeAsDouble() {
        return genome;
    }

    @Override
    public Double[] getGenomeCopy() {
        final Double[] gen = new Double[genome.length];
        for(int i = 0; i < genome.length; i++){
            gen[i] = genome[i];
        }
        
        return gen;
    }

    @Override
    public void setGenome(Double[] geno) {
        for(int i = 0; i < genome.length; i++){
            genome[i] = geno[i];
        }
    }
    
    void setGenome(final double[] geno){
        for(int i = 0; i < genome.length; i++){
            genome[i] = geno[i];
        }
    }
    
    int getNoOfCoords(){
        return genome.length;
    }
}
