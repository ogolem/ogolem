/**
Copyright (c) 2014, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive;

import java.util.Random;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericInitializer;

/**
 * An adapter from the generic to ParameterInit. One day remove..
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class ParameterInit implements GenericInitializer<Double,AdaptiveParameters>{
    
    private static final long serialVersionUID = (long) 20200429;
    private final Random random = new Random();
    private final GenericFitnessFunction<Double,AdaptiveParameters> fitness;
    private final double[] lower;
    private final double[] upper;
    
    ParameterInit(final GenericFitnessFunction<Double,AdaptiveParameters> fitness, final double[] minBorders,
            final double[] maxBorders){
        this.upper = maxBorders;
        this.lower = minBorders;
        this.fitness = fitness;
    }
    
    ParameterInit(final ParameterInit orig){
        this.upper = orig.upper.clone();
        this.lower = orig.lower.clone();
        this.fitness = orig.fitness.copy();
    }

    @Override
    public GenericInitializer<Double, AdaptiveParameters> copy() {
        return new ParameterInit(this);
    }

    @Override
    public AdaptiveParameters initialize(final AdaptiveParameters ref, final long futureID) {
        
        final AdaptiveParameters copy = new AdaptiveParameters(ref);
        copy.setID(futureID);
        
        final double[] params = copy.getAllParamters();
        for(int i = 0; i < params.length; i++){
            final double d = random.nextDouble();
            assert(upper[i] >= lower[i]);
            params[i] = d*(upper[i]-lower[i])+lower[i];
        }
        
        // evaluate fitness
        final AdaptiveParameters evaled = fitness.fitness(copy,false);
        evaled.setID(futureID);
        
        return evaled;
    }    
}
