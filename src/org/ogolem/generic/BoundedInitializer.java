/**
Copyright (c) 2013, J. M. Dieterich and B. Hartke
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
package org.ogolem.generic;

import org.ogolem.random.Lottery;

/**
 * A bounded initializer for continuous problems.
 * @author Johannes Dieterich
 * @version 2020-05-09
 */
public class BoundedInitializer<T extends Optimizable<Double>> implements GenericInitializer<Double,T> {
    
    private static final long serialVersionUID = (long) 20200509;
    private final GenericFitnessFunction<Double,T> fitness;
    private final Lottery r;
    private final double[] upper;
    private final double[] lower;
    
    public BoundedInitializer(final double[] low, final double[] up,
            final GenericFitnessFunction<Double,T> fitness){
        this.upper = up;
        this.lower = low;
        this.r = Lottery.getInstance();
        this.fitness = fitness;
    }
    
    public BoundedInitializer(final BoundedInitializer<T> orig){
        this.lower = orig.lower;
        this.upper = orig.upper;
        this.r = Lottery.getInstance();
        this.fitness = orig.fitness.copy();
    }

    @Override
    public BoundedInitializer<T> copy() {
        return new BoundedInitializer<>(this);
    }

    @Override
    public T initialize(final T ref, final long futureID) {
        
        @SuppressWarnings("unchecked")
        final T init = (T) ref.copy();
        
        final Double[] params = init.getGenomeCopy();
        final int dims = params.length;
        
        for(int i = 0; i < dims; i++){
            final double d = r.nextDouble();
            assert(upper[i] >= lower[i]);
            params[i] = d*(upper[i]-lower[i])+lower[i];
        }
        
        init.setGenome(params);
        
        // evaluate fitness
        final T evaled = fitness.fitness(init,false);
        evaled.setID(futureID);
        
        return evaled;
    }

    @Override
    public T initializeOnly(T ref, long futureID) {

        @SuppressWarnings("unchecked")
        final T init = (T) ref.copy();

        final Double[] params = init.getGenomeCopy();
        final int dims = params.length;

        for(int i = 0; i < dims; i++){
            final double d = r.nextDouble();
            assert(upper[i] >= lower[i]);
            params[i] = d*(upper[i]-lower[i])+lower[i];
        }

        init.setGenome(params);

        return init;
    }
}
