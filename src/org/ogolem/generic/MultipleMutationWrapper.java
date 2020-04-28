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
package org.ogolem.generic;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.random.Lottery;

/**
 * Wraps multiple mutations and assings execution probabilities to them.
 * @author Johannes Dieterich
 * @version 2014-03-28
 */
public class MultipleMutationWrapper<E,T extends Optimizable<E>> implements GenericMutation<E,T>{

    private static final long serialVersionUID = (long) 20140328;
    private final Lottery random = Lottery.getInstance();
    private final List<GenericMutation<E,T>> mutations;
    private final double[] probabilities;
    
    public MultipleMutationWrapper(final List<GenericMutation<E,T>> muts,
            final List<Double> probs){
        
        assert(muts.size() == probs.size());
        this.mutations = muts;
        this.probabilities = new double[muts.size()];
        double currProb = 0.0;
        int i = 0;
        for(final double prob : probs){
            currProb += prob;
            probabilities[i] = currProb;
            i++;
        }
        
        if(Math.abs(1.0-currProb) >= 1E-3){
            throw new RuntimeException("Probabilites for Mutations do not add up to 1.0 (100%) within a reasonable criterion.");
        }
    }
    
    public MultipleMutationWrapper(final MultipleMutationWrapper<E,T> orig){
        this.probabilities = orig.probabilities.clone();
        this.mutations = new ArrayList<>(orig.mutations.size());
        orig.mutations.forEach((mut) -> {
            this.mutations.add(mut.clone());
        });
    }
    
    @Override
    public MultipleMutationWrapper<E,T> clone() {
        return new MultipleMutationWrapper<>(this);
    }

    @Override
    public String getMyID() {
        String s = "MULTIPLE MUTATION WRAPPER\n";
        int i = 0;
        for(final GenericMutation<E,T> mutation : mutations){
            final double prev = (i==0) ? 0.0 : probabilities[i-1];
            s += (probabilities[i]-prev)*100 + "% of \t" + mutation.getMyID() + "\n\t";
            i++;
        }
        
        return s;
    }

    @Override
    public T mutate(final T individual) {
        
        assert(probabilities.length == mutations.size());
        final double r = random.nextDouble();
        for(int i = 0; i < probabilities.length; i++){
            if(r <= probabilities[i]){
                return mutations.get(i).mutate(individual);
            }
        }
        
        // numerical inaccuracies and stuff;
        return mutations.get(mutations.size()-1).mutate(individual);
    }
    
}
