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
 * A chain of mutation operators.
 * @author Johannes Dieterich
 * @version 2014-05-04
 */
public class GenericChainedMutation<E, T extends Optimizable<E>> implements GenericMutation<E,T> {

    private static final long serialVersionUID = (long) 2014054;
    private final List<GenericMutation<E,T>> mutations;
    private final List<Double> probs;
    private final Lottery r;
    
    public GenericChainedMutation(final List<GenericMutation<E,T>> muts,
            final List<Double> probs){
        assert(muts.size() == probs.size());
        this.mutations = muts;
        this.probs = probs;
        this.r = Lottery.getInstance();
    }
    
    public GenericChainedMutation(final GenericChainedMutation<E,T> orig){
        this.mutations = new ArrayList<>(orig.mutations.size());
        orig.mutations.forEach((mut) -> {
            mutations.add(mut.clone());
        });
        this.probs = new ArrayList<>(orig.probs.size());
        for(final double d : orig.probs){
            assert(d >= 0.0 && d <= 1.0);
            probs.add(d);
        }
        assert(mutations.size() == probs.size());
        this.r = Lottery.getInstance();
    }
    
    @Override
    public GenericMutation<E, T> clone() {
        return new GenericChainedMutation<>(this);
    }

    @Override
    public String getMyID() {
        String s = "CHAINED MUTATION\n";
        for(final GenericMutation<E,T> mut : mutations){
            s += mut.getMyID() + "\n";
        }
        
        return s;
    }

    @SuppressWarnings("unchecked")
    @Override
    public T mutate(final T orig) {
        
        assert(mutations.size() == probs.size());
        T mutated = (T) orig.clone();
        for(int i = 0; i < mutations.size(); i++){
            if(r.nextDouble() < probs.get(i)){
                mutated = mutations.get(i).mutate(mutated);
            }
        }
        
        return mutated;
    }
}
