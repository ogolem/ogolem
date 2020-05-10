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
package org.ogolem.generic;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.Lottery;

/**
 * A chain of crossover operators.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class GenericChainedXOver<E,T extends Optimizable<E>> implements GenericCrossover<E,T> {
    
    private static final long serialVersionUID = (long) 20200429;
    private final List<GenericCrossover<E,T>> xovers;
    private final List<Double> probs;
    private final Lottery r;
    
    public GenericChainedXOver(final List<GenericCrossover<E,T>> xovers,
            final List<Double> probs){
        assert(xovers.size() == probs.size());
        this.xovers = xovers;
        this.probs = probs;
        this.r = Lottery.getInstance();
    }
    
    public GenericChainedXOver(final GenericChainedXOver<E,T> orig){
        this.xovers = new ArrayList<>(orig.xovers.size());
        orig.xovers.forEach((xover) -> {
            xovers.add(xover.clone());
        });
        this.probs = new ArrayList<>(orig.probs.size());
        for(final double d : orig.probs){
            assert(d >= 0.0 && d <= 1.0);
            probs.add(d);
        }
        assert(xovers.size() == probs.size());
        this.r = Lottery.getInstance();
    }
    
    @Override
    public GenericCrossover<E, T> clone() {
        return new GenericChainedXOver<>(this);
    }

    @Override
    public String getMyID() {
        String s = "CHAINED CROSSOVER\n";
        for(final GenericCrossover<E,T> xover : xovers){
            s += xover.getMyID() + "\n";
        }
        
        return s;
    }

    @SuppressWarnings("unchecked")
    @Override
    public Tuple<T, T> crossover(final T mother, final T father, final long futureID) {
        
        assert(xovers.size() == probs.size());
        Tuple<T,T> children = new Tuple<>((T) mother.copy(), (T)father.copy());
        for(int  i = 0; i < xovers.size(); i++){
            if(r.nextDouble() < probs.get(i)){
                children = xovers.get(i).crossover(children.getObject1(), children.getObject2(), futureID);
            }
        }
        
        return children;
    }

    @Override
    public short hasPriority() {
        return -1;
    }
}
