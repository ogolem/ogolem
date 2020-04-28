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
import org.ogolem.helpers.Tuple;
import org.ogolem.random.Lottery;

/**
 * Wraps multiple X-overs and assings execution probabilities to them.
 * @author Johannes Dieterich
 * @version 2014-03-28
 */
public class MultipleXOverWrapper<E,T extends Optimizable<E>> implements GenericCrossover<E,T>{

    private static final long serialVersionUID = (long) 20140328;
    private final Lottery random = Lottery.getInstance();
    private final List<GenericCrossover<E,T>> xovers;
    private final double[] probabilities;
    private int lastXOver;
    
    public MultipleXOverWrapper(final List<GenericCrossover<E,T>> xovers,
            final List<Double> probs){
        
        assert(xovers.size() == probs.size());
        this.xovers = xovers;
        this.probabilities = new double[xovers.size()];
        double currProb = 0.0;
        int i = 0;
        for(final double prob : probs){
            currProb += prob;
            probabilities[i] = currProb;
            i++;
        }
        
        if(Math.abs(1.0-currProb) >= 1E-3){
            throw new RuntimeException("Probabilites for XOvers do not add up to 1.0 (100%) within a reasonable criterion.");
        }
    }
    
    public MultipleXOverWrapper(final MultipleXOverWrapper<E,T> orig){
        this.probabilities = orig.probabilities.clone();
        this.xovers = new ArrayList<>(orig.xovers.size());
        orig.xovers.forEach((mut) -> {
            this.xovers.add(mut.clone());
        });
    }
    
    @Override
    public MultipleXOverWrapper<E, T> clone() {
        return new MultipleXOverWrapper<>(this);
    }

    @Override
    public String getMyID() {
        String s = "MULTIPLE XOVER WRAPPER\n";
        int i = 0;
        for(final GenericCrossover<E,T> xover : xovers){
            final double prev = (i==0) ? 0.0 : probabilities[i-1];
            s += (probabilities[i]-prev)*100 + "% of \t" + xover.getMyID() + "\n\t";
            i++;
        }
        
        return s;
    }

    @Override
    public Tuple<T, T> crossover(final T mother, final T father, final long futureID) {
        
        assert(probabilities.length == xovers.size());
        final double r = random.nextDouble();
        for(int i = 0; i < probabilities.length; i++){
            if(r <= probabilities[i]){
                lastXOver = i;
                return xovers.get(i).crossover(mother, father, futureID);
            }
        }
        
        // numerical inaccuracies and stuff;
        return xovers.get(xovers.size()-1).crossover(mother, father, futureID);
    }

    @Override
    public short hasPriority() {
        return xovers.get(lastXOver).hasPriority();
    }
}
