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
import org.ogolem.random.Lottery;

/**
 * Multiple GLOBAL OPTIMIZATIONS wrapped nicely.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class GenericMultipleGlobOpt<E,T extends Optimizable<E>> implements GenericGlobalOptimization<E,T> {

    private static final long serialVersionUID = (long) 20200429;
    private final Lottery random = Lottery.getInstance();
    private final List<GenericGlobalOptimization<E,T>> globopts;
    private final double[] probabilities;
    
    public GenericMultipleGlobOpt(final List<GenericGlobalOptimization<E,T>> opts,
            final List<Double> probs){
        
        assert(opts.size() == probs.size());
        this.globopts = opts;
        this.probabilities = new double[opts.size()];
        double currProb = 0.0;
        int i = 0;
        for(final double prob : probs){
            currProb += prob;
            probabilities[i] = currProb;
            i++;
        }
        
        if(Math.abs(1.0-currProb) >= 1E-3){
            throw new RuntimeException("Probabilites for global optimizations do not add up to 1.0 (100%) within a reasonable criterion. Is: " + currProb);
        }
    }
    
    public GenericMultipleGlobOpt(final GenericMultipleGlobOpt<E,T> orig){
        this.probabilities = orig.probabilities.clone();
        this.globopts = new ArrayList<>(orig.globopts.size());
        orig.globopts.forEach((opt) -> {
            this.globopts.add(opt.copy());
        });
    }
    
    @Override
    public GenericGlobalOptimization<E, T> copy() {
        return new GenericMultipleGlobOpt<>(this);
    }

    @Override
    public String getMyID() {
        
        String s = "MULTIPLE GLOBAL OPTIMIZATIONS WRAPPER\n\n";
        int i = 0;
        for(final GenericGlobalOptimization<E,T> opt : globopts){
            final double prev = (i==0) ? 0.0 : probabilities[i-1];
            s += (probabilities[i]-prev)*100 + " % \t" + opt.getMyID() + "\n\n";
            i++;
        }
        
        return s;
    }

    @Override
    public T globalOptimization(final long futureID, final T mother, final T father) {
        
        assert(probabilities.length == globopts.size());
        final double r = random.nextDouble();
        for(int i = 0; i < probabilities.length; i++){
            if(r <= probabilities[i]){
                return globopts.get(i).globalOptimization(futureID, mother, father);
            }
        }
        
        // numerical inaccuracies and stuff;
        return globopts.get(globopts.size()-1).globalOptimization(futureID, mother, father);
    }
    
}
