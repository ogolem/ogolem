/**
Copyright (c) 2014-2015, J. M. Dieterich
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
package org.ogolem.locopt;

import org.ogolem.generic.*;
import java.util.ArrayList;
import java.util.List;

/**
 * A generic chain of fitness functions.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class GenericChainedLocalOpt<E,T extends ContinuousProblem<E>> implements GenericLocOpt<E,T> {
    
    private static final long serialVersionUID = (long) 20200429;

    private final List<GenericLocOpt<E,T>> functions;
    private final double cutoff;
    
    public GenericChainedLocalOpt(final List<GenericLocOpt<E,T>> fitFuncs, final double cutoff) throws Exception {
        if(fitFuncs.isEmpty()){throw new Exception("Wrong input for chained fitness function.");}
        this.functions = fitFuncs;
        this.cutoff = cutoff;
    }
    
    public GenericChainedLocalOpt(final GenericChainedLocalOpt<E,T> orig){
        if(orig.functions.isEmpty()){throw new RuntimeException("Wrong input for chained fitness function.");}
        this.functions = new ArrayList<>();
        orig.functions.forEach((function) -> {
            this.functions.add(function.copy());
        });
        this.cutoff = orig.cutoff;
    }
    
    @Override
    public GenericChainedLocalOpt<E, T> copy() {
        return new GenericChainedLocalOpt<>(this);
    }

    @Override
    public String getMyID() {
        String s = "CHAINED LOCAL OPTIMIZATION:\n\t";
        for (final GenericFitnessFunction<E, T> function : functions) {
            final String sx = function.getMyID();
            s += sx;
        }
        s += "\n\t cutoff: " + cutoff;
        
        return s;
    }

    @Override
    public T fitness(final T individual, final boolean forceOneEval) {
        
        if(forceOneEval){
            // just do the last
            return functions.get(functions.size()-1).fitness(individual, forceOneEval);
        }
        
        T ind = individual;
        for(final GenericLocOpt<E,T> func : functions){
            
            ind = func.fitness(ind, forceOneEval);
            if(ind.getFitness() >= cutoff){
                ind.setFitness(cutoff);
                return ind;
            }
        }
        
        return ind;
    }

    @Override
    public GenericFitnessBackend<E, T> getFitnessBackend() {
        return getBackend();
    }

    @Override
    public GenericBackend<E, T> getBackend() {
        
        // return the LAST backend in the list.
        return functions.get(functions.size()-1).getBackend();
    }
}
