/**
Copyright (c) 2013-2015, J. M. Dieterich
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

import org.ogolem.generic.stats.GenericDetailStatistics;

/**
 * An abstract implementation for a local optimization for continuous problems.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public abstract class GenericAbstractLocOpt<E,T extends ContinuousProblem<E>> implements GenericLocOpt<E,T>{
    
    private static final long serialVersionUID = (long) 20200429;
    protected final GenericBackend<E,T> back;
    
    public GenericAbstractLocOpt(final GenericBackend<E,T> back){
        this.back = back;
    }
    
    protected GenericAbstractLocOpt(final GenericAbstractLocOpt<E,T> orig){
        this.back = orig.back.copy();
    }
    
    @Override
    public abstract GenericAbstractLocOpt<E,T> copy();
    
    @Override
    public T fitness(final T individual, final boolean forceOneEval){     
        if(forceOneEval){
            GenericDetailStatistics.incrementFitnessEvals();
            return back.fitness(individual, forceOneEval);
        } else{
            GenericDetailStatistics.incrementLocalOpts();
            return optimize(individual);
        }
    }
    
    protected abstract T optimize(final T individual);
    
    @Override
    public GenericFitnessBackend<E,T> getFitnessBackend(){
        return back;
    }
    
    @Override
    public GenericBackend<E,T> getBackend(){
        return back;
    }
}
