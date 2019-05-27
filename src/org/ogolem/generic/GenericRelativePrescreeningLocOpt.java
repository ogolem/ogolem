/**
Copyright (c) 2013-2014, J. M. Dieterich
              2016, J. M. Dieterich and B. Hartke
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

import org.ogolem.generic.genericpool.GenericPool;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A prescreening version of the local optimization. Based off a maximum percentage of 
 * the current pool fitness spread. Will obviously only work in a) shared memory
 * environments or b) in environments where the pool is being propagated to the clients
 * (regularly).
 * @author Johannes Dieterich
 * @version 2016-01-14
 */
public class GenericRelativePrescreeningLocOpt <E,T extends ContinuousProblem<E>> extends GenericAbstractLocOpt<E,T>{
    
    private static final long serialVersionUID = (long) 20131124;
    private static final Logger l = LoggerFactory.getLogger(GenericRelativePrescreeningLocOpt.class);
    protected final GenericAbstractLocOpt<E,T> locopt;
    protected final GenericPool<E,T> pool;
    protected final double maxPerc;
    
    public GenericRelativePrescreeningLocOpt(final GenericAbstractLocOpt<E,T> locopt,
            final GenericBackend<E,T> back, final GenericPool<E,T> pool,
            final double maxPercentage){
        super(back);
        this.locopt = locopt;
        this.pool = pool;
        this.maxPerc = maxPercentage;
    }
    
    public GenericRelativePrescreeningLocOpt(final GenericRelativePrescreeningLocOpt<E,T> orig){
        super(orig);
        this.locopt = orig.locopt.clone();
        this.maxPerc = orig.maxPerc;
        this.pool = orig.pool; // is a singleton!
    }
    
    @Override
    public GenericRelativePrescreeningLocOpt<E,T> clone(){
        return new GenericRelativePrescreeningLocOpt<>(this);
    }
    
    @Override
    public String getMyID(){
        return "prescreening locopt using:\n\t" + locopt.getMyID();
    }
    
    @Override
    public T optimize(final T individual){
        
        // one shot
        @SuppressWarnings("unchecked")
        final double[] genome = back.getActiveCoordinates((T) individual.clone());
        final double e = back.fitness(genome, 0);
        
        final double currBest = pool.getFitnessOfIndividualAtPos(0);
        final double currLast = pool.getFitnessOfIndividualAtPos(pool.getCurrentPoolSize()-1);
        
        l.debug("Curr best: " + currBest + " curr last: " + currLast + " this: " + e);
        
        final double perc = (e-currBest)*100/(currLast-currBest);
        
        if(perc > maxPerc){
            l.debug("Individual does not survive prescreening: " + perc + " vs " + maxPerc);
            individual.setFitness(e);
            return individual;
        }
        
        l.debug("Individual survives prescreening, will go into optimization now.");
        
        return locopt.optimize(individual);
    }
    
}
