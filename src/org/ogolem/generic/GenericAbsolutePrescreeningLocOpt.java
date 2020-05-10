/**
Copyright (c) 2014, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A prescreening local optimization.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class GenericAbsolutePrescreeningLocOpt <E,T extends ContinuousProblem<E>> extends GenericAbstractLocOpt<E,T> {
    
    private static final long serialVersionUID = (long) 20200429;
    private static final Logger l = LoggerFactory.getLogger(GenericAbsolutePrescreeningLocOpt.class);
    protected final GenericLocOpt<E,T> locopt;
    protected final double maxFitness;

  public GenericAbsolutePrescreeningLocOpt(final GenericLocOpt<E,T> locopt,
            final GenericBackend<E,T> back, final double maxFitness){
        super(back);
        this.locopt = locopt;
        this.maxFitness = maxFitness;
    }
    
 public GenericAbsolutePrescreeningLocOpt(final GenericAbsolutePrescreeningLocOpt<E,T> orig){
        super(orig);
        this.locopt = orig.locopt.copy();
        this.maxFitness = orig.maxFitness;
    }
    
    @Override
    public GenericAbstractLocOpt<E, T> copy() {
        return new GenericAbsolutePrescreeningLocOpt<>(this);
    }

    @Override
    protected T optimize(final T individual) {
        
        // one shot
        @SuppressWarnings("unchecked")
        final double[] genome = back.getActiveCoordinates((T) individual.copy());
        final double e = back.fitness(genome, 0);
        
        l.debug("Fitness of one shot is " + e + " compared to " + maxFitness);
        
        if(e > maxFitness){
            l.debug("Discarding " + individual.getID() + " due to prescreeing.");
            individual.setFitness(e);
            return individual;
        }
        
        return locopt.fitness(individual, false);
    }

    @Override
    public String getMyID() {
        return "fitness prescreening locopt using:\n\t" + locopt.getMyID()
                + "\n\t with max fitness " + maxFitness;
    }
}
