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
package org.ogolem.locopt;

import contrib.newuoa.NEWUOAMethod;
import contrib.newuoa.NEWUOAOptimizer;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.generic.GenericAbstractGradFreeLocOpt;
import org.ogolem.generic.GenericFitnessBackend;
import org.ogolem.generic.Optimizable;
import org.ogolem.helpers.Tuple;

/**
 * A generic interface to the NEWUOA local optimization.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class NEWUOALocOpt<E,T extends ContinuousProblem<E>> extends GenericAbstractGradFreeLocOpt<E,T> {

    private static final long serialVersionUID = (long) 20200429;
    private final NEWUOAOptimizer opt;

    /**
     * Constructor
     * @param back the backend
     * @param rhoStart the start size of the trust radius
     * @param rhoEnd the converged size of the trust radius
     * @param maxIter the maximal number of iterations
     */
    public NEWUOALocOpt(final GenericFitnessBackend<E,T> back,
            final double rhoStart, final double rhoEnd, final int maxIter){
        super(back);
        this.opt = new NEWUOAOptimizer(rhoStart, rhoEnd, maxIter);
    }
    
    public NEWUOALocOpt(final NEWUOALocOpt<E,T> orig){
        super(orig);
        this.opt = new NEWUOAOptimizer(orig.opt);
    }
    
    @Override
    public NEWUOALocOpt<E,T> copy() {
        return new NEWUOALocOpt<>(this);
    }
    
    @Override
    public String getMyID() {
        return "NEWUOA LOCAL OPTIMIZATION using\n\t" + back.getMyID();
    }


    @Override
    protected T optimize(final T individual) {
        
        // first build the Adapter
        final int dims = back.numberOfActiveCoordinates(individual);
        @SuppressWarnings("unchecked")
        final T work = (T) individual.copy();
        final double[] guess = back.getActiveCoordinates(work);
        final NEWUOALocOpt.Adapter<E,T> adap = new NEWUOALocOpt.Adapter<>(back);
        
        // optimize
        final Tuple<Double,double[]> t = opt.doOptimize(dims,guess, adap);        

        @SuppressWarnings("unchecked")
        final T res = (T) individual.copy();
        res.setFitness(t.getObject1());
        back.updateActiveCoordinates(res, guess);
        
        return res;
    }

    private static class Adapter<E,T extends Optimizable<E>> implements NEWUOAMethod {
        
        private final GenericFitnessBackend<E,T> back;
        private int iter = 0;
        
        /**
         * Constructor
         * @param back backend, must be pre-seeded.
         */
        Adapter(final GenericFitnessBackend<E,T> back){
            this.back = back;
        }
        
        @Override
        public double computeObjectiveValue(final int n, final double[] point){
            iter++;
            return back.fitness(point, iter);
        }
    }
}
