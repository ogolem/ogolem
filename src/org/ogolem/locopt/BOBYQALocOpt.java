/**
Copyright (c) 2013-2014, J. M. Dieterich
                2015, J. M. Dieterich and B. Hartke
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

import contrib.bobyqa.AbstractBOBYQAMethod;
import contrib.bobyqa.BOBYQAOptimizer;
import org.ogolem.core.FixedValues;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.generic.GenericAbstractGradFreeLocOpt;
import org.ogolem.generic.GenericFitnessBackend;
import org.ogolem.generic.GenericFitnessBackend.BOUNDSTYPE;
import org.ogolem.generic.Optimizable;
import org.ogolem.helpers.Tuple;

/**
 * A generic BOBYQA interface.
 * @author Johannes Dieterich
 * @version 2014-03-28
 */
public class BOBYQALocOpt<E,T extends ContinuousProblem<E>> extends GenericAbstractGradFreeLocOpt<E,T> {
    
    private static final long serialVersionUID = (long) 20150328;
    private static final boolean DEBUG = false;
    private final int maxIter;
    private final double initialTrustRadius;
    private final double convergedTrustRadius;
    
    /**
     * Constructor
     * @param backend the backend
     * @param maxIter the maximum number of function evaluations
     * @param initialTrustRadius the initial trust radius
     * @param convergedTrustRadius the converged trust radius
     */
    public BOBYQALocOpt(final GenericFitnessBackend<E,T> backend,
            final int maxIter, final double initialTrustRadius,
            final double convergedTrustRadius){
        
        super(backend);
        this.maxIter = maxIter;
        this.initialTrustRadius = initialTrustRadius;
        this.convergedTrustRadius = convergedTrustRadius;
    }
    
    private BOBYQALocOpt(final BOBYQALocOpt<E,T> orig){
        super(orig);
        this.maxIter = orig.maxIter;
        this.initialTrustRadius = orig.initialTrustRadius;
        this.convergedTrustRadius = orig.convergedTrustRadius;
    }
    
    @Override
    public BOBYQALocOpt<E,T> clone() {
        return new BOBYQALocOpt<>(this);
    }

    @Override
    public String getMyID() {
        return "BOBYQA LOCAL OPTIMIZATION using:\n\t" + back.getMyID();
    }

    @Override
    protected T optimize(final T individual) {
        
        @SuppressWarnings("unchecked")
        final T clone = (T) individual.clone();
        
        // first build the Adapter
        final double[] optCoords = back.getActiveCoordinates(clone);
        
        // optimize
        final double[] g = optCoords.clone();
        final double[] normalized = new double[optCoords.length];

        // we don't care
        final BOUNDSTYPE boundsT = back.boundariesInRepresentation(individual);
        if(boundsT != BOUNDSTYPE.ALL){
            // print a warning
            System.err.println("WARNING: BOBYQALocOpt is a fully bound-constraint local optimization. The backend representation however is only bounded to " + boundsT.name() + " degree. Continuing...");
        }
        
        final double[][] borders = new double[2][optCoords.length];
        back.bestEstimateBoundaries(optCoords, borders[0], borders[1]);
        
        final BOBYQAOptimizer opt = new BOBYQAOptimizer(borders, 2*borders[0].length+1, initialTrustRadius, convergedTrustRadius, maxIter, true);
        final Adapter<E,T> adap = new Adapter<>(back, new double[optCoords.length], borders);
        
        // normalize
        adap.normalizeFromBounds(borders, g, normalized);
        
        final Tuple<Double,double[]> t = opt.doOptimize(normalized, adap);        
        
        // denormalize and build new object
        adap.denormalizeToBounds(borders, t.getObject2(), g);
        
        @SuppressWarnings("unchecked")
        final T res = (T) individual.clone();
        res.setFitness(t.getObject1());
        back.updateActiveCoordinates(res, g);
        
        return res;
    }
    
    private class Adapter<E,T extends Optimizable<E>> extends AbstractBOBYQAMethod{
        
        private final GenericFitnessBackend<E,T> backend;
        private final double[] point;
        private final double[][] bounds;
        private double lastFitness = FixedValues.NONCONVERGEDENERGY;
        int iter = 0;
        
        Adapter(final GenericFitnessBackend<E,T> back,
                final double[] point, final double[][] bounds){
            this.backend = back;
            this.point = point;
            this.bounds = bounds;
        }
        
        @Override
        public double computeObjectiveValue(final double[] normalized){
            
            // denormalize
            denormalizeToBounds(bounds, normalized, point);
                                    
            lastFitness = backend.fitness(point,iter);
            
            if(DEBUG){System.out.println("DEBUG: Backend evaluates fitness to " + lastFitness + " for " + iter);}
            iter++;
            
            return lastFitness;
        }
                
        @Override
        public boolean doesNormalize(){
            return true;
        }
    }
}
