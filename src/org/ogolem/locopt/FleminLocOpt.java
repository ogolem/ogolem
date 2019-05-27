/**
Copyright (c) 2013-2014, J. M. Dieterich
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

import numal.AP_linemin_method;
import numal.Analytic_problems;
import org.ogolem.core.FixedValues;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.generic.GenericAbstractLocOpt;
import org.ogolem.generic.GenericBackend;
import org.ogolem.helpers.Machine;
import org.ogolem.locopt.numalhelpers.NumalGradientWrapper;

/**
 * An adapter to Flemin from the NUMAL package.
 * @author Johannes Dieterich
 * @version 2014-12-15
 */
public class FleminLocOpt<E,T extends ContinuousProblem<E>> extends GenericAbstractLocOpt<E,T> {
    
    private static final long serialVersionUID = (long) 20141215;
    private final double machPrec;
    private final int maxIter;
    private final double convThresh;
    private final double convThreshGrad;

    /**
     * Constructor.
     * @param back the backend
     * @param convThresh the convergence threshold on the fitness
     * @param convThreshGrad the convergence threshold on the gradient
     * @param maxIter the maximal number of iterations
     */
    public FleminLocOpt(final GenericBackend<E,T> back, final double convThresh,
            final double convThreshGrad, final int maxIter){
        super(back);
        this.machPrec = Machine.calcMachinePrecision();
        this.convThresh = convThresh;
        this.convThreshGrad = convThreshGrad;
        this.maxIter = maxIter;
    }
    
    public FleminLocOpt(final FleminLocOpt<E,T> orig){
        super(orig);
        this.machPrec = orig.machPrec;
        this.maxIter = orig.maxIter;
        this.convThresh = orig.convThresh;
        this.convThreshGrad = orig.convThreshGrad;
    }
    
    @Override
    public FleminLocOpt<E,T> clone() {
        return new FleminLocOpt<>(this);
    }

    @Override
    public String getMyID() {
        return "FLEMIN LOCAL OPTIMIZATION: \n\t" + back.getMyID();
    }
    
    @Override
    protected T optimize(final T individual) {
        
        @SuppressWarnings("unchecked")
        final T work = (T) individual.clone();
        final int dims = back.numberOfActiveCoordinates(work);
        final double[] start = back.getActiveCoordinates(work).clone();
        
        // set the "method" for numal up
        final AP_linemin_method method = new NumalGradientWrapper<>(back, start);


        final double[] numalParams = new double[dims + 1];
        System.arraycopy(start, 0, numalParams, 1, dims);

        final double[] in = {
            machPrec,
            convThresh/dims,              // relative tolerance
            convThresh,                   // absolute tolerance
            0.0001,                       // line minimization control
            convThreshGrad,               // absolute tolerance for gradient
            -Double.MAX_VALUE+1,          // a lower bound for the function value
                                          // we know that -Double.MAX-1 is the absolute lowest the fitness can get
            1,                            // initialization of the metric
            (double) maxIter,             // number of function evals

        };

        final double[] out = new double[5];
        final double[] gradOut = new double[numalParams.length];
        final double[] hessianIn = new double[numalParams.length*numalParams.length];

        // call the numal optimization
        double fitness;
        try{
            fitness = Analytic_problems.flemin(dims, numalParams,
                    gradOut, hessianIn, method, in, out);
        } catch(Exception e){
            e.printStackTrace(System.err);

            individual.setFitness(FixedValues.NONCONVERGEDENERGY);
            return individual;
        }

        // put the optimized parameters back in
        System.arraycopy(numalParams, 1, start, 0, dims);
        
        @SuppressWarnings("unchecked")
        final T res = (T) individual.clone();
        res.setFitness(fitness);
        back.updateActiveCoordinates(res, start);
        
        return res;
    }
}
