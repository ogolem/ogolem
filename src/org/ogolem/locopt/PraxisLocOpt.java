/**
Copyright (c) 2014, J. M. Dieterich
              2016-2020, J. M. Dieterich and B. Hartke
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

import numal.AP_praxis_method;
import numal.Analytic_problems;
import org.ogolem.core.FixedValues;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.generic.GenericFitnessBackend;
import org.ogolem.generic.GenericAbstractGradFreeLocOpt;
import org.ogolem.helpers.Machine;
import org.ogolem.locopt.numalhelpers.NumalFitnessWrapper;

/**
 * An interface to NUMAL's gradient-free Praxis local optimization implementation.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class PraxisLocOpt <E,T extends ContinuousProblem<E>> extends GenericAbstractGradFreeLocOpt<E,T> {
    
    private static final long serialVersionUID = (long) 20200429;
    private final double machPrec;
    private final int maxIter;
    private final double convThresh;
    private final double maxStep;
    private final double maxScalingFactor;

    /**
     * Constructor.
     * @param back the backend
     * @param convThresh the convergence threshold on the fitness
     * @param maxStep the maximum step size
     * @param maxScaling the maximum scaling factor
     * @param maxIter the maximal number of iterations
     */
    public PraxisLocOpt(final GenericFitnessBackend<E,T> back, final double convThresh,
            final double maxStep, final double maxScaling, final int maxIter){
        super(back);
        this.machPrec = Machine.calcMachinePrecision();
        this.convThresh = convThresh;
        this.maxIter = maxIter;
        this.maxStep = maxStep;
        this.maxScalingFactor = maxScaling;
    }
    
    public PraxisLocOpt(final PraxisLocOpt<E,T> orig){
        super(orig);
        this.machPrec = orig.machPrec;
        this.maxIter = orig.maxIter;
        this.convThresh = orig.convThresh;
        this.maxStep = orig.maxStep;
        this.maxScalingFactor = orig.maxScalingFactor;
    }
    
    @Override
    public PraxisLocOpt<E,T> copy() {
        return new PraxisLocOpt<>(this);
    }

    @Override
    public String getMyID() {
        return "PRAXIS LOCAL OPTIMIZATION: \n\t" + back.getMyID();
    }
    
    @Override
    protected T optimize(final T individual) {
        
        @SuppressWarnings("unchecked")
        final T work = (T) individual.copy();
        final int dims = back.numberOfActiveCoordinates(work);
        final double[] start = back.getActiveCoordinates(work).clone();
        
        // set the "method" for numal up
        final AP_praxis_method method = new NumalFitnessWrapper<>(back, start);


        final double[] numalParams = new double[dims + 1];
        System.arraycopy(start, 0, numalParams, 1, dims);

        final double[] in = {
            machPrec,
            10*Math.sqrt(machPrec), // relative tolerance for the stepvector
            convThresh,             // absolute tolerance for the stepvector
            0.0,                    // not changed by numal
            0.0,                    // not changed by numal
            (double) maxIter,       // number of function evals
            maxStep,                // maximum step size
            maxScalingFactor,       // maximum scaling factor
            3.0,                    // terminates if no improvement has been made in X steps (this is our convergence check)
            0.0,                    // we assume the problem to be NOT ill-conditioned
        };

        final double[] out = new double[7];

        // call the numal optimization
        try{
            Analytic_problems.praxis(dims, numalParams, method, in, out);
        } catch(Exception e){
            e.printStackTrace(System.err);

            individual.setFitness(FixedValues.NONCONVERGEDENERGY);
            return individual;
        }
        
        // check for strange circumstances
        if(out[1] == 0.0){
            // everything OK
        } else if(out[1] == 1.0){
            // too many iterations, no real problem
            System.err.println("WARNING: Too many iterations in Praxis. ID is " + work.getID());
        } else if(out[1] == 2.0){
            // condition is bad
            System.err.println("WARNING: Condition is bad. Contact author for a switch. ;-)");
        }

        // put the optimized parameters back in
        System.arraycopy(numalParams, 1, start, 0, dims);

        // put the fitness in
        @SuppressWarnings("unchecked")
        final T res = (T) individual.copy();
        final double fitness = out[2];
        res.setFitness(fitness);
        back.updateActiveCoordinates(res, start);
        
        return res;
    }
    
}
