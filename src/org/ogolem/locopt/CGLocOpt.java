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

import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.NonLinearConjugateGradientOptimizer;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.ogolem.core.FixedValues;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.generic.GenericAbstractLocOpt;
import org.ogolem.generic.GenericBackend;
import org.ogolem.locopt.apachehelpers.MultivariateToEnergyProvider;
import org.ogolem.locopt.apachehelpers.MultivariateToGradientProvider;

/**
 * Interface to apache commons implementation of CG.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class CGLocOpt <E,T extends ContinuousProblem<E>> extends GenericAbstractLocOpt<E,T> {

    private static final long serialVersionUID = (long) 20200429;
    private static final boolean DEBUG = false;
    private final boolean useFletcherReevesUpdate;
    private final double absThresh;
    private final int maxIter;
    
    /**
     * Instantiation of the CG.
     * @param provider the backend provider
     * @param useFletcher whether Fletcher-Reeves or Polak-Ribiere updates should be done
     * @param absThresh absolute convergence threshold
     * @param maxIter maximum number of iterations
     */
    public CGLocOpt(final GenericBackend<E,T> provider, final boolean useFletcher,
            final double absThresh, final int maxIter){
        super(provider);
        this.useFletcherReevesUpdate = useFletcher;
        this.absThresh = absThresh;
        this.maxIter = maxIter;
    }
    
    public CGLocOpt(final CGLocOpt<E,T> orig){
        super(orig);
        this.useFletcherReevesUpdate = orig.useFletcherReevesUpdate;
        this.absThresh = orig.absThresh;
        this.maxIter = orig.maxIter;
    }
    
    @Override
    public CGLocOpt<E,T> copy() {
        return new CGLocOpt<>(this);
    }

    @Override
    public String getMyID() {
        return "conjugate gradient optimizer from APACHE COMMONS using:\n\t " + back.toString();
    }

    @Override
    protected T optimize(final T individual) {

        @SuppressWarnings("unchecked")
        final T work = (T) individual.copy();
        final MultivariateToEnergyProvider<E,T> provE = new MultivariateToEnergyProvider<>(back, work);
        final MultivariateToGradientProvider<E,T> provG = new MultivariateToGradientProvider<>(back, work);
        
        final NonLinearConjugateGradientOptimizer.Formula updater = (useFletcherReevesUpdate) ?  NonLinearConjugateGradientOptimizer.Formula.FLETCHER_REEVES
                : NonLinearConjugateGradientOptimizer.Formula.POLAK_RIBIERE;
        //XXX 
        final double relThreshold = absThresh;
        final ConvergenceChecker<PointValuePair> checker = new SimpleValueChecker(relThreshold,absThresh,maxIter);//SimpleConvChecker(absThresh,maxIter,work.getID());
        
        final NonLinearConjugateGradientOptimizer opter = new NonLinearConjugateGradientOptimizer(updater, checker);
        
        final ObjectiveFunctionGradient oG = new ObjectiveFunctionGradient(provG);
        final ObjectiveFunction oE = new ObjectiveFunction(provE);
        final GoalType goal = GoalType.MINIMIZE;
        final InitialGuess init = new InitialGuess(back.getActiveCoordinates(work));
        final MaxEval maxE = new MaxEval(maxIter);
        
        boolean wasConverged = true;
        PointValuePair result = null;
        try{
            result = opter.optimize(oG, oE, goal, init, maxE);
        } catch (TooManyEvaluationsException e){
            System.out.println("WARNING: Too many function evaluations in Apache's CG. Not converged! ID was " + individual.getID());
            
            // get the last status from within the interface. yes, it is ugly, but that is due to Apache NOT
            // returning it, as it should, IMHO.
            final double currE = provE.bestEnergy();
            final double[] currP = provE.bestPoint();

            if(DEBUG && !wasConverged){
                System.out.println("DEBUG: CG optimization did not converged.");
            }
            
            // apparently it likes to set this to null internally, nice one...
            if(currP == null){
                @SuppressWarnings("unchecked")
                final T res = (T) individual.copy();
                res.setFitness(FixedValues.NONCONVERGEDENERGY);
                return res;
            } else {
                @SuppressWarnings("unchecked")
                final T res = (T) individual.copy();
                res.setFitness(currE);
                back.updateActiveCoordinates(res, currP);
                
                return res;
            }
        }
        
        if(DEBUG && wasConverged){
            System.out.println("DEBUG: CG optimization was converged.");
        }
        
        final double e = result.getValue();
        final double[] p = result.getPoint();
                
        @SuppressWarnings("unchecked")
        final T res = (T) individual.copy();
        res.setFitness(e);
        back.updateActiveCoordinates(res, p);
        
        return res;
    }
    
    // this was my own, very primitive attempt at a convergence checker for the interface. it is still here for reference and testing
    /*private static class SimpleConvChecker implements ConvergenceChecker<PointValuePair>, Serializable {

        private static final long serialVersionUID = (long) 20140907;
        private static final boolean DEBUG = false;
        private final double thresh;
        private final double varThresh;
        private final int maxIter;
        private final long lID;
        private double[] currP;
        private double currE;
        
        SimpleConvChecker(final double thresh, final int maxIter, final long lID){//,
                //final double varThresh){
            this.thresh = thresh;
            this.varThresh = 1e42;//varThresh;
            this.maxIter = maxIter;
            this.lID = lID;
        }
        
        @Override
        public boolean converged(final int i, final PointValuePair pair1, final PointValuePair pair2) {
            
            
            final double e1 = pair1.getValue();
            final double e2 = pair2.getValue();
            final double[] p1 = pair1.getPoint();
            final double[] p2 = pair2.getPoint();
            
            currP = pair2.getPoint();
            currE = e2;
            
            if(i >= maxIter){
                System.out.println("WARNING: CG Local Optimization: maximum number of iterations reached for ID " + lID);
                return true;
            }
            
            if(DEBUG){System.out.println("Old energy: " + e1 + " new energy: " + e2);}
            
            final boolean convFitness = (Math.abs(e1-e2) < thresh);
            
            final boolean convPoints = true;
            /*for(int i){
            }*/
            
        /*    return (convFitness && convPoints);
        }
    }*/
}
