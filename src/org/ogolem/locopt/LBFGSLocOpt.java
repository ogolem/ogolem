/**
Copyright (c) 2012-2015, J. M. Dieterich
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

import contrib.lbfgs.LBFGS;
import org.ogolem.core.FixedValues;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.generic.GenericAbstractLocOpt;
import org.ogolem.generic.GenericBackend;
import org.ogolem.helpers.Machine;

/**
 * An interface to the slightly adjusted version of the L-BFGS implementation
 * from RISO.
 * @author Johannes Dieterich
 * @version 2015-01-07
 */
public class LBFGSLocOpt <E,T extends ContinuousProblem<E>> extends GenericAbstractLocOpt<E,T> {
    
    private static final long serialVersionUID = (long) 20141216;
    private static final boolean DEBUG = false;
    
    private final int maxSteps;
    private final int noOfCorr;
    private final double conv;
    private final double machPrec;
    private final int maxResets;
    private final double gtol;
    private final int lineIter;
    
    /**
     * Constructor.
     * @param back the backend
     * @param noCorrs number of corrections. Reasonable default: 7, maximum 15.
     * @param maxSteps the max number of optimization steps
     * @param conv the convergence criterion
     * @param gtol metric how expensive a gradient/energy is w.r.t. self time. Must be more than 1e-4, typical small value 0.1, typical value 0.9
     * @param lineIter the maximum number of line search steps per search
     * @param maxResets the number of resets the L-BFGS interface should attempt if it the L-BDFGS signals a problem w/ negative flag
     */
    public LBFGSLocOpt(final GenericBackend<E,T> back, final int noCorrs, final int maxSteps,
            final double conv, final double gtol, final int lineIter, final int maxResets){
        super(back);
        this.machPrec = Machine.calcMachinePrecision();
        this.noOfCorr = noCorrs;
        this.maxSteps = maxSteps;
        this.conv = conv;
        this.maxResets = maxResets;
        this.lineIter = lineIter;
        this.gtol = gtol;
    }
    
    public LBFGSLocOpt(final LBFGSLocOpt<E,T> orig){
        super(orig);
        this.conv = orig.conv;
        this.machPrec = orig.machPrec;
        this.maxSteps = orig.maxSteps;
        this.noOfCorr = orig.noOfCorr;
        this.maxResets = orig.maxResets;
        this.gtol = orig.gtol;
        this.lineIter = orig.lineIter;
    }
    
    @Override
    public LBFGSLocOpt<E, T> clone() {
        return new LBFGSLocOpt<>(this);
    }
    
    @Override
    public String getMyID() {
        return "L-BFGS LOCAL OPTIMIZATION using:\n\t" + back.getMyID();
    }

    @Override
    protected T optimize(final T individual) {
        
        @SuppressWarnings("unchecked")
        final T work = (T) individual.clone();
        final int dims = back.numberOfActiveCoordinates(work);
        final double[] p = back.getActiveCoordinates(work).clone();
        
        final double[] gradient = new double[dims];
        double lastE = back.gradient(p, gradient, 0);                
        
        if(DEBUG) System.out.println("DEBUG: Energy in LBFGS is " + lastE);

        // be silent!
        final int[] print = {-1, 0};
        final boolean beReallySilent = (maxResets > 0);
        final double[] diag = new double[dims];
        final int[] flags = {0, 0};

        // LBFGS iterations
        boolean converged = false;
        LBFGS lbfgs = new LBFGS(gtol, lineIter);
        int resets = 0;
        for (int iter = 0; iter < maxSteps; iter++) {

            try {
                lbfgs.lbfgs(dims, noOfCorr, p, lastE, gradient, false, diag, print, conv, machPrec, flags, beReallySilent);
            } catch (Exception e) {
                System.err.println("WARNING: Problem in LBFGS for ID " + work.getID() + ". Problem is: " + e.toString());
                if(DEBUG) System.out.println("DEBUG: Energy in LBFGS exit is " + lastE);
                break;
            }
            
            if(DEBUG){
                for(int i = 0; i < dims; i++){
                    if(p[i] != p[i]){
                        System.err.println("ERROR: Found NaN in LBFGS at " + i + " iteration " + iter);
                        lastE = FixedValues.NONCONVERGEDENERGY;
                        break;
                    }
                }
            }

            if (flags[0] < 0) {
                // there is a problem (likely linesearch)
                if(resets < maxResets && flags[0] == -1){
                    System.out.println("INFO: L-BFGS reset taking place for ID " + work.getID());
                    flags[0] = 0; // this is internally used to signal reset
                    flags[1] = 0;
                    back.resetToStable(p);
                    // evaluate one energy and gradient
                    lastE = back.gradient(p, gradient, iter);
                    //lbfgs = new LBFGS(gtol, lineIter); // JUST to make sure...
                    for(int i = 0; i < dims; i++){
                        diag[i] = 0.0;
                    }
                    resets++;
                    continue;
                }
                break;
            } else if (flags[0] == 0) {
                // converged
                converged = true;
                break;
            } else if (flags[0] == 1) {
                // evaluate gradient and energy
                lastE = back.gradient(p, gradient, iter);
                if(DEBUG) System.out.println("DEBUG: Energy in LBFGS loop is " + lastE);
            } else {
                System.err.println("ERROR: You should not end up here in LBFGS. Contact developer(s).");
            }
            
            if(DEBUG){
                for(int i = 0; i < gradient.length; i++){
                    if(gradient[i] != gradient[i]){
                        System.err.println("ERROR: Found NaN2 in LBFGS at " + i + " iteration " + iter);
                        
                        lastE = FixedValues.NONCONVERGEDENERGY;
                        break;
                    }
                }
            }
        }
        
        if(!converged){
            System.out.println("INFO: L-BFGS for individual " + work.getID() + " did not converge.");
        }

        @SuppressWarnings("unchecked")
        final T res = (T) individual.clone();
        res.setFitness(lastE);
        back.updateActiveCoordinates(res, p);
        
        return res;
    }    
}
