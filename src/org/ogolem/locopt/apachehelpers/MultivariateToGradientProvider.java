/**
Copyright (c) 2014, J. M. Dieterich
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
package org.ogolem.locopt.apachehelpers;

import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.generic.GenericBackend;

/**
 * Translates apache's multivariate to an gradient provider call.
 * @author Johannes Dieterich
 * @version 2014-12-15
 */
public class MultivariateToGradientProvider<E,T extends ContinuousProblem<E>> implements MultivariateVectorFunction {

    private static final long serialVersionUID = (long) 20140907;
    private static final boolean DEBUG = false;
    protected final GenericBackend<E,T> gradprov;
    protected T individual;
    protected int iter = 0;
    private double lastE = Double.MAX_VALUE;
    
    public MultivariateToGradientProvider(final GenericBackend<E,T> prov,
            final T individual) {
        this.gradprov = prov;
        this.individual = individual;
    }

    @Override
    public double[] value(final double[] doubles) throws IllegalArgumentException {
        
        final double[] gradient = new double[doubles.length];
        final double e = gradprov.gradient(doubles, gradient, iter);
        lastE = e;
        iter++;

        if(DEBUG){
            System.out.println("DEBUG: gradient is: XXX");
        }
        
        return gradient;
    }
    
    public double lastEnergy(){
        return lastE;
    }
}