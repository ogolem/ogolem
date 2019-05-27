/**
Copyright (c) 2014-2015, J. M. Dieterich
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

import java.io.Serializable;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.ogolem.generic.GenericFitnessBackend;
import org.ogolem.generic.Optimizable;

/**
 * Translates apache's multivariate to an energy provider call.
 * @author Johannes Dieterich
 * @version 2015-01-12
 */
public class MultivariateToEnergyProvider<E, T extends Optimizable<E>> implements MultivariateFunction, Serializable {
    
    private static final long serialVersionUID = (long) 20150112;
    private static final boolean DEBUG = false;
    protected final T individual;
    protected final GenericFitnessBackend<E,T> prov;
    protected int iter = 0;
    private double[] bestP = null;
    private double bestE = Double.MAX_VALUE;
    
    public MultivariateToEnergyProvider(final GenericFitnessBackend<E,T> prov,
            final T individual){
        this.prov = prov;
        this.individual = individual;
    }

    @Override
    public double value(final double[] doubles) {
        
        bestE = Double.MAX_VALUE;
        if(bestP == null){
            bestP = new double[doubles.length];
        }
        if(DEBUG){
            System.out.println("DEBUG: iter " + iter + " coordinates are:");
            for(final double d : doubles){
                System.out.println("\t" + d);
            }
        }
        
        final double e = prov.fitness(doubles, iter);
        if(e < bestE){
            bestE = e;
            System.arraycopy(doubles, 0, bestP, 0, doubles.length);
        }
        
        if(DEBUG){
            System.out.println("DEBUG: iter " + iter + " energy is: " + e);
        }
        iter++;
        
        return e;
    }
    
    public double bestEnergy(){
        return bestE;
    }
    
    public double[] bestPoint(){
        return bestP;
    }
}