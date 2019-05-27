/**
Copyright (c) 2014, J. M. Dieterich
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

import org.ogolem.generic.GenericBackend;

/**
 * A basic harmonic function, y = a*(x-x0)**2 + b
 * @author Johannes Dieterich
 * @version 2015-03-28
 */
public class HarmonicFunction implements GenericBackend<Double,BasicOptimizableType>{
    
    private static final long serialVersionUID = (long) 20150328;
    private final double x0;
    private final double a;
    private final double b;
    private final double boundUp;
    private final double boundLow;
    
    public HarmonicFunction(final double a, final double x0, final double b){
        this.a = a;
        this.x0 = x0;
        this.b = b;
        this.boundLow = Double.NaN;
        this.boundUp = Double.NaN;
    }
    
    public HarmonicFunction(final double a, final double x0, final double b,
            final double boundLow, final double boundUp){
        this.a = a;
        this.x0 = x0;
        this.b = b;
        this.boundLow = boundLow;
        this.boundUp = boundUp;
    }
    
    @Override
    public HarmonicFunction clone() {
        return new HarmonicFunction(a,x0,b,boundUp,boundLow);
    }

    @Override
    public String getMyID() {
        return "harmonic function test: a = " + a + "; x0 = " + x0 + "; b = " + b;
    }
    
    @Override
    public int numberOfActiveCoordinates(BasicOptimizableType individual) {
        return individual.getNoOfCoords();
    }

    @Override
    public double[] getActiveCoordinates(BasicOptimizableType individual) {
        return individual.getGenomeAsDouble().clone();
    }

    @Override
    public double gradient(double[] currCoords, double[] gradient, int iteration) {
        
        double fitness = 0.0;
        for(int i = 0; i < currCoords.length; i++){
            final double d = currCoords[i]-x0;
            fitness += a*d*d+b;
            gradient[i] = 2*a*(currCoords[i]-x0);
        }
        
        return fitness;
    }

    @Override
    public void updateActiveCoordinates(final BasicOptimizableType individual, double[] coordinates) {
        individual.setGenome(coordinates);
    }

    @Override
    public double fitness(double[] currCoords, int iteration) {
        
        double fitness = 0.0;
        for(int i = 0; i < currCoords.length; i++){
            final double d = currCoords[i]-x0;
            fitness += a*d*d+b;
        }
        
        return fitness;
    }

    @Override
    public BasicOptimizableType fitness(BasicOptimizableType individual, boolean forceOneEval) {
        final double[] actCoord = getActiveCoordinates((BasicOptimizableType) individual.clone());
        final double fitness = fitness(actCoord,1);
        
        individual.setFitness(fitness);
        
        return individual;
    }

    @Override
    public void resetToStable(final double[] coordinates) {
    }

    @Override
    public BOUNDSTYPE boundariesInRepresentation(final BasicOptimizableType individual) {
        return (boundLow != boundLow) ? BOUNDSTYPE.NONE : BOUNDSTYPE.ALL;
    }

    @Override
    public void bestEstimateBoundaries(final double[] currCoords, final double[] low, final double[] high) {
        
        if(boundLow != boundLow){
            // not really bounds....
            for(int i = 0; i < currCoords.length; i++){
                low[i] = currCoords[i] - 9999.9999;
                high[i] = currCoords[i] + 9999.9999;
            }
        } else {
            // use the bounds
            for(int i = 0; i < currCoords.length; i++){
                low[i] = boundLow;
                high[i] = boundUp;
            }
        }
    }
}
