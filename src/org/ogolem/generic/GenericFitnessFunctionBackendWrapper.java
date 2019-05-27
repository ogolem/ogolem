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
package org.ogolem.generic;

/**
 * A basic wrapper for interface translation
 * @author Johannes Dieterich
 * @version 2015-03-28
 */
public class GenericFitnessFunctionBackendWrapper<T extends Optimizable<Double>> implements GenericFitnessBackend<Double,T> {

    private static final long serialVersionUID = (long) 20141217;
    private final GenericFitnessFunction<Double,T> func;
    private final double[] lower;
    private final double[] upper;
    private T cache;
    
    /**
     * Constructor to wrap a fitness function into a backend.
     * @param function the fitness function to be wrapped
     * @param lowerBound the lower bound, should be null if not bounded
     * @param upperBound the upper bound, should be null if not bounded
     */
    public GenericFitnessFunctionBackendWrapper(final GenericFitnessFunction<Double,T> function, final double[] lowerBound, final double[] upperBound){
        this.func = function;
        this.lower = lowerBound;
        this.upper = upperBound;
    }
    
    public GenericFitnessFunctionBackendWrapper(final GenericFitnessFunctionBackendWrapper<T> orig){
        this.func = orig.func.clone();
        this.lower = (orig.lower == null) ? null : orig.lower.clone();
        this.upper = (orig.upper == null) ? null : orig.upper.clone();
    }
    
    @Override
    public GenericFitnessFunctionBackendWrapper<T> clone() {
        return new GenericFitnessFunctionBackendWrapper<>(this);
    }

    @Override
    public String getMyID() {
        return func.getMyID();
    }

    
    @Override
    public int numberOfActiveCoordinates(final T individual) {
        return individual.getGenomeCopy().length;
    }

    @Override
    public double[] getActiveCoordinates(final T individual) {
        this.cache = individual;
        return cache.getGenomeAsDouble();
    }

    @Override
    public void resetToStable(final double[] coordinates){
        // nothing
    }
    
    @Override
    public void updateActiveCoordinates(final T individual, final double[] coordinates) {
        final Double[] gen = new Double[coordinates.length];
        for(int x = 0; x < coordinates.length; x++){
            gen[x] = coordinates[x];
        }
        individual.setGenome(gen);
    }

    @Override
    public double fitness(final double[] currCoords, final int iteration) {
        updateActiveCoordinates(cache,currCoords);
        return func.fitness(cache, true).getFitness();
    }

    @Override
    public T fitness(T individual, boolean forceOneEval) {
        return func.fitness(individual, forceOneEval);
    }

    @Override
    public BOUNDSTYPE boundariesInRepresentation(final T individual) {
        if(this.lower == null || this.upper == null){
            return BOUNDSTYPE.NONE;
        } else {
            // XXX not too fine-grained control
            return BOUNDSTYPE.ALL;
        }
    }

    @Override
    public void bestEstimateBoundaries(final double[] currCoords, final double[] low, final double[] high) {
        
        if(lower != null && upper != null){
            // just copy
            if(low.length != high.length || low.length != currCoords.length || low.length != upper.length || low.length != lower.length){
                throw new RuntimeException("Something is off with the lengths of arrays here, sorry!");
            }
            System.arraycopy(lower, 0, low, 0, low.length);
            System.arraycopy(upper, 0, high, 0, low.length);
        }
        
        // oh well, this is a tough one
        System.err.println("WARNING: Entering uncharted land, no bounds are known for the generic fitness function wrapper but it is asked to provide some...");
        for(int i = 0; i < low.length; i++){
            low[i] = -(Double.MAX_VALUE-10);
            high[i] = Double.MAX_VALUE;
        }
    }
}