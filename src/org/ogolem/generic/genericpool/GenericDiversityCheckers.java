/**
Copyright (c) 2012,2014, J. M. Dieterich
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
package org.ogolem.generic.genericpool;

import org.ogolem.generic.Optimizable;

/**
 * A collection of generic diversity checkers
 * @author Johannes Dieterich
 * @version 2014-06-12
 */
public class GenericDiversityCheckers {
    
    public static class FitnessDiversityChecker<E,T extends Optimizable<E>> implements DiversityChecker<E,T>{
        
        private static final long serialVersionUID = (long) 20140612;
        private final double diversity;
        
        public FitnessDiversityChecker(final double fitnessDiv){
            assert(fitnessDiv >= 0.0);
            diversity = fitnessDiv;
        }
        
        @Override
        public boolean areDiverse(final GenericPoolEntry<E,T> individuum1, final GenericPoolEntry<E,T> individuum2){
            
            final double fit1 = individuum1.getFitness();
            final double fit2 = individuum2.getFitness();
            
            return (Math.abs(fit1-fit2) > diversity);
        }

        @Override
        public String getMyName() {
            return "standard fitness diversity, threshold " + diversity;
        }
    }
    
    /**
     * A fitness diversity checker based on the percentage difference w.r.t. to the
     * individual with smaller fitness.
     * @param <E> a type, e.g., Double
     * @param <T> an optimizable of type E
     */
    public static class PercentageFitnessDiversityChecker<E,T extends Optimizable<E>> implements DiversityChecker<E,T>{
        
        private static final long serialVersionUID = (long) 20140612;
        private final double diversity;
        
        /**
         * Constructor.
         * @param percentage percentage in the interval [0,1.0]
         */
        public PercentageFitnessDiversityChecker(final double percentage){
            assert(percentage >= 0.0);
            assert(percentage <= 1.0);
            diversity = percentage;
        }
        
        @Override
        public boolean areDiverse(final GenericPoolEntry<E,T> individuum1, final GenericPoolEntry<E,T> individuum2){
            
            final double fit1 = individuum1.getFitness();
            final double fit2 = individuum2.getFitness();
            
            return (Math.abs(fit1-fit2)/Math.abs(Math.min(fit1, fit2)) > diversity);
        }

        @Override
        public String getMyName() {
            return "percentage fitness diversity, percentage " + diversity*100 + "%";
        }
    }
    
    public static class ScaledFitnessDiversityChecker<E,T extends Optimizable<E>> implements DiversityChecker<E,T>{
        
        private static final long serialVersionUID = (long) 20131231;
        private final double diversity;
        
        public ScaledFitnessDiversityChecker(final double fitnessDiv){
            diversity = fitnessDiv;
        }
        
        @Override
        public boolean areDiverse(final GenericPoolEntry<E,T> individuum1, final GenericPoolEntry<E,T> individuum2){
            
            final double fit1 = individuum1.getFitness();
            final double fit2 = individuum2.getFitness();
            
            return (Math.abs(fit1-fit2) > diversity);
        }
    
    //private boolean checkScaledEnergyDiversity(Geometry geom1, Geometry geom2) {
    //    boolean bDiversity = false;
    //    //TODO scaled diversity check, adjusted sigmoid function again
    //    return bDiversity;
    //}

        @Override
        public String getMyName() {
            return "scaled fitness diversity checker, threshold " + diversity;
        }
    }
}
