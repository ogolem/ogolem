/**
Copyright (c) 2014, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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
 * A fitness-only generic backend to be used for gradient-free optimization
 * implementations. Despite not actually demanding the problem to be of continuous
 * nature yet, we strongly suspect this to be the case and may change this in
 * the future.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public interface GenericFitnessBackend<E,T extends Optimizable<E>> extends GenericFitnessFunction<E,T> {
    
    public enum BOUNDSTYPE{NONE, SOME, ALL};
    
    @Override
    public GenericFitnessBackend<E,T> copy();
    
    /**
     * Returns the number of active, i.e., to be optimized, coordinates in this
     * geometry in this representation.
     * @param individual the individual
     * @return the number of coordinates
     */
    int numberOfActiveCoordinates(final T individual);
    
    /**
     * Get all the active coordinates in this coordinate representation as an
     * array of doubles. This has no guarantees upon the order or anything or
     * magnitude of these values, these are implementation dependent. However,
     * there is one guarantee. The number of elements in the array will be the
     * same as the above returned number of active coordinates. Also may set
     * this individual as a cache to be used in subsequent fitness/gradient
     * evaluations! Hence, it is mutable!
     * @param individual the individual, also potentially used as a mutable cache
     * @return an array of doubles representing all relevant coordinates. do not change without cloning!
     */
    double[] getActiveCoordinates(final T individual);

    /**
     * Does this representation have a concept of boundaries for this individual.
     * @param individual the individual
     * @return whether all coordinates have bounds, or only some, or none
     */
    BOUNDSTYPE boundariesInRepresentation(final T individual);
    
    /**
     * Get best estimate boundaries in this representation. Obviously, the quality of this best guess can
     * only be trusted if a call to boundariesInRepresentation() returned ALL.
     * @param currCoords the current coordinates (from after a getActiveCoordinates call)
     * @param low the lower bounds. Please note, if boundariesInRepresentation has not returned ALL, this will be a best guess!
     * @param high the upper bounds. Please note, if boundariesInRepresentation has not returned ALL, this will be a best guess!
     */
    void bestEstimateBoundaries(final double[] currCoords, final double[] low, final double[] high);
    
    /**
     * Resets, if possible, the backend to some stable state, e.g., to deal with
     * instabilities in the local optimization. May manipulate said coordinate array
     * as if getActiveCoordinates was called!
     * @param coordinates the coordinates, may be changed on exit
     */
    void resetToStable(final double[] coordinates);
    
    /**
     * Update the active coordinates with a new set. The order must be preserved
     * with respect to the getActiveCoordinates returned one. Also, the array
     * must be of length number of active coordinates. Might use information as
     * preserved from the previously cached individual!
     * @param individual the individual
     * @param coordinates the coordinate values to be updated
     */
    void updateActiveCoordinates(final T individual, final double[] coordinates);
    
    /**
     * Calculates the fitness of this individual. getActiveCoordinates(T) must be called prior the first call to this routine!
     * @param currCoords the current coordinates w.r.t. the conditions described. Guaranteed to not be touched!
     * @param iteration a counter to make different invocations easily distinguishable
     * @return the fitness of this individual
     */
    double fitness(final double[] currCoords, final int iteration);
    
}
