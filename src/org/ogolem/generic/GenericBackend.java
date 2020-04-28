/**
Copyright (c) 2013-2020, J. M. Dieterich
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
 * A generic interface to be used for energy and gradient evaluations. Demands a
 * continuous problem.
 * @author Johannes Dieterich
 * @author Dominik Behrens
 * @version 2020-04-22
 */
public interface GenericBackend<E,T extends ContinuousProblem<E>> extends GenericFitnessBackend<E,T>{
    
    @Override
    public GenericBackend<E,T> clone();
        
    /**
     * Calculates the fitness and gradient of this individual. getActiveCoordinates(T) must be called prior the first call to this routine!
     * @param currCoords the current coordinates w.r.t. the conditions described. Guaranteed to not be touched!
     * @param gradient an array for the gradient of this individual, replaced upon return. Length must be greater or equal to number of active coordinates!
     * @param iteration a counter to make different invocations easily distinguishable
     * @return the fitness of this individual
     */
    double gradient(final double[] currCoords, final double[] gradient, final int iteration);

    /**
     * By default, GenericBackends do not provide history functionality (historyRememberPoint,
     * historyResetToBestPoint, historyGetBestPoint) and will return false here.
     * However, specific GenericBackends might implement these and will then return true.
     * @return True whenever the Backend supports history methods.
     */
    default boolean supportsHistory() {
        return false;
    }

    /**
     * If this GenericBackend supports history functionality, this method will remember the given values.
     * @param currGenome The promising current genome to be remembered.
     * @param gradient The gradient corresponding to the promising genome.
     * @param value The fitness value for the promising genome.
     */
    default void historyRememberPoint(final double[] currGenome, final double[] gradient, final double value) {
        // This is just a dummy method and will do nothing if no history is implemented.
    }

    /**
     * If this GenericBackend supports history functionality, this method will recall a remembered point and input
     * the remembered data back into the given double arrays. The backend might also do additional internal resets,
     * if necessary.
     * @param currGenome The incoming stub double array to be written to. It will be filled with the remembered
     *                   Genome,
     * @param gradient The incoming stub double array to be written to. It will be filled with the corresponding
     *                 remembered gradient.
     */
    default void historyResetToBestPoint(final double[] currGenome, final double[] gradient) {
        // This is just a dummy method and will do nothing if no history is implemented.
    }

    /**
     * If this GenericBackend supports history functionality, this method will impart the given individual with the
     * remembered genome information.
     * @param individual The individual stub that should be updated with the remembered Information.
     */
    default void historyGetBestPoint(T individual) {
        // This is just a dummy method and will do nothing if no history is implemented.
    }
}
