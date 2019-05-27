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
package org.ogolem.adaptive;

import org.ogolem.generic.GenericAbstractDarwin;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericMutation;
import org.ogolem.generic.GenericSanityCheck;
import org.ogolem.generic.IndividualWriter;

/**
 * A generified (or de-generified?) Darwin implementation based off GenericAbstractDarwin.
 * @author Johannes Dieterich
 * @version 2014-05-02
 */
public class GenericParameterDarwin extends GenericAbstractDarwin<Double,AdaptiveParameters> {
    
    private static final long serialVersionUID = (long) 20140502;

    GenericParameterDarwin(final GenericCrossover<Double,AdaptiveParameters> cross, final GenericMutation<Double,AdaptiveParameters> mut,
            final GenericSanityCheck<Double,AdaptiveParameters> sanity, final GenericFitnessFunction<Double,AdaptiveParameters> fitness,
            final IndividualWriter<AdaptiveParameters> writer, final double crossPoss, final double mutPoss,
            final boolean printBeforeFitness, final int noOfTries){
        super(cross, mut, sanity, fitness, writer, crossPoss, mutPoss,
            printBeforeFitness, noOfTries);
    }
    
    GenericParameterDarwin(final GenericParameterDarwin orig){
        super(orig);
    }
    
    @Override
    public GenericParameterDarwin clone() {
        return new GenericParameterDarwin(this);
    }

    @Override
    public String getMyID() {
        return "GENERIC PARAMETER DARWIN IMPLEMENTATION"
                + "\nusing: "
                + "\n\txover    " + xover.getMyID() + " which probability " + crossPoss*100 + "%"
                + "\n\tmutation " + mutation.getMyID() + " which probability " + mutPoss*100 + "%"
                + "\nfitness  " + fitness.getMyID();
    }

    @Override
    protected void postXOver(final AdaptiveParameters individual1, final AdaptiveParameters individual2, final long futureID) {
        // nothing
    }

    @Override
    protected void postMutation(final AdaptiveParameters individual) {
        // nothing
    }

    @Override
    protected void runAfterEachTry() {
        // does not need to be implemented
    }
}
