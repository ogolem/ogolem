/**
Copyright (c) 2013-2014, J. M. Dieterich
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
package org.ogolem.freqs;

import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.helpers.Tuple;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Provides a tiny helper function to decide if a structure is a local minimum
 * based on calculating its harmonic frequencies. Obviously, this is a rather
 * painful process and therefore we try to not forget any data we assembled during.
 * @author Johannes Dieterich
 * @version 2014-11-21
 */
public class MinimumDecider {
    
    private static final Logger log = LoggerFactory.getLogger(MinimumDecider.class);
    
    public static Tuple<Boolean,Frequencies> frequencyAnalysis(final CartesianCoordinates cartes,
            final BondInfo bonds, final FrequencyMethod method){
        
        assert(cartes != null);
        assert(bonds != null);
        assert(method != null);
        
        if(log.isDebugEnabled()){
            log.debug("Cartesian we are working on: ");
            final String[] sa = cartes.createPrintableCartesians();
            String s2 = "";
            for(final String s : sa){
                s2 += s + "\n";
            }
            log.debug(s2);
        }
        
        final Frequencies freqs = HarmonicFrequencyCalculator.calculateFrequencies(cartes, method, bonds);
        assert(freqs != null);
        log.debug("Finished calculating harmonic frequencies. Now returning...");
        
        return new Tuple<>(!freqs.anyImaginaryFrequency(), freqs);
    }
}
