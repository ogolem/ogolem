/**
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
package org.ogolem.generic.stats;

import java.util.HashMap;
import java.util.Map;

/**
 * A basic backend.
 * @author Johannes Dieterich
 * @version 2014-07-08
 */
public class BasicDetailedBackend implements DetailedStatistics {
    
    private long countTotTrials = 0l;
    private long countLocOpts = 0l;
    private long countFitEval = 0l;
    private long countGradEval = 0l;
    private long countSanDisc = 0l;
    private final Map<String,Long> countUnknown = new HashMap<>();
    
    @Override
    public synchronized void incrementTrials(){
        countTotTrials++;
    }
    
    @Override
    public long getTotalTrials(){
        return countTotTrials;
    }

    @Override
    public synchronized void incrementLocalOpts() {
        countLocOpts++;
    }

    @Override
    public long getTotalLocOpts() {
        return countLocOpts;
    }

    @Override
    public synchronized void incrementFitnessEvals() {
        countFitEval++;
    }

    @Override
    public long getTotalFitnessEvals() {
        return countFitEval;
    }
    
    @Override
    public synchronized void incrementGradientEvals(){
        countGradEval++;
    }
    
    @Override
    public long getTotalGradientEvals(){
        return countGradEval;
    }

    @Override
    public synchronized void incrementSanityDiscards() {
        countSanDisc++;
    }

    @Override
    public long getTotalSanityDiscards() {
        return countSanDisc;
    }

    @Override
    public synchronized void incrementUnknown(final String id) {
        countUnknown.put(id, 0l);
    }

    @Override
    public Map<String, Long> getAllUnknownCounters() {
        return countUnknown;
    }
}
