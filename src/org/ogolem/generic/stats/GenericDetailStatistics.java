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
package org.ogolem.generic.stats;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * A statistic object for detailed statistics. Please note: this MUST ensure
 * that calls are NO-OP if statistics are not enabled AND also note that this
 * obviously does not work when running in an MPP fashion (via either RMI or MPI).
 * @author Johannes Dieterich
 * @version 2014-07-09
 */
public class GenericDetailStatistics {
    
    private static boolean isEnabled = false;
    private static DetailedStatistics backend = new DummyDetailedStats();
    
    private GenericDetailStatistics(){}; // disallow instantiation
    
    public synchronized static void enableAllDetails(){
        backend = new BasicDetailedBackend();
        isEnabled = true;
    }
    
    public static void incrementTrials(){
        backend.incrementTrials();
    }
    
    public long getTotalTrials(){
        return backend.getTotalTrials();
    }
    
    public static void incrementLocalOpts(){
        backend.incrementLocalOpts();
    }
    
    public static long getTotalLocOpts(){
        return backend.getTotalLocOpts();
    }
    
    public static void incrementFitnessEvals(){
        backend.incrementFitnessEvals();
    }
    
    public static long getTotalFitnessEvals(){
        return backend.getTotalFitnessEvals();
    }
    
    public static void incrementGradientEvals(){
        backend.incrementGradientEvals();
    }
    
    public static long getTotalGradientEvals(){
        return backend.getTotalGradientEvals();
    }
    
    public static void incrementSanityDiscards(){
        backend.incrementSanityDiscards();
    }
    
    public static long getTotalSanityDiscards(){
        return backend.getTotalSanityDiscards();
    }
    
    public static void incrementUnknown(final String id){
        backend.incrementUnknown(id);
    }
    
    public static Map<String,Long> getAllUnknownCounters(){
        return backend.getAllUnknownCounters();
    }
    
    public static List<String> getOutput(){
        
        if(!isEnabled){
            final List<String> out = new LinkedList<>();
            out.add("#######################################");
            out.add("");
            out.add("Detailed statistics were not enabled.");
            out.add("");
            out.add("#######################################");
            return out;
        }
        
        final List<String> out = new LinkedList<>();
        out.add("");
        out.add("");
        out.add("#######################################");
        out.add("");
        out.add("");
        out.add("                 WARNING");
        out.add("IF YOU INTEND TO BASE ANYTHING FOR PUBLICATIONS");
        out.add("ON THESE COUNTERS REMEMBER THEY ARE APPROXIMATIVE.");
        out.add("WE STRIVE TO IMPROVE THEM AND THEY SHOULD BE SENSIBLE");
        out.add("FOR RELATIVE COMPARISONS (APPLE TO APPLE) BUT WE MAKE");
        out.add("NO GUARANTEES!");
        out.add("YOU HAVE BEEN WARNED!");
        out.add("");
        out.add("");
        out.add("#######################################");
        out.add("");
        out.add("Total number of trials:                " + backend.getTotalTrials());
        out.add("      translated to child individuals: " + 2*backend.getTotalTrials());
        out.add("Number of sanity based discards:       " + backend.getTotalSanityDiscards());
        out.add("Number of local optimizations:         " + backend.getTotalLocOpts());
        out.add("Number of fitness evaluations:         " + backend.getTotalFitnessEvals());
        out.add("Number of gradient evaluations:        " + backend.getTotalGradientEvals());
        
        final Set<Entry<String,Long>> customs = backend.getAllUnknownCounters().entrySet();
        for(final Entry<String,Long> entry : customs){
            out.add("Number of " + entry.getKey() + ": " + entry.getValue());
        }
        
        out.add("");
        out.add("#######################################");
        out.add("");
        out.add("");

        return out;
    }
}
