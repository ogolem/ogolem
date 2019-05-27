/**
Copyright (c) 2010-2013, J. M. Dieterich
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
package org.ogolem.core;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Any number of any global optimization algorithms may be used here.
 * @author Johannes Dieterich
 * @version 2013-11-20
 */
final class UNGlobOpt implements GlobalOptimization {

    private static final long serialVersionUID = (long) 20110606;
    private final int[] percentEnds;
    private final List<GlobalOptimization> globopts;
    private final Random random = new Random();

    UNGlobOpt(final List<GlobalOptimization> globopts, final List<Integer> percents,
            final GlobalConfig config){
        this.globopts = globopts;
        this.percentEnds = new int[globopts.size()];
        if(percents.size() != globopts.size()){
            throw new RuntimeException("Not enough percentages for all global optimizations set.");
        }
        int off = 0;
        for(int i = 0; i < percentEnds.length; i++){
            off += percents.get(i);
            percentEnds[i] = off;
        }
        if(off != 100){
            throw new RuntimeException("Percentages do not add up to 100. Get your math right! Total percents: " + off);
        }
    }

    UNGlobOpt(final UNGlobOpt orig){
        this.globopts = new ArrayList<>();
        for(final GlobalOptimization opt : orig.globopts){
            GlobalOptimization clone = opt.clone();
            this.globopts.add(clone);
        }
        
        this.percentEnds = orig.percentEnds.clone();
    }

    @Override
    public String getMyID(){
        String s = "unitednations:";
        for(final GlobalOptimization opt : globopts){
            s += opt.getMyID() + ";";
        }

        return s;
    }

    @Override
    public UNGlobOpt clone(){
        return new UNGlobOpt(this);
    }

    @Override
    public Geometry doTheGlobOpt(final long futureID, final Geometry geom1, final Geometry geom2) {

        final int point = random.nextInt(100)+1; // to get values 1-100 (inclusive)
        for(int i = 0; i < percentEnds.length; i++){
            if(point <= percentEnds[i]){
                return globopts.get(i).doTheGlobOpt(futureID, geom2, geom2);
            }
        }
        
        System.err.println("ERROR: Excecution flow reached end of UN globopt w/o dispatching to globopt.");
        throw new RuntimeException("Excecution flow reached end of UN globopt w/o dispatching to globopt.");
    }
}
