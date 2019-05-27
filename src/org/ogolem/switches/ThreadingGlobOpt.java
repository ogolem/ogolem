/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
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
package org.ogolem.switches;

import java.util.concurrent.*;
import java.util.List;
import org.ogolem.generic.genericpool.GenericPool;

/**
 * Does the globopt of the switches in a threading fashion.
 * @author Johannes Dieterich
 * @version 2010-08-03
 */
final class ThreadingGlobOpt {
    
    private final SwitchesConfig conf;
    private final GenericPool<Color,Switch> pool;

    private final int iThreads;

    private final int iOffset;

    private final int iIterations;

    ThreadingGlobOpt(SwitchesConfig globconf, int iNoOfThreads, GenericPool<Color,Switch> pool){
        this.conf = globconf;
        this.iThreads = iNoOfThreads;
        this.iOffset = SwitchesConfig.iPoolSize;
        this.iIterations = SwitchesConfig.iNoOfGlobIters;
        this.pool = pool;
    }

    void doAllGlobOpts(){

        final ExecutorService threadpool = Executors.newFixedThreadPool(iThreads);

        // do it for all
        for(int i = iOffset; i < (iIterations+iOffset); i++){
            //System.OUT.println("Task " + i + " fired off.");
            threadpool.submit(createGlobOptTask(pool,i,conf, Taboos.getReference()));
        }

        // shut the pool down
        threadpool.shutdown();

        // wait for all threads to finish
        try{
            threadpool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        } catch(InterruptedException e){
            System.err.println("Threadpool reached wallclock limit. This should really NEVER happen! " + e.toString());
        }
    }


  /**
   * Creates a simple Runnable doing the actual globopt.
   */
    private static Runnable createGlobOptTask(final GenericPool<Color,Switch> pool, final int position,
            final SwitchesConfig config, final Taboos taboos) {
        
        return new Runnable() {

            @Override
            public void run() {

                final List<Switch> vParents = pool.getParents();
                final SwitchesGlobOpt globopt = new SwitchesGlobOpt(config);
                final Switch sChild = globopt.doTheGlobOpt(position,
                        vParents.get(0), vParents.get(1));

                /*
                 * currently not used, remains here for the eventual creation
                 * of a history object
                 */
                boolean accepted;

                if (sChild != null){
                    // non-null'd switch returned, fine that is
                    
                    // switches do not support niching
                    accepted = pool.addIndividual(sChild, sChild.getFitness());

                    taboos.addTaboo(sChild);
                }
            }
        };
    }
}
