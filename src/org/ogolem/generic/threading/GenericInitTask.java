/**
Copyright (c) 2013, J. M. Dieterich
              2015-2016, J. M. Dieterich and B. Hartke
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
package org.ogolem.generic.threading;

import org.ogolem.core.GlobalConfig;
import org.ogolem.generic.GenericInitializer;
import org.ogolem.generic.Optimizable;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;
import org.ogolem.helpers.Tuple;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A generic implementation of an initialization task.
 * @author Johannes Dieterich
 * @version 2016-10-26
 */
public class GenericInitTask<E,T extends Optimizable<E>,V extends GenericInitializer<E,T>> implements TaskFactory<E,T,V>{
    
    private static final Logger l = LoggerFactory.getLogger(GenericInitTask.class);
    private final boolean DEBUG;
    
    public GenericInitTask(){
        DEBUG = (GlobalConfig.DEBUGLEVEL > 0);
    }
    
    @SuppressWarnings("unchecked")
    @Override
    public Runnable createTask(final GenericPool<E,T> pool, final GenericHistory<E,T> history,
            final V refStuff, final boolean useCache, final ObjectCache<V> cache,
            final boolean doNiching, final NicheComputer<E,T> nicheComp,
            final ObjectCache<NicheComputer<E,T>> nicheCompCache, final long taskID) {
    
        return new Runnable() {

            @Override
            public void run() {
                
                if(DEBUG){System.out.println("Starting generic init task for " + taskID);}
                
                if(useCache && !DEBUG){
                    l.debug("Trying to use cached object as helper for initializing...");
                    final Tuple<Boolean,V> stuff = cache.getUnusedEntry();
                    final V cached = stuff.getObject2();
                    runme(cached,pool,taskID);
                    l.debug("Returning helper for initializing to cache...");
                    stuff.setObject1(false);
                } else if(!useCache && !DEBUG){
                    l.debug("Trying to use new object as helper for initializing...");
                    final V helpers = (V) cache.getOriginalEntry().clone();
                    runme(helpers,pool,taskID);
                } else{
                    try{
                        l.debug("Trying to use new object as helper for initializing (try/catch)...");
                        final V helpers = (V) cache.getOriginalEntry().clone();
                        runme(helpers,pool,taskID);
                    } catch(Throwable t){
                        t.printStackTrace(System.err);
                    }
                }
            }
            
            private void runme(final V helper, final GenericPool<E,T> pool, final long taskID){
                
                l.debug("Starting initialization for " + taskID);
                final T individual = helper.initialize(pool.getExample(), taskID);
                assert(individual != null);

                l.debug("Initialization complete for " + taskID + " with fitness " + individual.getFitness());

                if(doNiching){
                    final Tuple<Boolean,NicheComputer<E,T>> comp = nicheCompCache.getUnusedEntry();
                    final Niche niche = comp.getObject2().computeNiche(individual);
                    comp.setObject1(false);
                    pool.addIndividual(individual, niche, individual.getFitness());
                } else {
                    pool.addIndividual(individual, individual.getFitness());
                }
                l.debug("Finished with initialization for " + taskID);
            }
        };
    }
}
