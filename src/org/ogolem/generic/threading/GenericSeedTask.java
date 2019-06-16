/**
Copyright (c) 2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import org.ogolem.core.GlobalConfig;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.IndividualReader;
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
 * @version 2015-05-26
 */
public class GenericSeedTask<E,T extends Optimizable<E>> implements TaskFactory<E,T,Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>>>{
    
    private static final Logger l = LoggerFactory.getLogger(GenericSeedTask.class);
    private final boolean DEBUG = (GlobalConfig.DEBUGLEVEL > 0);
    private static List<String> seedData = null;
    private static Iterator<String> seedIterator = null;
    
    public GenericSeedTask(final String seedPath) throws Exception {
        seedData = Collections.synchronizedList(new ArrayList<>());
        final File f = new File(seedPath);
        if (!f.exists() || !f.isDirectory()) {
            throw new RuntimeException("Seeding path must be a directory!");
        }
        for (final String singleSeed : f.list()) {
            seedData.add(seedPath + File.separator + singleSeed);
        }
        seedIterator = seedData.iterator();
    }
    
    public static int getNoOfSeeds(){
        return seedData.size();
    }
    
    @SuppressWarnings("unchecked")
    @Override
    public Runnable createTask(final GenericPool<E,T> pool, final GenericHistory<E,T> history,
            final Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>> refStuff, final boolean useCache,
            final ObjectCache<Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>>> cache, 
            final boolean doNiching, final NicheComputer<E,T> nicheComp,
            final ObjectCache<NicheComputer<E,T>> nicheCompCache, final long taskID) {
    
        return new Runnable() {

            @Override
            public void run() {
                
                if(DEBUG){System.out.println("Starting generic init task for " + taskID);}
                
                if(useCache && !DEBUG){
                    l.debug("Trying to use cached object as helper for initializing...");
                    final Tuple<Boolean,Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>>> stuff = cache.getUnusedEntry();
                    final Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>> cached = stuff.getObject2();
                    runme(cached,pool,taskID);
                    l.debug("Returning helper for initializing to cache...");
                    stuff.setObject1(false);
                } else if(!useCache && !DEBUG){
                    l.debug("Trying to use new object as helper for initializing...");
                    Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>> helpers = null;
                    try{
                        helpers = cache.getOriginalEntry().clone();
                    } catch(Exception e){
                        System.err.println("Exception in cloneing reader and fitness function.");
                        e.printStackTrace(System.err);
                        throw new RuntimeException(e);
                    }
                    assert(helpers != null);
                    runme(helpers,pool,taskID);
                } else{
                    try{
                        l.debug("Trying to use new object as helper for initializing (try/catch)...");
                        final Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>> helpers = cache.getOriginalEntry().clone();
                        runme(helpers,pool,taskID);
                    } catch(Throwable t){
                        t.printStackTrace(System.err);
                    }
                }
            }
            
            private void runme(final Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>> helper, final GenericPool<E,T> pool, final long taskID){
                
                l.debug("Starting seeding for " + taskID);
                
                l.debug("our helper... " + helper.toString());
                l.debug("our initer... " + helper.getObject1().toString());
                l.debug("our locopt... " + helper.getObject2().toString());
        
                String thisSeed;
                synchronized(seedData){
                    if(!seedIterator.hasNext()){return;}
                    thisSeed = seedIterator.next();
                }
                l.info("Seed " + thisSeed + " used for individual " + taskID);
                
                // read
                final T work = (T) pool.getExample().clone();
                try{
                    helper.getObject1().populateIndividualFromFile(work, thisSeed);
                } catch(Exception e){
                    System.err.println("Failure to populate individual from file " + thisSeed);
                    e.printStackTrace(System.err);
                }
                work.setID(taskID);
                assert(work != null);
                                
                final T x = helper.getObject2().fitness(work, false);
                  
                l.info("Seeding complete for " + taskID + " with fitness " + x.getFitness());

                if(doNiching){
                    final Tuple<Boolean,NicheComputer<E,T>> comp = nicheCompCache.getUnusedEntry();
                    final Niche niche = comp.getObject2().computeNiche(x);
                    comp.setObject1(false);
                    pool.addIndividualForced(x, niche, x.getFitness());
                } else {
                    pool.addIndividualForced(x, x.getFitness());
                }
                l.debug("Finished with seeding for " + taskID);
            }
        };
    }
}
