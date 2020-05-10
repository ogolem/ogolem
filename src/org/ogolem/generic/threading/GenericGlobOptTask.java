/**
Copyright (c) 2013, J. M. Dieterich
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
package org.ogolem.generic.threading;

import java.util.List;
import org.ogolem.core.GlobalConfig;
import org.ogolem.generic.GenericGlobalOptimization;
import org.ogolem.generic.Optimizable;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;
import org.ogolem.helpers.Tuple;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A generic implementation of a global optimization task.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class GenericGlobOptTask <E,T extends Optimizable<E>,V extends GenericGlobalOptimization<E,T>> implements TaskFactory<E,T,V>{
    
    private static final Logger l = LoggerFactory.getLogger(GenericInitTask.class);
    private final boolean DEBUG = (GlobalConfig.DEBUGLEVEL > 0);
    
    @SuppressWarnings("unchecked")
    @Override
    public Runnable createTask(final GenericPool<E,T> pool, final GenericHistory<E,T> history,
            final V refStuff, final boolean useCache, final ObjectCache<V> cache,
            final boolean doNiching, final NicheComputer<E,T> nicheComp,
            final ObjectCache<NicheComputer<E,T>> nicheCompCache, final long taskID) {
        
        return new Runnable() {
            
            @Override
            public void run() {
                
                // first figure out whether any globopt iteration is needed
                if(pool.acceptableFitnessReached()) {return;}
                
                if(useCache && !DEBUG){
                    try{
                    l.debug("Trying to use cached object as helper for global opt...");
                    final Tuple<Boolean,V> stuff = cache.getUnusedEntry();
                    final V cached = stuff.getObject2();
                    runme(cached,pool,history,taskID);
                    l.debug("Returning helper for global opt to cache...");
                    stuff.setObject1(false);
                    } catch(Throwable e){
                        e.printStackTrace(System.err);
                    }
                } else if(!useCache && !DEBUG){
                    l.debug("Trying to use new object as helper for global opt...");
                    final V helpers = (V) cache.getOriginalEntry().copy();
                    runme(helpers,pool,history,taskID);
                } else{
                    try{
                        l.debug("Trying to use new object as helper for global opt (try/catch)...");
                        final V helpers = (V) cache.getOriginalEntry().copy();
                        runme(helpers,pool,history,taskID);
                    } catch(Throwable t){
                        t.printStackTrace(System.err);
                    }
                }
            }
            
            private void runme(final V helper, final GenericPool<E,T> pool, final GenericHistory<E,T> history, final long taskID){
                
                l.debug("Starting global opt for " + taskID);
                final List<T> parents = pool.getParents();
                if(l.isDebugEnabled()){
                    l.debug("Took geometries " + parents.get(0).getID()
                            + " and " + parents.get(1).getID() + " out.");
                }
                
                final T child = helper.globalOptimization(taskID, parents.get(0), parents.get(1));
                
                if (child == null) {
                    history.addFamily(parents.get(0).getID(),
                            parents.get(1).getID(), (int) taskID, false,
                            true);
                    l.debug("Globopt individual " + taskID + " was null.");
                } else {
                    // non-null'd geometry returned, fine that is
                    assert(child.getID() == taskID);
                    child.setID(taskID); // yes, somewhat superfluous, but in case we run w/o assertation, better be safe than sorry
                    l.debug("Globopt individual " + taskID + " was not null and has fitness " + child.getFitness() + ".");
                    
                    final boolean hasChance = pool.hasChanceToBeAdded(child, child.getFitness());
                    if(!hasChance) {
                    	l.debug("Individual " + taskID + " has no chance of being added to pool.");
                        history.addFamily(parents.get(0).getID(),
                                parents.get(1).getID(), taskID, false,
                                false, child);
                        l.debug("Finished with global opt for " + taskID);
                        return;
                    }
                    
                    boolean accepted = false;
                    if(doNiching){
                        final Tuple<Boolean,NicheComputer<E,T>> comp = nicheCompCache.getUnusedEntry();
                        final Niche niche = comp.getObject2().computeNiche(child);
                        comp.setObject1(false);
                        accepted = pool.addIndividual(child, niche, child.getFitness());
                    }  else {
                        accepted = pool.addIndividual(child, child.getFitness());
                    }
                    l.debug("Was globopt individual " + taskID + " accepted? " + accepted);
                    history.addFamily(parents.get(0).getID(),
                            parents.get(1).getID(), taskID, accepted,
                            false, child);
                }
                
                l.debug("Finished with global opt for " + taskID);
            }
        };
    }
}
