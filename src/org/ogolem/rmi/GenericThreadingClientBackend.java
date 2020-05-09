/**
Copyright (c) 2016-2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.rmi;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.ogolem.generic.Configuration;
import org.ogolem.generic.Copyable;
import org.ogolem.generic.GenericGlobalOptimization;
import org.ogolem.generic.GenericInitializer;
import org.ogolem.generic.Optimizable;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.generichistory.GenericHistoryConfig;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.GenericPoolConfig;
import org.ogolem.generic.genericpool.GenericPoolEntry;
import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;
import org.ogolem.generic.threading.ObjectCache;
import org.ogolem.generic.threading.TaskFactory;
import org.ogolem.helpers.Tuple;
import org.ogolem.helpers.Tuple3D;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * The actual backend for the threading RMI client.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
class GenericThreadingClientBackend<E,T extends Optimizable<E>> {
    
    private static final Logger l = LoggerFactory.getLogger(GenericThreadingClientBackend.class);
    private static final boolean DEBUG = false;
    private static final short MAXCONTACTATTEMPTS = 20;
    
    private final Configuration<E,T> config;
    private final int threads;
    private final GenericPool<E,T> pool;
    private final GenericHistory<E,T> history;
    private final boolean doNiching;
    private final boolean useObjCache;
    private final NicheComputer<E,T> nicheComp;
    private final ObjectCache<NicheComputer<E,T>> nicheCompCache;
    private final RMICommunication<T> comm;
    private final int indsToMerge;
    private final int myID;
    private final int maxTasks;
    private final boolean doMaxStructsExchange;
    private final long sleepTime;
    
    private long currIDStart = 0;
    private long tasksGotten;
    private int noInits = 0;
    
    private int initTasks = 0;
    private int initChunksGotten = 0;
    private int chunksGotten = 0;
    
    private long initWorkTimeMS = 0l;
    private long initFetchTimeMS = 0l;
    private long optWorkTimeMS = 0l;
    private long optFetchTimeMS = 0l;
    private long waitTimeMS = 0l;
    
    private static boolean ISINIT = true;
    
    GenericThreadingClientBackend(final Configuration<E,T> config, final int threads,
            final int myID, final RMICommunication<T> comm, final int maxTasks,
            final boolean doMaxStructsExchange,
            final int indsToMerge, final long sleepTime) throws Exception {
        
        final GenericPoolConfig<E,T> poolConfig = config.getGenericPoolConfig();
        final T example = config.getExample();
        this.pool = GenericPool.getInstance(poolConfig, example);
        
        final GenericHistoryConfig hisConf = config.getGenericHistoryConfig();
        this.history = GenericHistory.getReference(hisConf);
        
        this.config = config;
        this.threads = threads;
        this.doNiching = (config.getNicheComputer() != null);
        this.nicheComp = config.getNicheComputer();
        if(doNiching){
            try{
                this.nicheCompCache = new ObjectCache<>(2*threads,nicheComp);
            } catch (Exception e){
                throw new RuntimeException(e);
            }
        } else {
            this.nicheCompCache = null;
        }
        this.comm = comm;
        this.indsToMerge = indsToMerge;
        this.myID = myID;
        this.doMaxStructsExchange = doMaxStructsExchange;
        this.sleepTime = sleepTime;
        this.useObjCache = (!DEBUG);
        this.maxTasks = maxTasks;
    }
    
    int getNoOfInitChunks(){
        return initChunksGotten;
    }
    
    int getNoOfGlobOptChunks(){
        return chunksGotten;
    }
    
    long getTimeInSInitFetch(){
        return initFetchTimeMS / 1000;
    }
    
    long getTimeInSInitWork(){
        return initWorkTimeMS / 1000;
    }
    
    long getTimeInSOptWork(){
        return optWorkTimeMS / 1000;
    }
    
    long getTimeInSOptFetch(){
        return optFetchTimeMS / 1000;
    }
    
    long getTimeInSWaited(){
        return waitTimeMS / 1000;
    }
    
    public void doAllTasks() throws Exception {
        
        // first do initialization tasks (setup task factory and cachers)
        final GenericInitializer<E,T> initializer = config.getInitializer();
        final TaskFactory<E,T,GenericInitializer<E,T>> initTaskFac = config.getInitFactory();
        final ObjectCache<GenericInitializer<E,T>> initCache = new ObjectCache<>(2*threads,initializer);
        
        final long sT = System.currentTimeMillis();
        List<Task<T>> initTasks = comm.getInitialTasks(maxTasks,myID);
        final long mT = System.currentTimeMillis();
        initFetchTimeMS += mT-sT;
        initChunksGotten++;
        
        do{
            noInits = initTasks.size();
            if(noInits == 0){
                // nothing to be done, kinda weird.
                System.out.println("WARNING: No tasks gotten for threading backend " + myID + " continuing but wondering...");
                break;
            }
            // kind of subptimal but we need to get the ID offset
            currIDStart = (int) initTasks.get(0).getDummyAnswer(myID).getResult().getID();
            
            if(DEBUG){System.out.println("DEBUG: Starting with ID " + currIDStart  + " we have " + initTasks.size() + " tasks to init.");}

            final long tOff = System.currentTimeMillis();
            doXXX(threads, initTaskFac, initCache, initTasks);
            final long tMid = System.currentTimeMillis();
            initWorkTimeMS += (tMid-tOff);
            
            final boolean doneYet = synchronizePools(true);
            if(!doneYet){
                initTasks = comm.getInitialTasks(maxTasks,myID);
                initChunksGotten++;
            } else {
                initTasks.clear(); // not strictly necessary, I think
            }
            final long tEnd = System.currentTimeMillis();
            initFetchTimeMS += (tEnd-tMid);
            
        } while(!initTasks.isEmpty());
        
        // merge pools, technically this should have happened before already. However, this does not hurt and synchronizes even further.
        synchronizePools(false);
        
        ISINIT = false;
        
        // then do globopt tasks (setup task factory and cachers)
        final GenericGlobalOptimization<E,T> globopt = config.getGlobalOptimization();
        final ObjectCache<GenericGlobalOptimization<E,T>> globCache = new ObjectCache<>(2*threads,globopt);
        TaskFactory<E,T,GenericGlobalOptimization<E,T>> globTasks = new GlobTaskFactory<>();

        while(!comm.isEverythingDone()){
            
            final long tOff = System.currentTimeMillis();
            doXXX(threads, globTasks, globCache, (int) currIDStart, (int) tasksGotten);
            final long tMid = System.currentTimeMillis();
            optWorkTimeMS += (tMid-tOff);
            
            final boolean done = synchronizePools(false);
            final long tEnd = System.currentTimeMillis();
            optFetchTimeMS += (tEnd-tMid);
            
            if(done) {break;}
        }
    }
    
    private boolean synchronizePools(final boolean isInInit) throws Exception {
        
        if (doMaxStructsExchange) {
            final List<T> myPool = new ArrayList<>();
            int c = 0;
            int offset = 0;
            // always synchronize our top, NEW individuals with the master            
            while(c < indsToMerge && pool.getCurrentPoolSize() > 0 && offset < pool.getCurrentPoolSize()){
                final T ind = pool.getIndividualAtPosition(offset); // we remove from the top (ignoring known IDs)
                final long id = ind.getID();
                if(id < currIDStart){
                    offset++; // increment the offset
                    continue;
                }
                myPool.add(ind);
                pool.removeIndividualAtPos(offset);
                c++;
            }
            final List<T> merged = comm.synchronizePool(myPool, myID, indsToMerge, currIDStart);
            if (merged.size() > indsToMerge) {
                throw new Exception("Getting more individuals back than I thought!");
            }
            merged.forEach((t) -> {
                // augment pool... we however will need to check if the IDs are the same!
                // (i.e., the main pool may give something back that actually originated from us
                // but then, we also want to get "older" IDs that originated from other proxies
                if (t.getID() < currIDStart) {
                    // as this is "expensive", only do this for the IDs that are older
                    if (doNiching) {
                        final Niche niche = nicheComp.computeNiche(t);
                        pool.addIndividualForcedUnsyncCheckID(t, niche, t.getFitness(), 5);
                    } else {
                        pool.addIndividualForcedUnsyncCheckID(t, null, t.getFitness(), 5);
                    }
                } else {
                    if (doNiching) {
                        final Niche niche = nicheComp.computeNiche(t);
                        pool.addIndividualForced(t, niche, t.getFitness());
                    } else {
                        pool.addIndividualForced(t, t.getFitness());
                    }
                }
            });
        } else {            
            final List<T> myPool = new ArrayList<>();
            for (final GenericPoolEntry<E, T> entry : pool) {
                final T ind = entry.getIndividual();
                final long id = ind.getID();
                if(id < currIDStart){
                    if(isInInit && DEBUG){
                        System.out.println("DEBUG: Skipping ID " + id + " vs " + currIDStart);
                    }
                    continue;
                }
                myPool.add(ind);
            }
            
            if(isInInit && noInits > myPool.size()){
                System.err.println("WARNING: Gotten " + noInits + " inits but returning only " + myPool.size() + " pool size is: " + pool.getCurrentPoolSize());
            }
            
            final List<T> merged = comm.synchronizePool(myPool, myID, pool.getPoolSize(), currIDStart);
            if(!doNiching){
                pool.replacePoolContent(merged);
            } else {
                final List<Niche> niches = new ArrayList<>();
                for (final T t : merged) {
                    // replace pool...
                    final Niche niche = nicheComp.computeNiche(t);
                    niches.add(niche);
                }
                pool.replacePoolContent(merged, niches);
            }
        }
        
        if(isInInit){
            return false;
        }
                    
        if (DEBUG) {
            System.out.println("DEBUG: Merging done. Getting new chunk of glob opt tasks.");
        }

        short contactAttempts = 0;
        boolean chunkNotAcquired = true;
        while(chunkNotAcquired){
                    
            final Tuple3D<RMICodes.JOBSTATE, Long, Integer> tup = comm.getGlobOptChunk(maxTasks, myID);
            if (null != tup.getObject1()) switch (tup.getObject1()) {
                case CONTINUE:
                    // all fine :-)
                    chunksGotten++;
                    this.currIDStart = tup.getObject2();
                    this.tasksGotten = tup.getObject3();
                    if (DEBUG) {
                        System.out.println("DEBUG: Chunk request (II): all fine. " + currIDStart + " " + tasksGotten);
                    }
                    return false;
                case FINISH:
                    // all exhausted
                    if (DEBUG) {
                        System.out.println("DEBUG: Chunk request (II): empty.");
                    }
                    return true;
                case WAITING:
                    // do wait
                    if (DEBUG) {
                        System.out.println("DEBUG: Chunk request (II): wait.");
                    }
                    break;
                default:
                    // unknown status for server-proxy communication
                    throw new Exception("Unknown status " + tup.getObject1() + " " + tup.getObject2() + "  " + tup.getObject3() + " for server-proxy communication.");
            }
            
            contactAttempts++;
            if(contactAttempts >= MAXCONTACTATTEMPTS){
                System.err.println("Tried " + MAXCONTACTATTEMPTS + " to contact server w/o success. Finishing up... This is a problem... Please debug...");
                return true;
            }
            
            // sleep
            waitTimeMS += sleepTime;
            Thread.sleep(sleepTime);
        }
        
        System.err.println("How did you end up here in ThreadingBackend?!");
        return false;
    }
    
    private <V extends Copyable> void doXXX(final int threads, final TaskFactory<E,T,V> tasker, final ObjectCache<V> cache,
            final List<Task<T>> initTasks){

        // start thread pool
        final ExecutorService threadpool = Executors.newFixedThreadPool(threads);
        
        // do it for all
        for(int i = 0; i < initTasks.size(); i++){

            final long id = initTasks.get(i).getDummyAnswer(myID).getResult().getID();
            l.debug("DEBUG: Starting task " + id);

            assert(tasker != null);
            assert(history != null);
            assert(cache != null);
            threadpool.submit(tasker.createTask(pool,history,cache.getOriginalEntry(),useObjCache,cache,
                    doNiching, nicheComp, nicheCompCache, id));
        }
        
        threadpool.shutdown();

        /*
         * wait for the last submitted task to be done before returning to get the next chunk
         */
        try{
            // this will block until there is an actual result.
            threadpool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e1) {
            l.error("Threadpool reached wallclock limit. This should really NEVER happen! ", e1);
            System.err.println("ERROR: Can't wait for the last task of an init set to finish. "
                    + e1.toString() + ". Continuing on own risk now.");
        } catch (CancellationException e2) {
            System.err.println("ERROR: Can't wait for the last task of an init set to finish. "
                    + "It got cancelled." + e2.toString() + ". Continuing now.");
        }
    }
    
    private <V extends Copyable> void doXXX(final int threads, final TaskFactory<E,T,V> tasker, final ObjectCache<V> cache,
            final int offset, final int chunkSize){

        assert(chunkSize > 0);
        
        // start thread pool
        final ExecutorService threadpool = Executors.newFixedThreadPool(threads);
        
        // do it for all
        for(int i = offset; i < (chunkSize + offset); i++){

            l.debug("DEBUG: Starting task " + i);

            assert(tasker != null);
            assert(history != null);
            assert(cache != null);
            threadpool.submit(tasker.createTask(pool,history,cache.getOriginalEntry(),useObjCache,cache,
                    doNiching, nicheComp, nicheCompCache, (long) i));
        }
        
        threadpool.shutdown();

        /*
         * wait for the last submitted task to be done before returning to get the next chunk
         */
        try{
            // this will block until there is an actual result.
            threadpool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e1) {
            l.error("Threadpool reached wallclock limit. This should really NEVER happen! ", e1);
            System.err.println("ERROR: Can't wait for the last task of a set to finish. "
                    + e1.toString() + ". Continuing on own risk now.");
        } catch (CancellationException e2) {
            System.err.println("ERROR: Can't wait for the last task of a set to finish. "
                    + "It got cancelled." + e2.toString() + ". Continuing now.");
        }
    } 
    
    
    private static class GlobTaskFactory<U, W extends Optimizable<U>> implements TaskFactory<U,W,GenericGlobalOptimization<U,W>> {

        @Override
        public Runnable createTask(final GenericPool<U, W> pool, final GenericHistory<U, W> history, final GenericGlobalOptimization<U,W> refStuff, final boolean useCache, final ObjectCache<GenericGlobalOptimization<U,W>> cache, 
                final boolean doNiching, final NicheComputer<U, W> nicheComp, final ObjectCache<NicheComputer<U, W>> nicheCompCache, final long taskID) {
            
            return new Runnable() {
            
            @Override
            public void run() {
                
                final boolean isInit = ISINIT;
                
                if(useCache && !DEBUG){
                    try{
                    l.debug("Trying to use cached object as helper for global opt...");
                    final Tuple<Boolean,GenericGlobalOptimization<U,W>> stuff = cache.getUnusedEntry();
                    final GenericGlobalOptimization<U,W> cached = stuff.getObject2();
                    runme(cached,pool,history,taskID, isInit);
                    l.debug("Returning helper for global opt to cache...");
                    stuff.setObject1(false);
                    } catch(Throwable e){
                        e.printStackTrace(System.err);
                    }
                } else if(!useCache && !DEBUG){
                    l.debug("Trying to use new object as helper for global opt...");
                    final GenericGlobalOptimization<U,W> helpers = cache.getOriginalEntry().copy();
                    runme(helpers,pool,history,taskID, isInit);
                } else{
                    try{
                        l.debug("Trying to use new object as helper for global opt (try/catch)...");
                        final GenericGlobalOptimization<U,W> helpers = cache.getOriginalEntry().copy();
                        runme(helpers,pool,history,taskID, isInit);
                    } catch(Throwable t){
                        t.printStackTrace(System.err);
                    }
                }
            }
            
            private void runme(final GenericGlobalOptimization<U,W> helper, final GenericPool<U,W> pool, final GenericHistory<U,W> history, final long taskID,
                    final boolean isInit){
                
                l.debug("Starting global opt for " + taskID);
                final List<W> parents = pool.getParents();
                if(l.isDebugEnabled()){
                    l.debug("Took geometries " + parents.get(0).getID()
                            + " and " + parents.get(1).getID() + " out.");
                }
                
                final W child = helper.globalOptimization(taskID, parents.get(0), parents.get(1));
                
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
                    boolean accepted = false;
                    if(doNiching){
                        final Tuple<Boolean,NicheComputer<U,W>> comp = nicheCompCache.getUnusedEntry();
                        final Niche niche = comp.getObject2().computeNiche(child);
                        comp.setObject1(false);
                        if(isInit){
                            pool.addIndividualForced(child, niche, child.getFitness());
                            accepted = true;
                        } else {
                            accepted = pool.addIndividual(child, niche, child.getFitness());
                        }
                    }  else {
                        if(isInit){
                            pool.addIndividualForced(child,  child.getFitness());
                            accepted = true;
                        } else {
                            accepted = pool.addIndividual(child, child.getFitness());
                        }
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
}
