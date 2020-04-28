/**
Copyright (c) 2013-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
              2017, J. M. Dieterich and B. Hartke
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
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import org.ogolem.generic.Optimizable;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.NicheComputer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A generic thread dispatcher.
 * @author Johannes Dieterich
 * @version 2017-03-18
 */
public class GenericThreadDispatcher<E,T extends Optimizable<E>,V extends Cloneable> {
    
    private static final Logger l = LoggerFactory.getLogger(GenericThreadDispatcher.class);
    private final int threads;
    private final int offset;
    private final int iterations;
    private final int subsToWait;
    private final GenericPool<E,T> pool;
    private final GenericHistory<E,T> history;
    private final boolean useObjCache;
    private final ObjectCache<V> cache;
    private final boolean doNiching;
    private final NicheComputer<E,T> nicheComp;
    private final ObjectCache<NicheComputer<E,T>> nicheCompCache;
    private final TaskFactory<E,T,V> tasker;
    
    public GenericThreadDispatcher(final int noOfThreads, final GenericPool<E,T> pool,
            final GenericHistory<E,T> hist, final ObjectCache<V> cache,
            final int offset, final int steps, final TaskFactory<E,T,V> tasker,
            final int subsToWait, final boolean doNiching, final NicheComputer<E,T> nicheComp){
        
        assert(tasker != null);
        assert(hist != null);
        assert(cache != null);
        assert(noOfThreads > 0);
        
        this.threads = noOfThreads;
        this.offset = offset;
        this.iterations = steps;
        this.pool = pool;
        this.history = hist;
        this.tasker = tasker;
        
        // setup the object cache to avoid re-instantiating things for EVERY task
        if(cache != null){
            this.useObjCache = true;
            this.cache = cache;
        } else{
            this.useObjCache = false;
            this.cache = null;
        }
        this.subsToWait = subsToWait;
        
        this.doNiching = doNiching;
        if(doNiching){
            try{
                this.nicheComp = nicheComp;
                this.nicheCompCache = new ObjectCache<>(2*noOfThreads,nicheComp);
            } catch (Exception e){
                throw new RuntimeException(e);
            }
        } else {
            this.nicheCompCache = null;
            this.nicheComp = null;
        }
    }

    public void doAllTasks(){

        final ExecutorService threadpool = Executors.newFixedThreadPool(threads);

        // do it for all
        int countToWait = 0;
        for(int i = offset; i < (iterations + offset); i++){

            l.debug("DEBUG: Starting task " + i);

            assert(tasker != null);
            assert(history != null);
            assert(cache != null);
            final Future<?> future = threadpool.submit(tasker.createTask(pool,history,cache.getOriginalEntry(),useObjCache,cache,
                    doNiching, nicheComp, nicheCompCache, (long) i));
            countToWait++;
            
            if(countToWait >= subsToWait){
                /*
                 * wait for the last submitted task to be done before submiting fresh ones
                 * this is sufficient for the targeted use (reduce footprint in heap)
                 */
                try{
                    // this will block until there is an actual result.
                    future.get();
                } catch(InterruptedException e1){
                    System.err.println("ERROR: Can't wait for the last task of a set to finish. "
                            + e1.toString() + ". Continuing on own risk now.");
                } catch(CancellationException e2){
                    System.err.println("ERROR: Can't wait for the last task of a set to finish. "
                            + "It got cancelled." + e2.toString() + ". Continuing now.");
                } catch(ExecutionException e3){
                    System.err.println("ERROR: The globopt task threw an exception. Providing stack trace now.");
                    e3.printStackTrace(System.err);
                    System.err.println("ERROR: Continuing on own risk now.");
                }
                countToWait = 0;
                /*
                 * make sure that we actually wanting to continue, by checking if a file named
                 * STOP exists in our running directory. If so break out, be done.
                 */
                final File f = new File("STOP");
                if(f.exists()){
                    // give some more output
                    l.info("Manual stopping intervention received, breaking loop and shutting down thread pool.");
                    break;
                }
            }
            
            if(pool.acceptableFitnessReached()){break;}
        }

        // shut the pool down
        threadpool.shutdown();

        // wait for all threads to finish
        try{
            threadpool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        } catch(InterruptedException e){
            l.error("Threadpool reached wallclock limit. This should really NEVER happen! ", e);
        }
    }    
}
