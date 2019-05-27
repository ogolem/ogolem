/**
Copyright (c) 2013-2014, J. M. Dieterich
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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.Output;
import org.ogolem.generic.Configuration;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericGlobalOptimization;
import org.ogolem.generic.GenericInitializer;
import org.ogolem.generic.IndividualReader;
import org.ogolem.generic.Optimizable;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.GenericPoolConfig;
import org.ogolem.generic.genericpool.GenericPoolEntry;
import org.ogolem.generic.IndividualWriter;
import org.ogolem.generic.genericpool.NicheComputer;
import org.ogolem.generic.stats.GenericDetailStatistics;
import org.ogolem.helpers.Fortune;
import org.ogolem.helpers.Tuple;
import org.ogolem.io.OutputPrimitives;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A fully generic, shared-memory global optimization manager.
 * @author Johannes Dieterich
 * @version 2016-10-26
 */
public class GenericOGOLEMOptimization<E,T extends Optimizable<E>> {
    
    private static final Logger l = LoggerFactory.getLogger(GenericOGOLEMOptimization.class);
    public static final int DEFAULTSUBSTOWAIT = 1000;
    private final Configuration<E,T> config;
    private final GenericPool<E,T> pool;
    private final GenericHistory<E,T> history;
    private final IndividualWriter<T> writer;
    private final IndividualReader<T> reader;
    private final NicheComputer<E,T> nicheComp;
    private final boolean doNiching;
    private final String outFile;
    private final String outFolder;
    private final int subsToWait;
    
    public GenericOGOLEMOptimization(final Configuration<E,T> config, 
            final GenericPool<E,T> pool, final GenericHistory<E,T> history, final String outFile,
            final String outFolder, final IndividualWriter<T> writer, final IndividualReader<T> reader,
            final int subsToWait){
        this.config = config;
        this.pool = pool;
        this.history = history;
        this.outFile = outFile;
        this.outFolder = outFolder;
        this.writer = writer;
        this.reader = reader;
        this.subsToWait = subsToWait;
        NicheComputer<E,T> tmp = null;
        try{
            tmp = config.getNicheComputer();
        } catch(Exception e){
            throw new RuntimeException(e);
        }
        this.nicheComp = tmp;
        this.doNiching = (nicheComp != null);
    }
    
    public T globOpt(final int noThreads){
        
        l.debug("Entering generic globopt. BRACE YOURSELF!");
        
        final long offsetTime = System.currentTimeMillis();
        
        final String[] header = Output.getHeader();
        final String[] generific = new String[]{
            "####################################################################",
            "####################################################################",
            "                  AND OUR PROBLEM OF THE DAY IS:",
            "                  " + pool.getExample().getClass().getCanonicalName(),
            "####################################################################",
            "####################################################################",
        };
        final List<String> confDat = config.getFormattedConfig();
        try{
            OutputPrimitives.writeOut(outFile, header, true);
            OutputPrimitives.writeOut(outFile, generific, true);
            OutputPrimitives.writeOut(outFile, confDat, true);
        } catch(IOException e){
            System.err.println("Failed to write header or config to output file. Ignoring...");
            e.printStackTrace(System.err);
        }

        GenericPoolConfig<E,T> poolConfig = null;
        final GenericInitializer<E,T> initializer;
        final GenericFitnessFunction<E,T> fitness;
        final GenericGlobalOptimization<E,T> globopt;
        ObjectCache<GenericInitializer<E,T>> initCache = null;
        ObjectCache<GenericGlobalOptimization<E,T>> globCache = null;
        ObjectCache<Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>>> fitFuncCache = null;
        TaskFactory<E,T,GenericInitializer<E,T>> initTasks = null;
        TaskFactory<E,T,GenericGlobalOptimization<E,T>> globTasks = null;
        int noGlobSteps = -1;
        boolean enableDetailedStats = false;
        try{
            poolConfig = config.getGenericPoolConfig();
            initializer = config.getInitializer();
            globopt = config.getGlobalOptimization();
            initTasks = config.getInitFactory();
            globTasks = config.getGlobalFactory();
            noGlobSteps = config.getNumberOfGlobalSteps();
            fitness = config.getFitnessFunction();
            initCache = new ObjectCache<>(noThreads,initializer);
            globCache = new ObjectCache<>(noThreads,globopt);
            final Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>> tup = new Tuple<>(reader,fitness);
            fitFuncCache = new ObjectCache<>(noThreads,tup);
            enableDetailedStats = config.wantsDetailedStats();
        } catch(final Exception e){
            System.err.println("Exception while trying to get stuff from config or initializing caches. Failing.");
            e.printStackTrace(System.err);
            System.exit(42);
        }
        
        assert(poolConfig != null);
        assert(initCache != null);
        assert(globCache != null);
        assert(fitFuncCache != null);
        assert(initTasks != null);
        assert(globTasks != null);
        assert(noGlobSteps >= 0);
        
        if(enableDetailedStats){
            GenericDetailStatistics.enableAllDetails();
        }
                
        final long initTime = System.currentTimeMillis();
        
        // first try to use seeds
        final String seedFolder = config.seedFolder();
        int noSeeded = 0;
        if(seedFolder != null && !seedFolder.isEmpty()){
            // apparently we do want to use them
            TaskFactory<E,T,Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>>> seedTasks = null;
            try{
                seedTasks = new GenericSeedTask<>(seedFolder);
            } catch(Exception e){
                System.err.println("Failure to setup seeding initialization. Exception was: " + e.getMessage());
                e.printStackTrace(System.err);
                System.exit(-20);
            }
            assert(seedTasks != null);
            final GenericThreadDispatcher<E,T,Tuple<IndividualReader<T>,GenericFitnessFunction<E,T>>> initSMPPool
                    = new GenericThreadDispatcher<>(noThreads, pool,history,fitFuncCache,0,Math.min(GenericSeedTask.getNoOfSeeds(),
                            poolConfig.getPoolSize()),seedTasks,subsToWait,doNiching,nicheComp);
            initSMPPool.doAllTasks();
            noSeeded = Math.min(GenericSeedTask.getNoOfSeeds(),poolConfig.getPoolSize());
        }
        
        final GenericThreadDispatcher<E,T,GenericInitializer<E,T>> initSMPPool =
                new GenericThreadDispatcher<>(noThreads, pool,history,initCache,noSeeded,poolConfig.getPoolSize()-noSeeded,initTasks,
                subsToWait,doNiching,nicheComp);
        initSMPPool.doAllTasks();
        final List<String> initFitnesses = new ArrayList<>();
        initFitnesses.add("#-----------------------------------------------------------");
        initFitnesses.add("#");
        initFitnesses.add("# Individuals in initial pool coming, fitness units unknown.");
        initFitnesses.add("#");
        initFitnesses.add("# pool position        individual id                 fitness");
        
        final List<String> poolCont = pool.getFormattedPool();
        for(final String s : poolCont){
            initFitnesses.add(s);
        }
        
        initFitnesses.add("#");
        initFitnesses.add("#-----------------------------------------------------------");
        try{
            OutputPrimitives.writeObjToBinFile(outFolder + File.separator + "postInitializationPool.bin", pool);
            OutputPrimitives.writeOut(outFile, initFitnesses, true);
            
            // print the individuals out as well
            int c = 0;
            for(final GenericPoolEntry<E,T> entry : pool){
                final T ind = entry.getIndividual();
                final String file = outFolder + File.separator + "initrank" + c + "individual" + ind.getID();
                writer.writeIndividual(ind, file);
                c++;
            }
        } catch (Exception e){
            System.err.println("Failure in writing initial pool fitnesses out. Ignoring.");
            e.printStackTrace(System.err);
        }
        
        final long poolFillTime = System.currentTimeMillis();
        
        // make sure that the pool is indeed as full as expected and otherwise throw a runtime exception
        final int currGenPoolSize = pool.getCurrentPoolSize();
        if(currGenPoolSize != pool.getPoolSize() && currGenPoolSize >= 2){
            System.out.println("INFO: Genetic pool should be " + pool.getPoolSize() + " but is "
                    + currGenPoolSize + ". This seems to be a problem where a lot of initial guesses collapse to the same solution.");
        } else if(currGenPoolSize != pool.getPoolSize() && currGenPoolSize < 2){
            throw new RuntimeException("INFO: Genetic pool should be " + pool.getPoolSize() + " but is "
                    + currGenPoolSize + ". After init, there must be at least two individuals in the pool. "
                    + "Most likely reason for this is a failure in the init. Rerun with higher debug level to "
                    + "see any exceptions/error that occur.");
        }
        
        final GenericThreadDispatcher<E,T,GenericGlobalOptimization<E,T>> globSMPPool =
                new GenericThreadDispatcher<>(noThreads, pool,history,globCache,poolConfig.getPoolSize(),noGlobSteps,globTasks,
                subsToWait,doNiching,nicheComp);
        globSMPPool.doAllTasks();
        
        history.flushRecords();
        history.writeTotalStats();
        
        final List<String> finalFitnesses = new ArrayList<>();
        finalFitnesses.add("#-----------------------------------------------------------");
        finalFitnesses.add("#");
        finalFitnesses.add("# Individuals in final pool coming, fitness units unknown.");
        finalFitnesses.add("#");
        finalFitnesses.add("# pool position        individual id                 fitness");
        
        final List<String> poolCont2 = pool.getFormattedPool();
        for(final String s : poolCont2){
            finalFitnesses.add(s);
        }
        
        finalFitnesses.add("#");
        finalFitnesses.add("#-----------------------------------------------------------");
        try{
            OutputPrimitives.writeObjToBinFile(outFolder + File.separator + "finalPool.bin", pool);
            OutputPrimitives.writeOut(outFile, finalFitnesses, true);
            
            // print the individuals out as well
            int c = 0;
            for(final GenericPoolEntry<E,T> entry : pool){
                final T ind = entry.getIndividual();
                final String file = outFolder + File.separator + "rank" + c + "individual" + ind.getID();
                writer.writeIndividual(ind,file);
                c++;
            }
        } catch (Exception e){
            System.err.println("Failure in writing initial pool fitnesses out. Ignoring.");
            e.printStackTrace(System.err);
        }
        
        final long globoptTime = System.currentTimeMillis();
        
        // timing stuff
        final String[] timing = new String[7];
        timing[0] = "Timing information for this generified globopt run coming:";
        timing[1] = "   startup time:        " + (initTime-offsetTime)/1000. + " s";
        timing[2] = "   pool filling:        " + (poolFillTime-initTime)/1000. + " s";
        timing[3] = "   global optimization: " + (globoptTime-poolFillTime)/1000. + " s";
        timing[4] = "   total time:          " + (globoptTime-initTime)/1000. + " s";
        timing[5] = "Thanks for running OGOLEM for optimization!";
        timing[6] = Fortune.randomFortune();
        
        try{
            OutputPrimitives.writeOut(outFile, timing, true);
        } catch(IOException e){
            System.err.println("Failed to write timing to output file. Ignoring...");
            e.printStackTrace(System.err);
        }
        
       if(enableDetailedStats){
           final List<String> detailedStats = GenericDetailStatistics.getOutput();
           try{
                OutputPrimitives.writeOut(outFile, detailedStats, true);
            } catch(IOException e){
                System.err.println("Failed to write detailed statistics to output file. Ignoring...");
                e.printStackTrace(System.err);
            }
       }
        
        // serialize the final pool
        try{
            OutputPrimitives.writeObjToBinFile(outFolder + File.separator + "pool.bin", pool);
        } catch(IOException e){
            System.err.println("Failed to write binary final pool to file. Ignoring...");
            e.printStackTrace(System.err);
        }
        
        // return the best individual
        return pool.getIndividualAtPosition(0);
    }
}
