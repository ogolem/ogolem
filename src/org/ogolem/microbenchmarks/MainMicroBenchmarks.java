/**
Copyright (c) 2019, J. M. Dieterich and B. Hartke
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
package org.ogolem.microbenchmarks;

import org.ogolem.helpers.StatisticUtils;
import org.ogolem.helpers.Tuple;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Designed to run microbenchmarks of some performance-critical functionalities
 * of ogolem.
 * Benchmarked functionalities should be micro - hence fast to benchmark and
 * fundamental.
 * @author Johannes Dieterich
 * @version 2019-12-29
 */
public class MainMicroBenchmarks {
    
    private static final Logger LOG = LoggerFactory.getLogger(MainMicroBenchmarks.class);
    
    /*
    number of macro iterations for performance benchmarking
    */
    private static final int NOPERFMACROITERATIONS = 3;
        
    /*
    number of micro iterations for performance benchmarking
    */
    private static final int NOPERFMICROITERATIONS = 1000;
        
    /*
    this should be a threshold large enough to ensure the JIT kicks in and
    optimizes the functions. Typical rule of thumb is 1000 function calls
    will do that - so play it safe
    */
    private static final int NOWARMUPITERATIONS = 10000;
    
    public static void run(final String[] args){
        
        // run aligning benchmark
        final AligningBench alignBench = new AligningBench();
        runOne(alignBench, 1000);
        
    }
    
    /**
     * Run one micro benchmark
     * @param bench the micro benchmark
     * @param microBenchMultiplier multiplier for the micro perf iterations. Should be selected to yield a macro iteration larger than 1000 ms.
     */
    private static void runOne(final SingleMicroBenchmark bench, final int microBenchMultiplier){
        
        final String name = bench.name();
        
        // first do the warmup iterations to give JIT/JVM some time to equillibrate
        LOG.debug("Benchmark " + name + " warming up...");
        double sideEffectAvoider = 0.0;
        try{
            for(int i = 0; i < NOWARMUPITERATIONS; i++){
                sideEffectAvoider += bench.runSingle();
            }
        } catch (Exception e){
            System.err.println("ERROR: Failure in warm up - should never happen!");
            e.printStackTrace(System.err);
        }
        LOG.debug("Benchmark " + name + " side effect avoider (warmup) " + sideEffectAvoider);
        
        // now do the real performance run
        final double[] times = new double[NOPERFMACROITERATIONS];
        LOG.debug("Benchmark " + name + " starting...");        
        
        for(int macro = 0; macro < NOPERFMACROITERATIONS; macro++){
            sideEffectAvoider = 0.0;
            final long tStart = System.currentTimeMillis();
            try{
                for(int i = 0; i < NOPERFMICROITERATIONS*microBenchMultiplier; i++){
                    sideEffectAvoider += bench.runSingle();
                }
            } catch (Exception e){
                System.err.println("ERROR: Failure in perf run - should never happen!");
                e.printStackTrace(System.err);
            }
            
            long tEnd = System.currentTimeMillis();
            LOG.debug("Benchmark " + name + " side effect avoider (perf run) " + sideEffectAvoider);
            times[macro] = (tEnd-tStart);
        }
        
        final Tuple<Double,Double> meanStdDev = StatisticUtils.meanAndStdDev(times);
        
        LOG.info("Benchmark " + name + " took " + Math.round(meanStdDev.getObject1()) + " +/- " + Math.round(meanStdDev.getObject2()) + " ms per macro iteration\n"
            + "                         or " + Math.round(1000.0/meanStdDev.getObject1()*microBenchMultiplier*NOPERFMICROITERATIONS) + " calls per second on average");
    }
}
