/**
Copyright (c) 2019-2020, J. M. Dieterich and B. Hartke
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
 * @version 2020-03-22
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
        
        // run advanced CD benchmark
        final AdvPairwiseCDBenchmark advCDBench = new AdvPairwiseCDBenchmark();
        runOne(advCDBench, 1000);
        
        // run advanced CD benchmark (check only)
        final AdvPairwiseCDCheckOnlyBenchmark advCDBench2 = new AdvPairwiseCDCheckOnlyBenchmark();
        runOne(advCDBench2, 1000);
        
        // run aligning benchmark
        final AligningBench alignBench = new AligningBench();
        runOne(alignBench, 1000);

        // run angle benchmark
        final AngleBench angleBench = new AngleBench();
        runOne(angleBench, 10000);
        
        // run angle 2 benchmark
        final AngleBench2 angleBench2 = new AngleBench2();
        runOne(angleBench2, 10000);
        
        // run dihedral benchmark
        final DihedralBench dihedralBench = new DihedralBench();
        runOne(dihedralBench, 10000);

	// run matrix multiplication benchmarks
        final MatMultBenchmark matMultBench = new MatMultBenchmark(3, 3, 3);
        runOne(matMultBench, 20000);
        final MatMultBenchmark matMultBench0 = new MatMultBenchmark(3, 55, 3);
        runOne(matMultBench0, 10000);
        final MatMultBenchmark matMultBench1 = new MatMultBenchmark(3, 1000, 3);
        runOne(matMultBench1, 1000);
        final MatMultBenchmark matMultBench2 = new MatMultBenchmark(256, 256, 256);
        runOne(matMultBench2, 1);

        // run specialized 3x3 matrix multiplication benchmarks
        final Mat3x3MultBenchmark mat3x3MultBench = new Mat3x3MultBenchmark(3);
        runOne(mat3x3MultBench, 20000);
        final Mat3x3MultBenchmark mat3x3MultBench0 = new Mat3x3MultBenchmark(55);
        runOne(mat3x3MultBench0, 10000);
        final Mat3x3MultBenchmark mat3x3MultBench1 = new Mat3x3MultBenchmark(1000);
        runOne(mat3x3MultBench1, 1000);

        // run a Norway packing mutation benchmark
        final NorwayPackingLJBench norwayLJ38Bench = new NorwayPackingLJBench(38);
        runOne(norwayLJ38Bench, 5);
        final NorwayPackingLJBench norwayLJ55Bench = new NorwayPackingLJBench(55);
        runOne(norwayLJ55Bench, 2);
        
        // run LJ benchmarks
        final LJFFEnergyBench ljEnergyBench = new LJFFEnergyBench();
        runOne(ljEnergyBench, 100);
        final LJFFGradientBench ljGradientBench = new LJFFGradientBench();
        runOne(ljGradientBench, 100);
        
        // run mixed LJ benchmarks
        final MixedLJFFEnergyBench mixedljEnergyBench = new MixedLJFFEnergyBench();
        runOne(mixedljEnergyBench, 100);
        final MixedLJFFGradientBench mixedljGradientBench = new MixedLJFFGradientBench();
        runOne(mixedljGradientBench, 100);
        
        // run adaptive LJ benchmarks
        final AdaptiveLJFFEnergyBench adaptiveljEnergyBench = new AdaptiveLJFFEnergyBench();
        runOne(adaptiveljEnergyBench, 100);
        final AdaptiveLJFFGradientBench adaptiveljGradientBench = new AdaptiveLJFFGradientBench();
        runOne(adaptiveljGradientBench, 100);

	// run TIP3P benchmarks
        final TIP3PEnergyBench tip3pEBench = new TIP3PEnergyBench();
        runOne(tip3pEBench, 10);
        final TIP3PGradientBench tip3pGBench = new TIP3PGradientBench();
        runOne(tip3pGBench, 10);
        
        // run TIP4P benchmarks
        final TIP4PEnergyBench tip4pEBench = new TIP4PEnergyBench();
        runOne(tip4pEBench, 10);
        final TIP4PGradientBench tip4pGBench = new TIP4PGradientBench();
        runOne(tip4pGBench, 10);

	// run small water TTM3F benchmark
        final WaterTTM3FSmallBenchmark ttm3fSmall = new WaterTTM3FSmallBenchmark();
        runOne(ttm3fSmall, 10);
        
        // run large water TTM3F benchmark
        final WaterTTM3FLargeBenchmark ttm3fLarge = new WaterTTM3FLargeBenchmark();
        runOne(ttm3fLarge, 10);
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
