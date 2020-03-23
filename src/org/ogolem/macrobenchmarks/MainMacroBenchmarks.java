/**
Copyright (c) 2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.macrobenchmarks;

import org.ogolem.io.InputPrimitives;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Designed to run macrobenchmarks, i.e., full global optimizations.
 * @author Johannes Dieterich
 * @version 2020-03-22
 */
public class MainMacroBenchmarks {
    
    private static final Logger LOG = LoggerFactory.getLogger(MainMacroBenchmarks.class);
            
    public static void run(final String[] args){
        
        if(args[0].equalsIgnoreCase("help")){
            System.out.println("This is the macro benchmarking functionality.");
            System.out.println("Required arguments:");
            System.out.println(" * which optimization submodule these benchmarks target (-cluster and -adaptive supported)");            
            System.out.println(" * the input file (csv format)");
            System.out.println(" * the number of threads to be used for each global optimization");
	    System.out.println("Note that the execution must be from a directory that contains the paths in the input csv file.\n I.e., with the default csv files, this must be executed from the ogolem main directory.");
            System.exit(0);
        }
        
        if(args.length < 3){
            throw new RuntimeException("Must specify all three mandatory arguments.");
        }
        
        
        BenchmarkRunner runner = null;
        if(args[0].equalsIgnoreCase("-cluster")){
            LOG.info("Benchmarking cluster structure optimization.");
            runner = new ClusterBenchmarkRunner();
        } else if(args[0].equalsIgnoreCase("-adaptive")){
            LOG.info("Benchmarking adaptive parameter optimization.");
            runner = new AdaptiveBenchmarkRunner();
        } else {
            throw new RuntimeException("Unsupported optimization submodule " + args[0]);
        }
        
        assert(runner != null);
        
        final String refFile = args[1];
        final int noThreads = Integer.parseInt(args[2]);
        
        String[] csvData = null;
        try {
            csvData = InputPrimitives.readFileIn(refFile);
        } catch(Exception e){
            throw new RuntimeException("Failure to read csv file in.", e);
        }
        
        assert(csvData != null);
        
        int errOut = 0;
        for(final String line : csvData){
            
            if(line.trim().startsWith("#") || line.trim().startsWith("//")){
                // comment line, ignore
                continue;
            }
            
            LOG.debug("Working on line: " + line);
            
            try {
                final boolean success = runner.runBenchmark(line, noThreads);
                if(!success){
                    LOG.error("Benchmark " + line + " was NOT a success!");
                    errOut++;
                }
            } catch (Exception e){
                errOut++;
                e.printStackTrace(System.err);
            }
        }
        
        if(errOut != 0){
            throw new RuntimeException("Failure to successfully complete some benchmarks. Number is " + errOut);
        }
    }
}
