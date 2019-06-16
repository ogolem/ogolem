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
package org.ogolem.generic;

import java.util.LinkedList;
import java.util.List;
import org.ogolem.generic.generichistory.GenericHistoryConfig;
import org.ogolem.generic.genericpool.GenericPoolConfig;
import org.ogolem.generic.threading.GenericGlobOptTask;
import org.ogolem.generic.threading.GenericInitTask;
import org.ogolem.generic.threading.TaskFactory;

/**
 * An abstract default config. Essentially, a bare-bones starting point if you
 * want to solve your own optimization problem using the ogolem kernels.
 * @author Johannes Dieterich
 * @version 2015-05-26
 * @param <E> the optimizable quantity, e.g., Double or Molecule
 * @param <T> the actual optimizable class
 */
public abstract class AbstractDefaultConfig<E,T extends Optimizable<E>> implements Configuration<E,T>{
    
    private static final long serialVersionUID = (long) 20140830;
    
    /**
     * where to find the seeding individuals. Also works as an identifier whether
     * to use seeding or not.
     */
    protected String seedingFolder = "N/A";
    
    /**
     * whether or not detailed stats are wanted at the expense of an additional
     * synchronization point (potential slowdown). only applicable for the SMP
     * case.
     */
    protected boolean enableDetailedStats = false;
    
    /**
     * the number of global steps to carry out in this optimization run.
     */
    protected int noOfGlobalSteps = 9900;

    /**
     * the generic history configuration. a default configuration.
     */
    protected GenericHistoryConfig hisConf = new GenericHistoryConfig();
    
    /**
     * the generic pool configuration. a default configuration.
     */
    protected GenericPoolConfig<E,T> poolConfig = new GenericPoolConfig<>();
    
    /**
     * an example object (fully setup!) for, e.g., initialization purposes
     */
    protected final T example;
    
    protected AbstractDefaultConfig(final T example){
        this.example = example;
    }
    
    @Override
    public List<String> getFormattedConfig() {
        
        final LinkedList<String> output = new LinkedList<>();
        output.add("############################################################");
        output.add("############################################################");
        output.add("");
        output.add("                     GLOBAL CONFIGURATION");
        output.add("");
        output.add("############################################################");
        output.add("############################################################");
        output.add("");
        output.add("GENERIC POOL CONFIGURATION:");
        final List<String> poolC = poolConfig.getMyConfig();
        poolC.forEach((s) -> {
            output.add(s);
        });
        output.add("");
        output.add("GENERIC HISTORY CONFIGURATION:");
        output.add(hisConf.getMyConfig());
        output.add("");
        output.add("############################################################");
        output.add("");
        output.add("INDIVIDUAL INITIALIZATION:");
        final String ini = (seedFolder() == null) ? "no seeding taking place" : "seeding from folder " + seedFolder();
        output.add(ini);
        output.add("");
        output.add("GLOBAL OPTIMIZATION:");
        output.add("number of global steps to finish " + noOfGlobalSteps);
        output.add("detailled statistics enabled (threading only) " + enableDetailedStats);        
        
        return output;
    }

    @Override
    public final String seedFolder() {
        return (seedingFolder == null || seedingFolder.equalsIgnoreCase("N/A")) ? null : seedingFolder;
    }

    @Override
    public abstract GenericFitnessFunction<E, T> getFitnessFunction() throws Exception;

    @Override
    public GenericHistoryConfig getGenericHistoryConfig() throws Exception {
        if(hisConf == null){throw new RuntimeException("Do not set the generic history configuration to null.");}
        if(hisConf.recordsToASCII <= 0 || hisConf.recordsToSerial <= 0){throw new RuntimeException("Records to serial and records to ASCII must be positive!");}
        if(hisConf.asciiAppend == null || hisConf.binOut == null){throw new RuntimeException("Both the ascii log and the binary log must not be NULL!");}
        
        return hisConf;
    }

    @Override
    public GenericPoolConfig<E, T> getGenericPoolConfig() throws Exception {
        if(poolConfig == null){throw new RuntimeException("Do not set the pool configuration to null.");}
        if(poolConfig.getPoolSize() < 2){throw new RuntimeException("The pool size must be two or more!");}
        if(poolConfig.doNiching() && poolConfig.getNicher() == null){
            throw new RuntimeException("Niching is enabled but no niche computer or nicher exists!");
        }
        
        return poolConfig;
    }

    @Override
    public int getNumberOfGlobalSteps() throws Exception {
        if(noOfGlobalSteps < 0){throw new RuntimeException("Number of global steps not allowed to be less than 0.");}
        return noOfGlobalSteps;
    }

    @Override
    public abstract GenericInitializer<E, T> getInitializer() throws Exception;

    @Override
    public <V extends GenericInitializer<E, T>> TaskFactory<E, T, V> getInitFactory() throws Exception {
        return new GenericInitTask<>();
    }

    @Override
    public abstract GenericGlobalOptimization<E, T> getGlobalOptimization() throws Exception;

    @Override
    public <V extends GenericGlobalOptimization<E, T>> TaskFactory<E, T, V> getGlobalFactory() throws Exception {
        return new GenericGlobOptTask<>();
    }

    @SuppressWarnings("unchecked")
    @Override
    public T getExample() throws Exception {
        return (T) example.clone();
    }

    @Override
    public abstract IndividualReader<T> getReader() throws Exception;

    @Override
    public abstract IndividualWriter<T> getWriter(final String folder) throws Exception;

    @Override
    public boolean wantsDetailedStats() throws Exception {
        return enableDetailedStats;
    }
}
