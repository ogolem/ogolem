/**
Copyright (c) 2012-2014, J. M. Dieterich
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
package org.ogolem.generic.genericpool;

import org.ogolem.generic.IndividualWriter;
import java.io.File;
import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;
import org.ogolem.generic.Optimizable;

/**
 * Configuration for the generic pool.
 * @author Johannes Dieterich
 * @version 2016-04-19
 */
public class GenericPoolConfig<E,T extends Optimizable<E>> implements Serializable {
    
    private static final long serialVersionUID = (long) 20150803;
    
    private boolean serializeAfterNewBest = false;
    private boolean writeEveryAdd = false;
    private boolean beSilent = false;
    private double acceptableFitness = Double.NEGATIVE_INFINITY;
    private int poolSize = 1000;
    private int addsToSerial = 1000;
    private int addsToStats = 10000;
    private String interBinFile = "IntermediateGenericPool.bin";
    
    // helper objects and (if applicable) configuration
    private DiversityChecker<E,T> diversity = null;
    private boolean doNiching = false;
    private Nicher<E,T> nicher = null;
    private GenericStatistics stats = null;
    private IndividualWriter<T> writer = null;
    private ParentSelector<E,T> selector = null;
    
    public GenericPoolConfig(){
    }
    
    public GenericPoolConfig(final String[] conf, final String outFolder) throws Exception {
        
        if(!conf[0].trim().equalsIgnoreCase("<GENERICPOOL>") || !conf[conf.length-1].trim().equalsIgnoreCase("</GENERICPOOL>")){
            throw new Exception("Generic pool configuration not properly embedded into <GENERICPOOL> tags.");
        }
        for(int i = 1; i < conf.length-1; i++){
            final String s = conf[i].trim();
            // remove comments
            if(s.startsWith("#") || s.startsWith("//")){continue;}
            
            // parse options
            if(s.startsWith("PoolSize=")){
                final String sx = s.substring(9).trim();
                final int po = Integer.parseInt(sx);
                if(po <= 0){
                    throw new Exception("Empty or negative pool not allowed.");
                }
                setPoolSize(po);
            } else if(s.equalsIgnoreCase("WriteEveryAdd")){
                writeEveryAdd = true;
            } else if(s.startsWith("AcceptableFitness=")){
                final String sx = s.substring(18).trim();
                final double da = Double.parseDouble(sx);
                setAcceptableFitness(da);
            } else if(s.startsWith("AddsToStats=")){
                final String sx = s.substring(12).trim();
                final int addsToStats = Integer.parseInt(sx);
                setAddsToStats(addsToStats);
            } else if(s.equalsIgnoreCase("BeSilent")){
                beSilent();
            } else{
                throw new Exception("Unknown option " + s);
            }
            
            //TODO more? what about adding parsing objects?
        }
        
        stats = new GenericStatistics(outFolder + File.separator + outFolder + ".log", addsToStats);
        diversity = new GenericDiversityCheckers.FitnessDiversityChecker<>(1E-5);
//        selector = new GenericParentSelectors.FitnessWeightedParentSelector<>(0.5,false); // bxh testing...
        selector = new GenericParentSelectors.FitnessWeightedParentSelector<>(0.5,0.5,false,false,0.0);
    }

    /*
     * GETTERS AND SETTERS
     */
    public double getAcceptableFitness() {
        return acceptableFitness;
    }

    public final void setAcceptableFitness(double acceptableFitness) {
        if(acceptableFitness != acceptableFitness) throw new IllegalArgumentException("Please specify a reasonable fitness");
        this.acceptableFitness = acceptableFitness;
    }

    public int getAddsToSerial() {
        return addsToSerial;
    }

    public void setAddsToSerial(final int addsToSerial) {
        if(addsToSerial <= 0){
            throw new IllegalArgumentException("AddsToSerial must be positive!");
        }
        this.addsToSerial = addsToSerial;
    }

    public int getAddsToStats() {
        return addsToStats;
    }

    public final void setAddsToStats(final int addsToStats) {
        if(addsToSerial <= 0){
            throw new IllegalArgumentException("AddsToStats must be positive!");
        }
        this.addsToStats = addsToStats;
    }

    public String getInterBinFile() {
        return interBinFile;
    }

    public void setInterBinFile(String interBinFile) {
        this.interBinFile = interBinFile;
    }

    public DiversityChecker<E,T> getDiversityChecker() {
        return diversity;
    }

    public void setDiversityChecker(DiversityChecker<E,T> diversity) {
        this.diversity = diversity;
    }

    public boolean doNiching() {
        return doNiching;
    }

    public void setDoNiching(boolean doNiching) {
        this.doNiching = doNiching;
    }

    public Nicher<E,T> getNicher() {
        return nicher;
    }

    public void setNicher(Nicher<E,T> nicher) {
        this.nicher = nicher;
    }

    public int getPoolSize() {
        return poolSize;
    }

    public final void setPoolSize(int poolSize) {
        this.poolSize = poolSize;
    }

    public boolean serializeAfterNewBest() {
        return serializeAfterNewBest;
    }

    public void setSerializeAfterNewBest(boolean serializeAfterNewBest) {
        this.serializeAfterNewBest = serializeAfterNewBest;
    }

    public boolean writeEveryAdd() {
        return writeEveryAdd;
    }

    public void setWriteEveryAdd(boolean writeEveryAdd) {
        this.writeEveryAdd = writeEveryAdd;
    }
    
    public GenericStatistics getStats() {
        return stats;
    }

    public void setStats(GenericStatistics stats) {
        this.stats = stats;
    }

    public IndividualWriter<T> getWriter() {
        return writer;
    }

    public void setWriter(IndividualWriter<T> writer) {
        this.writer = writer;
    }

    public ParentSelector<E,T> getSelector() {
        return selector;
    }

    public void setSelector(ParentSelector<E,T> selector) {
        this.selector = selector;
    }
    
    public final void beSilent(){
        this.beSilent = true;
    }
    
    public boolean shouldBeSilent(){
        return beSilent;
    }
    
    public List<String> getMyConfig(){
        
        final List<String> out = new LinkedList<>();
        out.add("##################################");
        out.add("GENERIC POOL CONFIGURATION");
        out.add("");
        out.add("\tpool size " + poolSize);
        out.add("\tniching enabled? " + doNiching);
        out.add("\tsilent enabled?" + beSilent);
        out.add("");
        out.add("##################################");
        // XXX extend
        
        return out;
    }
}
