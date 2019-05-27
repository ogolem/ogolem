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
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import org.ogolem.core.FixedValues;
import org.ogolem.generic.Optimizable;
import org.ogolem.io.OutputPrimitives;

/**
 * A generic genetic pool.
 * @author Johannes Dieterich
 * @version 2016-01-14
 */
public class GenericPool<E,T extends Optimizable<E>> implements Serializable, Iterable<GenericPoolEntry<E,T>> {
    
    private static final long serialVersionUID = (long) 2050803;
    private static final boolean DEBUG = false;
    
    // the pool
    private final List<GenericPoolEntry<E,T>> geneticPool;
    
    // pool configuration
    private final boolean serializeAfterNewBest;
    private final boolean writeEveryAdd;
    private final boolean beSilent;
    private final double acceptableFitness;
    private final int poolSize;
    private final int addsToSerial;
    private final int addsToStats;
    private final String interBinFile;
    
    private final T ref;
    
    private int countAddsSerial = 0;
    private int countAddsStats = 0;
    
    // helper objects
    private final transient DiversityChecker<E,T> diversity;
    private final boolean doNiching;
    private final transient Nicher<E,T> nicher;
    
    private final transient GenericStatistics stats;
    private final transient IndividualWriter<T> writer;
    
    private final transient ParentSelector<E,T> selector;
    
    /**
     * This is specific for the singleton design pattern, since the
     * constructor is private, we need to construct the object ourselves.
     * Please note that due to Java's type erasure, this is untyped!
     */
    @SuppressWarnings("rawtypes")
    private static GenericPool instance;
    
    private GenericPool(final GenericPoolConfig<E,T> conf, final T reference){
        serializeAfterNewBest = conf.serializeAfterNewBest();
        writeEveryAdd = conf.writeEveryAdd();
        beSilent = conf.shouldBeSilent();
        acceptableFitness = conf.getAcceptableFitness();
        poolSize = conf.getPoolSize();
        addsToSerial = conf.getAddsToSerial();
        addsToStats = conf.getAddsToStats();
        interBinFile = conf.getInterBinFile();
    
        // helper objects and (if applicable) configuration
        diversity = conf.getDiversityChecker();
        doNiching = conf.doNiching();
        if(DEBUG && doNiching){ System.out.println("DEBUG: We are running with niching enabled.");}
        nicher = conf.getNicher();
        stats = conf.getStats();
        if(beSilent){stats.silence();}
        writer = conf.getWriter();
        selector = conf.getSelector();
        ref = reference;
        
        // the genetic pool
        geneticPool = Collections.synchronizedList(new ArrayList<GenericPoolEntry<E,T>>(poolSize));    
    }
    
    
    /**
     * If needed lazy instantiation of the object.
     * @param config The configuration
     * @param ref A reference object of type T.
     * @return An instance of this GenericPool.
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    public synchronized static <E,T extends Optimizable<E>> GenericPool<E,T> getInstance(final GenericPoolConfig<E,T> config,
            final T ref) {
        if (instance == null){ instance = new GenericPool(config,ref);}
        return instance;
    }
    
    @SuppressWarnings({"rawtypes", "unchecked"})
    public synchronized static <E,T extends Optimizable<E>> GenericPool<E,T> getInstance() {
        if (instance == null){System.err.println("ERROR: Asking for instance of pool before we had a chance to instantiate it!"); System.exit(22);}
        return instance;
    }
    
    public int getPoolSize(){
        return poolSize;
    }
    
    public int getCurrentPoolSize(){
        return Math.min(poolSize, geneticPool.size());
    }
    
    public Iterator<GenericPoolEntry<E,T>> getPoolIterator(){
        return geneticPool.iterator();
    }
    
    public List<GenericPoolEntry<E,T>> getAllIndividuals(){
        return Collections.unmodifiableList(geneticPool);
    }
    
    @Override
    public Iterator<GenericPoolEntry<E,T>> iterator(){
        return geneticPool.iterator();
    }
    
    public GenericPoolEntry<E,T> getEntryAtPosition(final int position){
        
        if(position >= geneticPool.size()){
            System.err.println("Pool has " + geneticPool.size() + " entries, requested entry " + position + " which does not work. Returning null.");
            return null;
        }

        return geneticPool.get(position);
    }
    
    public T getIndividualAtPosition(final int position){
        
        if(position >= geneticPool.size()){throw new RuntimeException("Pool has " + geneticPool.size() + " entries, requested entry " + position + " which does not work.");}
        
        return geneticPool.get(position).getIndividual();
    }
    
    public double getFitnessOfIndividualAtPos(final int position){
        
        if(position >= geneticPool.size()){
            System.err.println("Pool has " + geneticPool.size() + " entries, requested entry " + position + " which does not work. Returning NONCONVERGEDENERGY.");
            return FixedValues.NONCONVERGEDENERGY;
        }
        
        return geneticPool.get(position).getFitness();
    }
    
    public Niche getNicheOfIndividualAtPos(final int position){
        
        if(position >= geneticPool.size()){throw new RuntimeException("Pool has " + geneticPool.size() + " entries, requested entry " + position + " which does not work.");}
        
        return geneticPool.get(position).getNiche();
    }
    
    public List<String> getFormattedPool(){
        
        final List<String> output = new ArrayList<>();
        int pos = 0;
        for(final GenericPoolEntry<E,T> entry : geneticPool){
            final T ind = entry.getIndividual();
            final long id = ind.getID();
            final double fit = ind.getFitness();
            final String s = String.format(Locale.US, "%6d   %10d  %18.10f",pos,id,fit);
            output.add(s);
            pos++;
        }
        
        return output;
    }
    
    public void removeIndividualAtPos(final int position){
        if(doNiching){
            final Niche n = geneticPool.get(position).getNiche();
            nicher.delete(n);
        }
        geneticPool.remove(position);
    }
    
    public void emptyPool(){
        if(doNiching){
            for(int pos = 0; pos < geneticPool.size(); pos++){
                final GenericPoolEntry<E,T> entry = geneticPool.get(pos);
                nicher.delete(entry.getNiche());
            }
        }
        geneticPool.clear();
    }
    
    public synchronized boolean addIndividualForced(final T individual, final double fitness){
        return addIndividualForcedUnsync(individual,null,fitness);
    }
    
    public synchronized boolean addIndividualForced(final T individual, final Niche niche, final double fitness){
        return addIndividualForcedUnsync(individual,niche,fitness);
    }
    
    /**
     * The unsychronized version of addIndividualForced. HANDLE WITH CARE! MUST ONLY BE USED IF ACCESS IS
     * SYNCHRONIZED EXPLICITLY IN THE CALLING CODE!
     * @param individual the individual. Must not be null.
     * @param niche the niche. Can be null if no niching is employed.
     * @param fitness the fitness. Must be a finite number.
     * @return true if the individual was accepted to the pool, false otherwise.
     */
    public boolean addIndividualForcedUnsync(final T individual, final Niche niche, final double fitness){
        
        assert(individual != null);
        assert(!Double.isInfinite(fitness));
        assert(!Double.isNaN(fitness));
        
        if(DEBUG){System.out.println("DEBUG: Adding individual " + individual.getID() + " to the pool (forced).");}
        
        if(writeEveryAdd) {writer.writeIndividual(individual);}
        
        final int currentSize = geneticPool.size();
        for(int pos = 0; pos < currentSize; pos++){
            final GenericPoolEntry<E,T> entry = geneticPool.get(pos);
            final double posFit = entry.getFitness();
            if(fitness < posFit){
                final GenericPoolEntry<E,T> newEntry = new GenericPoolEntry<>(individual, fitness, niche);
                geneticPool.add(pos, newEntry);
                if(doNiching && niche != null){
                    nicher.report(niche);
                }
                // ensure pool size
                ensureSize(poolSize);
                // notify...
                microManage(newEntry, pos, true);
                return true;
            }
        }
        
        if(currentSize < poolSize){
            final GenericPoolEntry<E,T> newEntry = new GenericPoolEntry<>(individual, fitness, niche);
            geneticPool.add(newEntry);
            if (doNiching && niche != null) {
                nicher.report(niche);
            }
            microManage(newEntry, currentSize, true);
            return true;
        }
        
        return false;
    }
    
    /**
     * The unsychronized version of addIndividualForced that also checks if the ID is *known* within a configurable threshold of IDs.
     * HANDLE WITH CARE! MUST ONLY BE USED IF ACCESS IS SYNCHRONIZED EXPLICITLY IN THE CALLING CODE!
     * @param individual the individual. Must not be null.
     * @param niche the niche. Can be null if no niching is employed.
     * @param fitness the fitness. Must be a finite number.
     * @param checkBackForth how many IDs to check backwards and forwards if this ID is known already. 5 seems to be a fair number.
     * @return true if the individual was accepted to the pool, false otherwise.
     */
    public boolean addIndividualForcedUnsyncCheckID(final T individual, final Niche niche, final double fitness, final int checkBackForth){
        
        assert(individual != null);
        assert(!Double.isInfinite(fitness));
        assert(!Double.isNaN(fitness));
        
        if(DEBUG){System.out.println("DEBUG: Adding individual " + individual.getID() + " to the pool (forced).");}
        
        if(writeEveryAdd) {writer.writeIndividual(individual);}
        
        final int currentSize = geneticPool.size();
        for(int pos = 0; pos < currentSize; pos++){
            final GenericPoolEntry<E,T> entry = geneticPool.get(pos);
            final double posFit = entry.getFitness();
            if(fitness < posFit){
                /*
                 * now check back and forth if this ID is known
                 */
                final long myID = individual.getID();
                
                int c = 0;
                int checkPos = pos;
                while(c < checkBackForth && checkPos >= 0){
                    final long thisID = geneticPool.get(checkPos).getIndividual().getID();
                    if(thisID == myID){
                        // known
                        return false;
                    }
                    c++;
                    checkPos--;
                }
                
                c = 1;
                checkPos = pos+1;
                while(c < checkBackForth && checkPos < currentSize){
                    final long thisID = geneticPool.get(checkPos).getIndividual().getID();
                    if(thisID == myID){
                        // known
                        return false;
                    }
                    c++;
                    checkPos++;
                }
                
                final GenericPoolEntry<E,T> newEntry = new GenericPoolEntry<>(individual, fitness, niche);
                geneticPool.add(pos, newEntry);
                if(doNiching && niche != null){
                    nicher.report(niche);
                }
                // ensure pool size
                ensureSize(poolSize);
                // notify...
                microManage(newEntry, pos, true);
                return true;
            }
        }
        
        if(currentSize < poolSize){
            final GenericPoolEntry<E,T> newEntry = new GenericPoolEntry<>(individual, fitness, niche);
            geneticPool.add(newEntry);
            if (doNiching && niche != null) {
                nicher.report(niche);
            }
            microManage(newEntry, currentSize, true);
            return true;
        }
        
        return false;
    }
    
    public synchronized void replacePoolContent(final List<T> newIndividuals){
        
        unSyncReplacePoolContent(newIndividuals, null);
    }
    
    public synchronized void replacePoolContent(final List<T> newIndividuals, final List<Niche> niches){
        unSyncReplacePoolContent(newIndividuals,niches);
    }
    
    public void unSyncReplacePoolContent(final List<T> newIndividuals, final List<Niche> niches){
        
        if(DEBUG){System.out.println("DEBUG: Replacing all individuals in the pool (forced).");}
        if(newIndividuals.size() < poolSize && newIndividuals.size() < geneticPool.size()){
            System.out.println("INFO: Replacing existing pool with less than the specified pool size and less than the previous number of individuals. " + poolSize + " " + newIndividuals.size() + " " + geneticPool.size());
        }
        
        emptyPool(); // remove what we had
        // XXX probably, we should also inform the nicher of the purge here
        
        // now simply replace
        double lastFitness = Double.NEGATIVE_INFINITY;
        for(int i = 0; i < newIndividuals.size(); i++){
            final T individual = newIndividuals.get(i);
            final double fitness = individual.getFitness();
            
            if(fitness < lastFitness){throw new RuntimeException("Input not ordered: Last fitness " + lastFitness + ", new fitness " + fitness);}
            lastFitness = fitness;
            
            final GenericPoolEntry<E, T> newEntry = (niches == null) ? new GenericPoolEntry<>(individual, fitness, null)
                    : new GenericPoolEntry<>(individual, fitness, niches.get(i));
            geneticPool.add(newEntry);
            if (doNiching && niches != null) {
                nicher.report(niches.get(i));
            }
        }
        
        // only now ensure size. do not micromanage. We do not need this here.
        ensureSize(poolSize);
    }
    
    public synchronized boolean addIndividual(final T individual, final double fitness){
        return addIndividualUnsync(individual, null, fitness);
    }
    
    public synchronized boolean addIndividual(final T individual, final Niche niche, final double fitness){
        return addIndividualUnsync(individual, niche, fitness);
    }

    /**
     * The unsychronized version of addIndividual. HANDLE WITH CARE! MUST ONLY BE USED IF ACCESS IS
     * SYNCHRONIZED EXPLICITLY IN THE CALLING CODE!
     * @param individual the individual. Must not be null.
     * @param niche the niche. Can be null if no niching is employed.
     * @param fitness the fitness. Must be a finite number.
     * @return true if the individual was accepted to the pool, false otherwise.
     */
    public boolean addIndividualUnsync(final T individual, final Niche niche, final double fitness){
        
        assert(individual != null);
        assert(!Double.isInfinite(fitness));
        assert(!Double.isNaN(fitness));
        
        if(writeEveryAdd) {writer.writeIndividual(individual);}
        
        // only check this if the pool is fully filled
        final int currentSize = geneticPool.size();
        if(currentSize >= poolSize && fitness >= geneticPool.get(geneticPool.size()-1).getFitness()) {
            stats.registerIndividualNotAdded(individual.getID());
            return false;
        } // fitness out of range
        
        
        if(doNiching && niche != null) {
            
            final GenericPoolEntry<E,T> newEntry = new GenericPoolEntry<>(individual, fitness, niche);
            
            boolean bottomOfNiche = true;
            // get all the individuals of the same niche and check the diversity with them
            boolean firstPlaceSeen = false;
            int removePos = -1;
            for(int pos = 0; pos < geneticPool.size(); pos++){
                
                final Niche other = getNicheOfIndividualAtPos(pos);
                if(other.comp(niche)){
                    if(!firstPlaceSeen && geneticPool.get(pos).getFitness() > fitness){
                        // new best individual in this niche: accept anyways
                        final GenericPoolEntry<E,T> thisEntry = geneticPool.get(pos);
                        final boolean areDiverse = diversity.areDiverse(newEntry, thisEntry);
                        if(!areDiverse){
                            // mark previous individual for removal as it was worse and is not diverse enough!
                            removePos = pos;
                        }
                        break;
                    }
                    
                    // let's check the diversity
                    final GenericPoolEntry<E,T> thisEntry = geneticPool.get(pos);
                    if(bottomOfNiche){
                        bottomOfNiche = false;
                        if(fitness <= thisEntry.getFitness()) break; // accept new best entry w/o diversity-check!
                    }
                    final boolean areDiverse = diversity.areDiverse(newEntry, thisEntry);
                    if(!areDiverse){
                        // tough luck
                        stats.registerIndividualNotAdded(individual.getID());
                        return false;
                    }
                    
                    firstPlaceSeen = true;
                }
            }
            
            // if we end up here, the individual indeed is diverse w.r.t. all the other individuals in the same
            // niche or is the new best individual in this niche, now just add :-)
            for(int pos = 0; pos < geneticPool.size(); pos++){
                final GenericPoolEntry<E,T> entry = geneticPool.get(pos);
                final double posFit = entry.getFitness();
                if(fitness < posFit){
                
                    final GenericPoolEntry<E,T> addEntry = newEntry.clone();
                    geneticPool.add(pos, addEntry);
                    
                    nicher.report(niche);
                    boolean removed = false;
                    if(removePos < 0){                        
                        removed = nicher.cleanUp(this);
                    } else {
                        if(pos <= removePos){removePos++;}
                        
                        nicher.delete(geneticPool.get(removePos).getNiche());
                        geneticPool.remove(removePos);
                    }
                    if(DEBUG) {
                        System.out.println("DEBUG: This fitness " + fitness + " compared to " + posFit);
                        System.out.println("DEBUG: Did Nicher remove something? " + removed);
                        System.out.println("DEBUG: Worst individual: " + geneticPool.get(geneticPool.size()-1).getFitness());
                        if(removed){
                            for(int i = 0; i < geneticPool.size(); i++){
                                System.out.println(" " + i + " fitness " + geneticPool.get(i).getFitness() + " niche " + geneticPool.get(i).getNiche().getID());
                            }
                        }
                    }

                    // ensure pool size
                    ensureSize(poolSize);
                    // notify...
                    microManage(addEntry, pos, false);          
                    return true;
                }
            }
            
            // if the pool is not entirely filled, add in the end
            if(currentSize < poolSize){
                geneticPool.add(newEntry);
                microManage(newEntry, currentSize, true);
                nicher.report(niche);
                return true;
            }
        } else {
            // not niching, simpler case
            for(int pos = 0; pos < geneticPool.size(); pos++){
                final GenericPoolEntry<E,T> entry = geneticPool.get(pos);
                final double posFit = entry.getFitness();
                if(fitness < posFit){
                
                    final GenericPoolEntry<E,T> newEntry = new GenericPoolEntry<>(individual, fitness, null);
                    final boolean divFront = (pos == 0) ? true : diversity.areDiverse(newEntry, geneticPool.get(pos-1));
                    if(!divFront){
                        // not enough diverstiy to the front -> reject
                        stats.registerIndividualNotAdded(individual.getID());
                        return false;
                    }
                    
                    final boolean divBack = diversity.areDiverse(newEntry, entry);
                    if(!divBack){
                        // not enough diversity compared to the previous pos individual -> "replace" the other individual, thereby keeping the pool size constant
                        final GenericPoolEntry<E,T> addEntry = newEntry.clone();
                        geneticPool.set(pos, addEntry);
                        // ensure pool size
                        ensureSize(poolSize);
                        // notify...
                        microManage(addEntry, pos, false);          
                        return true;
                    }
                    
                    // both diversity to the front AND to the back. just add in between.
                    final GenericPoolEntry<E,T> addEntry = newEntry.clone();
                    geneticPool.add(pos, addEntry);
                    
                    // ensure pool size
                    ensureSize(poolSize);
                    // notify...
                    microManage(addEntry, pos, false);          
                    return true;
                }
            }
            
            // if the pool is not entirely filled, add in the end
            if(currentSize < poolSize){
                final GenericPoolEntry<E,T> newEntry = new GenericPoolEntry<>(individual, fitness, niche);
                geneticPool.add(newEntry);
                microManage(newEntry, currentSize, true);
                return true;
            }
        }
        
        stats.registerIndividualNotAdded(individual.getID());
        return false;
    }
    
    public List<T> getParents(){
        return selector.getParents(this);
    }
    
    public double[] getAllFitnesses(){
        
        final double[] fs = new double[geneticPool.size()];
        for(int i = 0; i < geneticPool.size(); i++){
            fs[i] = geneticPool.get(i).getFitness();
        }
        
        return fs;
    }
    
    public boolean acceptableFitnessReached(){
        // I could construct some really weird situation in which this could fail because of a race on geneticPool
        // short of locking everything inside here, there is not much we can do about this. I am willing to take a
        // risk
        return (geneticPool.isEmpty()) ? false : (geneticPool.get(0).getFitness() <= acceptableFitness);
    }
    
    public T getExample(){
        return ref;
    }
    
    public List<String> myConfiguration(){
        
        final List<String> conf = new ArrayList<>();
        conf.add(" size of genetic pool " + poolSize);
        conf.add(" serialize after new best " + serializeAfterNewBest);
        conf.add(" write every add " + writeEveryAdd);
        conf.add(" additions to serialization " + addsToSerial);
        conf.add(" additions to statistics " + addsToStats);
        conf.add(" intermediate pool file " + interBinFile);
        conf.add(" target fitness for pool " + acceptableFitness);
        conf.add(" are we niching? " + doNiching);
        
        return conf;
    }
    
    private void microManage(final GenericPoolEntry<E,T> newEntry, final int posAdded, final boolean forced){
        
        assert(newEntry != null);
        assert(posAdded >= 0);
        assert(newEntry.getIndividual() != null);
        
        countAddsSerial++;
        if(!forced) {countAddsStats++;}
        stats.individualAddedToPool(newEntry.getIndividual().getID(), posAdded, newEntry.getFitness());
        if((posAdded == 0 && serializeAfterNewBest) || countAddsSerial >= addsToSerial){
            serializeMe();
            countAddsSerial = 0;
        }
        if(countAddsStats >= addsToStats && doNiching){
            if(DEBUG) {System.out.println("DEBUG: Presenting some niching statistics...");}
            nicher.calculateNichePopulation(true, this);
            if(DEBUG) {System.out.println("DEBUG: Done presenting niching statistics.");}
            countAddsStats = 0;
        }
    }
    
    void serializeMe(){

        if(beSilent){return;}
        
        try {
            OutputPrimitives.writeObjToBinFile(interBinFile, this);
        } catch (Exception e) {
            // some error occured during the serialization attempt, try again
            e.printStackTrace(System.err);
            try {
                OutputPrimitives.writeObjToBinFile(interBinFile, this);
            } catch (Exception e2) {
                // ok, it failed a second time. too bad...
                e2.printStackTrace(System.err);
                System.err.println("Double failure in trying to write the genetic pool out. Continuing without guarantee now!");
            }

        }
    }
    
    private void ensureSize(final int allowedSize){
        while(geneticPool.size() > allowedSize){
            if(doNiching){
                // report decrement in niche population
                nicher.delete(geneticPool.get(geneticPool.size()-1).getNiche());
            }
            geneticPool.remove(geneticPool.size()-1);
        }
    }
}
