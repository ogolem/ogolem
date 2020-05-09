/**
Copyright (c) 2013-2014, J. M. Dieterich
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
package org.ogolem.generic;

import org.ogolem.core.GlobalConfig;
import org.ogolem.generic.stats.GenericDetailStatistics;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.Lottery;
import java.util.ArrayList;
import java.util.List;

/**
 * A generic, abstract implementation of a classical GA global optimization.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public abstract class GenericAbstractDarwin<E,T extends Optimizable<E>> implements GenericDarwin<E,T> {
    
    private static final long serialVersionUID = (long) 20200429;
    protected final GenericCrossover<E,T> xover;
    protected final GenericMutation<E,T> mutation;
    protected final GenericSanityCheck<E,T> sanitizer;
    protected final GenericFitnessFunction<E,T> fitness;
    protected final IndividualWriter<T> writer;
    protected final Lottery r;
    protected final boolean DEBUG;
    protected final boolean printBeforeFitness;
    protected final int noOfTries;
    protected double crossPoss;
    protected double mutPoss;
    
    public GenericAbstractDarwin(final GenericCrossover<E,T> cross, final GenericMutation<E,T> mut,
            final GenericSanityCheck<E,T> sanity, final GenericFitnessFunction<E,T> fitness,
            final IndividualWriter<T> writer, final double crossPoss, final double mutPoss,
            final boolean printBeforeFitness, final int noOfTries){
        this.DEBUG = (GlobalConfig.DEBUGLEVEL > 0);
        this.r = Lottery.getInstance();
        this.xover = cross;
        this.mutation = mut;
        this.sanitizer = sanity;
        this.fitness = fitness;
        this.writer = writer;
        this.crossPoss = crossPoss;
        this.mutPoss = mutPoss;
        this.printBeforeFitness = printBeforeFitness;
        this.noOfTries = noOfTries;
    }
    
    public GenericAbstractDarwin(final GenericAbstractDarwin<E,T> orig){
        this.xover = orig.xover.clone();
        this.mutation = orig.mutation.clone();
        this.sanitizer = orig.sanitizer.clone();
        this.fitness = orig.fitness.copy();
        this.writer = orig.writer.clone();
        this.crossPoss = orig.crossPoss;
        this.mutPoss = orig.mutPoss;
        this.printBeforeFitness = orig.printBeforeFitness;
        this.noOfTries = orig.noOfTries;
        
        this.DEBUG = (GlobalConfig.DEBUGLEVEL > 0);
        this.r = Lottery.getInstance();
    }
    
    @Override
    public abstract GenericAbstractDarwin<E,T> copy();
    
    @Override
    public abstract String getMyID();
    
    @Override
    @SuppressWarnings("unchecked")
    public T globalOptimization(final long futureID, final T mother, final T father){
        
        @SuppressWarnings("unchecked")
        final List<T> gs = new ArrayList<>(2);
        int off = 0;
        for (int tryc = 0; tryc < noOfTries; tryc++) {
            
            GenericDetailStatistics.incrementTrials();
            
            T child1;
            T child2;
            
            // should we even cross?
            final double dc = r.nextDouble();
            if(dc <= crossPoss){
                final Tuple<T,T> childGeoms = cross(mother, father, futureID);
                if(childGeoms.getObject1() == null) continue; // crossing signals problem
                child1 = childGeoms.getObject1();
                child2 = childGeoms.getObject2();
            } else {
                child1 = (T) mother.copy();
                child2 = (T) father.copy();
            }
            
            // set numbers
            if(child1 != null) child1.setID(futureID);
            if(child2 != null) child2.setID(futureID);
            
            // should we mutate?
            final double dm = r.nextDouble();
            if(dc > crossPoss || dm <= mutPoss){
                // either because mutation should be done or because we did no crossover
                if(child1 != null) child1 = mutate(child1);
                if(child2 != null) child2 = mutate(child2);
            }
                        
            // sanity check (and if applicable) fitness function evaluation
            if(child1 != null){
                if(DEBUG){System.out.println("DEBUG: Child 1 for " + futureID + " was not null.");}
                
                final boolean sane = sanitizer.isSane(child1);
                if(sane){
                    if(DEBUG){System.out.println("DEBUG: CHILD1 for " + futureID  + " was sane.");}
                    if(printBeforeFitness){
                        System.out.println("INFO: CHILD1 " + child1.getID() + " BEFORE FITNESS COMING.");
                        writer.writeIndividual(child1);
                    }
                    child1 = fitness.fitness(child1, false);
                } else{
                    GenericDetailStatistics.incrementSanityDiscards();
                    if(DEBUG){System.out.println("DEBUG: Child 1 insanity found for " + futureID);}
                    child1 = null;
                }
            }
            
            if(child2 != null){
                if(DEBUG){System.out.println("DEBUG: Child 2 for " + futureID + " was not null.");}
                
                final boolean sane = sanitizer.isSane(child2);
                if(sane){
                    if(DEBUG){System.out.println("DEBUG: CHILD2 for " + futureID  + " was sane.");}
                    if(printBeforeFitness){
                        System.out.println("INFO: CHILD2 " + child2.getID() + " BEFORE FITNESS COMING.");
                        writer.writeIndividual(child2);
                    }
                    child2 = fitness.fitness(child2, false);
                } else{
                    GenericDetailStatistics.incrementSanityDiscards();
                    if(DEBUG){System.out.println("DEBUG: Child 2 insanity found for" + futureID);}
                    child2 = null;
                }
            }
            
            // set numbers
            if(child1 != null){
                child1.setID(futureID);
                child1.setFatherID(father.getID());
                child1.setMotherID(mother.getID());
            }
            if(child2 != null){
                child2.setID(futureID);
                child2.setFatherID(father.getID());
                child2.setMotherID(mother.getID());
            }
            
            // lets see...
            boolean whereStart;
            final short prio = xover.hasPriority();
            if(prio >= 0){
                whereStart = (prio == 0);
            } else{
                whereStart = whereToStart();
            }
            if(whereStart){
                // first check one then two
                if(child1 != null && off < 2){
                    gs.add(child1);
                    off++;
                }
                if(child2 != null && off < 2){
                    gs.add(child2);
                    off++;
                }
            } else{
                // first check two then one
                if(child2 != null && off < 2){
                    gs.add(child2);
                    off++;
                }
                if(child1 != null && off < 2){
                    gs.add(child1);
                    off++;
                }
            }
            
            // all filled in
            if(off >= 2) break;
            
            // any last wishes?
            runAfterEachTry();
        }
        
        if(gs.isEmpty()) return null;
        else if(gs.size() == 1) return gs.get(0);
        else {
        	final T candidate = (gs.get(0).getFitness() <= gs.get(1).getFitness()) ? gs.get(0) : gs.get(1);
            return candidate;
        }
    }
    
    @Override
    public T mutate(final T start){
        return mutation.mutate(start);
    }
	
    @Override
    public Tuple<T,T> cross(final T mother, final T father, final long futureID){
        return xover.crossover(mother, father, futureID);
    }
    
    protected abstract void postXOver(final T individual1, final T individual2, final long futureID);
    
    protected abstract void postMutation(T individual);
    
    protected abstract void runAfterEachTry();
    
    protected boolean whereToStart(){
        return r.nextBoolean();
    }
}
