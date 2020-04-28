/**
Copyright (c) 2012-2014, J. M. Dieterich
              2016, J. M. Dieterich and B. Hartke
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

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import org.ogolem.generic.Optimizable;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * A collection of some simple parent selection algorithms.
 * @author Johannes Dieterich
 * @author Bernd Hartke
 * @version 2016-01-14
 */
public class GenericParentSelectors {
    
    public static <E,T extends Optimizable<E>> ParentSelector<E,T> buildSelector(final String input) throws Exception {
        
        if(input.startsWith("random:")){
            final String s = input.substring(7).trim();
            final String[] tokens = s.trim().split("\\,");
            
            double stepAt1 = 1.0;
            double stepAt2 = 1.0;
            for(final String token : tokens){
                if(token.equalsIgnoreCase("stepAt=")){
                    stepAt1 = Double.parseDouble(token.substring(7).trim());
                    stepAt2 = stepAt1;
                }
                if(token.equalsIgnoreCase("stepAt1=")){
                    stepAt1 = Double.parseDouble(token.substring(8).trim());
                }
                if(token.equalsIgnoreCase("stepAt2=")){
                    stepAt2 = Double.parseDouble(token.substring(8).trim());
                }
            }
            
            return new RandomParentSelector<>(stepAt1,stepAt2);
            
        } else if(input.startsWith("fitnessrankbased:")){
            final String s = input.substring(17).trim();
            final String[] tokens = s.trim().split("\\,");
            
            boolean bothFitness = false;
            boolean perNiche = false;
            double gaussWidth1 = 0.5;
            double gaussWidth2 = 0.5;
            double sameNiche = 0.0;
            for(final String token : tokens){
                if(token.equalsIgnoreCase("bothfitness")){
                    bothFitness = true;
                } else if(token.startsWith("perniche")){
                    perNiche = true;
                } else if(token.startsWith("sameniche=")){
                    sameNiche = Double.parseDouble(token.substring(10).trim());
                } else if(token.startsWith("gausswidth=")){
                    gaussWidth1 = Double.parseDouble(token.substring(11).trim());
                    gaussWidth2 = gaussWidth1;
                } else if(token.startsWith("gausswidth1=")){
                    gaussWidth1 = Double.parseDouble(token.substring(12).trim());
                } else if(token.startsWith("gausswidth2=")){
                    gaussWidth2 = Double.parseDouble(token.substring(12).trim());
                } else{
                    throw new RuntimeException("Illegal modifier " + token + " for fitness rank based selector.");
                }
            }
            
            return new FitnessWeightedParentSelector<>(gaussWidth1,gaussWidth2,bothFitness,perNiche,sameNiche);
            
        } else if(input.startsWith("fitnessvaluebased:")){
            final String s = input.substring(18).trim();
            final String[] tokens = s.trim().split("\\,");
            
            boolean useGaussShape = false;
            boolean bothFitness = false;
            double endProb = 0.1;
            for(final String token : tokens){
                if(token.equalsIgnoreCase("bothfitness")){
                    bothFitness = true;
                } else if(token.equalsIgnoreCase("usegausshape")){
                    useGaussShape = true;
                } else if(token.startsWith("tailprob=")){
                    endProb = Double.parseDouble(token.substring(9).trim());
                } else{
                    throw new RuntimeException("Illegal modifier " + token + " for fitness value based selector.");
                }
            }
            
            return new RealFitnessWeightedParentSelector<>(useGaussShape,bothFitness,endProb);
            
        } else {
            throw new RuntimeException("Illegal parent selector " + input + ".");
        }
    }
    
    public static class RandomParentSelector<E,T extends Optimizable<E>> implements ParentSelector<E,T> {
        
        private static final long serialVersionUID = (long) 20120215;
        private final double stepAt1;
        private final double stepAt2;
        
        public RandomParentSelector(final double stepAt1, final double stepAt2){
            this.stepAt1 = stepAt1;
            this.stepAt2 = stepAt2;
        }
        
        @SuppressWarnings("unchecked")
        @Override
        public List<T> getParents(final GenericPool<E,T> pool){
            
            assert(0.0 <= stepAt1 && stepAt1 <= 1.0);
            assert(0.0 <= stepAt2 && stepAt2 <= 1.0);
            final int poolSize = pool.getCurrentPoolSize();
            final Lottery r = Lottery.getInstance();
            // with stepAt=1.0, this really is purely random over the whole pool;
            // with stepAt<1.0, this corresponds to Karl Heinz Hoffmann's favorite selector: a step function
            final int index1 = r.nextInt((int) Math.round(stepAt1 * poolSize));
            final int index2 = r.nextInt((int) Math.round(stepAt2 * poolSize));
            
            final List<T> parents = new LinkedList<>();
            parents.add((T) pool.getIndividualAtPosition(index1).clone());
            parents.add((T) pool.getIndividualAtPosition(index2).clone());
            
            return parents;
        }
    }
    
    public static class FitnessWeightedParentSelector<E,T extends Optimizable<E>> implements ParentSelector<E,T> {
        
        private static final long serialVersionUID = (long) 20120215;
        private final double gaussShaper1;
        private final double gaussShaper2;
        private final boolean bothFitness;
        private boolean perNiche;
        private double sameNiche;
        
        public FitnessWeightedParentSelector(final double gaussShaper1, final double gaussShaper2,
                final boolean bothFitness, boolean perNiche, final double sameNiche){
            this.gaussShaper1 = gaussShaper1;
            this.gaussShaper2 = gaussShaper2;
            this.bothFitness = bothFitness;
            this.perNiche = perNiche;
            this.sameNiche = sameNiche;
        }
        
        @SuppressWarnings("unchecked")
        @Override
        public List<T> getParents(final GenericPool<E,T> pool){
            
            final int poolSize = pool.getCurrentPoolSize();
            final Lottery r = Lottery.getInstance();
            final int binSize1, binSize2;
            int niche1 = 0;
            int niche2 = 0;
            
            final boolean haveNiches = (pool.getNicheOfIndividualAtPos(0) != null);
            if(perNiche && !haveNiches){ // logically inconsistent, fix silently
                perNiche = false;
                sameNiche = 0.0;
            }
            List<String> listNicheNames = new ArrayList<>();
            List<List<Integer>> listNicheMemberPositions = new ArrayList<>();
            if(perNiche){
                // necessary (costly!) preparation: get needed niche info
                for(int i = 0; i < poolSize; i++){
                    final Niche niche = pool.getNicheOfIndividualAtPos(i);
                    if(niche == null){
                        System.err.println("ERROR: No niche info for individual at pool position " + i + ". Returning position 0 twice as parent.");
                        final List<T> parents = new LinkedList<>();
                        parents.add((T) pool.getIndividualAtPosition(0).clone());
                        parents.add((T) pool.getIndividualAtPosition(0).clone());
                        
                        return parents;
                    }
                    final String myNicheName = niche.getID();
                    if(listNicheNames.contains(myNicheName)){
                        listNicheMemberPositions.get(listNicheNames.indexOf(myNicheName)).add(i);
                    }else{
                        listNicheNames.add(myNicheName);
                        List<Integer> myMemberPositions = new ArrayList<>();
                        myMemberPositions.add(i);
                        listNicheMemberPositions.add(myMemberPositions);
                    }
                }
                // now choose two of the niches randomly (for mother and father),
                // the same or different, depending on sameNiche
                final int numberOfNiches = listNicheNames.size();
                if(r.nextDouble() <= sameNiche){
                    niche1 = r.nextInt(numberOfNiches);
                    niche2 = niche1;
                    binSize1 = listNicheMemberPositions.get(niche1).size();
                    binSize2 = binSize1;
                }else{
                    niche1 = r.nextInt(numberOfNiches);
                    niche2 = r.nextInt(numberOfNiches);
                    binSize1 = listNicheMemberPositions.get(niche1).size();
                    binSize2 = listNicheMemberPositions.get(niche2).size();
                }
            } else{
                binSize1 = poolSize;
                binSize2 = binSize1;
            }
            
            final int index1;
            if(bothFitness){
                final double gaussR = RandomUtils.halfgaussDouble(0.0, 1.0, gaussShaper1);
                // the bin where our mother is
                index1 = (int) Math.floor(gaussR * binSize1);
            } else{
                // the mother is completely random.... guys.... ;-)
                index1 = r.nextInt(binSize1);
            }
            
            // the father is always fitness-based
            final double gaussR = RandomUtils.halfgaussDouble(0.0, 1.0, gaussShaper2);
            
            // the bin where our father is
            final int index2 = (int) Math.floor(gaussR * binSize2);
            
            assert(index1 < binSize1);
            assert(index2 < binSize2);
            
            final List<T> parents = new LinkedList<>();
            if(perNiche){
                parents.add((T) pool.getIndividualAtPosition(listNicheMemberPositions.get(niche1).get(index1)).clone());
                parents.add((T) pool.getIndividualAtPosition(listNicheMemberPositions.get(niche2).get(index2)).clone());
            }else{
                parents.add((T) pool.getIndividualAtPosition(index1).clone());
                parents.add((T) pool.getIndividualAtPosition(index2).clone());
            }
            
            return parents;
        }
    }
    
    public static class RealFitnessWeightedParentSelector<E,T extends Optimizable<E>> implements ParentSelector<E,T> {
        
        private static final long serialVersionUID = (long) 20140713;
        private final boolean bothFitness;
        private final double endProb;
        private final boolean useGaussShape;
        
        public RealFitnessWeightedParentSelector(final boolean useGaussShape, final boolean bothFitness,
                final double endProb){
            assert(endProb >= 0.0 && endProb <= 1.0);
            this.useGaussShape = useGaussShape;
            this.bothFitness = bothFitness;
            this.endProb = endProb;
        }
        
        @SuppressWarnings("unchecked")
        @Override
        public List<T> getParents(final GenericPool<E,T> pool){
            
            // WARNING: as this does not lock the pool, it is therefore intrinsically inaccurate (due to structures moving)
            
            final int poolSize = pool.getCurrentPoolSize();
            
            // figure out the top and bottom of the pool
            final double[] allFits = pool.getAllFitnesses();
            assert(allFits.length == poolSize);
            final double bestFitness = allFits[0];
            final double worstFitness = allFits[poolSize-1];
            final double diffFitness = worstFitness - bestFitness;
            
            final double m = -(1.0-endProb)/diffFitness;
            
            final double[] probs = new double[poolSize];
            double sumProbs = 0.0;
            // we use a basic linear dependence upon fitness
            for(int i = 0; i < poolSize; i++){
                final double thisFitDiff = allFits[i]-bestFitness; // this is a positive number
                assert(thisFitDiff >= -1E-5);
                // we define the first to get 1.0, the last to get endProb
                final double prob = m*thisFitDiff+1.0;
                sumProbs += prob;
                probs[i] = prob;
                assert(prob <= 1.0001 && prob >= (endProb-1E5)) : "Prob " + prob + " end " + endProb;
            }
            
            // normalize and turn to endsets
            for(int i = 0; i < poolSize; i++){
                final double d = probs[i] / sumProbs;
                probs[i] = (i == 0) ? d : probs[i-1]+d;
            }
            probs[poolSize-1] = 1.0001; // to get rid of numerical inaccuracies
            
            
            final Lottery r = Lottery.getInstance();
            
            int index1 = -1;
            if(bothFitness){
                final double rnd = (useGaussShape) ? r.nextDouble() : RandomUtils.gaussDouble(0.1, 1.0);
                for(int i = 0; i < poolSize; i++){if(rnd < probs[i]){index1 = i; break;}}
            } else{
                // the mother is completely random.... guys.... ;-)
                index1 = r.nextInt(poolSize);
            }
            
            // the father is always fitness-based
            final double rnd = (useGaussShape) ? r.nextDouble() : RandomUtils.gaussDouble(0.1, 1.0);
            int index2 = -1;
            for(int i = 0; i < poolSize; i++){if(rnd < probs[i]){index2 = i; break;}}
                        
            assert(index1 < poolSize && index1 >= 0);
            assert(index2 < poolSize && index2 >= 0);
            
            final List<T> parents = new LinkedList<>();
            parents.add((T) pool.getIndividualAtPosition(index1).clone());
            parents.add((T) pool.getIndividualAtPosition(index2).clone());
            
            return parents;
        }
    }
}
