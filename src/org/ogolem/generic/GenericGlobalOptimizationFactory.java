/**
Copyright (c) 2014, J. M. Dieterich
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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.helpers.Tuple;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A factory for the generic global optimization algorithms.
 * @author Johannes Dieterich
 * @version 2020-02-12
 */
public abstract class GenericGlobalOptimizationFactory<E, T extends Optimizable<E>> implements Serializable {
    
    private static final long serialVersionUID = (long) 20200212;
    private static final Logger LOG = LoggerFactory.getLogger(GenericGlobalOptimizationFactory.class);
    
    protected final long totalSteps;
    
    protected GenericGlobalOptimizationFactory(final long totalGlobSteps){
        this.totalSteps = totalGlobSteps;
    }
    
    public abstract GenericGlobalOptimization<E,T> translateToGlobOpt(final String globOptString) throws Exception;

    public final Tuple<List<Double>,List<GenericCrossover<E,T>>> getXoversAndProbs(final String xOverString) throws Exception {
        
        final String[] sa = tokenizeFirstLevel(xOverString);
        final List<GenericCrossover<E, T>> xovers = new ArrayList<>();
        final List<Double> probs = new ArrayList<>();
        for (final String s : sa) {
            final String[] sa2 = s.split("\\%");
            final double prob = Double.parseDouble(sa2[0]) / 100;
            final GenericCrossover<E, T> xover = getXOver(sa2[1]);
            probs.add(prob);
            xovers.add(xover);
        }
        
        return new Tuple<>(probs,xovers);
    }
    
    public final Tuple<List<Double>,List<GenericMutation<E,T>>> getMutationsAndProbs(final String mutString) throws Exception {
        
        final String[] sa = tokenizeFirstLevel(mutString);
        final List<GenericMutation<E, T>> mutations = new ArrayList<>();
        final List<Double> probs = new ArrayList<>();
        for (final String s : sa) {
            final String[] sa2 = s.split("\\%");
            final double prob = Double.parseDouble(sa2[0]) / 100;
            final GenericMutation<E,T> mutation = getMutation(sa2[1]);
            probs.add(prob);
            mutations.add(mutation);
        }
        
        return new Tuple<>(probs,mutations);
    }
    
    public final GenericCrossover<E,T> getXOver(final String xOverString) throws Exception {
        
        /*
         * first try the specialized one
         */
        final GenericCrossover<E,T> special = specializedXOver(xOverString);
        if(special != null){
            return special;
        }
        
        LOG.debug("No special crossover found, trying generic ones.");
        
        if(xOverString.equalsIgnoreCase("nocrossover") || xOverString.equalsIgnoreCase("nocrossover:")
                || xOverString.equalsIgnoreCase("noxover") || xOverString.equalsIgnoreCase("noxover:")){
            return new NoCrossover<>();
        } else if(xOverString.startsWith("multiple:")){
            final Tuple<List<Double>,List<GenericCrossover<E,T>>> stuff = getXoversAndProbs(xOverString.substring(9));
            
            return new MultipleXOverWrapper<>(stuff.getObject2(),stuff.getObject1());
        } else if(xOverString.startsWith("chained:")){
            final Tuple<List<Double>,List<GenericCrossover<E,T>>> stuff = getXoversAndProbs(xOverString.substring(8));
            
            return new GenericChainedXOver<>(stuff.getObject2(),stuff.getObject1());
        } else if(xOverString.startsWith("mutationasxover:")){
            final String s = xOverString.substring(16).trim();
            final GenericMutation<E,T> mut = getMutation(s);
            
            return new GenericMutationAsXOver<>(mut);
        } else if(xOverString.startsWith("changingxover:")){
            final String s = xOverString.substring(14).trim();
            
            final String[] dat = s.split("\\]");
            final List<GenericCrossover<E,T>> xovers = new ArrayList<>();
            final List<Double> percs = new ArrayList<>();
            for(final String datum : dat){
                final String[] spl = datum.trim().split("\\%\\[");
                final double perc = Double.parseDouble(spl[0].trim())/100;
                percs.add(perc);
                final GenericCrossover<E,T> xover = getXOver(spl[1].trim());
                xovers.add(xover);
            }
            
            return new GenericStepChangingXOver<>(totalSteps, xovers, percs);
        } else if(xOverString.startsWith("portugal:")){
            
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(9));
            int noCuts = 1;
            for (final String token : tokens) {
                if (token.startsWith("nocuts=")) {
                    noCuts = integerToken("nocuts=",token);
                    if(noCuts < 0){throw new RuntimeException("Number of cuts must be positive!");}
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (portugal).");
                }
            }

            return new GenericPortugalCrossover<>(noCuts);
        } else if(xOverString.startsWith("germany:")){
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(8));
            double gaussWidth = 0.3;
            for (final String token : tokens) {
                if (token.startsWith("gausswidth=")) {
                    gaussWidth = doubleToken("gausswidth=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (germany).");
                }
            }

            return new GenericGermanyCrossover<>(gaussWidth);
        }
        
        throw new IllegalArgumentException("Unknown crossover: " + xOverString + " both generic and special");
    }
    
    public final GenericMutation<E,T> getMutation(final String mutString) throws Exception {
        
        /*
         * first try the specialized one
         */
        final GenericMutation<E,T> special = specializedMutation(mutString);
        if(special != null){
            return special;
        }
        
        LOG.debug("No special mutation found for " + mutString +  ", trying generic ones.");
        
        if(mutString.equalsIgnoreCase("nomutation") || mutString.equalsIgnoreCase("nomutation:")){
            return new NoMutation<>();
        } else if(mutString.startsWith("multiple:")){
            final Tuple<List<Double>,List<GenericMutation<E,T>>> stuff = getMutationsAndProbs(mutString.substring(9));
            
            return new MultipleMutationWrapper<>(stuff.getObject2(),stuff.getObject1());
        } else if(mutString.startsWith("chained:")){
            final Tuple<List<Double>,List<GenericMutation<E,T>>> stuff = getMutationsAndProbs(mutString.substring(8));
            
            return new GenericChainedMutation<>(stuff.getObject2(),stuff.getObject1());
        }  else if(mutString.startsWith("changingmut:")){
            final String s = mutString.substring(12).trim();
            
            final String[] dat = s.split("\\]");
            final List<GenericMutation<E,T>> muts = new ArrayList<>();
            final List<Double> percs = new ArrayList<>();
            for(final String datum : dat){
                final String[] spl = datum.trim().split("\\%\\[");
                final double perc = Double.parseDouble(spl[0].trim())/100;
                percs.add(perc);
                final GenericMutation<E,T> mut = getMutation(spl[1].trim());
                muts.add(mut);
            }
            
            return new GenericStepChangingMut<>(totalSteps, muts, percs);
        }
        
        
        throw new IllegalArgumentException("Unknown mutation: " + mutString + " both generic and special");
    }
    
    protected GenericCrossover<E,T> specializedXOver(final String xOverString) throws Exception {
        return null;
    }
    
    protected GenericMutation<E,T> specializedMutation(final String mutString) throws Exception {
        return null;
    }
    
    /*
     * HELPER FUNCTIONS FOR OUR GRAMMAR
     */
    protected final String[] tokenizeFirstLevel(final String s){
        
        if(s.trim().isEmpty()){
            return new String[0];
        }
        
        final String[] tokens = s.trim().split("\\|");
        for(int i = 0; i < tokens.length; i++){
            tokens[i] = tokens[i].trim();
        }
        
        return tokens;
    }
    
    protected final String[] tokenizeSecondLevel(final String s){
        
        if(s.trim().isEmpty()){
            return new String[0];
        }
        
        final String[] tokens = s.trim().split("\\;");
        for(int i = 0; i < tokens.length; i++){
            tokens[i] = tokens[i].trim();
        }
        
        return tokens;
    }
    
    protected final String[] tokenizeThirdLevel(final String s){
        
        if(s.trim().isEmpty()){
            return new String[0];
        }
        
        final String[] tokens = s.trim().split("\\,");
        for(int i = 0; i < tokens.length; i++){
            tokens[i] = tokens[i].trim();
        }
        
        return tokens;
    }
    
    protected final String[] tokenizeFourthLevel(final String s){
        
        if(s.trim().isEmpty()){
            return new String[0];
        }
        
        final String[] tokens = s.trim().split("\\/");
        for(int i = 0; i < tokens.length; i++){
            tokens[i] = tokens[i].trim();
        }
        
        return tokens;
    }
    
    protected final String stringToken(final String key, final String token){
        return token.substring(key.length()).trim();
    }
    
    protected final boolean booleanToken(final String key, final String token){
        return Boolean.parseBoolean(token.substring(key.length()).trim());
    }
    
    protected final int integerToken(final String key, final String token){
        return Integer.parseInt(token.substring(key.length()).trim());
    }
    
    protected final double doubleToken(final String key, final String token){
        return Double.parseDouble(token.substring(key.length()).trim());
    }
}