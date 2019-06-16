/**
Copyright (c) 2014, J. M. Dieterich
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
package org.ogolem.adaptive;

import org.ogolem.generic.BoundedGenericMutation;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericGlobalOptimization;
import org.ogolem.generic.GenericGlobalOptimizationFactory;
import org.ogolem.generic.GenericMutation;
import org.ogolem.generic.GenericSanityCheck;
import org.ogolem.generic.IndividualWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Specialized parameter factory.
 * @author Johannes Dieterich
 * @version 2014-08-11
 */
public class ParamGlobOptFactory extends GenericGlobalOptimizationFactory<Double, AdaptiveParameters>{

    private static final long serialVersionUID = (long) 20140811;
    private static final Logger log = LoggerFactory.getLogger(ParamGlobOptFactory.class);
    
    private final double[] upper;
    private final double[] lower;
    private final GenericFitnessFunction<Double,AdaptiveParameters> fitness;
    private final GenericSanityCheck<Double,AdaptiveParameters> sanity;
    private final IndividualWriter<AdaptiveParameters> writer;
    private final int noTries;
    private final double crossPoss;
    private final double mutPoss;
    private final boolean printBeforeFitness;
    
    
    public ParamGlobOptFactory(final GenericSanityCheck<Double,AdaptiveParameters> sanity,
            final GenericFitnessFunction<Double,AdaptiveParameters> fitness,
            final IndividualWriter<AdaptiveParameters> writer, final double crossPoss, final double mutPoss,
            final boolean printBeforeFitness, final int noOfTries, final double[] lowerBorder,
            final double[] upperBorder, final long totalNoGlobSteps) {
        
        super(totalNoGlobSteps);
        this.crossPoss = crossPoss;
        this.fitness = fitness;
        this.noTries = noOfTries;
        this.sanity = sanity;
        this.printBeforeFitness = printBeforeFitness;
        this.writer = writer;
        this.mutPoss = mutPoss;
        this.lower = lowerBorder;
        this.upper = upperBorder;
    }
    
    @Override
    public GenericGlobalOptimization<Double, AdaptiveParameters> translateToGlobOpt(String globOptString) throws Exception {
        
        
            // syntax: parameters{xover() mutation()}
            final String w = globOptString.trim();
            if(!w.startsWith("parameters{")){
                throw new RuntimeException("Syntax: parameters{xover()mutation()} not met!"
                        + " Starts with " + w);
            }
            final String w2 = w.substring(11, w.indexOf("}")).trim();

            if(!w2.startsWith("xover(")){
                throw new RuntimeException("Syntax: parameters{xover()mutation()} not met!");
            }
            final String xStr = w2.substring(6,w2.indexOf(")"));
            final String mStr = w2.substring(w2.lastIndexOf("(")+1,w2.lastIndexOf(")"));
                
            try{
                final GenericCrossover<Double,AdaptiveParameters> x = getXOver(xStr);
                final GenericMutation<Double,AdaptiveParameters> m = getMutation(mStr);
            
                final GenericParameterDarwin globOpt = new GenericParameterDarwin(x,m,sanity,fitness,
                    writer, crossPoss, mutPoss, printBeforeFitness, noTries);
            
                return globOpt;
                
            } catch(Exception e){    
                throw new RuntimeException("No global optimization engine found for input: " + globOptString + "."
                    + " Tried both specialized and generic ones, this is FATAL!",e);
            }
    }
        
    @Override
    protected GenericMutation<Double,AdaptiveParameters> specializedMutation(final String mutString) throws Exception {
        log.info("Working on mutation string " + mutString);
        
        if(mutString.startsWith("germany:")){
            
            int mode = BoundedGenericMutation.SINGLEMUT;
            final String[] tokens = tokenizeThirdLevel(mutString.substring(8));
            for(final String token : tokens){
                if(token.startsWith("mode=")){
                    final String s = token.substring(5).trim();
                    if(s.equalsIgnoreCase("one")){
                        mode = BoundedGenericMutation.SINGLEMUT;
                    } else if(s.equalsIgnoreCase("all")){
                        mode = BoundedGenericMutation.MULTIMUT;
                    } else{
                        throw new RuntimeException("Unknown mode " + mode + " in mutation mode for germany.");
                    }
                } else {
                    throw new RuntimeException("Unknown option for specialized mutation germany " + token);
                }
            }
            
            return new BoundedGenericMutation<>(mode,lower,upper);
        }
        
        return null;
    }
}
