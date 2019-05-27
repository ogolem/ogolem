/**
Copyright (c) 2012-2013, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
              2018, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive.genericfitness;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.math3.util.FastMath;
import org.ogolem.adaptive.AdaptiveParameters;
import static org.ogolem.adaptive.FixedValues.NONCONVERGEDGRADIENT;
import org.ogolem.adaptive.ParameterGradient;
import static org.ogolem.core.FixedValues.NONCONVERGEDENERGY;
import org.ogolem.properties.Property;

/**
 * A generic fitness function for parametrization purposes.
 * @author Johannes Dieterich
 * @version 2018-01-11
 */
public class GenericFitnessFunction implements Cloneable, Serializable {
    
    private static final long serialVersionUID = (long) 20180111;
    private static final boolean EXACTPOW = false;
    private final List<GenericFitnessTerm<?>> terms;
    private final List<Double> termWeights;
    private final BatchedPropertyCalculator batcher;
    private final boolean printContributions;
    
    @SuppressWarnings("rawtypes")
    public GenericFitnessFunction(final BatchedPropertyCalculator batcher, final List<GenericFitnessTerm<?>> terms, final List<Double> termWeights, final boolean printContributions){
        this.terms = terms;
        this.termWeights = termWeights;
        assert(terms.size() == termWeights.size());
        assert(batcher != null);
        this.batcher = batcher;        
        for(final GenericFitnessTerm<?> term : terms){
            if(term instanceof SerialGenericFitnessTerm){
                final SerialGenericFitnessTerm t = (SerialGenericFitnessTerm) term;
                final PropertyCalculator<? extends Property,? extends ReferenceInputData<?>> calc = t.getCalculatorReference();
                if(calc instanceof PseudoPropertyCalculator){
                    final PseudoPropertyCalculator pseudo = (PseudoPropertyCalculator) calc;
                    pseudo.setBatcherReference(batcher);
                }
            }
        }
        this.printContributions = printContributions;
    }

    /**
     * Copy constructor returning a shallow copy.
     * @param orig 
     */
    @SuppressWarnings("rawtypes")
    private GenericFitnessFunction(final GenericFitnessFunction orig){
        this.termWeights = orig.termWeights;
        this.batcher = orig.batcher.clone();
        this.terms = new ArrayList<>();        
        for(final GenericFitnessTerm<?> termOrig : orig.terms){
            final GenericFitnessTerm<?> term = termOrig.clone();
            terms.add(term);
            if(term instanceof SerialGenericFitnessTerm){
                final SerialGenericFitnessTerm t = (SerialGenericFitnessTerm) term;
                final PropertyCalculator<? extends Property,? extends ReferenceInputData<?>> calc = t.getCalculatorReference();
                if(calc instanceof PseudoPropertyCalculator){
                    final PseudoPropertyCalculator pseudo = (PseudoPropertyCalculator) calc;
                    pseudo.setBatcherReference(batcher);
                }
            }
        }
        this.printContributions = orig.printContributions;
    }
    
    @Override
    public GenericFitnessFunction clone(){
        return new GenericFitnessFunction(this);
    }
    
    public double evaluateFitness(final AdaptiveParameters parameters){
        
        // recall the batcher!
        batcher.recalcForNewParameters(parameters);
        
        if(printContributions){System.out.println("CONTRIBUTIONS:");}
        
        double d = 0.0;
        for(int i = 0; i < terms.size(); i++){
            final GenericFitnessTerm<?> term = terms.get(i);
            final double weight = termWeights.get(i);
            if(printContributions){
                System.out.println("Next fitness contribution:");
                System.out.println("\tWeight is: " + weight);
            }
            d += weight * term.calculateFitnessForProperty(parameters);
            if(printContributions){System.out.println("-----------------------------------------");}
        }
        
        d = makeSensible(d);
        
        return d;
    }
    
    public ParameterGradient evaluateGradient(final AdaptiveParameters parameters){
        
        // recall the batcher! XXX not yet implemented
        // batcher.recalcForNewParameters(parameters);
        
        if(printContributions){System.out.println("CONTRIBUTIONS:");}
        
        double totFit = 0.0;
        final double[] totGrad = new double[parameters.getNumberOfParamters()];
        final double[] tmpGrad = new double[parameters.getNumberOfParamters()];
        for(int i = 0; i < terms.size(); i++){
            for(int j = 0; j < tmpGrad.length; j++){tmpGrad[j] = 0.0;}
            final GenericFitnessTerm<?> term = terms.get(i);
            final double weight = termWeights.get(i);
            if(printContributions){
                System.out.println("Next fitness contribution:");
                System.out.println("\tWeight is: " + weight);
            }
            totFit += weight*term.calculateGradientForProperty(parameters, tmpGrad);
            for(int j = 0; j < tmpGrad.length; j++){totGrad[j] += weight*tmpGrad[j];}
            if(printContributions){System.out.println("-----------------------------------------");}
        }
        
        final ParameterGradient total = new ParameterGradient();
        total.setGradientTotal(totGrad);
        total.setTotalEnergy(totFit);
        makeSensible(total);

        return total;
    }
    
    static double makeSensible(double fit){
        if(Double.isInfinite(fit) || Double.isNaN(fit) || Double.compare(Double.MAX_VALUE,fit) <= 0){ return NONCONVERGEDENERGY;}
        else{ return fit;}
    }
    
    static void makeSensible(final double[] g){
        
        for(int i = 0; i < g.length; i++){
            if(Double.isInfinite(g[i]) || Double.isNaN(g[i]) || Double.compare(Double.MAX_VALUE,g[i]) <= 0){
                g[i] = NONCONVERGEDGRADIENT;
            }
        }
    }
    
    static void makeSensible(final ParameterGradient grad){

        final double[] g = grad.getTotalGradient();
        for(int i = 0; i < g.length;i++){
            if(Double.isInfinite(g[i]) || Double.isNaN(g[i]) || Double.compare(Double.MAX_VALUE,g[i]) <= 0){
                g[i] = NONCONVERGEDGRADIENT;
            }
        }
    }
    
    static double polynomialPenaltyFunction(final double optimal,
            final double actual, final double increase, final double pow,
            final boolean normalize){

        final double x = Math.abs(actual-optimal)*100/Math.abs(optimal);
        final double d = (EXACTPOW) ? Math.pow(x,pow)*increase : FastMath.pow(x,pow)*increase;

        if(Double.isInfinite(d)|| Double.isNaN(d)){
            // make sensible
            return NONCONVERGEDENERGY;
        }

        return d;
    }
    
    /**
     * @param optimal
     * @param actual
     * @param increase
     * @param pow
     * @param normalize whether we normalize with getValue(). Must only be true if getValue for this Property works!
     * @return the penalty
     */
    static double polynomialPenaltyFunction(final Property optimal,
            final Property actual, final double increase, final double pow, final boolean normalize){

        final double x = (normalize) ? optimal.absoluteDifference(actual)*100/Math.abs(optimal.getValue()) : optimal.absoluteDifference(actual);
        final double d = (EXACTPOW) ? Math.pow(x,pow)*increase : FastMath.pow(x,pow)*increase;

        if(Double.isInfinite(d)|| Double.isNaN(d)){
            // make sensible
            return NONCONVERGEDENERGY;
        }

        return d;
    }
}
