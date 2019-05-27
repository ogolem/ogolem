/**
Copyright (c) 2012-2014, J. M. Dieterich
              2015-2018, J. M. Dieterich and B. Hartke
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

import static java.lang.Math.abs;
import java.util.List;
import org.ogolem.adaptive.AdaptiveParameters;
import static org.ogolem.core.FixedValues.NONCONVERGEDENERGY;
import org.ogolem.properties.Property;

/**
 * A generic, property based fitness term.
 * @author Johannes Dieterich
 * @version 2018-01-11
 */
public class SerialGenericFitnessTerm<T extends Property, V extends ReferenceInputData<T>> implements GenericFitnessTerm<T> {
    
    private static final long serialVersionUID = (long) 20180111;
    private static final boolean DEBUG = false;
    private final List<GenericReferencePoint<T,V>> referenceValues;
    private final PropertyCalculator<T,V> calculator;
    private final boolean useMaxAllowedDiffs;
    private final boolean useOnlyMaxAllowedDiffs;
    private final boolean normalizePenality;
    private final boolean exactComp;
    private final boolean refToFirst;
    private final double increaseConstant;
    private final double penaltyPow;
    private final boolean printContribution;
    
    public SerialGenericFitnessTerm(final List<GenericReferencePoint<T,V>> referenceValues, final PropertyCalculator<T,V> calc,
            final FitnessTermConfig<T> config, final boolean printContribution){
        assert(referenceValues != null);
        assert(config != null);
        assert(calc != null);
        this.referenceValues = referenceValues;
        this.calculator = calc;
        this.useMaxAllowedDiffs = config.useMaxAllowedDiffs();
        this.useOnlyMaxAllowedDiffs = config.useOnlyMaxAllowedDiffs();
        if(useOnlyMaxAllowedDiffs && !useMaxAllowedDiffs) {System.err.println("ERROR: Using useOnlyMaxAllowedDiffs and NOT using useMaxAllowedDiffs is useless!");}
        this.normalizePenality = config.normalizePenality();
        this.exactComp = config.doExactEnergyComp();
        this.refToFirst = config.doReferenceToFirst();
        this.increaseConstant = config.getPenaltyPrefactor();
        this.penaltyPow = config.getPenaltyPower();
        this.printContribution = printContribution;
    }
    
    /**
     * Copy constructor providing a SHALLOW copy.
     * @param orig 
     */
    public SerialGenericFitnessTerm(final SerialGenericFitnessTerm<T,V> orig){
        assert(orig != null);
        this.referenceValues = orig.referenceValues;
        this.calculator = orig.calculator.clone();
        this.useMaxAllowedDiffs = orig.useMaxAllowedDiffs;
        this.useOnlyMaxAllowedDiffs = orig.useOnlyMaxAllowedDiffs;
        this.normalizePenality = orig.normalizePenality;
        this.exactComp = orig.exactComp;
        this.refToFirst = orig.refToFirst;
        this.increaseConstant = orig.increaseConstant;
        this.penaltyPow = orig.penaltyPow;
        assert(referenceValues != null);
        this.printContribution = orig.printContribution;
    }
    
    @Override
    public SerialGenericFitnessTerm<T,V> clone(){
        return new SerialGenericFitnessTerm<>(this);
    }
    
    @Override
    public double calculateFitnessForProperty(final AdaptiveParameters p){
        
        assert(p != null);
        
        T prop1 = null;
        // only used if(!refToFirst)
        int offset = -1;
        
        double fitness = 0.0;
        int geomC = -1;
        DataLoop: for(final GenericReferencePoint<T,V> point : referenceValues){
            
            geomC++;
            
            if(DEBUG) {System.out.println("DEBUG: Before loop " + geomC + " fitness is " + fitness + ".");}
            
            final V inpData = point.getReferenceInputData();
            final T refProp = point.getReferenceProperty();
            final double refWeight = point.getRefWeight();
            final double refMaxDiff = point.getMaxAllowedDiff();
            
            
            final T prop = calculator.calculateProperty(p, inpData);
            if(DEBUG) {System.out.println("DEBUG: Property for reference " + geomC + " is " + prop.printableProperty());}
            final boolean insane = prop.makeSensible();
            
            if(printContribution){
                System.out.println("\tProperty is:              " + prop.name());
                System.out.println("\tValue is:                 " + prop.printableProperty());
                System.out.println("\tReference property value: " + refProp.printableProperty());
            }
            
            if(refToFirst){
                if(geomC == 0){
                    prop1 = prop;
                    continue DataLoop;
                }
                
                // compute differences
                assert(prop1 != null);
                final double diffFit = prop.signedDifference(prop1);
                final T zeroProp = referenceValues.get(0).getReferenceProperty();
                final double diffRef = refProp.signedDifference(zeroProp);
                final double delta = abs(diffFit-diffRef);

                if (useMaxAllowedDiffs && delta >= refMaxDiff) {
                    final double fx = GenericFitnessFunction.polynomialPenaltyFunction(diffRef, diffFit, increaseConstant, penaltyPow, normalizePenality);
                    if(printContribution){
                        System.out.println("\t I Added " + fx + " b/c larger than max allowed diff.");
                    }
                    fitness += fx;
                }
                if(!useOnlyMaxAllowedDiffs) {
                    final double fx = abs(refWeight * delta);
                    if(printContribution){
                        System.out.println("\t II Added " + fx + " for difference of properties.");
                    }
                    fitness += fx;
                }
            } else {
                
                if (insane) {
                    fitness += NONCONVERGEDENERGY;
                    if(printContribution){
                        System.out.println("\t Added " + NONCONVERGEDENERGY + " b/c insane.");
                    }
                    continue DataLoop;
                } else if (!insane && prop1 == null) {
                    offset = geomC;
                    prop1 = prop;
                    if (!exactComp) {
                        // we do not need to add the difference, it is considered to be the zero point
                    } else {
                        final double delta = refProp.absoluteDifference(prop);
                        if (useMaxAllowedDiffs && delta >= refMaxDiff) {
                            final double fx = GenericFitnessFunction.polynomialPenaltyFunction(refProp, prop, increaseConstant, penaltyPow, normalizePenality);
                            if(printContribution){
                                System.out.println("\t III Added " + fx + " b/c larger than max allowed diff.");
                            }
                            fitness += fx;
                        }
                        if(!useOnlyMaxAllowedDiffs) {
                            final double fx = refWeight * delta;
                            if(printContribution){
                                System.out.println("\t IV Added " + fx + " for difference of properties. Weight employed is " + refWeight);
                            }
                            fitness += fx;
                        }
                    }
                    continue DataLoop;
                } else if (!insane && prop1 != null) {
                    if (!exactComp) {
                        // so now we create the deviation to the reference value relative to the first realistic and sum it up
                        // just compare the shapes, not the absolute values
                        final double delta1 = prop.signedDifference(prop1);
                        final double delta2 = refProp.signedDifference(referenceValues.get(offset).getReferenceProperty());
                        if (useMaxAllowedDiffs && (abs(delta1-delta2) >= refMaxDiff)) {
                            final double fx = GenericFitnessFunction.polynomialPenaltyFunction(delta2, delta1, increaseConstant, penaltyPow, normalizePenality);
                            if(printContribution){
                                System.out.println("\t V Added " + fx + " b/c larger than max allowed diff.");
                            }
                            fitness += fx;
                        }
                        if(!useOnlyMaxAllowedDiffs) {
                            final double fx = abs(refWeight * delta1-delta2);
                            if(printContribution){
                                System.out.println("\t VI Added " + fx + " for difference of properties.");
                            }
                            fitness += fx;
                        }
                    } else {
                        final double delta = refProp.absoluteDifference(prop);
                        if (useMaxAllowedDiffs && delta >= refMaxDiff) {
                            final double fx = GenericFitnessFunction.polynomialPenaltyFunction(refProp, prop, increaseConstant, penaltyPow, normalizePenality);
                            if(printContribution){
                                System.out.println("\t VII Added " + fx + " b/c larger than max allowed diff.");
                            }
                            fitness += fx;
                        }
                        // we do an exact computation
                        if(!useOnlyMaxAllowedDiffs) {
                            final double fx = refWeight * delta;
                            if(printContribution){
                                System.out.println("\t VIII Added " + fx + " for difference of properties.");
                            }
                            fitness += fx;
                        }
                    }
                    continue DataLoop;
                } else {
                    System.err.println("ERROR: You really shouldn't end up here in the generic fitness term calculation, contact the author.");
                    fitness += NONCONVERGEDENERGY;
                    continue;
                }
            }
        }

        return fitness;
    }
    
    @Override
    public double calculateGradientForProperty(final AdaptiveParameters p, final double[] gradient){
        
        assert(p != null);
        assert(gradient != null);
        
        T prop1 = null;
        double[] grad1 = null;
        // only used if(!refToFirst)
        int offset = -1;
        
        final double[] g = new double[p.getNumberOfParamters()];
        double fitness = 0.0;
        int geomC = -1;
        DataLoop: for(final GenericReferencePoint<T,V> point : referenceValues){
            
            geomC++;
            
            final V inpData = point.getReferenceInputData();
            final T refProp = point.getReferenceProperty();
            final double refWeight = point.getRefWeight();
            final double refMaxDiff = point.getMaxAllowedDiff();
            
            final T prop = calculator.calculatePropertyGradient(p, inpData, g);
            final boolean insane = prop.makeSensible();
            GenericFitnessFunction.makeSensible(g);
            
            if(printContribution){
                System.out.println("\tProperty is:              " + prop.name());
                System.out.println("\tValue is:                 " + prop.printableProperty());
                System.out.println("\tReference property value: " + refProp.printableProperty());
            }
            
            if(refToFirst){
                if(geomC == 0){
                    prop1 = prop;
                    grad1 = g.clone();
                    continue;
                }
                
                // compute differences
                assert(prop1 != null);
                assert(grad1 != null);
                final double diffFit = prop.signedDifference(prop1);
                final T zeroProp = referenceValues.get(0).getReferenceProperty();
                final double diffRef = refProp.signedDifference(zeroProp);
                final double delta = abs(diffFit-diffRef);

                if (useMaxAllowedDiffs && delta >= refMaxDiff) {
                    final double fx = GenericFitnessFunction.polynomialPenaltyFunction(diffRef, diffFit, increaseConstant, penaltyPow, normalizePenality);
                    if(printContribution){
                        System.out.println("\t I Added " + fx + " b/c larger than max allowed diff.");
                    }
                    fitness += fx;
                }
                if(!useOnlyMaxAllowedDiffs) {
                    final double fx = abs(refWeight * delta);
                    if(printContribution){
                        System.out.println("\t II Added " + fx + " for difference of properties.");
                    }
                    fitness += fx;
                }
                // adjust and add the gradient
                for(int i = 0; i < g.length; i++) {gradient[i] += refWeight*(g[i] - grad1[i]);}
            } else {
                
                if (insane) {
                    fitness += NONCONVERGEDENERGY;
                    if(printContribution){
                        System.out.println("\t Added " + NONCONVERGEDENERGY + " b/c insane.");
                    }
                    continue DataLoop;
                } else if (!insane && prop1 == null) {
                    offset = geomC;
                    prop1 = prop;
                    if (!exactComp) {
                        // we do not need to add the difference, it is considered to be the zero point
                    } else {
                        final double delta = refProp.absoluteDifference(prop);
                        if (useMaxAllowedDiffs && delta >= refMaxDiff) {
                            final double fx = GenericFitnessFunction.polynomialPenaltyFunction(refProp, prop, increaseConstant, penaltyPow, normalizePenality);
                            if(printContribution){
                                System.out.println("\t III Added " + fx + " b/c larger than max allowed diff.");
                            }
                            fitness += fx;
                        }
                        if(!useOnlyMaxAllowedDiffs) {
                            final double fx = refWeight * delta;
                            if(printContribution){
                                System.out.println("\t IV Added " + fx + " for difference of properties. Weight employed is " + refWeight);
                            }
                            fitness += fx;
                        }
                    }
                    for(int i = 0; i < g.length; i++) {gradient[i] += refWeight*g[i];} //XXX think about this...
                    continue DataLoop;
                } else if (!insane && prop1 != null) {
                    if (!exactComp) {
                        // so now we create the deviation to the reference value relative to the first realistic and sum it up
                        // just compare the shapes, not the absolute values
                        final double delta1 = prop.signedDifference(prop1);
                        final double delta2 = refProp.signedDifference(referenceValues.get(offset).getReferenceProperty());
                        if (useMaxAllowedDiffs && (abs(delta1-delta2) >= refMaxDiff)) {
                            final double fx = GenericFitnessFunction.polynomialPenaltyFunction(delta2, delta1, increaseConstant, penaltyPow, normalizePenality);
                            if(printContribution){
                                System.out.println("\t V Added " + fx + " b/c larger than max allowed diff.");
                            }
                            fitness += fx;
                        }
                        if(!useOnlyMaxAllowedDiffs) {
                            final double fx = abs(refWeight * delta1-delta2);
                            if(printContribution){
                                System.out.println("\t VI Added " + fx + " for difference of properties.");
                            }
                            fitness += fx;
                        }
                    } else {
                        final double delta = refProp.absoluteDifference(prop);
                        if (useMaxAllowedDiffs && delta >= refMaxDiff) {
                            final double fx = GenericFitnessFunction.polynomialPenaltyFunction(refProp, prop, increaseConstant, penaltyPow, normalizePenality);
                            if(printContribution){
                                System.out.println("\t VII Added " + fx + " b/c larger than max allowed diff.");
                            }
                            fitness += fx;
                        }
                        // we do an exact computation
                        if(!useOnlyMaxAllowedDiffs) {
                            final double fx = refWeight * delta;
                            if(printContribution){
                                System.out.println("\t VIII Added " + fx + " for difference of properties.");
                            }
                            fitness += fx;
                        }
                    }
                    for(int i = 0; i < g.length; i++) {gradient[i] += refWeight*g[i];}
                    continue DataLoop;
                } else {
                    System.err.println("ERROR: You really shouldn't end up here in the generic fitness term calculation, contact the author.");
                    fitness += NONCONVERGEDENERGY;
                    for(int i = 0; i < g.length; i++) {gradient[i] += refWeight*g[i];}
                    continue DataLoop;
                }
            }
            
        }
                
        return fitness;
    }
    
    PropertyCalculator<T,V> getCalculatorReference(){
        return calculator;
    }
}
