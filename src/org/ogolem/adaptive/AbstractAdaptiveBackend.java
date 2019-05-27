/**
Copyright (c) 2015, 2017, J. M. Dieterich and B. Hartke
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

import java.util.ArrayList;
import java.util.List;
import org.ogolem.adaptive.genericfitness.GenericReferencePoint;
import org.ogolem.adaptive.genericfitness.PropertyCalculator;
import org.ogolem.adaptive.genericfitness.ReferenceInputData;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.properties.BulkModulus;
import org.ogolem.properties.CellVolume;
import org.ogolem.properties.DeltaGauge;
import org.ogolem.properties.Energy;
import org.ogolem.properties.Forces;
import org.ogolem.properties.Property;

/**
 * An abstract adaptivable.
 * @author Johannes Dieterich
 * @version 2017-09-14
 */
public abstract class AbstractAdaptiveBackend implements AdaptiveBackend {

    private static final long serialVersionUID = (long) 20151028;
    
    AbstractAdaptiveBackend(){
    }
    
    @Override
    public abstract AbstractAdaptiveBackend clone();

    @Override
    public abstract double energyOfStructWithParams(final CartesianCoordinates cartes, final AdaptiveParameters params, final int geomID, final BondInfo bonds);

    @Override
    public abstract double gradientOfStructWithParams(final CartesianCoordinates cartes, final AdaptiveParameters params, final int geomID, final BondInfo bonds, final double[] grad);

    @Override
    public abstract double[][] minMaxBordersForParams(final AdaptiveParameters params);

    @Override
    public abstract AdaptiveParameters createInitialParameterStub(final ArrayList<CartesianCoordinates> refCartes, final String sMethod);

    @Override
    public <T extends Property, V extends ReferenceInputData<T>> PropertyCalculator<T,V> getCalculatorForProperty(final T property, final V data) {
        System.out.println("INFO: Default implementation of getCalculatorForProperty from abstract super class AbstractAdaptiveBackend called.");
        return null;
    }
    
    @SuppressWarnings("unchecked")
    @Override
    public List<? extends Property> runAllPropertyCalcs(final AdaptiveParameters params, final List<GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>>> referencePoints){
        
        final List<Property> allProps = new ArrayList<>();
        
        for(final GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>> refPoint : referencePoints){
            
            // figure type of property out
            final Property p = refPoint.getReferenceProperty();
            if(p instanceof Energy){
                
                final Energy typedP = (Energy) p;
                final ReferenceInputData<Energy> data = (ReferenceInputData<Energy>) refPoint.getReferenceInputData();
                
                // get the property calculator
                final PropertyCalculator<Energy,ReferenceInputData<Energy>> calc = getCalculatorForProperty(typedP,data);
                
                // calculate property
                final Energy calcP = calc.calculateProperty(params, data);
            
                // add
                allProps.add(calcP);

            } else if(p instanceof Forces){
                
                final Forces typedP = (Forces) p;
                final ReferenceInputData<Forces> data = (ReferenceInputData<Forces>) refPoint.getReferenceInputData();
                
                // get the property calculator
                final PropertyCalculator<Forces,ReferenceInputData<Forces>> calc = getCalculatorForProperty(typedP,data);
                
                // calculate property
                final Forces calcP = calc.calculateProperty(params, data);
            
                // add
                allProps.add(calcP);
                
            } else if(p instanceof BulkModulus){
                
                final BulkModulus typedP = (BulkModulus) p;
                final ReferenceInputData<BulkModulus> data = (ReferenceInputData<BulkModulus>) refPoint.getReferenceInputData();
                
                // get the property calculator
                final PropertyCalculator<BulkModulus,ReferenceInputData<BulkModulus>> calc = getCalculatorForProperty(typedP,data);
                
                // calculate property
                final BulkModulus calcP = calc.calculateProperty(params, data);
            
                // add
                allProps.add(calcP);
                
            } else if(p instanceof CellVolume){
                
                final CellVolume typedP = (CellVolume) p;
                final ReferenceInputData<CellVolume> data = (ReferenceInputData<CellVolume>) refPoint.getReferenceInputData();
                
                // get the property calculator
                final PropertyCalculator<CellVolume,ReferenceInputData<CellVolume>> calc = getCalculatorForProperty(typedP,data);
                
                // calculate property
                final CellVolume calcP = calc.calculateProperty(params, data);
            
                // add
                allProps.add(calcP);
                
            } else if(p instanceof DeltaGauge){
                
                final DeltaGauge typedP = (DeltaGauge) p;
                final ReferenceInputData<DeltaGauge> data = (ReferenceInputData<DeltaGauge>) refPoint.getReferenceInputData();
                
                // get the property calculator
                final PropertyCalculator<DeltaGauge,ReferenceInputData<DeltaGauge>> calc = getCalculatorForProperty(typedP,data);
                
                // calculate property
                final DeltaGauge calcP = calc.calculateProperty(params, data);
            
                // add
                allProps.add(calcP);
                
            } else {
                // error
                throw new RuntimeException("Unknown property type " + p.printableProperty() + ".");
            }
        }
        
        return allProps;
    }
}
