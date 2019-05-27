/**
Copyright (c) 2015, J. M. Dieterich and B. Hartke
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
import org.ogolem.adaptive.Adaptivable;
import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.properties.Property;

/**
 * A cache for batched property calculations.
 * @author Johannes Dieterich
 * @version 2015-10-29
 */
public class SerialBatchedPropertyCalculator implements BatchedPropertyCalculator {
    
    private static final long serialVersionUID = (long) 20151028;
    
    private final Adaptivable adaptivable;
    private final List<List<? extends Property>> calculatedProperties;
    private final List<PropertyBatch> batches;
    
    private long uniqueParamID;
    
    public SerialBatchedPropertyCalculator(final Adaptivable adaptivable, final List<PropertyBatch> refBatches){
        this.uniqueParamID = 0;
        this.calculatedProperties = new ArrayList<>();
        this.batches = refBatches;
        this.adaptivable = adaptivable;
    }
    
    private SerialBatchedPropertyCalculator(final SerialBatchedPropertyCalculator orig){
        this.uniqueParamID = 0;
        this.calculatedProperties = new ArrayList<>();
        this.adaptivable = orig.adaptivable.clone();
        this.batches = orig.batches; // no clone needed
    }
    
    @Override
    public SerialBatchedPropertyCalculator clone(){
        return new SerialBatchedPropertyCalculator(this);
    }
    
    @Override
    public void recalcForNewParameters(final AdaptiveParameters params){
        this.uniqueParamID = params.getUniqueID();
        
        // do all the needed calculations
        calculatedProperties.clear();
        
        for(final PropertyBatch batch : batches){
            final List<? extends Property> props = adaptivable.runAllPropertyCalcs(params, batch.getBatch());
            final int batchID = batch.getBatchID();
            calculatedProperties.add(batchID, props);
        }
    }
    
    @SuppressWarnings("unchecked")
    @Override
    public <T extends Property,V extends ReferenceInputData<T>> T obtainProperty(final AdaptiveParameters params, final int refID, final T propertyType){
        
        final long myUniqueID = params.getUniqueID();
        if(myUniqueID != uniqueParamID){
            throw new RuntimeException("No recalculation of properties has taken place prior to asking for property " + propertyType.name()
            + " Have " + uniqueParamID + " should be " + myUniqueID + " for " + this.hashCode());
        }
        
        // simply get the proper entry in the list
        final List<? extends Property> allPropsForPoint = calculatedProperties.get(refID);
        if(allPropsForPoint == null){
            throw new RuntimeException("No properties calculated for reference point " + refID + ". Logic error.");
        }
        
        // loop over them and figure the correct property out
        for(int i = 0; i < allPropsForPoint.size(); i++){
            final Property p = allPropsForPoint.get(i);
            if(p.name().equalsIgnoreCase(propertyType.name())){
                // correct one found
                return (T) p;
            }
        }
        
        System.err.println("No such property " + propertyType.name() + " precalculated. Returning null.");
        return null;
    }
    
    public static class PropertyBatch implements Serializable, Cloneable {
        
        private static final long serialVersionUID = (long) 20151028;
        
        private final List<GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>>> batch;
        private final int batchID;
        
        public PropertyBatch(final List<GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>>> refBatch,
                final int batchID){
            
            assert(!refBatch.isEmpty());
            
            this.batch = refBatch;
            this.batchID = batchID;
            
            // check!
            for(final GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>> point : batch){
                final int id = point.getReferenceID();
                if(id != batchID){
                    throw new RuntimeException("Wrong batch ID specified. Must be ID of *all* reference input data.");
                }
            }
        }
        
        private PropertyBatch(final PropertyBatch orig){
            this.batch = orig.batch; // no clone needed
            this.batchID = orig.batchID;
        }
        
        @Override
        public PropertyBatch clone(){
            return new PropertyBatch(this);
        }
        
        int getBatchID(){
            return batchID;
        }
        
        List<GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>>> getBatch(){
            return batch;
        }
    }
}
