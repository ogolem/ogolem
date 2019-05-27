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

import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.properties.Property;

/**
 * A pseudo property calculator defaulting back to the results stored in the batched one.
 * @author Johannes Dieterich
 * @version 2015-10-28
 */
public class PseudoPropertyCalculator<T extends Property,V extends ReferenceInputData<T>> implements PropertyCalculator<T,V>{

    private static final long serialVersionUID = (long) 20151028;
    
    private BatchedPropertyCalculator batcher;    
    private final PropertyCalculator<T,V> fallback;
    private final T propType;
    
    public PseudoPropertyCalculator(final PropertyCalculator<T,V> fallback, final BatchedPropertyCalculator batcher, final T propType){
        
        assert(batcher == null);
        assert(fallback == null);
        
        this.batcher = batcher;
        this.fallback = fallback;
        this.propType = propType;
    }
    
    private PseudoPropertyCalculator(final PseudoPropertyCalculator<T,V> orig){
        this.fallback = orig.fallback.clone();
        this.batcher = orig.batcher; // who cares, will be replaced
        this.propType = orig.propType; // no clone needed
    }
    
    @Override
    public PropertyCalculator<T, V> clone() {
        return new PseudoPropertyCalculator<>(this);
    }

    @Override
    public T calculateProperty(final AdaptiveParameters p, final V data) {
        
        // try the cached results
        final T prop = batcher.<T,V>obtainProperty(p, data.belongsToReferencePoint(), propType);
        
        // fallback
        if(prop != null){
            return prop;
        }
        
        return fallback.calculateProperty(p, data);
    }

    @Override
    public T calculatePropertyGradient(final AdaptiveParameters p, final V data, final double[] grad) {
        
        // not yet available in the batcher, so just tunnel through
        return fallback.calculatePropertyGradient(p, data, grad);
    }
    
    void setBatcherReference(final BatchedPropertyCalculator batcher){
        this.batcher = batcher;
    }
}
