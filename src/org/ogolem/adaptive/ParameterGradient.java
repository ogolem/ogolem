/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010     , J. M. Dieterich
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

import java.util.Iterator;
import java.util.List;
import org.ogolem.core.AbstractGradient;

/**
 * A wrapper class around the Gradient class from the core package. Existing
 * to reduce internal dependencies.
 * @author Johannes Dieterich
 * @version 2010-09-17
 */
public class ParameterGradient extends AbstractGradient{

    /**
     * Constructs a single gradient from a list of gradients.
     * @param gradients
     * @param dim The dimension of the resulting gradient.
     */
    public ParameterGradient(final List<ParameterGradient> gradients,
            final List<Double> gradientWeights, final int dim){

        final Iterator<ParameterGradient> itGrads = gradients.iterator();
        final Iterator<Double> itWeights = gradientWeights.iterator();

        this.gradientValues = new double[dim];
        this.functionValue = 0.0;

        // this will be rather slow, but that is how it is, unfortunately
        while(itGrads.hasNext() && itWeights.hasNext()){
            final ParameterGradient pgTemp = itGrads.next();

            final double tempWeight = itWeights.next();
            final double[] gradTemp = pgTemp.gradientValues;

            functionValue += tempWeight*Math.abs(pgTemp.functionValue);
            
            for(int i = 0; i < gradTemp.length; i++){
                gradientValues[i] += tempWeight*gradTemp[i];
            }
        }
    }

    /**
     * An empty legacy constructor.
     */
    public ParameterGradient(){
        super();
    }

    ParameterGradient(int iNoOfGradientRows){
        super(iNoOfGradientRows);
    }

    private static final boolean isParamGradient = true;

    boolean isParameterGradient(){
        return isParamGradient;
    }
    
    public void setTotalEnergy(final double energy){
        setFunctionValue(energy);
    }
    
    public double getTotalEnergy(){
        return getFunctionValue();
    }
    
    public double[] getTotalGradient(){
        return getGradient();
    }
    
    public void setGradientTotal(final double[] gradient){
        setGradient(gradient);
    }
}
