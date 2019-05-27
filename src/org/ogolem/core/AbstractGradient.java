/**
Copyright (c) 2012, J. M. Dieterich
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
package org.ogolem.core;

import java.util.List;

/**
 * An abstract gradient.
 * @author Johannes Dieterich
 * @version 2012-03-09
 */
public class AbstractGradient {
    
    protected double[] gradientValues;
    protected double functionValue;
    
    public AbstractGradient(){
    }
    
    public AbstractGradient(final int dim){
        this.gradientValues = new double[dim];
    }
    
    public AbstractGradient(final double val, final double[] grad){
        this.functionValue = val;
        this.gradientValues = grad;
    }
    
    public AbstractGradient(final List<AbstractGradient> grads, final int dim){
        gradientValues = new double[dim];
        for(final AbstractGradient grad : grads){
            functionValue += grad.functionValue;
            final double[] g = grad.gradientValues;
            for(int i = 0; i < dim; i++){
                gradientValues[i] += g[i];
            }
        }
    }
    
    public void setGradient(final double[] grad){
        this.gradientValues = grad;
    }
    
    public double[] getGradient(){
        return gradientValues;
    }
    
    public void setFunctionValue(final double val){
        functionValue = val;
    }
    
    public double getFunctionValue(){
        return functionValue;
    }
}
