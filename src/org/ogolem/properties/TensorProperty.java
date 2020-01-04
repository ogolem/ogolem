/**
Copyright (c) 2017, J. M. Dieterich and B. Hartke
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
package org.ogolem.properties;

import org.ogolem.core.FixedValues;

/**
 * The base class for a tensor property
 * @author Johannes Dieterich
 * @version 2017-12-15
 */
public abstract class TensorProperty implements Property{
    
    private static final long serialVersionUID = (long) 20171215;
    
    protected final boolean normDifferences;
    protected final double[][][] data;
    
    protected TensorProperty(final double[][][] data, final boolean normDifferences){
        this.normDifferences = normDifferences;
        this.data = data;
    }
    
    protected TensorProperty(final TensorProperty orig){
        this.normDifferences = orig.normDifferences;
        if(orig.data != null){
        
            this.data = new double[orig.data.length][][];
            for(int i = 0; i < orig.data.length; i++){
                this.data[i] = orig.data[i].clone();
            }
        } else {
            this.data = null;
        }
    }
    
    @Override
    public abstract TensorProperty clone();
    
    /**
     * Return the absolute sum of all matrix elements. If this is not wished overriding is necessary.
     * @return the pseudo-vector norm
     */
    @Override
    public double getValue(){
        
        if(data == null){return FixedValues.NONCONVERGEDENERGY;}
        
        double sum = 0.0;
        for(int i = 0; i < data.length; i++){
            for(int j = 0; j < data[i].length; j++){
                for(int k = 0; k < data[i][j].length; k++){
                    sum += Math.abs(data[i][j][k]);
                }
            }
        }
        
        return sum;
    }
    
    @Override
    public double signedDifference(Property p){
        
        if (!ensureCorrectProperty(p)) {throw new IllegalArgumentException("Property should be an instance of " + name());}

        final TensorProperty tp = (TensorProperty) p;
        
        if(this.data == null || tp.data == null){return FixedValues.NONCONVERGEDENERGY;}
        
        if(tp.data.length != this.data.length){throw new RuntimeException("Tensor properties differ in lengths I: " + tp.data.length + " vs " + this.data.length);}
        if(tp.data[0].length != this.data[0].length){throw new RuntimeException("Tensor properties differ in lengths II: " + tp.data[0].length + " vs " + this.data[0].length);}
        if(tp.data[0][0].length != this.data[0][0].length){throw new RuntimeException("Tensor properties differ in lengths III: " + tp.data[0][0].length + " vs " + this.data[0].length);}
        
        double diff = 0.0;
        for (int i = 0; i < data.length; i++) {
            if(data[0].length != data[i].length || data[i].length != tp.data[i].length){throw new IllegalArgumentException("First matrix dimension must be the same! Are " + data[i].length + " vs " + tp.data[i].length);}
            for (int j = 0; j < data[0].length; j++) {
                if(data[0][0].length != data[0][i].length || data[0][i].length != tp.data[0][i].length){throw new IllegalArgumentException("Second matrix dimension must be the same! Are " + data[i][j].length + " vs " + tp.data[i][j].length);}
                for (int k = 0; k < data[0][0].length; k++) {
                    diff += (this.data[i][j][k] - tp.data[i][j][k]); // XXX this is suboptimal, but for the time being I have no better definition in my head
                }
            }
        }
        
        if(this.normDifferences){
            diff /= (data.length * data[0].length);
        }

        return diff;
    }
    
    @Override
    public double absoluteDifference(Property p){
        
        if (!ensureCorrectProperty(p)) {throw new IllegalArgumentException("Property should be an instance of " + name());}

        final TensorProperty tp = (TensorProperty) p;
        
        if(this.data == null || tp.data == null){return FixedValues.NONCONVERGEDENERGY;}
        
        if(tp.data.length != this.data.length){throw new RuntimeException("Tensor properties differ in lengths I: " + tp.data.length + " vs " + this.data.length);}
        if(tp.data[0].length != this.data[0].length){throw new RuntimeException("Tensor properties differ in lengths II: " + tp.data[0].length + " vs " + this.data[0].length);}
        if(tp.data[0][0].length != this.data[0][0].length){throw new RuntimeException("Tensor properties differ in lengths III: " + tp.data[0][0].length + " vs " + this.data[0].length);}
        
        double diff = 0.0;
        for (int i = 0; i < data.length; i++) {
            if(data[0].length != data[i].length || data[i].length != tp.data[i].length){throw new IllegalArgumentException("First matrix dimension must be the same! Are " + data[i].length + " vs " + tp.data[i].length);}
            for (int j = 0; j < data[0].length; j++) {
                if(data[0][0].length != data[0][i].length || data[0][i].length != tp.data[0][i].length){throw new IllegalArgumentException("Second matrix dimension must be the same! Are " + data[i][j].length + " vs " + tp.data[i][j].length);}
                for (int k = 0; k < data[0][0].length; k++) {
                    diff += Math.abs(this.data[i][j][k] - tp.data[i][j][k]);
                }
            }
        }
        
        if(this.normDifferences){
            diff /= (data.length * data[0].length * data[0][0].length);
        }

        return diff;
    }
    
    protected abstract boolean ensureCorrectProperty(final Property p);
}
