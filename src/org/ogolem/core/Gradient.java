/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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
 * This is a value object containing an (however derived) energy gradient for
 * a set of cartesian coordinates. We disallow extension due to a couple of dirty
 * tricks in the code.
 * @author Johannes Dieterich
 * @version 2015-03-27
 */
public final class Gradient extends AbstractGradient {

    /*
     * Fields
     */
    private double[][] local3DGradient; // kind of a cache, so to speak

    /*
     * Constructors
     */
    public Gradient(){
        // default constructor
        this.functionValue = FixedValues.NONCONVERGEDENERGY;
    }
    
    public Gradient(final int dim1, final int dim2){
        this.functionValue = FixedValues.NONCONVERGEDENERGY;
        this.local3DGradient = new double[dim1][dim2];
    }

    public Gradient(final List<Gradient> gradients, final int[] dims){

        this.local3DGradient = new double[dims[0]][dims[1]];
        this.functionValue = 0.0;

        for(final Gradient gradient : gradients){
            final double[][] tempGrad = gradient.local3DGradient;
            for(int i = 0; i < dims[0]; i++){
                for(int j = 0; j < dims[1]; j++){
                    local3DGradient[i][j] += tempGrad[i][j];
                }
            }

            functionValue += gradient.functionValue;
        }
    }

    /*
     * Getters and Setters
     */

    public double[][] getTotalGradient() {
        return local3DGradient;
    }

    public double getTotalEnergy() {
        return getFunctionValue();
    }

    public void setGradientTotal(final double[][] grad) {
        local3DGradient = grad;
    }
    
    @Override
    public void setGradient(final double[] gradient){
        throw new RuntimeException("Method setGradient(double[]) must not be used on the two-dimensional gradient.");
    }

    public void setTotalEnergy(final double energy) {
        setFunctionValue(energy);
    }

    /*
     * Methods
     */
    @Override
    public double[] getGradient() {
        
        assert(local3DGradient.length == 3); // make sure we are working with a cartesian gradient
        assert(gradientValues == null); // need to ensure this in order to avoid information loss!
        
        final int arrayLength = local3DGradient[0].length;
        gradientValues = new double[arrayLength * 3];
        
        System.arraycopy(local3DGradient[0], 0, gradientValues, 0, arrayLength);
        System.arraycopy(local3DGradient[1], 0, gradientValues, arrayLength, arrayLength);
        System.arraycopy(local3DGradient[2], 0, gradientValues, arrayLength * 2, arrayLength);
        
        return gradientValues;
    }
    
    /**
     * Get the gradient data.
     * @param gradient a scratch object, overwritten on exit with gradient values, first all x, then all y, then all z values.
     */
    public void getGradientData(final double[] gradient){
        
        assert(local3DGradient.length == 3); // make sure we are working with a cartesian gradient
        assert(gradientValues == null); // need to ensure this in order to avoid information loss!
        
        final int arrayLength = local3DGradient[0].length;
        
        System.arraycopy(local3DGradient[0], 0, gradient, 0, arrayLength);
        System.arraycopy(local3DGradient[1], 0, gradient, arrayLength, arrayLength);
        System.arraycopy(local3DGradient[2], 0, gradient, arrayLength * 2, arrayLength);
    }
    
    public void markProblem(){
        this.functionValue = FixedValues.NONCONVERGEDENERGY;
        for(int i = 0; i < local3DGradient.length; i++){
            for(int j = 0; j < local3DGradient[0].length; j++){
                local3DGradient[i][j] = FixedValues.NONCONVERGEDGRADIENT;
            }
        }
    }
    
    public void zeroGradient(){
        if(local3DGradient == null){System.err.println("WARNING: Attempting zeroing of non-existing gradient storage!"); return;}
        for(int i = 0; i < local3DGradient.length; i++){
            for(int j = 0; j < local3DGradient[0].length; j++){
                local3DGradient[i][j] = 0.0;
            }
        }
    }

    void newZeroGradient(final int noOfAtoms) {
        
        functionValue = 0.0;
        local3DGradient = new double[3][noOfAtoms];
    }
    
    public void copyDataIn(final Gradient otherGrad){
        
        this.functionValue = otherGrad.functionValue;
        
        if(local3DGradient == null){
            local3DGradient = new double[otherGrad.local3DGradient.length][otherGrad.local3DGradient[0].length];
        }
        
        if(local3DGradient.length != otherGrad.local3DGradient.length){
            throw new RuntimeException("Gradient dimensions mismatch (first dimension) " + local3DGradient.length + " vs "
                + otherGrad.local3DGradient.length);
        }
        
        for(int i = 0; i < local3DGradient.length; i++){
            // check if the sizes match
            if(otherGrad.local3DGradient[i].length != local3DGradient[i].length){
                throw new RuntimeException("Gradient dimensions mismatch (second dimension) " + local3DGradient[i].length + " vs "
                + otherGrad.local3DGradient[i].length);
            }
            System.arraycopy(otherGrad.local3DGradient[i], 0, local3DGradient[i], 0, local3DGradient[i].length);
        }
    }
}
