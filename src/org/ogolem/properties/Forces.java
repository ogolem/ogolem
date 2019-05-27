/**
Copyright (c) 2015 - 2017, J. M. Dieterich and B. Hartke
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
 * A forces property.
 * @author Johannes Dieterich
 * @version 2017-10-13
 */
public class Forces implements Property {
    
    private static final long serialVersionUID = (long) 20150727;

    public static final double DEFAULTFORCE = FixedValues.NONCONVERGEDGRADIENT;
    
    private final double[][] forces;
    
    public Forces(final double[][] forces){
        this.forces = forces;
    }
    
    private Forces(final Forces orig){
        
        if(orig.forces != null){
        
            this.forces = new double[orig.forces.length][];
            for(int i = 0; i < orig.forces.length; i++){
                this.forces[i] = orig.forces[i].clone();
            }
        } else {
            this.forces = null;
        }
    }
    
    @Override
    public Forces clone() {
        return new Forces(this);
    }

    /**
     * Will return the absolute sum of all the gradient components.
     * @return absolute sum of all gradient components.
     */
    @Override
    public double getValue() {
        
        if(forces == null){return FixedValues.NONCONVERGEDENERGY;}
        
        double sum = 0.0;
        for(int i = 0; i < forces.length; i++){
            for(int j = 0; j < forces[i].length; j++){
                sum += Math.abs(forces[i][j]);
            }
        }
        
        return sum;
    }

    /**
     * This will ALSO return absolute differences!
     * @param p the other property. Must be an instance of Forces.
     * @return the ABSOLUTE difference
     */
    @Override
    public double signedDifference(final Property p) {
        return absoluteDifference(p);
    }

    @Override
    public double absoluteDifference(final Property p) {
        if(!(p instanceof Forces)) {throw new IllegalArgumentException("Property should be an instance of Forces!");}
        else{
            final Forces other = (Forces) p;
            
            if(forces == null || other.forces == null){return FixedValues.NONCONVERGEDENERGY;}
            
            final double[][] oForces = other.forces;
            if(oForces.length != forces.length){throw new IllegalArgumentException("First force dimension must be the same! Are " + forces.length + " vs " + oForces.length);}
            double sumDiff = 0.0;
            for(int i = 0; i < forces.length; i++){
                if(oForces[i].length != forces[i].length){throw new IllegalArgumentException("Second force dimension " + i + " must be the same! Are " + forces[i].length + " vs " + oForces[i].length);}
                for(int j = 0; j < forces[i].length; j++){
                    final double diff = Math.abs(forces[i][j] - oForces[i][j]);
                    sumDiff += diff;
                }
            }
            
            // "norm"
            sumDiff /= forces.length;
            
            return sumDiff;
        }
    }

    @Override
    public boolean makeSensible() {
        
        if(forces == null){return false;}
        
        boolean wasTouched = false;
        for(int i = 0; i < forces.length; i++){
            for(int j = 0; j < forces[i].length; j++){
                if(Double.isInfinite(forces[i][j]) || Double.isNaN(forces[i][j]) || forces[i][j] > FixedValues.NONCONVERGEDGRADIENT){
                    forces[i][j] = FixedValues.NONCONVERGEDGRADIENT;
                    wasTouched = true;
                }
            }
        }
        
        return wasTouched;
    }

    @Override
    public String printableProperty() {
        
        if(forces == null){return "NULL'D FORCES";}
        
        String s = "";
        for(int i = 0; i < forces.length; i++){
            s += "Row " + i + ":";
            for(int j = 0; j < forces[i].length; j++){
                s += forces[i][j] + "\t";
            }
            
            s += "\n";
        }
        
        return s;
    }

    @Override
    public String name() {
        return "FORCES";
    }
    
    public static Forces getDefaultForces(final int noAtoms){
        
        final double[][] forceVals = new double[3][noAtoms];
        for(int coord = 0; coord < 3; coord++){
            for(int at = 0; at < noAtoms; at++){
                forceVals[coord][at] = DEFAULTFORCE;
            }
        }
        
        final Forces forces = new Forces(forceVals);
        
        return forces;
    }
}
