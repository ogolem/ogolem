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
package org.ogolem.properties;

import org.ogolem.core.FixedValues;

/**
 * A density on a grid (Cartesian, equidistantly spaced) property.
 * @author Johannes Dieterich
 * @version 2017-09-25
 */
public class Density implements Property {
    
    private static final long serialVersionUID = (long) 20170925;

    public static final double DEFAULTDENSITY = 0.0;
    
    private final double[][][] density;
    
    public Density(final double[][][] density){
        this.density = density;
    }
    
    private Density(final Density orig){
        
        if(orig.density != null){        
            this.density = new double[orig.density.length][orig.density[0].length][];
            for(int i = 0; i < orig.density.length; i++){
                for(int j = 0; j < orig.density[0].length; j++){
                    this.density[i][j] = orig.density[i][j].clone();
                }
            }
        } else {
            this.density = null;
        }
    }
    
    @Override
    public Density clone() {
        return new Density(this);
    }

    /**
     * Will return the integral over the density, except for negative densities!
     * @return integral over the density.
     */
    @Override
    public double getValue() {
        
        if(density == null){return FixedValues.NONCONVERGEDENERGY;}
        
        double sum = 0.0;
        for(int i = 0; i < density.length; i++){
            for(int j = 0; j < density[i].length; j++){
                for(int k = 0; k < density[i][j].length; k++){
                    if(density[i][j][k] < 0.0) continue;
                    sum += density[i][j][k];
                }
            }
        }
        
        return sum;
    }

    /**
     * This will ALSO return absolute differences!
     * @param p the other property. Must be an instance of Density.
     * @return the ABSOLUTE difference
     */
    @Override
    public double signedDifference(final Property p) {
        return absoluteDifference(p);
    }

    @Override
    public double absoluteDifference(final Property p) {
        if(!(p instanceof Density)) {throw new IllegalArgumentException("Property should be an instance of Density!");}
        else{                        
            final Density other = (Density) p;
            
            if(this.density == null || other.density == null){
                return FixedValues.NONCONVERGEDENERGY;
            }
            
            final double[][][] oDensity = other.density;
            if(oDensity.length != density.length){throw new IllegalArgumentException("First density dimension must be the same! Are " + density.length + " vs " + oDensity.length);}
            double sumDiff = 0.0;
            for(int i = 0; i < density.length; i++){
                if(oDensity[i].length != density[i].length){throw new IllegalArgumentException("Second density dimension in " + i + " must be the same! Are " + density[i].length + " vs " + oDensity[i].length);}
                for(int j = 0; j < density[i].length; j++){
                    for(int k = 0; k < density[i][j].length; k++){
                        if(oDensity[i][j].length != density[i][j].length){
                            System.out.println("INFO: Full dimensions other: " + oDensity.length + "/" + oDensity[i].length + "/" + oDensity[i][j].length);
                            System.out.println("INFO: Full dimensions this: " + density.length + "/" + density[i].length + "/" + density[i][j].length);
                            throw new IllegalArgumentException("Third density dimension in " + i + "/" + j + " must be the same! Are " + density[i][j].length + " vs " + oDensity[i][j].length);
                        }
                        if(density[i][j][k] < 0.0 || oDensity[i][j][k] < 0.0) continue; // flagged to be ignored.
                        final double diff = Math.abs(density[i][j][k] - oDensity[i][j][k]);
                        sumDiff += diff;
                    }
                }
            }
            
            // "norm"
            sumDiff /= density.length;
            
            return sumDiff;
        }
    }

    @Override
    public boolean makeSensible() {
        
        if(density == null){return false;}
        
        boolean wasTouched = false;
        for(int i = 0; i < density.length; i++){
            for(int j = 0; j < density[i].length; j++){
                for(int k = 0; k < density[i][j].length; k++){
                    if(Double.isInfinite(density[i][j][k]) || Double.isNaN(density[i][j][k]) || density[i][j][k] > FixedValues.NONCONVERGEDGRADIENT){
                        density[i][j][k] = 0.0;
                        wasTouched = true;
                    }
                }
            }
        }
        
        return wasTouched;
    }

    @Override
    public String printableProperty() {
        
        if(density == null){return "NULL'D DENSITY";}
        
        String s = "";
        for(int i = 0; i < density.length; i++){
            for(int j = 0; j < density[i].length; j++){
                for(int k = 0; k < density[i][j].length; k++){
                    s += density[i][j][k] + "\t";
                }
                
                s += "\n";
            }
            
            s += "\n";
        }
        
        return s;
    }

    @Override
    public String name() {
        return "DENSITY";
    }
    
    public static Density getDefaultDensity(final int gridX, final int gridY, final int gridZ){
        
        final double[][][] densityVals = new double[gridX][gridY][gridZ];
        for(int x = 0; x < gridX; x++){
            for(int y = 0; y < gridY; y++){
                for(int z = 0; z < gridZ; z++){
                    densityVals[x][y][z] = 42.0;
                }
            }
        }
        
        final Density density = new Density(densityVals);
        
        return density;
    }
    
    public int getGridDimX(){
        
        if(density == null){return 0;}
        
        assert(density.length > 0);
        return density.length;
    }
    
    public int getGridDimY(){
        
        if(density == null){return 0;}
        
        assert(density[0].length > 0);
        return density[0].length;
    }
    
    public int getGridDimZ(){
        
        if(density == null){return 0;}
        
        assert(density[0][0].length > 0);
        return density[0][0].length;
    }
}
