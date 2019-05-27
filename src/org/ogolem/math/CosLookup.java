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
package org.ogolem.math;

/**
 * A lookup table for the cos function.
 * @author Johannes Dieterich
 * @version 2012-06-12
 */
public class CosLookup extends AbstractLookup {
    
    private static final long serialVersionUID = (long) 20120612;
    
    private static final double twoPi = 2.0*Math.PI;
    
    public CosLookup(final int entries){
        super(entries, 0.0, Math.PI);
    }
    
    public CosLookup(final CosLookup orig){
        super(orig);
    }
    
    @Override
    public CosLookup clone(){
        return new CosLookup(this);
    }
    
    /**
     * A lookup function for the cos including linear interpolation.
     * @param x The angle.
     * @return The looked up and interpolated cos(x).
     */
    public double cosInter(final double x){
        
        double y = x % (twoPi);
        if(y < 0.0) y = -y;
        
        final boolean b = y < Math.PI;      
        final double z = (b) ? y : y - Math.PI;
        final double poi = z * disEntries;
        final int point = (int) poi;
        // jump in
        final double uncorr = (b) ? table[point] : -table[point];
        final double uncorr2 = (point == ENTRIESDECR) ? uncorr : table[point+1];
        final double uncorrT = (b) ? uncorr2 : -uncorr2;
        final double rest = poi-point;
        
        return (uncorr + rest*(uncorrT-uncorr));
    }
    
    /**
     * A lookup function for the cos without interpolation.
     * @param x The angle.
     * @return The looked up cos(x).
     */
    public double cosNonInter(final double x){
        
        double y = x % (twoPi);
        if(y < 0.0) y = -y;
        
        final boolean b = y < Math.PI;
        final double z = (b) ? y : y - Math.PI;
        final double poi = z * disEntries;
        final int point = (int) poi;
        // jump in
        final double uncorr = (b) ? table[point] : -table[point];
        return uncorr;
    }
    
    @Override
    protected double func(final double x){
        return Math.cos(x);
    }
}
