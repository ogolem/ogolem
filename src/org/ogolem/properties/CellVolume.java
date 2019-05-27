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
package org.ogolem.properties;

/**
 * An (equilibrium) cell volume property.
 * @author Johannes Dieterich
 * @version 2015-09-21
 */
public class CellVolume implements Property {
    
    private static final long serialVersionUID = (long) 20150921;

    private double cellVolume;
    
    public CellVolume(final double cellVolume){
        assert(!Double.isInfinite(cellVolume));
        assert(!Double.isNaN(cellVolume));
        assert(cellVolume >= 0.0);
        this.cellVolume = cellVolume;
    }
    
    private CellVolume(final CellVolume orig){
        this.cellVolume = orig.cellVolume;
    }
    
    @Override
    public CellVolume clone() {
        return new CellVolume(this);
    }

    @Override
    public double getValue() {
        return cellVolume;
    }

    @Override
    public double signedDifference(final Property p) {
        if(!(p instanceof CellVolume)) {throw new IllegalArgumentException("Property should be an instance of CellVolume!");}
        return (cellVolume - p.getValue());
    }

    @Override
    public double absoluteDifference(final Property p) {
        if(!(p instanceof CellVolume)) {throw new IllegalArgumentException("Property should be an instance of CellVolume!");}
        return Math.abs(cellVolume - p.getValue());
    }

    @Override
    public boolean makeSensible() {
        if(Double.isInfinite(cellVolume) || Double.isNaN(cellVolume) || cellVolume < 0.0){
            cellVolume = 0.0; // not too sane...
            return true;
        }
        
        return false;
    }

    @Override
    public String printableProperty() {
        return "" + cellVolume;
    }

    @Override
    public String name() {
        return "CELL VOLUME";
    }
}
