/**
Copyright (c) 2011-2012, J. M. Dieterich
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
 * A lookup table for the exp function.
 * @author Johannes Dieterich
 * @version 2012-06-12
 */
public class ExpLookup extends AbstractLookup {
       
    private static final long serialVersionUID = (long) 20120612;
    
    public ExpLookup(final int entries, final double start, final double end){
        super(entries, start, end);
    }
    
    public ExpLookup(final ExpLookup orig){
        super(orig);
    }
    
    @Override
    public ExpLookup clone(){
        return new ExpLookup(this);
    }
    
    /**
     * A lookup function for the exp including linear interpolation.
     * @param x Must be in the interval [end,start]
     * @return The looked up and interpolated exp(x). Canonical exp(x) if x is outside the interval.
     */
    public double expInter(final double x){
        return funcInter(x);
    }
    
    /**
     * A lookup function for the exp without interpolation.
     * @param x Must be in the interval [start,end]
     * @return The looked up exp(x). Canonical exp(x) if x is outside the interval.
     */
    public double expNonInter(final double x){
        return funcNonInter(x);
    }
    
    @Override
    protected double func(final double x){
        return Math.exp(x);
    }
}
