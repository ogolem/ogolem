/**
Copyright (c) 2011-2014, J. M. Dieterich
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
 * A lookup table for the any function f(x).
 * @author Johannes Dieterich
 * @version 2014-12-21
 */
public class GenericLookup extends AbstractLookup {
       
    private static final long serialVersionUID = (long) 20141221;
    private final LookupedFunction func;
    
    public GenericLookup(final LookupedFunction func, final int entries, final double start, final double end){
        super(entries, start, end, true);
        assert(func != null);
        this.func = func;
    }
    
    public GenericLookup(final GenericLookup orig){
        super(orig);
        assert(orig.func != null);
        this.func = orig.func.clone();
    }
    
    @Override
    public GenericLookup clone(){
        return new GenericLookup(this);
    }
    
    /**
     * A lookup function for the function including linear interpolation.
     * @param x Must be in the interval [end,start]
     * @return The looked up and interpolated func(x). Canonical func(x) if x is outside the interval.
     */
    public double inter(final double x){
        return funcInter(x);
    }
    
    /**
     * A lookup function for the function without interpolation.
     * @param x Must be in the interval [start,end]
     * @return The looked up func(x). Canonical func(x) if x is outside the interval.
     */
    public double nonInter(final double x){
        return funcNonInter(x);
    }
    
    @Override
    protected double func(final double x){
        return func.func(x);
    }
}
