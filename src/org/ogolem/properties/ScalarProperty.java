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

/**
 * The base class for a scalar property
 * @author Johannes Dieterich
 * @version 2017-12-15
 */
public abstract class ScalarProperty implements Property {
    
    private static final long serialVersionUID = (long) 20171215;
    
    protected double scalar;
    
    protected ScalarProperty (final double scalar){
        assert(!Double.isInfinite(scalar));
        assert(!Double.isNaN(scalar));
        this.scalar = scalar;
    }
    
    @Override
    public abstract ScalarProperty clone();
    
    @Override
    public double getValue(){
        return scalar;
    }
    
    @Override
    public double signedDifference(final Property p){
        
        if(!ensureCorrectProperty(p)){throw new IllegalArgumentException("Property should be an instance of " + name());}
        return (this.getValue() - p.getValue());
    }
    
    @Override
    public double absoluteDifference(final Property p){
        if(!ensureCorrectProperty(p)){throw new IllegalArgumentException("Property should be an instance of " + name());}
        return Math.abs(this.getValue() - p.getValue());
    }
    
    protected abstract boolean ensureCorrectProperty(final Property p);
}
