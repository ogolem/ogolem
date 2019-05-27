/**
Copyright (c) 2012-2014, J. M. Dieterich
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
package org.ogolem.properties;

import org.ogolem.core.FixedValues;

/**
 *
 * @author Johannes Dieterich
 * @version 2015-03-03
 */
public class Stability implements Property {
    
    private static final long serialVersionUID = (long) 20130104;
    private final boolean stable;
    
    public Stability(final boolean stable){
        this.stable = stable;
    }
    
    @Override
    public Stability clone(){
        return new Stability(stable);
    }
    
    @Override
    public double getValue(){
        return (stable) ? 0.0 : FixedValues.NONCONVERGEDENERGY;
    }
    
    @Override
    public double signedDifference(Property p){
        if(!(p instanceof Stability)) {throw new IllegalArgumentException("Property should be an instance of Stability!");}
        return (getValue() - p.getValue());
    }
    
    @Override
    public double absoluteDifference(Property p){
        if(!(p instanceof Stability)) {throw new IllegalArgumentException("Property should be an instance of Stability!");}
        final Stability sp = (Stability) p;
        if(stable == sp.stable) {return 0.0;}
        else {return FixedValues.NONCONVERGEDENERGY;}
    }
    
    @Override
    public boolean makeSensible(){
        // per definition, either is sensible
        return false;
    }
    
    @Override
    public String printableProperty(){
        return "" + stable;
    }

    @Override
    public String name() {
        return "STABILITY";
    }
}
