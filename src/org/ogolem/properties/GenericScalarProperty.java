/**
Copyright (c) 2017-2018, J. M. Dieterich and B. Hartke
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
 * A generic scalar property, i.e., a bare bones implementation.
 * @author Johannes Dieterich
 * @version 2018-01-02
 */
public class GenericScalarProperty extends ScalarProperty {

    private static final long serialVersionUID = (long) 20171215;
    
    private final long id;
    
    public GenericScalarProperty(final double data, final long id){
        super(data);
        this.id = id;
    }
    
    public GenericScalarProperty(final GenericScalarProperty orig){
        super(orig.scalar);
        this.id = orig.id;
    }
    
    @Override
    public GenericScalarProperty clone() {
        return new GenericScalarProperty(this);
    }

    @Override
    protected boolean ensureCorrectProperty(Property p) {
        if(!(p instanceof GenericScalarProperty)) {return false;}
        
        final GenericScalarProperty gp = (GenericScalarProperty) p;
        return (gp.id == id);
    }

    @Override
    public boolean makeSensible() {
        if(Double.isInfinite(this.getValue()) || Double.isNaN(this.getValue())){
            this.scalar = 0.0;
            return true;
        }
        return false;
    }

    @Override
    public String printableProperty() {
        return "" + this.getValue();
    }

    @Override
    public String name() {
        return "GENERICSCALAR" + id;
    }
    
}
