/**
Copyright (c) 2013, J. M. Dieterich
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

import java.util.ArrayList;
import java.util.List;

/**
 * An abstract stub for the bond info.
 * @author Johannes Dieterich
 * @version 2013-10-03
 */
abstract class AbstractBondInfo implements BondInfo {
    
    private static final long serialVersionUID = (long) 20131003;
    protected final int noAtoms;
    
    AbstractBondInfo(final int noAtoms){
        this.noAtoms = noAtoms;
    }
    
    AbstractBondInfo(final AbstractBondInfo orig){
        this.noAtoms = orig.noAtoms;
    }
    
    @Override
    public abstract BondInfo clone();
    
    @Override
    public List<String> translateToInput(){
        
        final List<String> li = new ArrayList<>();
        for(int at1 = 0; at1 < noAtoms-1; at1++){
            for(int at2 = at1+1; at2 < noAtoms; at2++){
                final short b = bondType(at1, at2);
                if(b == BondInfo.NOBOND){continue;}
                
                // print out...
                li.add(at1 + "\t" + at2 + "\t" + b);
            }
        }
        return li;
    }
}
