/**
Copyright (c) 2013, J. M. Dieterich
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
package org.ogolem.core;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;

/**
 * Chains multiple niche computations together.
 * @author Johannes Dieterich
 * @version 2015-05-26
 */
class ChainedNicheComp implements NicheComputer<Molecule,Geometry> {
    
    private static final long serialVersionUID = (long) 20140107;
    private final List<NicheComputer<Molecule,Geometry>> nicheComps;
    
    ChainedNicheComp(final List<NicheComputer<Molecule,Geometry>> nicheComps){
        this.nicheComps = nicheComps;
    }
    
    private ChainedNicheComp(final ChainedNicheComp orig){
        this.nicheComps = new ArrayList<>();
        for(final NicheComputer<Molecule,Geometry> nicher : orig.nicheComps){
            nicheComps.add(nicher.clone());
        }
    }
    
    @Override
    public ChainedNicheComp clone(){
        return new ChainedNicheComp(this);
    }
    
    @Override
    public Niche computeNiche(final Geometry g){
        
        String nicheString = "";
        for(final NicheComputer<Molecule,Geometry> comp : nicheComps){
            final Niche n = comp.computeNiche(g);
            nicheString += n.getID();
        }
        
        return new Niche(nicheString);
    }
}
