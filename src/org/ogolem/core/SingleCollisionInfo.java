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
package org.ogolem.core;

import java.util.ArrayList;
import java.util.List;

/**
 * A collision info only able to store a single collision.
 * @author Johannes Dieterich
 * @version 2015-07-23
 */
class SingleCollisionInfo extends AbstractCollisionInfo {

    private static final long serialVersionUID = (long) 20150720;
    
    private int atom1;
    private int atom2;
    private double strength;
    
    @Override
    public boolean reportCollision(final int atom1, final int atom2, final double strength) {
        
        if(noCollisions > 0){
            System.err.println("Previous collision already stored in SingleCollisionInfo.");
            return false;
        }
        
        this.atom1 = atom1;
        this.atom2 = atom2;
        this.strength = strength;
        noCollisions++;
        
        return true;
    }

    @Override
    public List<Collision> getCollisions() {
        
        if(noCollisions == 0){
            return new ArrayList<>();
        }
        
        final Collision coll = new Collision(atom1,atom2,strength);
        final List<Collision> colls = new ArrayList<>();
        colls.add(coll);
        
        return colls;
    } 

    @Override
    protected void cleanState() {
        atom1 = -1;
        atom2 = -1;
        strength = -1.0;
    }
}
