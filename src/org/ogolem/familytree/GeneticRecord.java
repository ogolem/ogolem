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
package org.ogolem.familytree;

import java.io.Serializable;

/**
 * A single genetic record. Only for package internal use!
 * @author Johannes Dieterich
 * @version 2015-07-10
 */
class GeneticRecord implements Serializable {
    
    private static final long serialVersionUID = (long) 20150710;
    
    private final int id;
    private final int father;
    private final int mother;
    private final boolean wasNull;
    private final boolean wasAccepted;
    private boolean inFinalPool = false;
    private int noChildren = 0;
    
    GeneticRecord(final int id, final int mother, final int father, final boolean wasAccepted, final boolean wasNull){
        this.id = id;
        this.mother = mother;
        this.father = father;
        this.wasAccepted = wasAccepted;
        this.wasNull = wasNull;
    }
    
    void incrementChildren(){
        noChildren++;
    }
    
    int id() {
        return id;
    }

    int father() {
        return father;
    }

    int mother() {
        return mother;
    }

    boolean wasNull() {
        return wasNull;
    }

    boolean wasAccepted() {
        return wasAccepted;
    }

    int noChildren() {
        return noChildren;
    }
    
    boolean isInFinalPool() {
        return inFinalPool;
    }

    void setInFinalPool(final boolean inFinalPool) {
        this.inFinalPool = inFinalPool;
    }
}