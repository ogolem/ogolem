/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2012, J. M. Dieterich
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
package org.ogolem.generic.generichistory;

import java.io.Serializable;
import org.ogolem.generic.Optimizable;

/**
 * This is a single genetic record wrapping all needed information.
 * @author Johannes Dieterich
 * @version 2015-07-15
 */
class GeneticRecord<E,T extends Optimizable<E>> implements Serializable {
    
    //TODO this is just a stub
    // fill it with life aka the grid

    private static final long serialVersionUID = (long) 20150715;

    private final long id;
    private final long mother;
    private final long father;
    private final boolean wasAccepted;
    private final boolean wasNull;

    public GeneticRecord(final long id, final long mother, final long father,
            final boolean wasAccepted, final boolean wasNull){
        this.id = id;
        this.father = father;
        this.mother = mother;
        this.wasAccepted = wasAccepted;
        this.wasNull = wasNull;
    }

    /*
     * Getters and setters
     */
    public long getID(){
        return id;
    }

    
    String formattedHistoryLine(){
        return ("\t" + id + "\t" + mother + "\t" + father + "\t" + wasAccepted + "\t\t\t" + wasNull);
    }
}
