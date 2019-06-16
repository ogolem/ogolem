/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2012, J. M. Dieterich
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
package org.ogolem.switches;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Lists all taboos and can be dynamically extended and, of course, can be asked
 * whether a certain switch is already known and therefore a taboo.
 * @author Johannes Dieterich
 * @version 2012-12-01
 */
public final class Taboos implements Serializable{

    private static final long serialVersionUID = (long) 20121201;

    private final List<String> knownColors;
    
    private static final Taboos nonos = new Taboos();

    /**
     * Constructs the taboo list with big enough of a list.
     */
    private Taboos(){
        this.knownColors = Collections.synchronizedList(new ArrayList<>(2*SwitchesConfig.iNoOfGlobIters
                + SwitchesConfig.iPoolSize));
    }

    /**
     * Since it is a singleton, this hands over a reference to it.
     * @return a reference to the taboos
     */
    public static Taboos getReference(){
        return nonos;
    }

    /**
     * Tells wether a taboo exists for this switch.
     * @param sw
     * @return true if this switch is known
     */
    boolean isThisKnown(final Switch sw){
        // ask the switch what its colorcode is
        final String sColorCode = sw.myColorCode();

        return knownColors.contains(sColorCode);
    }

    /**
     * Adds another taboo to the list.
     * @param sw
     */
    public void addTaboo(final Switch sw){
        final String sColorCode = sw.myColorCode();
        knownColors.add(sColorCode);
    }
}
