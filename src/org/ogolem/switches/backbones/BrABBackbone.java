/**
Copyright (c) 2010, J. M. Dieterich and B. Hartke
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
package org.ogolem.switches.backbones;

import java.io.Serializable;
import java.util.ArrayList;
import org.ogolem.core.CartesianCoordinates;

/**
 * The bridged azobenzene backbone.
 * @author Johannes Dieterich
 * @version 2012-06-24
 */
public final class BrABBackbone implements Serializable{
//TODO remains to be seen how to implement this in detail....

    private static final long serialVersionUID = (long) 20100622;

    private final ArrayList<CartesianCoordinates> alCartes;

    private final int[][] iaConnects;

    private final ArrayList<boolean[][]> alBonds;

    BrABBackbone(final CartesianCoordinates cartesiansCis, final CartesianCoordinates cartesiansTrans,
            final int[] connectsCis, final int[] connectsTrans, final boolean[][] bondsCis, final boolean[][] bondsTrans){
        this.alCartes = new ArrayList<>(2);
        // put cis and trans form in
        alCartes.add(cartesiansCis);
        alCartes.add(cartesiansTrans);

        this.iaConnects = new int[2][];
        iaConnects[0] = connectsCis;
        iaConnects[1] = connectsTrans;

        this.alBonds = new ArrayList<>(2);
        // put cis and trans form in
        alBonds.add(bondsCis);
        alBonds.add(bondsTrans);
    }

    BrABBackbone(final BrABBackbone refBackbone){
        this.alCartes = new ArrayList<>(2);
        // put cis and trans in
        alCartes.add(refBackbone.getCartesCopy(true));
        alCartes.add(refBackbone.getCartesCopy(false));

        this.iaConnects = new int[2][];
        iaConnects[0] = refBackbone.getConnectsCopy(true);
        iaConnects[1] = refBackbone.getConnectsCopy(false);

        this.alBonds = new ArrayList<>(2);
        // put cis and trans in
        alBonds.add(refBackbone.getBondsCopy(true));
        alBonds.add(refBackbone.getBondsCopy(false));
    }

    CartesianCoordinates getCartesCopy(final boolean bCis){
        if(bCis){
            // return the first
            return new CartesianCoordinates(alCartes.get(0));
        } else{
            // return the second
            return new CartesianCoordinates(alCartes.get(1));
        }
    }

    int[] getConnectsCopy(final boolean bCis){
        if(bCis){
            return iaConnects[0].clone();
        } else{
            return iaConnects[1].clone();
        }
    }

    boolean[][] getBondsCopy(final boolean bCis){
        if(bCis){
            // return the first
            return alBonds.get(0).clone();
        } else{
            // return the second
            return alBonds.get(1).clone();
        }
    }
}
