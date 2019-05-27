/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
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
package org.ogolem.switches.sidechains;

import org.ogolem.core.ZMatrix;
import org.ogolem.switches.SidechainInter;

/**
 * Represents single a hydrogen atom.
 * @author Johannes Dieterich
 * @version 2010-01-21
 */
public final class HydrogenSide implements SidechainInter{

    private static final long serialVersionUID = (long) 20091218;

    private final ZMatrix zmat;

    private final boolean[][] baBonds;

    public HydrogenSide(){
        zmat = new ZMatrix(2);

        final String[] saAtoms = {
            "XX",
            "H"
        };
        zmat.setAllAtomNames(saAtoms);

        /*
         * we do not need dihedral connects, angle connects and bond connects since the standard
         * zeros are just fine
         */
        final double[] daBondLengths = {
            0.0,
            2.101
        };
        zmat.setAllBondLengths(daBondLengths);

        final int[] iaConnects = {
            0,
            0
        };
        zmat.setAllBondConnects(iaConnects);

        /*
         * we (of course) also do not need any kind of angle and dihedral values
         */

        // now the bonds
        baBonds = new boolean[2][2];
        baBonds[0][0] = true;
        baBonds[0][1] = true;
        baBonds[1][0] = true;
        baBonds[1][1] = true;
    }

    @Override
    public ZMatrix returnZMatrixCopy(){
        return new ZMatrix(zmat);
    }

    @Override
    public boolean[][] returnBondingCopy(){
        boolean[][] baCopy = new boolean[baBonds.length][];

        for(int i = 0; i < baBonds.length; i++){
            baCopy[i] = baBonds[i].clone();
        }

        return baCopy;
    }
}
