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
 * Represents a hydroxyl sidechain.
 * @author Johannes Dieterich
 * @version 2010-01-21
 */
public class HydroxylSide implements SidechainInter{

    private static final long serialVersionUID = (long) 20091218;

    private final ZMatrix zmat;

    private final boolean[][] baBonds;

    public HydroxylSide(){
        zmat = new ZMatrix(3);

        final String[] saAtoms = {
            "XX",
            "O",
            "H"
        };
        zmat.setAllAtomNames(saAtoms);

        final int[] iaBondConnects = {
            0,
            0,
            1,
        };
        zmat.setAllBondConnects(iaBondConnects);

        /*
         * we do not need dihedral connects and angle connectssince the standard
         * zeros are just fine
         */
        final double[] daBondLengths = {
            0.0,
            1.7896,
            1.7896
        };
        zmat.setAllBondLengths(daBondLengths);

        final double[] daBondAngles = {
            0.0,
            0.0,
            105.0*Math.PI/180.0
        };
        zmat.setAllBondAngles(daBondAngles);

        /*
         * we (of course) also do not need any kind of dihedral values
         */

        // now the bonding
        baBonds = new boolean[3][3];
        baBonds[0][0] = true;
        baBonds[0][1] = true;
        baBonds[0][2] = false;
        baBonds[1][0] = true;
        baBonds[1][1] = true;
        baBonds[1][2] = true;
        baBonds[2][0] = false;
        baBonds[2][1] = true;
        baBonds[2][2] = true;
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
