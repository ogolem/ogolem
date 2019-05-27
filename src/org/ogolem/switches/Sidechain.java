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
package org.ogolem.switches;

import org.ogolem.core.ZMatrix;
import org.ogolem.switches.sidechains.*;

/**
 * A decorator to all the sidechains packed together in the sidechains subpackage.
 * @author Johannes Dieterich
 * @version 2010-01-21
 */
final class Sidechain implements SidechainInter{

    private static final long serialVersionUID = (long) 20100121;

    private final SidechainInter sidechain;

    Sidechain(final Color color){

        final int iColor = color.getThisColor();
        switch(iColor){
            case 0:
                this.sidechain = new HydrogenSide();
                break;
            case 1:
                this.sidechain = new HydroxylSide();
                break;
            case 2:
                this.sidechain = new MethylSide();
                break;
            case 3:
                this.sidechain = new ThiolSide();
                break;
            case 4:
                this.sidechain = new AminoSide();
                break;
            case 5:
                this.sidechain = new AldehydeSide();
                break;
            case 6:
                this.sidechain = new ChlorideSide();
                break;
            case 7:
                this.sidechain = new CarboxylSide();
                break;
            case 8:
                this.sidechain = new NitroSide();
                break;
            default:
                System.err.println("WARNING: I do not understand the color " + iColor + ". Using hydrogen now.");
                this.sidechain = new HydrogenSide();
        }
    }

    @Override
    public ZMatrix returnZMatrixCopy(){
        return sidechain.returnZMatrixCopy();
    }

    @Override
    public boolean[][] returnBondingCopy(){
        return sidechain.returnBondingCopy();
    }
}
