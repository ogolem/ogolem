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

import org.ogolem.core.CartesianCoordinates;

/**
 * A decorator to all possible local optimization engines.
 * @author Johannes Dieterich
 * @version 2009-10-21
 */
final class LocalOpt implements LocalOptimization{

    private final LocalOptimization locopt;

    LocalOpt(final SwitchesConfig swConfig, final int iWhichLocOpt){
        switch(iWhichLocOpt){
            case 0:
                locopt = new TinkerLocOpt(0, swConfig.dBlowBondsFac, swConfig.bDebug);
                break;
            case 1:
                locopt = new TinkerLocOpt(1, swConfig.dBlowBondsFac, swConfig.bDebug);
                break;
            case 2:
                locopt = new TinkerLocOpt(2, swConfig.dBlowBondsFac, swConfig.bDebug);
                break;
            case 100:
                locopt = new OpenBabelLocOpt(iWhichLocOpt - 100, swConfig);
                break;
            case 101:
                locopt = new OpenBabelLocOpt(iWhichLocOpt - 100, swConfig);
                break;
            case 102:
                locopt = new OpenBabelLocOpt(iWhichLocOpt - 100, swConfig);
                break;
            case 103:
                locopt = new OpenBabelLocOpt(iWhichLocOpt - 100, swConfig);
                break;
            case 200:
                locopt = new MopacLocOpt(0, swConfig.iMaxIterLocOpt);
                break;
            case 201:
                locopt = new MopacLocOpt(1, swConfig.iMaxIterLocOpt);
                break;
            case 202:
                locopt = new MopacLocOpt(2, swConfig.iMaxIterLocOpt);
                break;
            case 203:
                locopt = new MopacLocOpt(3, swConfig.iMaxIterLocOpt);
                break;
            case 204:
                locopt = new MopacLocOpt(4, swConfig.iMaxIterLocOpt);
                break;
            default:
                System.err.println("WARNING: No choice " + iWhichLocOpt +
                        " available for the local optimization. Using Tinker now.");
                locopt = new TinkerLocOpt(0, swConfig.dBlowBondsFac, swConfig.bDebug);
        }
    }

    @Override
    public CartesianCoordinates locOptThis(final CartesianCoordinates cartesStart, final int id){
        return locopt.locOptThis(cartesStart, id);
    }
}
