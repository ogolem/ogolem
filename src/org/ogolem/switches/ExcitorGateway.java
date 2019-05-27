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
 * A decorator for the excitor classes.
 * @author Johannes Dieterich
 * @version 2010-06-22
 */
final class ExcitorGateway implements Excitor{

    private final Excitor excitor;

    ExcitorGateway(final SwitchesConfig swConfig, final int iWhichExcitor){

        switch(iWhichExcitor){
            case 0:
                excitor = new GugaCI(swConfig, 0, false);
                break;
            case 1:
                excitor = new GugaCI(swConfig, 1, false);
                break;
            case 2:
                excitor = new GugaCI(swConfig, 2, false);
                break;
            case 3:
                excitor = new GugaCI(swConfig, 3, false);
                break;
            case 4:
                excitor = new GugaCI(swConfig, 4, false);
                break;
            case 5:
                excitor = new GugaCI(swConfig, 5, false);
                break;
            case 6:
                excitor = new GugaCI(swConfig, 6, false);
                break;
            case 7:
                excitor = new GugaCI(swConfig, 7, false);
                break;
            case 8:
                excitor = new GugaCI(swConfig, 8, false);
                break;
            case 9:
                excitor = new GugaCI(swConfig, 9, false);
                break;
            case 10:
                excitor = new GugaCI(swConfig, 10, false);
                break;
            case 11:
                excitor = new GugaCI(swConfig, 11, false);
                break;
            case 100:
                excitor = new MopacExcitor(swConfig, "AM1", false, false);
                break;
            case 101:
                excitor = new MopacExcitor(swConfig, "PM3", false, false);
                break;
            case 102:
                excitor = new MopacExcitor(swConfig, "PM5", false, false);
                break;
            case 103:
                excitor = new MopacExcitor(swConfig, "PM6", false, false);
                break;
            case 104:
                excitor = new MopacExcitor(swConfig, "AM1", true, true);
                break;
            case 105:
                excitor = new MopacExcitor(swConfig, "AM1", true, false);
                break;
            case 106:
                excitor = new MopacExcitor(swConfig, "PM3", true, false);
                break;
            case 107:
                excitor = new MopacExcitor(swConfig, "PM5", true, false);
                break;
            case 108:
                excitor = new MopacExcitor(swConfig, "PM6", true, false);
                break;
            case 150:
                excitor = new OrcaExcitor(swConfig, "BP", "def2-TZVP");
                break;
            case 151:
                excitor = new OrcaExcitor(swConfig, "BLYP", "def2-TZVP");
                break;
            default:
                System.err.println("WARNING: No excitor representing " + iWhichExcitor +
                        " available. Using GugaCI.");
                excitor = new GugaCI(swConfig, 0, false);
        }
    }

    @Override
    public double s0s1TransitionEnergy(final CartesianCoordinates cartes, final int iID){
        return excitor.s0s1TransitionEnergy(cartes, iID);
    }
}
