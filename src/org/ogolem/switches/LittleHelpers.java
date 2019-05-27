/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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

import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CoordTranslation;

/**
 * A collection of small helping functions not having much in common.
 * @author Johannes Dieterich
 * @version 2013-09-23
 */
final class LittleHelpers {

    /**
     * Tunnels the request through to the CoordinateTranslation of the core package.
     * @param cartes
     * @param dBlowFactor
     * @return
     */
    static BondInfo bondingInfo(final CartesianCoordinates cartes,
            final double dBlowFactor){
        return CoordTranslation.checkForBonds(cartes, dBlowFactor);
    }

    /**
     * Checks whether a set of four given cartesian coordinates has the proper
     * dihedral as specified within a specified threshhold.
     * @param dOptimalValue
     * @param dAllowDiff
     * @param daRespCoords
     * @return
     */
    static boolean checkDihedral(final double dOptimalValue, final double dAllowDiff,
            final double[][] daRespCoords){

        final double[] daCoordsOne = daRespCoords[0];

        final double[] daCoordsTwo = daRespCoords[1];

        final double[] daCoordsThree = daRespCoords[2];

        final double[] daCoordsFour = daRespCoords[3];

        double dDihedral = CoordTranslation.calcDihedral(daCoordsOne,
                daCoordsTwo, daCoordsThree, daCoordsFour);

        if(dDihedral >= 360.0 * Math.PI/180.0){
            dDihedral = dDihedral - Math.PI*2.0;
        }

        if(Math.abs(dDihedral - dOptimalValue) <= dAllowDiff){
            return true;
        } else {
            return false;
        }
    }
}
