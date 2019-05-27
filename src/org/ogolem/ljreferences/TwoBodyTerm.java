/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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
package org.ogolem.ljreferences;

import java.text.DecimalFormat;

/**
 *
 * @author Johannes Dieterich
 * @version 2014-09-06
 */
final class TwoBodyTerm {

    private static final boolean bPrintGeometry = true;

    private final org.ogolem.adaptive.AdaptiveLJFF ljff;

    TwoBodyTerm(final org.ogolem.adaptive.AdaptiveParameters params){
        ljff = new org.ogolem.adaptive.AdaptiveLJFF(false, false, 6, 16, 2, params, false,1.2,20.0,false,false);
    }

    double evaluateTwoBodyTerms(final String sAtom1, final String sAtom2,
            final String sAtom3, final double dDist12, final double dDist13,
            final double dAngle123){

        final double dRealAngle = dAngle123 * Math.PI / 180.0;

        org.ogolem.core.ZMatrix zmat = new org.ogolem.core.ZMatrix(3);
        final String[] saAts = {sAtom1, sAtom2, sAtom3};
        zmat.setAllAtomNames(saAts);

        final double[] daDists = new double[3];
        daDists[1] = dDist12 * org.ogolem.core.Constants.ANGTOBOHR;
        daDists[2] = dDist13 * org.ogolem.core.Constants.ANGTOBOHR;
        zmat.setAllBondLengths(daDists);

        final int[] iaBondConnects = {0, 0, 0};
        zmat.setAllBondConnects(iaBondConnects);

        final double[] daAngles = {0.0, 0.0, dRealAngle};
        zmat.setAllBondAngles(daAngles);

        final int[] iaAngleConnects = {0, 0, 1};
        zmat.setAllAnglesConnects(iaAngleConnects);

        final double[][] xxx = null;
        final org.ogolem.core.CartesianCoordinates cartes = zmat.translateToCartesianAndAlign(xxx);

        final double[] daXYZ1D = cartes.getAll1DCartes();

        // evaluate the two body terms
        final double dEnergyOf2BodyTerms = ljff.energyCalculation((long) -1000,
                -1, daXYZ1D, cartes.getAllAtomTypes(), cartes.getAllAtomNumbers(), cartes.getAllAtomsPerMol(), new double[3], 3, cartes.getAllCharges(),
                cartes.getAllSpins(), null);

        if(bPrintGeometry){

            final String[] saGeom = cartes.createPrintableCartesians();

            final DecimalFormat distFormat = new DecimalFormat("0.0");

            final String sFileName = "refgeoms" + System.getProperty("file.separator") + sAtom1 + sAtom2 + sAtom3 + "_"
                    + distFormat.format(dDist12) + "_" + distFormat.format(dDist13)
                    + "_" + distFormat.format(dAngle123) + ".xyz";
            try{
                LittleHelpers.writeNow(sFileName, saGeom);
            } catch(Exception e){
                System.err.println("ERROR: Couldn't write geometry.");
                e.printStackTrace(System.err);
            }
        }

        return dEnergyOf2BodyTerms;
    }

}
