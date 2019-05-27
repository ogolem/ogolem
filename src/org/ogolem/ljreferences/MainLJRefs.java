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
package org.ogolem.ljreferences;

import java.text.DecimalFormat;
import org.ogolem.core.FixedValues;

/**
 *
 * @author Johannes Dieterich
 * @version 2010-08-22
 */
public class MainLJRefs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        final boolean bWriteOrca = false;
        final boolean bReadMolpro = true;
        final boolean bWriteData = true;
        final boolean bDoDoubles = false;
        final boolean bDoTriples = true;
        final boolean bWriteOGOLEMInput = true;

        // for the doubles
        final float fDistStart = 1.0f;
        final float fDistEnd = 10.0f;
        final float fDistIncr = 0.2f;

        // for the triples
        final float fAngleStart = 20.0f;
        final float fAngleEnd = 180.0f;
        final float fAngleIncr = 20.0f;
        final float fDist12 = 3.5f;
        final float fDist13 = 4.0f;

        // for the 2BodyEvaluation
        final boolean bEvaluate2BodyTerms = true;
        final String sParameterFile = "adaptive-ljff.param";
        org.ogolem.adaptive.AdaptiveParameters params = null;

        final String[] saRareGases = {"He", "Ne", "Ar", "Kr", "Xe"};

        double[] daRefs = null;
        final DecimalFormat energyFormat = new DecimalFormat("0.00000000");
        final DecimalFormat distFormat = new DecimalFormat("0.0");

        if(bWriteOrca){
            Orca.writeAllOrcaInput();
        } else if(bReadMolpro){

            // get the references
            daRefs = Molpro.readMolproOutputForSingles(saRareGases);

            if (bDoDoubles) {
                // screen all pairs from 1.0 to 10.0 Angstrom in 0.2 Angstrom steps
                for (int i = 0; i < saRareGases.length; i++) {
                    final String sFirstGas = saRareGases[i];

                    for (int j = i; j < saRareGases.length; j++) {

                        final String sSecGas = saRareGases[j];

                        if (bWriteData) {
                            System.out.println("# distance     energy        normalized energy (and somehow rounded)");
                            System.out.println("# references: E(" + sFirstGas + ") = " + energyFormat.format(daRefs[i])
                                    + "    E(" + sSecGas + ") = " + energyFormat.format(daRefs[j]));
                        }

                        // now over the distances
                        float fCurrDist = fDistStart;
                        while (fCurrDist <= fDistEnd) {

                            final double dEnergy = Molpro.readMolproOutputForPair(sFirstGas, sSecGas, fCurrDist);

                            if (dEnergy >= FixedValues.NONCONVERGEDENERGY) {
                                // increment distance
                                fCurrDist += fDistIncr;
                                continue;
                            }

                            final double dNormEnergy = dEnergy - daRefs[i] - daRefs[j];

                            if (bWriteData) {
                                System.out.println(" " + distFormat.format(fCurrDist) + "   " + energyFormat.format(dEnergy)
                                        + "   " + energyFormat.format(dNormEnergy));
                            }

                            // increment distance
                            fCurrDist += fDistIncr;
                        }
                    }
                }
            }

            if (bDoTriples) {

                if(bEvaluate2BodyTerms){
                    // read parameters in
                    try{
                        final String[] saParamData = LittleHelpers.readFileIn(sParameterFile);
                        params = new org.ogolem.adaptive.AdaptiveParameters(saParamData, -1);
                    } catch(Exception e){
                        System.err.println("ERROR: Can't read parameters in!");
                        e.printStackTrace();
                    }
                }

                for(int i = 0; i < saRareGases.length; i++){
                    final String sAtom1 = saRareGases[i];
                    
                    for(int j = 0; j < saRareGases.length; j++){
                        final String sAtom2 = saRareGases[j];

                        for(int k = 0; k < saRareGases.length; k++){
                            final String sAtom3 = saRareGases[k];

                            if (bWriteData) {
                                System.out.println("# angle     energy        normalized energy (and somehow rounded)");
                                System.out.println("# references: E(" + sAtom1 + ") = " + energyFormat.format(daRefs[i])
                                    + "    E(" + sAtom2 + ") = " + energyFormat.format(daRefs[j])
                                    + "    E(" + sAtom3 + ") = " + energyFormat.format(daRefs[k]));
                            }

                            // now over the angles
                            float fCurrAngle = fAngleStart;
                            while(fCurrAngle <= fAngleEnd){

                                final double dEnergy = Molpro.readMolproOutputForTriple(sAtom1, sAtom2, sAtom3, fCurrAngle);

                                if (dEnergy >= FixedValues.NONCONVERGEDENERGY) {
                                    // increment angle
                                    fCurrAngle += fAngleIncr;
                                    continue;
                                }

                                final double dNormEnergy = dEnergy - daRefs[i] - daRefs[j] - daRefs[k];

                                double d2BodyEnergy = Double.NEGATIVE_INFINITY;
                                if(bEvaluate2BodyTerms){
                                    final TwoBodyTerm twoBodyTerm = new TwoBodyTerm(params);
                                    d2BodyEnergy = twoBodyTerm.evaluateTwoBodyTerms(sAtom1,
                                            sAtom2, sAtom3, (double) fDist12, (double) fDist13, (double) fCurrAngle);
                                }

                                if (bWriteData && !bEvaluate2BodyTerms) {
                                    System.out.println(" " + distFormat.format(fCurrAngle) + "   " + energyFormat.format(dEnergy)
                                        + "   " + energyFormat.format(dNormEnergy));
                                } else if(bWriteData && bEvaluate2BodyTerms && bWriteOGOLEMInput){
                                    System.out.println(" " + distFormat.format(fCurrAngle) + "   " + energyFormat.format(dEnergy)
                                        + "   " + energyFormat.format(dNormEnergy)
                                        + "   " + energyFormat.format(d2BodyEnergy)
                                        + "   " + energyFormat.format(dNormEnergy-d2BodyEnergy));

                                    // write OGOLEM input ;-)
                                    System.err.println("<REFERENCE>");
                                    System.err.println("<PATH>");
                                    System.err.println("refgeoms"+ System.getProperty("file.separator")
                                            + sAtom1 + sAtom2 + sAtom3 + "_" + distFormat.format(fDist12) + "_" + distFormat.format(fDist13)
                                            + "_" + distFormat.format(fCurrAngle) + ".xyz");
                                    System.err.println("</PATH>");
                                    System.err.println("<ENERGY>");
                                    System.err.println(" " + dNormEnergy);
                                    System.err.println("</ENERGY>");
                                    System.err.println("</REFERENCE>");
                                }

                                // increment angle
                                fCurrAngle += fAngleIncr;
                            }
                        }
                    }
                }
            }
        }
    }
}
