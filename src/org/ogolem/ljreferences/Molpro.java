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
final class Molpro {

    static double readMolproOutputForTriple(final String sFirstAtom, final String sSecAtom,
            final String sThirdAtom, final float fAngle123){
        
        final DecimalFormat formAngles = new DecimalFormat("0");
        final String sCurrAngle = formAngles.format(fAngle123);
        final String sOutputFile = sFirstAtom + "_" + sSecAtom + "_" + sThirdAtom
                + "_" + sCurrAngle + ".out";

        // read it in
        String[] saCont;
        try{
            saCont = LittleHelpers.readFileIn(sOutputFile);
        } catch(Exception e){
            e.printStackTrace();
            return FixedValues.NONCONVERGEDENERGY;
        }

        int iMemLine = -1;
        for (int i = 0; i < saCont.length; i++) {
            String s = saCont[i].trim();
            if (s.startsWith("Variable memory released")) {
                iMemLine = i;
            }
        }

        if (iMemLine == -1) {
            System.err.println("#WARNING: in finding variable memory released tag.");
            return FixedValues.NONCONVERGEDENERGY;
        } else {
            String sEnergyLine = saCont[iMemLine - 2].trim();
            String sEnergy = sEnergyLine.substring(0,sEnergyLine.indexOf(" ")).trim();
            try {
                double dEnergy = Double.parseDouble(sEnergy);
                return dEnergy;
            } catch (Exception e) {
                e.printStackTrace();
                return FixedValues.NONCONVERGEDENERGY;
            }
        }
    }

    static double readMolproOutputForPair(final String sFirstAtom, final String sSecAtom,
            final float fDist){
        
        final DecimalFormat formDists = new DecimalFormat("0.0");
        final String sCurrDist = formDists.format(fDist);
        final String sOutputFile = sFirstAtom + "_" + sSecAtom + "_" + sCurrDist + ".out";

        // read it in
        String[] saCont;
        try{
            saCont = LittleHelpers.readFileIn(sOutputFile);
        } catch(Exception e){
            e.printStackTrace();
            return FixedValues.NONCONVERGEDENERGY;
        }

        int iMemLine = -1;
        for (int i = 0; i < saCont.length; i++) {
            String s = saCont[i].trim();
            if (s.startsWith("Variable memory released")) {
                iMemLine = i;
            }
        }

        if (iMemLine == -1) {
            System.err.println("#WARNING: in finding variable memory released tag.");
            return FixedValues.NONCONVERGEDENERGY;
        } else {
            String sEnergyLine = saCont[iMemLine - 2].trim();
            String sEnergy = sEnergyLine.substring(0,sEnergyLine.indexOf(" ")).trim();
            try {
                double dEnergy = Double.parseDouble(sEnergy);
                return dEnergy;
            } catch (Exception e) {
                e.printStackTrace();
                return FixedValues.NONCONVERGEDENERGY;
            }
        }
    }

    static double[] readMolproOutputForSingles(final String[] saAtoms){

        final double[] daRefs = new double[saAtoms.length];

        for(int iAtom = 0; iAtom < saAtoms.length; iAtom++){
            final String sAt = saAtoms[iAtom].toLowerCase();

            String[] saCont;

            try{
                saCont = LittleHelpers.readFileIn(sAt + ".out");
            } catch(Exception e){
                e.printStackTrace();
                continue;
            }

            int iMemLine = -1;
            for(int i = 0; i < saCont.length; i++){
                String s = saCont[i].trim();
                if(s.startsWith("Variable memory released")){
                    iMemLine = i;
                }
            }

            if(iMemLine == -1){
                System.err.println("All of this is pointless!!!! Error in finding variable memory released tag.");
            } else{
                String sEnergyLine = saCont[iMemLine - 2].trim();
                String sEnergy = sEnergyLine.substring(0,sEnergyLine.indexOf(" ")).trim();
                try{
                    double dEnergy = Double.parseDouble(sEnergy);
                    daRefs[iAtom] = dEnergy;
                } catch(Exception e){
                    e.printStackTrace();
                }
            }
        }

        return daRefs;
    }
}
