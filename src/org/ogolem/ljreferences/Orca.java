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
package org.ogolem.ljreferences;

import java.util.LinkedList;
import java.util.List;
import java.text.DecimalFormat;

/**
 *
 * @author Johannes Dieterich
 * @version 2012-06-24
 */
final class Orca {

    static void writeAllOrcaInput(){
        final float fDistStart = 1.0f;
        final float fDistEnd = 10.0f;
        final float fDistIncr = 0.2f;

        final String sMethodLine = "!RHF MP2 CCSD(T) Extrapolate(4/5,aug-cc) VeryTightSCF SemiDirect";
        final String sMemLine = "%MaxCore 2000";

        //final String[] saRareGases = {"He", "Ne", "Ar", "Kr", "Xe"};
        // neither aug-cc nor cc nor aug-ano provide Xe basis sets
        final String[] saRareGases = {"He", "Ne", "Ar", "Kr"};

        final DecimalFormat formDists = new DecimalFormat("0.0");

        // screen all pairs from 1.0 to 10.0 Angstrom in 0.2 Angstrom steps
        for(int i = 0; i < saRareGases.length; i++){
            final String sFirstGas = saRareGases[i];

            for(int j = i; j < saRareGases.length; j++){
                final String sSecGas = saRareGases[j];

                // now over the distances
                float fCurrDist = fDistStart;
                while(fCurrDist <= fDistEnd){

                    // the distance as a string
                    final String sCurrDist = formDists.format(fCurrDist);

                    // create input
                    final String sInpFile = sFirstGas + "_" + sSecGas + "_" + sCurrDist + "_autogen.inp";

                    final List<String> lInp = new LinkedList<>();

                    lInp.add(sMethodLine);
                    lInp.add(sMemLine);
                    lInp.add("%basis");
                    lInp.add("  Aux auto");
                    lInp.add("end");
                    lInp.add("*xyz 0 1");

                    lInp.add(sFirstGas + "  0.0 0.0 0.0");
                    lInp.add(sSecGas + " " + sCurrDist + " 0.0 0.0");

                    lInp.add("*");

                    // write it out
                    final String[] saInp = new String[lInp.size()];
                    for(int x = 0; x < saInp.length; x++){
                        saInp[x] = lInp.get(x);
                    }

                    try{
                        LittleHelpers.writeNow(sInpFile, saInp);
                    } catch(Exception e){
                        e.printStackTrace(System.err);
                    }

                    // increment distance
                    fCurrDist += fDistIncr;
                }
            }
        }
    }
}
