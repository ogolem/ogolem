/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive;

import org.ogolem.core.CastException;
import static org.ogolem.core.Constants.EVTOHARTREE;
import org.ogolem.core.InitIOException;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.InquiryPrimitives;

/**
 * More advanced input features.
 * @author Johannes Dieterich
 * @version 2013-01-05
 */
public class Input{

    static boolean doesFileExist(final String sFile){
        return InquiryPrimitives.doesFileExist(sFile);
    }

    static String[] ReadFile(String sFile) throws InitIOException{
        try{
            String[] saFileContent = InputPrimitives.readFileIn(sFile);
            return saFileContent;
        } catch(Exception e){
            throw new InitIOException(e);
        }
    }

    static double ReadEnergyMopacOutput(final String sMopacOutput, final boolean bAzobenzeneMode)
            throws InitIOException, CastException{

        // first read the file in
        String[] saContent = ReadFile(sMopacOutput);

        double dEnergy = FixedValues.NONCONVERGEDENERGY;

        for (int i = saContent.length - 1; i >= 0; i--) {
            if (saContent[i].contains("TOTAL ENERGY")) {
                String sTemp = saContent[i].trim();
                // the 26 is really mopac output specific.
                sTemp = sTemp.substring(26);
                sTemp = sTemp.trim();
                // get rid of the EV in the end
                sTemp = sTemp.substring(0, sTemp.indexOf(" "));
                try {
                    dEnergy = Double.parseDouble(sTemp) * EVTOHARTREE;
                } catch (Exception e) {
                    System.err.println("WARNING: Problems casting the energy of the mopac output.");
                    throw new CastException(e);
                }
            }
        }

        if(dEnergy > FixedValues.NONCONVERGEDENERGY){
            throw new CastException("WARNING: Something went wrong in reading MOPACs energy.");
        }


        return dEnergy;
    }

    static double[][] readBordersIn(final String sBorderPath, final int iNoOfParams)
            throws InitIOException, CastException{

        // first read the file in: will already throw an exception of not existing
        final String[] cont = ReadFile(sBorderPath);

        // then check if it is a valid file
        if(!cont[0].trim().equalsIgnoreCase("###OGOLEMAUX###")){
            throw new InitIOException("Your border input is lacking the correct format.");
        }

        final double[][] borders = new double[2][iNoOfParams];
        for(int i = 1; i < iNoOfParams +1; i++){
            // read min and max in
            final String[] sa = cont[i].trim().split("\\s+");

            // try to cast the two
            try{
                final double min = Double.parseDouble(sa[0].trim());
                final double max = Double.parseDouble(sa[1].trim());
                borders[0][i-1] = min;
                borders[1][i-1] = max;
            } catch(Exception e){
                throw new CastException("Error in casting borders.",e);
            }
        }

        return borders;
    }

    @SuppressWarnings("unchecked")
    public static GenericPool<Double,AdaptiveParameters> readParamPoolIn(final String sPath) throws Exception {
        Object oObj = null;
        try {
            oObj = InputPrimitives.readBinInput(sPath);
        } catch (Exception e) {
            throw e;
        }
        // randomly chosen size of 10 (not too big), just works (TM)
        GenericPool<Double,AdaptiveParameters> pool;
        try {
            pool = (GenericPool<Double,AdaptiveParameters>) oObj;
        } catch (Exception e) {
            throw e;
        }
        return pool;
    }
}
