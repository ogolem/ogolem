/**
Copyright (c) 2010, J. M. Dieterich and B. Hartke
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

import java.io.File;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.StreamGobbler;

/**
 * Using Orca for TD-DFT calculations.
 * @author Johannes Dieterich
 * @version 2010-07-06
 */
final class OrcaExcitor implements Excitor {

    private static final boolean bDebug = false;

    private final String sWhichMethod;
    
    private final String sBasisSet;

    // hard-coded ATM
    private final int iNoOfRoots = 8;

    // these are not used at the moment
    private static final int iCharge=0;
    private static final int iSpin=1;

    OrcaExcitor(final SwitchesConfig config, final String method,
            final String basis){
        this.sWhichMethod = method;
        this.sBasisSet = basis;
    }

    @Override
    public double s0s1TransitionEnergy(final CartesianCoordinates cartes, final int iID){

        // file stuff
        final String sFolderName = "orca" + iID;
        final String sInputFile = sFolderName + System.getProperty("file.separator")
                + "orca" + iID + ".inp";
        final String sInput = "orca" + iID  +".inp";
        final String sXYZFile = "orca" + iID + ".xyz";

        // set input up and write it
        final String[] saInput = setInputUp(cartes, sXYZFile);
        try{
            Output.createAFolder(sFolderName);
            Output.printMiscToFile(sInputFile, saInput);
        } catch(Exception e){
            System.err.println("ERROR: Couldn't do I/O properly. " + e.toString());
            return FixedValues.UNEXCITABLEENERGY;
        }

        // call orca
        Process proc = null;
        String[] saOutput;
        try{
            final Runtime rt = Runtime.getRuntime();
            final String sCommand = "orca  " + sInput;
            final File dir = new File(sFolderName);

            // inherit the environment
            proc = rt.exec(sCommand,null,dir);

            // any error message?
            final StreamGobbler errorGobbler = new
                StreamGobbler(proc.getErrorStream(), "ERROR");

            // any output?
            final StreamGobbler outputGobbler = new
                StreamGobbler(proc.getInputStream(), "OUTPUT");

            // kick them off
            errorGobbler.start();
            outputGobbler.start();

            // any error???
            final int iExitValue = proc.waitFor();
            if(iExitValue != 0){
                System.err.println("WARNING: Orca returns non-zero return value (local optimization).");
                return FixedValues.UNEXCITABLEENERGY;
            } else{
                // orca should(!) have completed normally...
                saOutput = outputGobbler.getData();
            }
        } catch(Exception e){
            System.err.println("WARNING: Orca has a problem (local optimization)." + e.toString());
            return FixedValues.UNEXCITABLEENERGY;
        } catch(Error err){
            System.err.println("WARNING: Orca has a problem (local optimization)." + err.toString());
            err.printStackTrace(System.err);
            return FixedValues.UNEXCITABLEENERGY;
        }

        // read the output
        double dEnergy = FixedValues.UNEXCITABLEENERGY;
        for(int i = saOutput.length - 1; i >=0; i--){
            if(saOutput[i].startsWith("STATE  1")){
                final String[] sa = saOutput[i].split("\\s+");
                try{
                    dEnergy = Double.parseDouble(sa[3]);
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't parse excitation energy. " + e.toString());
                    dEnergy = FixedValues.UNEXCITABLEENERGY;
                }
            }
        }

        // clean up: delete directory
        if(!bDebug){
            try{
                SwitchesInput.removeFolder(sFolderName);
            }catch(Exception e){
                System.err.println("WARNING: Problem cleaning orca files up. " + e.toString());
            }
        }

        return dEnergy;
    }

    private String[] setInputUp(final CartesianCoordinates cartes,
            final String sXYZFile){

        final int iNoOfAtoms = cartes.getNoOfAtoms();
        final String[] saInput = new String[iNoOfAtoms + 9];

        // first the optimization job
	String sMeth;
	if(sWhichMethod.equalsIgnoreCase("bp")){
            sMeth = "RI " + sWhichMethod;
	} else{
            sMeth = sWhichMethod;
	}
        saInput[0] = "!RKS " + sMeth + " " + sBasisSet + " TightSCF";
        saInput[1] = "!Opt";
        saInput[2] = "*xyz " + iCharge + " " + iSpin;

        // coordinates
        final String[] saPrintCoords = cartes.createPrintableCartesians();
        for(int i = 0; i < iNoOfAtoms; i++){
            saInput[i+3] = saPrintCoords[i+2];
        }
        saInput[iNoOfAtoms + 3] = "*";

        // new job: TD-DFT
        saInput[iNoOfAtoms + 4] = "$new_job";
        saInput[iNoOfAtoms + 5] = saInput[0];
        saInput[iNoOfAtoms + 6] = "%tddft nroots " + iNoOfRoots;
        saInput[iNoOfAtoms + 7] = "       end";
        saInput[iNoOfAtoms + 8] = "* xyzfile " + iCharge + " " + iSpin + " "
                + sXYZFile;

        return saInput;
    }
}
