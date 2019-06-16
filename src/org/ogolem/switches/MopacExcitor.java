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
import org.ogolem.core.InitIOException;
import org.ogolem.core.StreamGobbler;

/**
 * Uses MOPAC to calculate vertical excitation energies
 * @author Johannes Dieterich
 * @version 2010-04-13
 */
final class MopacExcitor implements Excitor{

    private static final boolean bDebug = false;

    final String sWhichMethod;
    final boolean bUsePersico;
    final boolean bUseExternal;
    final int iNoOfAccElectrons;
    final int iNoOfMOs;

    // these are not used at the moment
    private static final int iCharge=0;
    private static final int iSpin=0;

    MopacExcitor(final SwitchesConfig config, final String method,
            final boolean persico, final boolean useexternal){
        this.sWhichMethod = method;
        this.bUsePersico = persico;
        this.bUseExternal = useexternal;
        // this obviously does NOT allow half occ MOs
        this.iNoOfAccElectrons = config.iActiveOcc*2;
        this.iNoOfMOs = config.iActiveNonOcc + config.iActiveOcc;
    }

    @Override
    public double s0s1TransitionEnergy(final CartesianCoordinates cartes, final int iID){

        // file stuff
        final String sMOPACInput = "mopac" + iID + ".dat";
        final String sMOPACOutput = "mopac" + iID + ".out";
        final String sMOPACDatBasis = "mopac" + iID;
        final String sMOPACBasis = "mopac" + iID + "exc";

        final String[] saAtoms = cartes.getAllAtomTypes();

        if(saAtoms.length != cartes.getAllXYZCoordsCopy()[0].length){
            System.err.println("ERROR: Atom types and coordinates not of same dimension.");
            return FixedValues.UNEXCITABLEENERGY;
        }

        /*
         * create a folder for this run
         */
        try{
            Output.createAFolder(sMOPACBasis);
        } catch(Exception e){
            System.err.println("ERROR: Can't creat folder. " + e.toString());
            return FixedValues.UNEXCITABLEENERGY;
        }

        /*
         * write the output aka input for the to be called program
         */
        try{
            Output.writeMopacExcInput(sMOPACBasis, sMOPACInput, sWhichMethod,
                    cartes.getAllXYZCoordsCopy(), saAtoms, iCharge, iSpin, 2000,
                    bUsePersico, bUseExternal, iNoOfAccElectrons, iNoOfMOs);
        } catch(InitIOException e){
            System.err.println("WARNING: Problem in writing geometry for MOPAC input (S0/S1 calculation)." + e.toString());
            return FixedValues.UNEXCITABLEENERGY;
        }

        /*
         * call mopac
         */
        Process proc = null;
        try{
            final Runtime rt = Runtime.getRuntime();
            final String sCommand = "mopac  " + sMOPACDatBasis;
            final String[] saEnvp = new String[0];
            final File dir = new File(sMOPACBasis);

            proc = rt.exec(sCommand,saEnvp,dir);

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
                System.err.println("WARNING: Mopac returns non-zero return value (local optimization).");
                return FixedValues.UNEXCITABLEENERGY;
            } else{
                // mopac should(!) have completed normally...
            }
        } catch(Exception e){
            System.err.println("Mopac has a problem (local optimization)." + e.toString());
            return FixedValues.UNEXCITABLEENERGY;
        }


        /*
         * read mopacs output
         */
        double dS0S1Energy = FixedValues.UNEXCITABLEENERGY;
        try{
            dS0S1Energy = SwitchesInput.readExcMopacOutput(
                    sMOPACBasis + System.getProperty("file.separator") + sMOPACOutput);
        } catch(Exception e){
            System.err.println("WARNING: Problem in reading the output of mopac for exciting.");
            e.printStackTrace(System.err);
            return  FixedValues.UNEXCITABLEENERGY;
        }

        /*
         * clean up
         */
        if(!bDebug){
            try{
                SwitchesInput.removeFolder(sMOPACBasis);
            }catch(Exception e){
                System.err.println("Problem cleaning mopac files up. " + e.toString());
            }
        }

        /*
         * return it
         */

        return dS0S1Energy;
    }
}
