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
package org.ogolem.switches;

import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.InitIOException;
import org.ogolem.core.StreamGobbler;

/**
 * Call MOPAC to locally optimize a switch.
 * @author Johannes Dieterich
 * @version 2014-09-06
 */
final class MopacLocOpt implements LocalOptimization{

    private static final boolean bDebug = false;

    private final int iWhichBackend;

    private final int iNoOfCycles;

    /**
     * The constructor for usage as a local optimizing engine.
     */
    MopacLocOpt(final int iWhichMethod, final int iMaxIterLocOpt){
        this.iWhichBackend = iWhichMethod;
        this.iNoOfCycles = iMaxIterLocOpt;
    }

    @Override
    public CartesianCoordinates locOptThis(final CartesianCoordinates cartes, int iID){

        final String sMopacInput = "mopac" + iID + ".dat";
        final String sMopacBasis = "mopac" + iID;

        final float[] faCharges = cartes.getAllCharges();
        final short[] iaSpins = cartes.getAllSpins();

        /*
         * write the output aka input for the to be called program
         */
        try{
            WriteOutput(sMopacInput,cartes.getAllXYZCoordsCopy(),cartes.getAllAtomTypes(), cartes.getTotalCharge(),
                    cartes.getTotalSpin());
        } catch(InitIOException e){
            System.err.println("ERROR: Problem in writing geometry for mopac input (local optimization). " + e.toString());
        }

        /*
         * call mopac
         */
        try{
            final Runtime rt = Runtime.getRuntime();
            final Process proc = rt.exec("mopac " + sMopacBasis);

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
                return null;
            } else{
                // mopac should(!) have completed normally...
            }


        } catch(Exception e){
            System.err.println("Mopac has a problem (local optimization)." + e.toString());
            return null;
        }


        /*
         * read mopacs output
         */
        final String sMopacOutput = sMopacInput.substring(0,sMopacInput.indexOf(".")) + ".out";

        CartesianCoordinates newCartes = new CartesianCoordinates(cartes);

        try{
            newCartes = SwitchesInput.readXYZMopacOutput(sMopacOutput, cartes.getNoOfAtoms(), cartes.getNoOfMolecules(), cartes.getAllAtomsPerMol());
        } catch(Exception e){
            System.err.println("WARNING: Problem in reading the output of mopac." + e.toString());
        }

        /*
         * clean up
         */
        if(!bDebug){
            try{
                org.ogolem.core.Input.RemoveMopacFiles(sMopacInput);
            }catch(Exception e){
                System.err.println("Problem cleaning mopac files up." + e.toString());
            }
        }

        /*
         * return it
         */
        newCartes.setAllCharges(faCharges);
        newCartes.setAllSpins(iaSpins);

        return newCartes;
    }


    private void WriteOutput(String sMopacInput, double[][] daXYZ, String[] saAtoms, final int iTotalCharge,
            final int iTotalSpin) throws InitIOException{

        final int iSpin = iTotalSpin;

        switch(iWhichBackend){
            case 0:
                try{
                    Output.writeMopacInput(sMopacInput, "mndo", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles);
                }catch(InitIOException e){
                    throw e;
                };
                break;
            case 1:
                try{
                    Output.writeMopacInput(sMopacInput, "am1", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles);
                }catch(InitIOException e){
                    throw e;
                };
                break;
            case 2:
                try{
                    Output.writeMopacInput(sMopacInput, "pm3", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles);
                }catch(InitIOException e){
                    throw e;
                };
                break;
            case 3:
                try {
                    Output.writeMopacInput(sMopacInput, "pm5", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles);
                } catch (InitIOException e) {
                    throw e;
                };
                break;
            case 4:
                try {
                    Output.writeMopacInput(sMopacInput, "pm6", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles);
                } catch (InitIOException e) {
                    throw e;
                };
                break;
            default:
                try{
                    Output.writeMopacInput(sMopacInput, "am1", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles);
                }catch(InitIOException e){
                    throw e;
                };
                break;
        }
    }
}
