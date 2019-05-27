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

import java.io.File;
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.InitIOException;
import org.ogolem.core.StreamGobbler;

/**
 * Uses MNDOs GugaCI for calculating excitation energies.
 * @author Johannes Dieterich
 * @version 2009-11-23
 */
final class GugaCI implements Excitor{

    private final int iWhichBackend;

    private final boolean bCOSMOWater;

    // these are not used at the moment
    private static final int iCharge=0;
    private static final int iSpin=0;

    /*
     * GUGA-CI related options
     */
    private final int iActiveOcc;
    private final int iActiveNonOcc;
    private final int iOccPi;
    private final int iNonOccPi;
    private final int iOrbDef;
    private final int iNoOfRefOcc;
    private final int iDefOfRefOcc;
    private final int iMaxExcitLevel;
    private final int iNoOfLowestStates;
    private final int iWantedCIState;

    GugaCI(final SwitchesConfig swConfig, final int iWhichMethod, final boolean bWaterSolv){
        this.iWhichBackend = iWhichMethod;
        this.bCOSMOWater = bWaterSolv;
        this.iActiveOcc = swConfig.iActiveOcc;
        this.iActiveNonOcc = swConfig.iActiveNonOcc;
        this.iOccPi = swConfig.iOccPi;
        this.iNonOccPi = swConfig.iNonOccPi;
        this.iOrbDef = swConfig.iOrbDef;
        this.iNoOfRefOcc = swConfig.iNoOfRefOcc;
        this.iDefOfRefOcc = swConfig.iDefOfRefOcc;
        this.iMaxExcitLevel = swConfig.iMaxExcitLevel;
        this.iNoOfLowestStates = swConfig.iNoOfLowestStates;
        this.iWantedCIState = swConfig.iWantedCiState;
    }

    @Override
    public double s0s1TransitionEnergy(final CartesianCoordinates cartes, final int iID){

        final String sMNDOInput = "mndo" + iID + ".inp";
        final String sMNDOOutput = "mndo" + iID + ".out";
        final String sMNDOBasis = "mndo" + iID;

        // get the atomic numbers together
        int[] iaCurrAtNos = new int[cartes.getAllAtomTypes().length];
        final String[] saAtoms = cartes.getAllAtomTypes();

        if(saAtoms.length != cartes.getAllXYZCoordsCopy()[0].length){
            System.err.println("ERROR: Atom types and coordinates not of same dimension.");
            return FixedValues.UNEXCITABLEENERGY;
        }

        for (int i = 0; i < saAtoms.length; i++) {
            iaCurrAtNos[i] = AtomicProperties.giveAtomicNumber(saAtoms[i]);
        }


        /*
         * create a folder for this run
         */
        try{
            Output.createAFolder(sMNDOBasis);
        } catch(Exception e){
            System.err.println("ERROR: Can't creat folder. " + e.toString());
            return FixedValues.UNEXCITABLEENERGY;
        }

        /*
         * write the output aka input for the to be called program
         */
        try{
            Output.writeMNDOInput(sMNDOBasis, sMNDOInput, iWhichBackend, cartes.getAllXYZCoordsCopy(), iaCurrAtNos,
                    iCharge, iSpin, bCOSMOWater, iActiveOcc, iActiveNonOcc, iOccPi, iNonOccPi, iOrbDef,
                    iNoOfRefOcc, iDefOfRefOcc, iMaxExcitLevel, iNoOfLowestStates, iWantedCIState);
        } catch(InitIOException e){
            System.err.println("WARNING: Problem in writing geometry for MNDO input (S0/S1 calculation)." + e.toString());
            return FixedValues.UNEXCITABLEENERGY;
        }

        /*
         * call mndo
         */
        Process proc = null;
        try{
            final Runtime rt = Runtime.getRuntime();
            final String sCommand = "mndo.sh " + sMNDOInput + " " + sMNDOOutput;
            final String[] saEnvp = new String[0];
            final File dir = new File(sMNDOBasis);

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
                System.err.println("WARNING: MNDO returns non-zero return value (local optimization).");
                return FixedValues.UNEXCITABLEENERGY;
            } else{
                // mndo should(!) have completed normally...
            }


        } catch(Exception e){
            System.err.println("WARNING: MNDO has a problem (local optimization)." + e.toString());
            return FixedValues.UNEXCITABLEENERGY;
        }

        /*
         * read mndo's output in
         */
        String[] saOutput;
        try{
            saOutput = SwitchesInput.readFileIn(sMNDOBasis + System.getProperty("file.separator") + sMNDOOutput);
        } catch(Exception e){
            System.err.println("WARNING: Couldn't read in MNDO output." + e.toString());
            return FixedValues.UNEXCITABLEENERGY;
        }


        /*
         * translate mndo's output
         */
        double dS0S1Energy = 0.0;
        try{
            dS0S1Energy = readEnergyIn(saOutput);
        } catch(Exception e){
            System.err.println("WARNING: Problem in translating the output of MNDO." + e);
            return FixedValues.UNEXCITABLEENERGY;
        }

        /*
         * clean up
         */
        try{
            SwitchesInput.removeFolder(sMNDOBasis);
        }catch(Exception e){
            System.err.println("WARNING: Problem cleaning MNDO files up." + e.toString());
        }

        /*
         * return it
         */
        return dS0S1Energy;
    }

    private static double readEnergyIn(final String[] saOutput){

        double dEnergy = FixedValues.UNEXCITABLEENERGY;

        // first check whether we even were successful
        for(int iLine = saOutput.length - 50; iLine < saOutput.length; iLine++){
            // should be within the last 50 lines
            final String sLine = saOutput[iLine].trim();
            if(sLine.equalsIgnoreCase("UNSUCCESSFUL GEOMETRY OPTIMIZATION.")){
                return FixedValues.UNEXCITABLEENERGY;
            }
        }

        for(int iLine = saOutput.length -1; iLine > -1; iLine--){
            final String sLine = saOutput[iLine].trim();
            if(sLine.startsWith("State  2,")){
                // correct line found
                String sEnergy = sLine.substring(28).trim();
                sEnergy = sEnergy.substring(0,sEnergy.indexOf("eV,")).trim();
                // now we should have the proper value
                try{
                    dEnergy = Double.parseDouble(sEnergy) * org.ogolem.core.Constants.EVTOHARTREE;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't parse energy difference in MNDO run. " + e.toString());
                    return FixedValues.UNEXCITABLEENERGY;
                }
            }
        }

        return dEnergy;
    }
}
