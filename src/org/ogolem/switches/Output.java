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
package org.ogolem.switches;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.text.DecimalFormat;
import static org.ogolem.core.Constants.*;
import org.ogolem.core.InitIOException;
import org.ogolem.core.SerialException;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.io.OutputPrimitives;

/**
 * All output related functions.
 * @author Johannes Dieterich
 * @version 2012-11-05
 */
final class Output{

    static void printMiscToFile(final String sFileName, final String[] saFileContent) throws IOException{
        WriteNow(sFileName, saFileContent, false);
    }

    static void createAFolder(final String sFolderName) throws IOException{
        final File f = new File(sFolderName);

        // try to figure whether the file exists
        if(f.exists()){
            return;
        }

        boolean bSuccess = f.mkdirs();
        if (bSuccess) {
            // Folder successfully and freshly created, no further action needed.
        } else {
            // Folder exists already
            System.err.println("Failure in creating needed folder " + sFolderName + ". Folder exists.");
        }
    }

    static void writeMNDOInput(final String sFolder, final String sInputFile,
            final int iWhichMethod, final double[][] daXYZ, final int[] iaAtomicNos,
            final int iCharge, final int iSpin, final boolean bCOSMOWater, final int iActiveOcc,
            final int iActiveNonOcc, final int iOccPi, final int iNonOccPi, final int iOrbDef,
            final int iNoOfRefOcc, final int iDefOfRefOcc, final int iMaxExcitLevel,
            final int iNoOfLowestStates, final int iWantedCIState)
            throws InitIOException{

        // method ID to MNDO method ID
        int iMethodID;
        switch(iWhichMethod){
            case 0:
                // MNDO/d
                iMethodID = -10;
                break;
            case 1:
                // OM3
                iMethodID = -8;
                break;
            case 2:
                // PM3
                iMethodID = -7;
                break;
            case 3:
                // OM2
                iMethodID = -6;
                break;
            case 4:
                // OM1
                iMethodID = -5;
                break;
            case 5:
                // AM1
                iMethodID = -2;
                break;
            case 6:
                // MNDOC
                iMethodID = -1;
                break;
            case 7:
                // MNDO
                iMethodID = 0;
                break;
            case 8:
                // MINDO/3
                iMethodID = 1;
                break;
            case 9:
                // CNDO/2
                iMethodID = 2;
                break;
            case 10:
                // SCC-DFTB
                iMethodID = 5;
                break;
            case 11:
                // SCC-DFTB w/ Jorgensens correction
                iMethodID = 6;
                break;
            default:
                iMethodID = -8;
                System.err.println("ERROR: MNDO method to method translation. " +
                        "Contact the author. Using OM3 now.");
        }

        final String[] saInput = new String[iaAtomicNos.length + 7];
        saInput[0] = "iop=" + iMethodID +" igeom=1 jop=0 iform=1 mplib=0 inrefd=2 kci=5 +";
        /*
         * small keyword explanantion (since they are rather cryptic):
         * igeom=1 : use cartesian coordinates
         * jop=0   : optimization for energy minimum
         * iform=1 : fortran free format input (are they serious?!)
         * mplib=0 : use one core only. it doesn't really follow that, hope remains though.
         * inrefd=2: detailed output
         * kci=5   : correlation treatment using GUGA-CI (wouldn't be the GugaCICaller otherwise, right? ;-) )
         */
        
        saInput[1] = "ici1=" + iActiveOcc + " ici2=" + iActiveNonOcc + " ioutci=2 movo=" + iOrbDef ;
        saInput[1] += " nciref=" + iNoOfRefOcc + " mciref=" + iDefOfRefOcc + " levexc=" + iMaxExcitLevel + " +";
        /*
         * ioutci=2: printing flag
         * movo=   : definition of orbitals involved in the active CI space
         * nciref= : number of reference occupations
         * mciref= : definition of reference occupations
         * levexc= : maximum excitation level relative to any of the reference configurations
         */

        saInput[2] = "iroot=" + iNoOfLowestStates + " lroot=" + iWantedCIState + " multci=0 jci1=" + iOccPi + " cidir=1 jci2=" + iNonOccPi;
        /*
         * iroot=  : total number of lowest CI states computed
         * lroot=  : if(>0): number of wanted CI state, if(<0): state with the biggest norm
         * multci= : spin multiplicity of CI states
         * jci1=   : total number of occ pi MOs
         * cidir=  : direct CI features (1: keep in mem)
         * jci2=   : total number of non occ pi MOs
         */

        // water solvation using the COSMO modell
        if(bCOSMOWater){
            saInput[2] += " icosmo=1 +";
        } else{
            // no solvation
            saInput[2] += " +";
        }

        /*
         * MNDO's definition of spin is (unfortunately...) slightly different
         * from normal.
         * 1 means two single occupied orbitals forming in total a singlet
         * 0 means a normal singlet state
         * the rest is as usual
         * WE DO NOT AT THE MOMENT USE THE SPECIAL SINGLET FEATURE!
         */
        int iMult;
        if(iSpin == 0){
            iMult = 0;
        } else {
            iMult = iSpin + 1;
        }
        saInput[3] = "kharge=" + iCharge +" imult=" + iMult + " icore=1024";
        saInput[4] = "Automatically created by OGOLEM";
        saInput[5] = " ";

        // coordinates
        final DecimalFormat form = new DecimalFormat("0.00000000");
        for(int i = 0; i < iaAtomicNos.length; i++){
            saInput[i+6] = " " + iaAtomicNos[i] + "   " + form.format(daXYZ[0][i] * BOHRTOANG)
                    + " 1 " + form.format(daXYZ[1][i] * BOHRTOANG) + " 1 "
                    + form.format(daXYZ[2][i]*BOHRTOANG) + " 1 ";
        }

        // and an empty line in the end is important...
        saInput[saInput.length -1] = "";

        // write it out
        final String sMNDOInp = sFolder + System.getProperty("file.separator") + sInputFile;
        try{
            WriteNow(sMNDOInp, saInput, false);
        } catch(Exception e){
            throw new InitIOException(e);
        }
    }

    private static void WriteNow(final String sOutputPath, final String[] saToWrite,
            final boolean bNoOverride) throws IOException {
        BufferedWriter buffwriter = null;

        try {
            buffwriter = new BufferedWriter(new FileWriter(sOutputPath, bNoOverride));
            for (final String sCont: saToWrite) {
                buffwriter.write(sCont);
                buffwriter.write(System.getProperty("line.separator"));
            }
        } catch (IOException e) {
            throw e;
        } finally {
            if (buffwriter != null) {
                try {
                    buffwriter.close();
                } catch (IOException e) {
                    throw e;
                }
            }
        }
    }

    static void writeMopacInput(final String sMopacInput, final String sMethod, final double[][] daXYZ,
            final String[] saAtoms, final int iTotalCharge, final int iTotalSpin, final int iNoOfCycles)
            throws InitIOException {

        // case switching for the spins
        String sSpin;

        // the multiplicity (since we operate with full spin units, just 1s+1)
        final int iMultiplicity = 1 * iTotalSpin + 1;

        switch (iMultiplicity) {
            // it is a singlet, but we omit the keyword since it is default
            case 1:
                sSpin = "";
                break;

            case 2:
                sSpin = "DOUBLET";
                break;

            case 3:
                sSpin = "TRIPLET";
                break;

            case 4:
                sSpin = "QUARTET";
                break;

            case 5:
                sSpin = "QUINTET";
                break;

            case 6:
                sSpin = "SEXTET";
                break;

            case 7:
                sSpin = "SEPTET";
                break;

            case 8:
                sSpin = "OCTET";
                break;

            case 9:
                sSpin = "NONET";
                break;

            default:
                sSpin = "";
                break;
        }

        final int iNoOfAtoms = saAtoms.length;
        final String[] saOutput = new String[iNoOfAtoms + 3];
        // first line: the method and all that
        saOutput[0] = "XYZ NOLOG GEO-OK T=100H CYCLES=" + iNoOfCycles + " " + sMethod + " charge=" + iTotalCharge + " " + sSpin;

        // the next line: a tag
        saOutput[1] = "CREATED BY OGOLEM";

        // third line: empty
        saOutput[2] = "";

        // the geometry in the wonderful MOPAC xyz format...
        for (int i = 3; i < iNoOfAtoms + 3; i++) {
            saOutput[i] = saAtoms[i - 3] + "\t" + (daXYZ[0][i - 3] * BOHRTOANG) + " 1"
                    + "\t" + (daXYZ[1][i - 3] * BOHRTOANG) + " 1"
                    + "\t" + (daXYZ[2][i - 3] * BOHRTOANG) + " 1";
        }

        //write it actually
        try {
            WriteNow(sMopacInput, saOutput, false);
        } catch (IOException e) {
            throw new InitIOException(e);
        }
    }

    static void writeMopacExcInput(final String sMopacDir, final String sMopacInput, final String sMethod,
            final double[][] daXYZ, final String[] saAtoms, final int iTotalCharge, final int iTotalSpin,
            final int iNoOfCycles, final boolean bUsePersico, final boolean bUseExternal,
            final int iNoOfAccElectrons, final int iNoOfMOs)
            throws InitIOException {
         // case switching for the spins
        String sSpin;

        // the multiplicity (since we operate with full spin units, just 1s+1)
        final int iMultiplicity = 1 * iTotalSpin + 1;

        switch (iMultiplicity) {
            // it is a singlet, but we omit the keyword since it is default
            case 1:
                sSpin = "SINGLET";
                break;

            case 2:
                sSpin = "DOUBLET";
                break;

            case 3:
                sSpin = "TRIPLET";
                break;

            case 4:
                sSpin = "QUARTET";
                break;

            case 5:
                sSpin = "QUINTET";
                break;

            case 6:
                sSpin = "SEXTET";
                break;

            case 7:
                sSpin = "SEPTET";
                break;

            case 8:
                sSpin = "OCTET";
                break;

            case 9:
                sSpin = "NONET";
                break;

            default:
                sSpin = "";
                break;
        }

        // since we want to do excitation business here, we need to populate the OPEN(,) keyword
        String sOpen = "OPEN(" + iNoOfAccElectrons +"," + iNoOfMOs + ")";

        String sAdditional="";
        if(bUsePersico && bUseExternal){

            // first copy the parameter file
            final String sParamFile = "bestpar";
            final String sParamEndFile = sMopacDir + System.getProperty("file.separator") + "bestpar";
            try{
                OutputPrimitives.copyFile(sParamFile, sParamEndFile);
            } catch(Exception e){
                throw new InitIOException("Couldn't copy parameters file bestpars.",e);
            }

            sAdditional += " EXTERNAL=bestpar FLOCC=0.10 MICROS=94 ";
        } else if(bUsePersico && !bUseExternal){
            sAdditional += " FLOCC=0.10 MICROS=94 ";
        }

        final int iNoOfAtoms = saAtoms.length;

        String[] saOutput;
        
        if(!bUsePersico){
            saOutput = new String[iNoOfAtoms + 4];
        } else{
            saOutput = new String[iNoOfAtoms + 104];
        }

        // first line: the method and all that
        saOutput[0] = sMethod + " VECTORS " + sOpen + " PREC "+ sSpin + " ROOT=1 GEO-OK +";
        saOutput[1] = sAdditional + " PULAY MICROS=94 MECI=20";

        // the next line: a tag
        saOutput[2] = "CREATED BY OGOLEM";

        // third line: empty
        saOutput[3] = "";

        // the geometry in the wonderful MOPAC xyz format...
        int iCounter = 4;
        for (int i = 4; i < iNoOfAtoms + 4; i++) {
            saOutput[i] = saAtoms[i - 4] + "\t" + (daXYZ[0][i - 4] * BOHRTOANG) + " 1"
                    + "\t" + (daXYZ[1][i - 4] * BOHRTOANG) + " 1"
                    + "\t" + (daXYZ[2][i - 4] * BOHRTOANG) + " 1";
            iCounter++;
        }

        if(bUsePersico){
            saOutput[iCounter   ] = "";
            saOutput[iCounter+1 ] = "OCCUP";
            //TODO consider whether this shouldn't be made better
            saOutput[iCounter+2 ] = "2.0 2.0 2.0 2.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0";
            saOutput[iCounter+3 ] = "";
            saOutput[iCounter+4 ] = "MICROS";
            saOutput[iCounter+5 ] = "11111110000001111111000000";
            saOutput[iCounter+6 ] = "11111101000001111110100000";
            saOutput[iCounter+7 ] = "11111011000001111110100000";
            saOutput[iCounter+8 ] = "11110111000001111110100000";
            saOutput[iCounter+9 ] = "11111101000001111101100000";
            saOutput[iCounter+10] = "11111011000001111101100000";
            saOutput[iCounter+11] = "11110111000001111101100000";
            saOutput[iCounter+12] = "11111101000001111011100000";
            saOutput[iCounter+13] = "11111011000001111011100000";
            saOutput[iCounter+14] = "11110111000001111011100000";
            saOutput[iCounter+15] = "11111101000001111111000000";
            saOutput[iCounter+16] = "11111100100001111111000000";
            saOutput[iCounter+17] = "11111100010001111111000000";
            saOutput[iCounter+18] = "11111100001001111111000000";
            saOutput[iCounter+19] = "11111100000101111111000000";
            saOutput[iCounter+20] = "11111100000011111111000000";
            saOutput[iCounter+21] = "11111011000001111111000000";
            saOutput[iCounter+22] = "11111010100001111111000000";
            saOutput[iCounter+23] = "11111010010001111111000000";
            saOutput[iCounter+24] = "11111010001001111111000000";
            saOutput[iCounter+25] = "11111010000101111111000000";
            saOutput[iCounter+26] = "11111010000011111111000000";
            saOutput[iCounter+27] = "11110111000001111111000000";
            saOutput[iCounter+28] = "11110110100001111111000000";
            saOutput[iCounter+29] = "11110110010001111111000000";
            saOutput[iCounter+30] = "11110110001001111111000000";
            saOutput[iCounter+31] = "11110110000101111111000000";
            saOutput[iCounter+32] = "11110110000011111111000000";
            saOutput[iCounter+33] = "11101111000001111111000000";
            saOutput[iCounter+34] = "11101110100001111111000000";
            saOutput[iCounter+35] = "11101110010001111111000000";
            saOutput[iCounter+36] = "11101110001001111111000000";
            saOutput[iCounter+37] = "11101110000101111111000000";
            saOutput[iCounter+38] = "11101110000011111111000000";
            saOutput[iCounter+39] = "11011111000001111111000000";
            saOutput[iCounter+40] = "11011110100001111111000000";
            saOutput[iCounter+41] = "11011110010001111111000000";
            saOutput[iCounter+42] = "11011110001001111111000000";
            saOutput[iCounter+43] = "11011110000101111111000000";
            saOutput[iCounter+44] = "11011110000011111111000000";
            saOutput[iCounter+45] = "10111111000001111111000000";
            saOutput[iCounter+46] = "10111110100001111111000000";
            saOutput[iCounter+47] = "10111110010001111111000000";
            saOutput[iCounter+48] = "10111110001001111111000000";
            saOutput[iCounter+49] = "10111110000101111111000000";
            saOutput[iCounter+50] = "10111110000011111111000000";
            saOutput[iCounter+51] = "11111110000001111110100000";
            saOutput[iCounter+52] = "11111110000001111110010000";
            saOutput[iCounter+53] = "11111110000001111110001000";
            saOutput[iCounter+54] = "11111110000001111110000100";
            saOutput[iCounter+55] = "11111110000001111110000010";
            saOutput[iCounter+56] = "11111110000001111110000001";
            saOutput[iCounter+57] = "11111110000001111101100000";
            saOutput[iCounter+58] = "11111110000001111101010000";
            saOutput[iCounter+59] = "11111110000001111101001000";
            saOutput[iCounter+60] = "11111110000001111101000100";
            saOutput[iCounter+61] = "11111110000001111101000010";
            saOutput[iCounter+62] = "11111110000001111101000001";
            saOutput[iCounter+63] = "11111110000001111011100000";
            saOutput[iCounter+64] = "11111110000001111011010000";
            saOutput[iCounter+65] = "11111110000001111011001000";
            saOutput[iCounter+66] = "11111110000001111011000100";
            saOutput[iCounter+67] = "11111110000001111011000010";
            saOutput[iCounter+68] = "11111110000001111011000001";
            saOutput[iCounter+69] = "11111110000001110111100000";
            saOutput[iCounter+70] = "11111110000001110111010000";
            saOutput[iCounter+71] = "11111110000001110111001000";
            saOutput[iCounter+72] = "11111110000001110111000100";
            saOutput[iCounter+73] = "11111110000001110111000010";
            saOutput[iCounter+74] = "11111110000001110111000001";
            saOutput[iCounter+75] = "11111110000001101111100000";
            saOutput[iCounter+76] = "11111110000001101111010000";
            saOutput[iCounter+77] = "11111110000001101111001000";
            saOutput[iCounter+78] = "11111110000001101111000100";
            saOutput[iCounter+79] = "11111110000001101111000010";
            saOutput[iCounter+80] = "11111110000001101111000001";
            saOutput[iCounter+81] = "11111110000001011111100000";
            saOutput[iCounter+82] = "11111110000001011111010000";
            saOutput[iCounter+83] = "11111110000001011111001000";
            saOutput[iCounter+84] = "11111110000001011111000100";
            saOutput[iCounter+85] = "11111110000001011111000010";
            saOutput[iCounter+86] = "11111110000001011111000001";
            saOutput[iCounter+87] = "01111111000001111111000000";
            saOutput[iCounter+88] = "01111110100001111111000000";
            saOutput[iCounter+89] = "01111110010001111111000000";
            saOutput[iCounter+90] = "01111110001001111111000000";
            saOutput[iCounter+91] = "01111110000101111111000000";
            saOutput[iCounter+92] = "01111110000011111111000000";
            saOutput[iCounter+93] = "11111110000000111111100000";
            saOutput[iCounter+94] = "11111110000000111111010000";
            saOutput[iCounter+95] = "11111110000000111111001000";
            saOutput[iCounter+96] = "11111110000000111111000100";
            saOutput[iCounter+97] = "11111110000000111111000010";
            saOutput[iCounter+98] = "11111110000000111111000001";
            saOutput[iCounter+99] = "";
        }
        

        //write it actually
        try {
            WriteNow(sMopacDir + System.getProperty("file.separator") + sMopacInput, saOutput, false);
        } catch (IOException e) {
            throw new InitIOException(e);
        }
    }

    static void writeBinaryPool(final GenericPool<Color,Switch> pool, final boolean bTemp) throws SerialException{

        final String sBinOutputPath = (bTemp) ? "snapswitchpool.bin" : "switchpool.bin";

        ObjectOutputStream outStream = null;
        try {
            outStream = new ObjectOutputStream(new FileOutputStream(sBinOutputPath));
            outStream.writeObject(pool);
        } catch (IOException e) {
            throw new SerialException("Error occured during binary writing!", e);
        } finally {
            if (outStream != null) {
                try {
                    outStream.close();
                } catch (IOException e) {
                    throw new SerialException("Error occured during binary writing!", e);
                }
            }
        }
    }
}
