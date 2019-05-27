/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
              2017, J. M. Dieterich and B. Hartke
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

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Iterator;
import java.util.LinkedList;
import org.ogolem.core.CartesianCoordinates;
import static org.ogolem.core.Constants.*;
import org.ogolem.core.ZMatrix;
import org.ogolem.io.OutputPrimitives;

/**
 * The set of output functions.
 * @author Johannes Dieterich
 * @version 2017-11-16
 */
public class Output{

    static void writeAdaptOrcaInput(String sOrcaInput, CartesianCoordinates cartes, boolean bLocOpt,
            boolean bGradient, boolean bEnergy, AdaptiveParameters params)
            throws IOException{

        String sMethod = params.getForWhichMethod();

        LinkedList<String> llInput = new LinkedList<>();

        // first the general thing
        llInput.add("!" + sMethod + " TightSCF");

        if(bLocOpt){
            llInput.add("!Opt");
        } else if(bGradient){
            // remains to be seen how we code it in orca
            //TODO once the work in orca is done
        } else if(bEnergy){
            // by default done, no extra input required
        }

        // explicitly specify the parameters
        String[] saForAtoms = params.getForWhichAtoms();
        for(int i = 0; i < saForAtoms.length; i++){
            
        }
        //TODO all the orca functionality

        // add the geometry as xyz (remember charge and spin)
        llInput.add("*xyz " + cartes.getTotalCharge() + " " + cartes.getTotalSpin());

        String[] saAtoms = cartes.getAllAtomTypes();
        double[][] daXYZ = cartes.getAllXYZCoordsCopy();

        for(int i = 0; i < saAtoms.length; i++){
            llInput.add(saAtoms[i] + "\t" + daXYZ[0][i] * BOHRTOANG + "\t"
                    + daXYZ[1][i] * BOHRTOANG + "\t" + daXYZ[2][i] * BOHRTOANG);
        }
        llInput.add("*");
        
        // copy the data
        String[] saInput = new String[llInput.size()];
        for(int i = 0; i < llInput.size(); i++){
            saInput[i] = llInput.get(i);
        }

        // actually write it
        PrintMiscToFile(sOrcaInput, saInput);
    }

    static void PrintMiscToFile (String sFile, String[] saToBePrinted) throws IOException {
        OutputPrimitives.writeOut(sFile, saToBePrinted, true);
    }

    static void writeAdaptMopacInput (final String sMopacFolder, final String sMopacInput, final CartesianCoordinates cartes,
            final boolean bLocOpt, final boolean bEnergy, final AdaptiveParameters params,
            final boolean bAzoBenzene, final String sMopacParams) throws IOException{

        final LinkedList<String> llInput = new LinkedList<>();
        final LinkedList<String> llParamInp = new LinkedList<>();

        final int iNoOfAtoms = cartes.getNoOfAtoms();
        final String[] saAtoms = cartes.getAllAtomTypes();
        final double[][] daXYZ = cartes.getAllXYZCoordsCopy();

        final double[] daParams = params.getAllParamters();
        final String[] saParamAtoms = params.getForWhichAtoms();

        final int iTotalSpin = cartes.getTotalSpin();

        // first create the mopac folder
        OutputPrimitives.createAFolder(sMopacFolder);

        if(bAzoBenzene){
            /*
             * we are doing azobenzene specific input. this sucks a little bit
             * but it seems to be the only practicable solution, IMHO
             */

            /*
             * ATTENTION: This is highly azobenzene specific and can be considered being rather ugly.
             * But it really is the only possibility to accomplish the (non-targeted!) use with
             * excited states having the same spin but still being different
             * So we assume SINGLET as given and consider the total "spin" just as an identifier
             * to WHICH PES we look at
             * on top, we assume that the reference geometries were given as zmatrices, so that the
             * zmatrix field is populated in the reference cartesian coordinates
             */
            final int iState = iTotalSpin;

            llInput.add("AM1 VECTORS OPEN(14,13) PREC FLOCC=0.1 GEO-OK SINGLET ROOT="+ iState + " +");
            llInput.add("ITRY=5000 PULAY EXTERNAL=" + sMopacParams + " MICROS=94 1SCF");
            llInput.add("input automatically created by OGOLEM");
            llInput.add("");

            final ZMatrix zmat = cartes.getZMatrices()[0];
            // the coordinates
            final double[] daBonds = zmat.getAllBondLengths();
            final double[] daAngles = zmat.getAllBondAngles();
            final double[] daDihedrals = zmat.getAllDihedrals();
            final int[] iaBondConn = zmat.getAllBondConnects();
            final int[] iaAngleConn = zmat.getAllAnglesConnects();
            final int[] iaDihedralConn = zmat.getAllDihedralConnects();

            final double dIntToAng = 180 / Math.PI;

            final DecimalFormat mop = new DecimalFormat("0.0000");

            for(int i = 0; i < iNoOfAtoms; i++){
                // ID bond 0 angle 0 dihedral 0 bondconn angleconn dihedralconn
                String sTemp = saAtoms[i];
                String sTemp2 = mop.format(daBonds[i] * BOHRTOANG);
                while(sTemp2.length() < 12){
                    sTemp2 = " " + sTemp2;
                }
                sTemp = sTemp + sTemp2 + "   0   ";
                sTemp2 = mop.format(daAngles[i] * dIntToAng);
                while(sTemp2.length() < 12){
                    sTemp2 = " " + sTemp2;
                }
                sTemp = sTemp + sTemp2 + "   0   ";
                sTemp2 = mop.format(daDihedrals[i] * dIntToAng);
                while(sTemp2.length() < 12){
                    sTemp2 = " " + sTemp2;
                }
                sTemp = sTemp + sTemp2 + "   0   ";
                if(i == 0){
                    sTemp = sTemp + "  0    0    0";
                } else if (i == 1){
                    sTemp = sTemp + "  " + (iaBondConn[i]+1)  + "    0    0";
                } else if (i == 2){
                    sTemp = sTemp + "  " + (iaBondConn[i]+1)  + "    "+ (iaAngleConn[i]+1) +"    0";
                } else{
                    sTemp = sTemp + (iaBondConn[i] + 1) + "   " + (iaAngleConn[i] + 1) + "   ";
                    sTemp = sTemp + (iaDihedralConn[i] + 1);
                }
                llInput.add(sTemp);
            }

            // more azobenzene specific stuff
            llInput.add("");
            llInput.add("OCCUP");
            llInput.add("2.0 2.0 2.0 2.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0");
            llInput.add("");
            llInput.add("MICROS");
            llInput.add("11111110000001111111000000");
            llInput.add("11111101000001111110100000");
            llInput.add("11111011000001111110100000");
            llInput.add("11110111000001111110100000");
            llInput.add("11111101000001111101100000");
            llInput.add("11111011000001111101100000");
            llInput.add("11110111000001111101100000");
            llInput.add("11111101000001111011100000");
            llInput.add("11111011000001111011100000");
            llInput.add("11110111000001111011100000");
            llInput.add("11111101000001111111000000");
            llInput.add("11111100100001111111000000");
            llInput.add("11111100010001111111000000");
            llInput.add("11111100001001111111000000");
            llInput.add("11111100000101111111000000");
            llInput.add("11111100000011111111000000");
            llInput.add("11111011000001111111000000");
            llInput.add("11111010100001111111000000");
            llInput.add("11111010010001111111000000");
            llInput.add("11111010001001111111000000");
            llInput.add("11111010000101111111000000");
            llInput.add("11111010000011111111000000");
            llInput.add("11110111000001111111000000");
            llInput.add("11110110100001111111000000");
            llInput.add("11110110010001111111000000");
            llInput.add("11110110001001111111000000");
            llInput.add("11110110000101111111000000");
            llInput.add("11110110000011111111000000");
            llInput.add("11101111000001111111000000");
            llInput.add("11101110100001111111000000");
            llInput.add("11101110010001111111000000");
            llInput.add("11101110001001111111000000");
            llInput.add("11101110000101111111000000");
            llInput.add("11101110000011111111000000");
            llInput.add("11011111000001111111000000");
            llInput.add("11011110100001111111000000");
            llInput.add("11011110010001111111000000");
            llInput.add("11011110001001111111000000");
            llInput.add("11011110000101111111000000");
            llInput.add("11011110000011111111000000");
            llInput.add("10111111000001111111000000");
            llInput.add("10111110100001111111000000");
            llInput.add("10111110010001111111000000");
            llInput.add("10111110001001111111000000");
            llInput.add("10111110000101111111000000");
            llInput.add("10111110000011111111000000");
            llInput.add("11111110000001111110100000");
            llInput.add("11111110000001111110010000");
            llInput.add("11111110000001111110001000");
            llInput.add("11111110000001111110000100");
            llInput.add("11111110000001111110000010");
            llInput.add("11111110000001111110000001");
            llInput.add("11111110000001111101100000");
            llInput.add("11111110000001111101010000");
            llInput.add("11111110000001111101001000");
            llInput.add("11111110000001111101000100");
            llInput.add("11111110000001111101000010");
            llInput.add("11111110000001111101000001");
            llInput.add("11111110000001111011100000");
            llInput.add("11111110000001111011010000");
            llInput.add("11111110000001111011001000");
            llInput.add("11111110000001111011000100");
            llInput.add("11111110000001111011000010");
            llInput.add("11111110000001111011000001");
            llInput.add("11111110000001110111100000");
            llInput.add("11111110000001110111010000");
            llInput.add("11111110000001110111001000");
            llInput.add("11111110000001110111000100");
            llInput.add("11111110000001110111000010");
            llInput.add("11111110000001110111000001");
            llInput.add("11111110000001101111100000");
            llInput.add("11111110000001101111010000");
            llInput.add("11111110000001101111001000");
            llInput.add("11111110000001101111000100");
            llInput.add("11111110000001101111000010");
            llInput.add("11111110000001101111000001");
            llInput.add("11111110000001011111100000");
            llInput.add("11111110000001011111010000");
            llInput.add("11111110000001011111001000");
            llInput.add("11111110000001011111000100");
            llInput.add("11111110000001011111000010");
            llInput.add("11111110000001011111000001");
            llInput.add("01111111000001111111000000");
            llInput.add("01111110100001111111000000");
            llInput.add("01111110010001111111000000");
            llInput.add("01111110001001111111000000");
            llInput.add("01111110000101111111000000");
            llInput.add("01111110000011111111000000");
            llInput.add("11111110000000111111100000");
            llInput.add("11111110000000111111010000");
            llInput.add("11111110000000111111001000");
            llInput.add("11111110000000111111000100");
            llInput.add("11111110000000111111000010");
            llInput.add("11111110000000111111000001");
            

            // now the parameter file
            int iCounter = 0;
            for(int i = 0; i < saParamAtoms.length; i++){
                if(saParamAtoms[i].equalsIgnoreCase("H")){
                    // no params ;-)
                } else if(saParamAtoms[i].equalsIgnoreCase("C")){
                    llParamInp.add("USS\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("UPP\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("BETST\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("BETPT\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("BETAS\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("BETAP\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("ZS\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("ZP\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("ALP\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GSS\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GSP\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GPP\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GP2\t C\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("HSP\t C\t" + daParams[iCounter]);
                    iCounter++;
                } else if(saParamAtoms[i].equalsIgnoreCase("N")){
                    llParamInp.add("USS\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("UPP\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("BETAS\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("BETAP\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GSS\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GSP\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GPP\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GP2\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("HSP\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("ZS\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("ZP\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("ALP\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN11\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN12\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN13\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN21\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN22\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN23\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN31\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN32\t N\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN33\t N\t" + daParams[iCounter]);
                    iCounter++;
                } else {
                    throw new IOException("ERROR: We have a not known element in the azobenzene input.");
                }
            }
        } else {
            
            final int iCharge = cartes.getTotalCharge();

            final int iMult = cartes.getTotalSpin() + 1;
            String sSpin;

            switch (iMult) {
                case 1:
                    // it is a singlet, but we omit the keyword since it is default
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
                    System.err.println("WARNING: We cannot translate the spin of the" +
                            " molecule to mopac input. Using SINGLET.");
                    sSpin = "";
                    break;
            }


            if(bEnergy){
                llInput.add("XYZ NOLOG GEO-OK " + params.getForWhichMethod() + " 1SCF charge=" + iCharge
                + " " + sSpin + " +");
                llInput.add("EXTERNAL=params");
            } else if(bLocOpt){
                llInput.add("XYZ NOLOG GEO-OK " + params.getForWhichMethod() + " charge=" + iCharge
                + " " + sSpin + " +");
                llInput.add("EXTERNAL=params");
            } else {
                throw new IOException("ERROR: Failure in writing MOPAC input. You shouldn't be here. Notify the author.");
            }

            // coordinates, independent of the method
            llInput.add("automatically generated by OGOLEM");
            llInput.add("");
            for(int i = 0; i < iNoOfAtoms; i++){
                llInput.add("" + saAtoms[i] + "\t" +(daXYZ[0][i]*BOHRTOANG)
                        + " 1" + "\t" + (daXYZ[1][i]*BOHRTOANG) + " 1"
                        + "\t" + (daXYZ[2][i]*BOHRTOANG) + " 1");
            }

            // parameters
            int iCounter = 0;
            for(int i = 0; i < saParamAtoms.length; i++){
                if(saParamAtoms[i].equalsIgnoreCase("H")){
                    // hydrogen is special (less parameters)
                    llParamInp.add("USS\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("BETAS\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("ZS\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("ALP\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GSS\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN11\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN12\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN13\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN21\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN22\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN23\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN31\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN32\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN33\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                } else {
                    // all others are the same
                    llParamInp.add("USS\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("UPP\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("BETAS\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("BETAP\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GSS\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GSP\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GPP\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("GP2\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("HSP\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("ZS\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("ZP\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("ALP\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN11\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN12\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN13\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN21\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN22\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN23\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN31\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN32\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                    llParamInp.add("FN33\t" + saParamAtoms[i] + "\t" + daParams[iCounter]);
                    iCounter++;
                }
            }
        }

        // copy the data
        String[] saInput = new String[llInput.size()];
        Iterator<String> itInput = llInput.iterator();

        for(int i = 0; i < saInput.length; i++){
            saInput[i] = itInput.next();
        }

        String[] saParamInput = new String[llParamInp.size()];
        Iterator<String> itParamInput = llParamInp.iterator();

        for(int i = 0; i < saParamInput.length; i++){
            saParamInput[i] = itParamInput.next();
        }

        // write it out
        final String sWhereTo = sMopacFolder + System.getProperty("file.separator");
        OutputPrimitives.writeOut(sWhereTo + sMopacInput, saInput, true);
        OutputPrimitives.writeOut(sWhereTo + "params", saParamInput, true);
    }
}
