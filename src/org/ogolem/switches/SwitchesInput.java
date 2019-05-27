/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.InitIOException;
import org.ogolem.core.SerialException;
import org.ogolem.generic.genericpool.GenericPool;
// get physical constants at hand
import static org.ogolem.core.Constants.*;

/**
 * Handles all input for the switches part.
 * @author Johannes Dieterich
 * @version 2013-09-23
 */
public final class SwitchesInput{

    public static SwitchesConfig readConfigIn(final String sInputFile)
            throws CastException, InitIOException, IOException{

        final SwitchesConfig config = new SwitchesConfig();

        if(!sInputFile.endsWith(".sogo")){
            throw new InitIOException("Input file does not have the correct ending (.sogo).");
        }

        // read it in
        final String[] saFileContent = readFileIn(sInputFile);

        // check the first line for the correct start
        if(!saFileContent[0].equalsIgnoreCase("###OGOLEMSWITCHES###")){
            throw new InitIOException("Input file does not start with ###OGOLEMSWITCHES###.");
        }

        // now check for keywords
        for(String sCurrentLine:saFileContent){

            sCurrentLine = sCurrentLine.trim();

            if(sCurrentLine.equalsIgnoreCase("###OGOLEMSWITCHES###")){
                continue;
            } else if(sCurrentLine.startsWith("#") || sCurrentLine.startsWith("//")){
                // sort comments out
                continue;
            } else if(sCurrentLine.startsWith("PoolSize=")){
                final String sTemp = sCurrentLine.substring(9).trim();
                int iPoolSize;
                try{
                    iPoolSize = Integer.parseInt(sTemp);
                    SwitchesConfig.iPoolSize = iPoolSize;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast integer choice for PoolSize, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("NoOfGlobIters=")){
                final String sTemp = sCurrentLine.substring(14).trim();
                int iNoOfIters;
                try{
                    iNoOfIters = Integer.parseInt(sTemp);
                    SwitchesConfig.iNoOfGlobIters = iNoOfIters;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast integer choice for NoOfGlobIters, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("StepsToSerialization=")){
                final String sTemp = sCurrentLine.substring(21).trim();
                int iStepsToSerial;
                try{
                    iStepsToSerial = Integer.parseInt(sTemp);
                    SwitchesConfig.iSwitchesToSerial = iStepsToSerial;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast integer choice for StepsToSerialization, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("MaxIterLocOpt=")){
                final String sTemp = sCurrentLine.substring(14).trim();
                try{
                    final int iNoOfIters = Integer.parseInt(sTemp);
                    config.iMaxIterLocOpt = iNoOfIters;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast integer choice for MaxIterLocOpt, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("PhonyDebug=")){
                final String sTemp = sCurrentLine.substring(11).trim();
                boolean bDebug;
                try{
                    bDebug = Boolean.parseBoolean(sTemp);
                    config.bDebug = bDebug;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast boolean choice for PhonyDebug, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("MoreMutation=")){
                final String sTemp = sCurrentLine.substring(13).trim();
                boolean bMoreMut;
                try{
                    bMoreMut = Boolean.parseBoolean(sTemp);
                    config.bMoreMutation = bMoreMut;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast boolean choice for MoreMutation, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("AnyDiversityCheck=")){
                final String sTemp = sCurrentLine.substring(18).trim();
                boolean bAnyDivCheck;
                try{
                    bAnyDivCheck = Boolean.parseBoolean(sTemp);
                    SwitchesConfig.bAnyDivCheck = bAnyDivCheck;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast boolean choice for AnyDiversityCheck, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("BlowBondDetect=")){
                final String sTemp = sCurrentLine.substring(15).trim();
                double dBlow;
                try{
                    dBlow = Double.parseDouble(sTemp);
                    config.dBlowBondsFac = dBlow;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for BlowBondDetect, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("FitnessDiversity=")){
                final String sTemp = sCurrentLine.substring(17).trim();
                double dFitDiv;
                try{
                    dFitDiv = Double.parseDouble(sTemp);
                    SwitchesConfig.dFitnessDiversity = dFitDiv;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for FitnessDiversity, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("SigmoidFactorA=")){
                final String sTemp = sCurrentLine.substring(15).trim();
                double dSigFacA;
                try{
                    dSigFacA = Double.parseDouble(sTemp);
                    config.dSigFactorA = dSigFacA;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for SigmoidFactorA, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("SigmoidFactorB=")){
                final String sTemp = sCurrentLine.substring(15).trim();
                double dSigFacB;
                try{
                    dSigFacB = Double.parseDouble(sTemp);
                    config.dSigFactorB = dSigFacB;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for SigmoidFactorB, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("DihedralThreshCis=")){
                final String sTemp = sCurrentLine.substring(18).trim();
                double dThreshDihedral;
                try{
                    dThreshDihedral = Double.parseDouble(sTemp);
                    config.dThreshDihedralCis = dThreshDihedral * Math.PI / 180.0;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for DihedralThreshCis, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("DihedralThreshTrans=")){
                final String sTemp = sCurrentLine.substring(20).trim();
                double dThreshDihedral;
                try{
                    dThreshDihedral = Double.parseDouble(sTemp);
                    config.dThreshDihedralTrans = dThreshDihedral * Math.PI / 180.0;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for DihedralThreshTrans, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("OptimalCisDihedral=")){
                final String sTemp = sCurrentLine.substring(19).trim();
                try{
                    final double dOptCis = Double.parseDouble(sTemp) * Math.PI / 180.0;
                    config.dOptimalCisVal = dOptCis;
                }catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for OptimalCisDihedral, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("OptimalTransDihedral=")){
                final String sTemp = sCurrentLine.substring(21).trim();
                try{
                    final double dOptTrans = Double.parseDouble(sTemp) * Math.PI / 180.0;
                    config.dOptimalTransVal = dOptTrans;
                }catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for OptimalTransDihedral, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("ExactEnergyFactor=")){
                final String sTemp = sCurrentLine.substring(18).trim();
                double dExactEFac;
                try{
                    dExactEFac = Double.parseDouble(sTemp);
                    config.dExactEnergyFactor = dExactEFac;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for ExactEnergyFactor, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("ExactEnergyCis=")){
                final String sTemp = sCurrentLine.substring(15).trim();
                double dExactECis;
                try{
                    dExactECis = recognizeEnergyUnit(sTemp);
                    config.dExactEnergyCis = dExactECis;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for ExactEnergyCis, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("ExactEnergyTrans=")){
                final String sTemp = sCurrentLine.substring(17).trim();
                double dExactETrans;
                try{
                    dExactETrans = recognizeEnergyUnit(sTemp);
                    config.dExactEnergyTrans = dExactETrans;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for ExactEnergyTrans, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("MinimumEnergyDiff=")){
                final String sTemp = sCurrentLine.substring(18).trim();
                try{
                    final double dMinE = recognizeEnergyUnit(sTemp);
                    config.dMinEnergyDiff = dMinE;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for MinimumEnergyDiff, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("LocOptAlgo=")){
                final String sTemp = sCurrentLine.substring(11).trim();
                final int iWhichLocOpt = SwitchesConfig.parseIntLocOptFromString(sTemp);
                config.iWhichLocAlgo = iWhichLocOpt;
            } else if(sCurrentLine.startsWith("ExcitorMethod=")){
                final String sTemp = sCurrentLine.substring(14).trim();
                final int iWhichExcitor = SwitchesConfig.parseIntExcitorFromString(sTemp);
                config.iWhichExcitor = iWhichExcitor;
            } else if(sCurrentLine.startsWith("NumberOfActiveOcc=")){
                final String sTemp = sCurrentLine.substring(18).trim();
                try{
                    final int iNoOfAcOcc = Integer.parseInt(sTemp);
                    config.iActiveOcc = iNoOfAcOcc;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for NumberOfActiveOcc, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("NumberOfNonActiveOcc=")){
                final String sTemp = sCurrentLine.substring(21).trim();
                try{
                    final int iNoOfNonAcOcc = Integer.parseInt(sTemp);
                    config.iActiveNonOcc = iNoOfNonAcOcc;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for NumberOfNonActiveOcc, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("NumberOfOccPi=")){
                final String sTemp = sCurrentLine.substring(14).trim();
                try{
                    final int iNoOfOccPi = Integer.parseInt(sTemp);
                    config.iOccPi = iNoOfOccPi;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for NumberOfOccPi, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("NumberOfNonOccPi=")){
                final String sTemp = sCurrentLine.substring(17).trim();
                try{
                    final int iNoOfNonOccPi = Integer.parseInt(sTemp);
                    config.iNonOccPi = iNoOfNonOccPi;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for NumberOfNonOccPi, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("GUGADefOfOrbs=")){
                final String sTemp = sCurrentLine.substring(14).trim();
                try{
                    final int iDefOfOrbs = Integer.parseInt(sTemp);
                    // check whether the choice makes sense
                    if(iDefOfOrbs != 0 && iDefOfOrbs != 1 && iDefOfOrbs != -1 && iDefOfOrbs != -2
                            && iDefOfOrbs != -3 && iDefOfOrbs != -4){
                        // this is for sure not valid
                        System.err.println("ERROR: Your choice for GUGADefOfOrbs makes no sense. We use the default.");
                    } else {
                        // valid choice
                        config.iOrbDef = iDefOfOrbs;
                    }
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for GUGADefOfOrbs, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("NumberOfRefOcc=")){
                final String sTemp = sCurrentLine.substring(15).trim();
                try{
                    final int iNoOfRefOcc = Integer.parseInt(sTemp);
                    config.iNoOfRefOcc = iNoOfRefOcc;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for NumberOfRefOcc, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("DefOfRefOcc=")){
                final String sTemp = sCurrentLine.substring(12).trim();
                try{
                    final int iDefOfRefOcc = Integer.parseInt(sTemp);
                    config.iDefOfRefOcc = iDefOfRefOcc;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for DefOfRefOcc, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("MaxExcitationLevel=")){
                final String sTemp = sCurrentLine.substring(19).trim();
                try{
                    final int iMaxExcitLevel = Integer.parseInt(sTemp);
                    config.iMaxExcitLevel = iMaxExcitLevel;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for MaxExcitationLevel, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("NumberOfLowestStates=")){
                final String sTemp = sCurrentLine.substring(21).trim();
                try{
                    final int iNoOfLowestStates = Integer.parseInt(sTemp);
                    config.iNoOfLowestStates = iNoOfLowestStates;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for NumberOfLowestStates, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("WantedCIState=")){
                final String sTemp = sCurrentLine.substring(14).trim();
                try{
                    final int iWantedCIState = Integer.parseInt(sTemp);
                    config.iWantedCiState = iWantedCIState;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't cast double choice for WantedCIState, using default. " +
                            e.toString());
                }
            } else if(sCurrentLine.startsWith("Backbone=")){
                // read the backend in
                final String sFileBase = sCurrentLine.substring(9).trim();
                final String sFileCis = sFileBase + "cis.xyz";
                final String sFileTrans = sFileBase + "trans.xyz";

                System.out.println("INFO: We assume this to be the base of the backbone names, so files "
                        + sFileCis + " and " + sFileTrans + " should be present.");

                final CartesianCoordinates backCartesCis = org.ogolem.core.Input.readCartesFromFile(sFileCis);
                final CartesianCoordinates backCartesTrans = org.ogolem.core.Input.readCartesFromFile(sFileTrans);

                /*
                 * now we still need to figure the connects out. those are the dummies denoted with
                 * XX.
                 */
                final String[] saAtomsCis = backCartesCis.getAllAtomTypes();

                final LinkedList<Integer> llConnectsCis = new LinkedList<>();
                for(int at = 0; at < saAtomsCis.length; at++){
                    final String sCurrAt= saAtomsCis[at];
                    if(sCurrAt.equalsIgnoreCase("XX")){
                        llConnectsCis.add(at);
                    }
                }

                final String[] saAtomsTrans = backCartesTrans.getAllAtomTypes();

                final LinkedList<Integer> llConnectsTrans = new LinkedList<>();
                for(int at = 0; at < saAtomsTrans.length; at++){
                    final String sCurrAt= saAtomsTrans[at];
                    if(sCurrAt.equalsIgnoreCase("XX")){
                        llConnectsTrans.add(at);
                    }
                }


                // check whether we even have any connects
                if(llConnectsCis.size() == 0 || llConnectsTrans.size() == 0){
                    throw new CastException(" WE DO NOT HAVE ANY CONNECTS!");
                }

                // check whether both are same length
                if(llConnectsCis.size() != llConnectsTrans.size()){
                    throw new CastException("WE DO NOT HAVE THE SAME NUMBER OF CONNECTS!");
                }

                final int[] iaConnectsCis = new int[llConnectsCis.size()];
                final int[] iaConnectsTrans = new int[llConnectsCis.size()];

                for(int i = 0; i < iaConnectsCis.length; i++){
                    iaConnectsCis[i] = llConnectsCis.get(i);
                    iaConnectsTrans[i] = llConnectsTrans.get(i);
                }

                // create the bonding info, again: sufficient for one of the two
                final BondInfo baBackBondsCis = LittleHelpers.bondingInfo(backCartesCis, config.dBlowBondsFac);
                final BondInfo baBackBondsTrans = LittleHelpers.bondingInfo(backCartesTrans, config.dBlowBondsFac);
     
                // create the actual backbone object
                final Backbone back = new Backbone(backCartesCis, backCartesTrans,
                        iaConnectsCis, iaConnectsTrans, baBackBondsCis,
                        baBackBondsTrans);

                // put it in the config
                SwitchesConfig.backbone = back;
            } else if(sCurrentLine.startsWith("AtomsForDihedralCis=")){
                final int[] iaAtsForDih = new int[4];
                final String sTemp = sCurrentLine.substring(20).trim();

                // now try to separate them using ;
                final String[] saAtsForDih = sTemp.split(";");
                try{
                    for(int i = 0 ; i < 4; i++){
                        // parse the atoms
                        iaAtsForDih[i] = Integer.parseInt(saAtsForDih[i]);
                    }
                    config.iaWhichAtsForDihedralCis = iaAtsForDih;
                } catch(Exception e){
                    throw new CastException("Couldn't parse which atoms for dihedral check of cis to use.",e);
                }
            } else if(sCurrentLine.startsWith("AtomsForDihedralTrans=")){
                final int[] iaAtsForDih = new int[4];
                final String sTemp = sCurrentLine.substring(22).trim();

                // now try to separate them using ;
                final String[] saAtsForDih = sTemp.split(";");
                try{
                    for(int i = 0 ; i < 4; i++){
                        // parse the atoms
                        iaAtsForDih[i] = Integer.parseInt(saAtsForDih[i]);
                    }
                    config.iaWhichAtsForDihedralTrans = iaAtsForDih;
                } catch(Exception e){
                    throw new CastException("Couldn't parse which atoms for dihedral check of trans to use.",e);
                }
            } else if(sCurrentLine.startsWith("Colors=")){
                final String sTemp = sCurrentLine.substring(7).trim();

                // first get all possible colors
                final ColorPalette palette = ColorPalette.getReference();
                final String[] saPossCols = palette.getPossibleColors();

                final ArrayList<String> alPossCols = new ArrayList<>();
                alPossCols.addAll(Arrays.asList(saPossCols));

                // the sTemp is a seperated list of input options
                final String[] saInputColors = sTemp.split(":");

                for(int i = 0; i < saInputColors.length; i++){
                    for(int j = 0; j < saPossCols.length; j++){
                        if(saInputColors[i].equalsIgnoreCase( "no" + saPossCols[j])){
                            // remove it from the arraylist
                            alPossCols.remove(saPossCols[j]);
                            break;
                        }

                        if(j == saPossCols.length -1){
                            System.err.println("WARNING: No such color " + saInputColors[i] + ".");
                        }
                    }
                }

                // now initialize the color palette
                final Iterator<String> itPossCols = alPossCols.iterator();
                final ArrayList<Color> alColors = new ArrayList<>(alPossCols.size());

                while(itPossCols.hasNext()){
                    alColors.add(new Color(itPossCols.next()));
                }

                palette.initializeColors(alColors);

            } else {
                System.err.println("WARNING: No idea what to do with " + sCurrentLine + " ignoring it.");
            }
        }

        /*
         * check for whether or not the ColorPalette has been initialized
         */
        final ColorPalette palette = ColorPalette.getReference();
        if(palette.getRandomColor() == null){
            System.out.println("INFO: ColorPalette was not explicitly initialized, doing that now.");
            final String[] saColors = palette.getPossibleColors();

            final ArrayList<Color> alColors = new ArrayList<>();
            for(int i = 0; i < saColors.length; i++){
                final Color color = new Color(saColors[i]);
                alColors.add(color);
            }

            palette.initializeColors(alColors);
        }

        /*
         * check whether or not we have a backbone
         */
        if(SwitchesConfig.backbone == null){
            throw new CastException("No backbone specified.");
        }

        /*
         * check whether or not the atoms for the dihedral measurement are
         * populated
         */
        if(config.iaWhichAtsForDihedralCis == null || config.iaWhichAtsForDihedralTrans == null){
            throw new CastException("No atoms for dihedral measurement specified.");
        }

        return config;
    }

    /**
     * Checks if a file exists. Is synchronized out of safety since that might be helpful to erase
     * a potential race-condition in the reaction (e.g. file creation) to this request.
     * @param sFileName
     * @return if this file exists
     */
    synchronized static boolean isFilePresents(final String sFileName){
        final File f = new File(sFileName);

        return f.exists();
    }

    static CartesianCoordinates readXYZTinker(String sCartesianFile, int iNoOfAtoms, int iNoOfMolecules, int[] iaAtomsPerMol) throws InitIOException, CastException{
        CartesianCoordinates cartes = new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaAtomsPerMol);
        String[] saTinkerXYZOut;
        try{
            saTinkerXYZOut = readFileIn(sCartesianFile);
        } catch(Exception e){
            throw new InitIOException(e);
        }

        // the chopping and casting games
        String sTemp;
        String sTemp2;
        int iIndex;
        double[] daCoords = new double[3];
        for(int i = 0; i < iNoOfAtoms; i++){
            sTemp = saTinkerXYZOut[i + 1];
            sTemp = sTemp.trim();
            //chop the first column off
            sTemp = sTemp.substring(sTemp.indexOf(" ")).trim();

            sTemp = sTemp.substring(sTemp.indexOf(" "));
            sTemp = sTemp.trim();
            // the coordinates
            for(int j = 0; j < 3; j++){
                iIndex = sTemp.indexOf(" ");
                sTemp2 = sTemp.substring(0, iIndex);
                try{
                    daCoords[j] = Double.parseDouble(sTemp2) * ANGTOBOHR;
                }catch(Exception e){
                    throw new CastException(e);
                }
                sTemp = sTemp.substring(iIndex);
                sTemp = sTemp.trim();
            }
            cartes.setXYZCoordinatesOfAtom(daCoords, i);
        }

        return cartes;
    }

    static void removeTinkerFiles(final String sTinkerBasis) throws InitIOException{
        // just remove the correct files
        removeFile(sTinkerBasis + ".xyz");
        removeFile(sTinkerBasis + ".xyz_2");
        removeFile(sTinkerBasis + ".key");
    }

    static void removeFolder(String sFolder) throws InitIOException{
        final File f = new File(sFolder);
        final File[] files = f.listFiles();
        for (int j = 0; j < files.length; j++) {
            if(files[j].isDirectory()) {
                // recursive call... :-)
                removeFolder(files[j].getAbsolutePath());
            }
            else {
                final boolean bReturn = files[j].delete();
                if(!bReturn){
                    System.err.println("ERROR: File couldn't be deleted. " + files[j]);
                }
            }
        }
        final boolean bDeleteFolderOK = f.delete();
        if (bDeleteFolderOK) {
            //Nothing, fine.
        } else {
            // no success. for whatever reason.
            throw new InitIOException("No success in deleting directory " + sFolder + " .");
        }
    }

    static String[] readFileIn(final String sInputPath) throws IOException{
        return org.ogolem.io.InputPrimitives.readFileIn(sInputPath);
    }

    static double readExcMopacOutput(final String sMopacOutput) throws Exception{

        /*
         * get the string array
         */
        String[] saData;
        try{
            saData = readFileIn(sMopacOutput);
        } catch(Exception e){
            throw new InitIOException("Error in reading mopac's output file.",e);
        }

        int iLineNumbStates = -1;
        for(int i = 0; i < saData.length; i++){
            final String s = saData[i].trim();
            if(s.startsWith("STATE") && !s.contentEquals("STATE VECTORS")){
                iLineNumbStates = i+3;
                break;
            } else if(s.startsWith("THE GEOMETRY MAY NOT BE COMPLETELY OPTIMIZED")){
                throw new InitIOException("MOPAC reports the geometry to be non-optimized.");
            }
        }

        if(iLineNumbStates == -1){
            // no such line
            throw new InitIOException("No state line in MOPAC output.");
        }

        // which state is ours?
        boolean bCorrStateFound = false;
        boolean bCorrExcFound = false;
        String sStateMult = "SINGLET";
        double dRelEnergy = 0.0;
        double dS0S1Energy = FixedValues.UNEXCITABLEENERGY;
        for(int i = iLineNumbStates; i < saData.length; i++){
            final String s = saData[i].trim();
            final String sWord1 = s.substring(0,s.indexOf(" "));
            boolean bFirstStateFound = false;
            if(!bCorrExcFound){
                if(sWord1.endsWith("+")){
                    bFirstStateFound = true;
                } else if(sWord1.endsWith("*")){
                    bFirstStateFound = true;
                }
            }
            if(bFirstStateFound){
                final String[] tokens = s.split("\\s+");
                // parse its relative energy
                try{
                    dRelEnergy = Double.parseDouble(tokens[2]);
                } catch(Exception e){
                    throw e;
                }

                // parse its multiplicity
                sStateMult = tokens[4];
                
                bCorrStateFound = true;
                
            } else if(bCorrStateFound && !bCorrExcFound){
                // split into parts
                final String[] tokens = s.split("\\s+");

                // check whether it has the correct multiplicity
                if(tokens[4].equalsIgnoreCase(sStateMult)){

                    // parse the energy
                    try{
                        final double dState2 = Double.parseDouble(tokens[2]);
                        dS0S1Energy = dState2 - dRelEnergy;
                        dS0S1Energy *= EVTOHARTREE;
                    } catch(Exception e){
                        throw e;
                    }

                    bCorrExcFound = true;
                    // we can break the loop
                    break;
                }
            }
        }

        if(!bCorrStateFound || !bCorrExcFound){
            throw new InitIOException("Couldn't parse states and excitation energies.");
        }

        return dS0S1Energy;
    }

    static CartesianCoordinates readXYZMopacOutput(String sMopacOutput, int iNoOfAtoms, int iNoOfMolecules, int[] iaNoAtsPerMol) throws InitIOException, CastException{
        /*
         * get the string array
         */
        String[] saData;
        try{
            saData = readFileIn(sMopacOutput);
        } catch(Exception e){
            throw new InitIOException("Error in reading mopac's output file.",e);
        }

        /*
         * get the starting point for the geometry
         */
        CartesianCoordinates cartesians = new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaNoAtsPerMol);
        int iGeometryStart = 0;
        for(int i = 0; i < saData.length; i++){
            // search for the last occurence of "CARTESIAN COORDINATES"
            if(saData[i].contains("CARTESIAN COORDINATES")){
                iGeometryStart = i+4;
            } else{
                // nothing
            }
        }

        // get the actual information...
        String sTemp;
        String sTemp2;
        String sTempAtom;
        double[] daTempCoord = new double[3];
        for(int i= iGeometryStart; i < iNoOfAtoms+iGeometryStart; i++){
            sTemp = saData[i];
            sTemp = sTemp.trim();
            // get rid of the atom id since it is redundant and therefore useless
            sTemp = sTemp.substring(sTemp.indexOf(" "));
            sTemp = sTemp.trim();
            sTempAtom = sTemp.substring(0, sTemp.indexOf(" "));
            cartesians.setAtom(sTempAtom, i-iGeometryStart);
            sTemp = sTemp.substring(sTemp.indexOf(" "));
            sTemp = sTemp.trim();

            // now the actual set of coordinates

            // x
            sTemp2 = sTemp.substring(0, sTemp.indexOf(" "));
            try {
                daTempCoord[0] = Double.parseDouble(sTemp2) * ANGTOBOHR;
            } catch (Exception e) {
                throw new CastException("Failure in mopac coordinate casting.", e);
            }
            sTemp = sTemp.substring(sTemp.indexOf(" "));

            // y
            sTemp = sTemp.trim();
            sTemp2 = sTemp.substring(0, sTemp.indexOf(" "));
            try {
                daTempCoord[1] = Double.parseDouble(sTemp2) * ANGTOBOHR;
            } catch (Exception e) {
                throw new CastException("Failure in mopac coordinate casting.", e);
            }
            sTemp = sTemp.substring(sTemp.indexOf(" "));

            // z
            sTemp2 = sTemp.trim();
            try {
                daTempCoord[2] = Double.parseDouble(sTemp2) * ANGTOBOHR;
            } catch (Exception e) {
                throw new CastException("Failure in mopac coordinate casting.", e);
            }

            cartesians.setXYZCoordinatesOfAtom(daTempCoord, i-iGeometryStart);
        }

        // get the line of the energy
        int iEnergyLine = 0;
        for(int i = 0; i < saData.length; i++){
            if(saData[i].contains("TOTAL ENERGY")){
                iEnergyLine = i;
            } else{
                // nothing
            }
        }

        // actually get to the energy
        sTemp = saData[iEnergyLine].trim();
        // the 26 is really mopac output specific.
        sTemp = sTemp.substring(26);
        sTemp = sTemp.trim();
        // get rid of the EV in the end
        sTemp = sTemp.substring(0,sTemp.indexOf(" "));
        double dEnergy = 0.0;
        try{
            dEnergy = Double.parseDouble(sTemp) * EVTOHARTREE;
        }catch(Exception e){
            System.err.println("Problems casting the energy of the mopac output.");
            throw new CastException(e);
        }

        // set the energy
        cartesians.setEnergy(dEnergy);

        return cartesians;
    }

    static void removeFile(final String sToFilePath) throws InitIOException{
        final File f = new File(sToFilePath);
        final boolean bSuccess = f.delete();
        if(bSuccess == false){
           throw new InitIOException("Couldn't remove file.");
        }
    }

    private static double recognizeEnergyUnit(final String sEnergy) throws Exception{

        if(sEnergy.endsWith("nm")){
            final String sTemp = sEnergy.substring(0, sEnergy.length() - 3).trim();
            final double dEnergy = 1.0/Double.parseDouble(sTemp) * NMTOHARTREE;
            return dEnergy;
        } else if(sEnergy.endsWith("hartree")){
            final String sTemp = sEnergy.substring(0, sEnergy.length() - 8).trim();
            final double dEnergy = Double.parseDouble(sTemp);
            return dEnergy;
        } else if(sEnergy.endsWith("kj/mol") ){
            final String sTemp = sEnergy.substring(0, sEnergy.length() - 7).trim();
            final double dEnergy = Double.parseDouble(sTemp) * KJTOHARTREE;
            return dEnergy;
        } else if(sEnergy.endsWith("kcal/mol")){
            final String sTemp = sEnergy.substring(0, sEnergy.length() - 9).trim();
            final double dEnergy = Double.parseDouble(sTemp) * KCALTOHARTREE;
            return dEnergy;
        } else{
            System.err.println("WARNING: No known extension for the energy found. Assuming hartree now.");
            final double dEnergy = Double.parseDouble(sEnergy);
            return dEnergy;
        }
    }

    @SuppressWarnings("unchecked")
    public static GenericPool<Color,Switch> readSwitchPool(final String sPath) throws Exception{
        
        Object oObj = null;

        try {
            oObj = ReadBinInput(sPath);
        } catch (Exception e) {
            throw e;
        }

        GenericPool<Color,Switch> pool;
        try {
            pool = (GenericPool<Color,Switch>) oObj;
        } catch (Exception e) {
            throw e;
        }

        return pool;
    }

    private static Object ReadBinInput(final String sToBinPath) throws InitIOException, SerialException {

        Object obj = null;
        ObjectInputStream objectStream = null;

        try {
            objectStream = new ObjectInputStream(new FileInputStream(sToBinPath));
            obj = objectStream.readObject();
            objectStream.close();
        } catch (IOException e) {
            throw new InitIOException("Error occured during binary reading!", e);
        } catch (ClassNotFoundException e) {
            throw new SerialException("Couldn't cast to object. This IS strange!", e);
        } finally {
            if (objectStream != null) {
                try {
                    objectStream.close();
                } catch (IOException e) {
                    throw new SerialException("Couldn't close file.", e);
                }
            }
        }

        return obj;
    }
}
