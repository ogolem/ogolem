/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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
package org.ogolem.core;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.List;
import static org.ogolem.core.Constants.BOHRTOANG;
import org.ogolem.io.InquiryPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Calls CP2K for local optimization.
 * @author Johannes Dieterich
 * @version 2015-04-21
 */
final class CP2KCaller extends AbstractLocOpt{

    // the ID
    private static final long serialVersionUID = (long) 20140421;

    static enum FUNCTIONAL{B3LYP,BLYP,BP,HCTH120,NONE,NO_SHORTCUT,OLYP,PADE,PBE,PBE0,TPSS};
    
    private static final int iCellIncrement = 10;

    private final FUNCTIONAL whichFunctional;

    private final int iMaxIter;

    private final String[][] saAuxInfos;

    /**
     * Constructor for use as a local optimization backend.
     * @param globconf
     */
    CP2KCaller(GlobalConfig globconf, FUNCTIONAL whichFunctional){
        super(globconf);
        this.whichFunctional = whichFunctional;
        this.iMaxIter = globconf.maxIterLocOpt;
        
        final String sWhichAuxFile = globconf.outputFolder + "-cp2k.aux";
        String[] saAuxInput;

        try {
            saAuxInput = Input.ReadFile(sWhichAuxFile);
        } catch (Exception e) {
            System.err.println("Error in reading the mandatory aux-file!" + e.toString());
            saAuxInput = null;
        }
        if (saAuxInput == null || !saAuxInput[0].equalsIgnoreCase("###OGOLEMAUX###")) {
            this.saAuxInfos = null;
        } else {
            // read in, which basis sets and potentials are to be taken
            String[][] saTempInfo = new String[saAuxInput.length-1][3];
            try{
                for(int i = 1; i < saAuxInput.length; i++){
                    // always three strings: for which element which basisset and which potential
                    int iIndex = saAuxInput[i].indexOf(" ");
                    saTempInfo[i-1][0] = saAuxInput[i].substring(0,iIndex).trim();
                    String sTemp = saAuxInput[i].substring(iIndex).trim();
                    iIndex = sTemp.indexOf(" ");
                    saTempInfo[i-1][1] = sTemp.substring(0,iIndex).trim();
                    saTempInfo[i-1][2] = sTemp.substring(iIndex).trim();
                }
            } catch(Exception e){
                // we just catch anything that comes along
                System.err.println("ERROR: Can't read in aux file properly. " + e.toString());
                saTempInfo = null;
            }
            saAuxInfos = saTempInfo;
        }
    }

    private CP2KCaller(final CP2KCaller orig){
        super(orig);
        // shallow copy should be enough
        this.whichFunctional = orig.whichFunctional;
        // again, a shallow copy should be enough
        this.saAuxInfos = orig.saAuxInfos.clone();
        this.iMaxIter = orig.iMaxIter;
    }

    @Override
    public CP2KCaller clone(){
        return new CP2KCaller(this);
    }

    @Override
    public String myIDandMethod(){
        return "CP2K: " + whichFunctional;
    }
    
    @Override
    public String myID(){
        return "CP2K: " + whichFunctional;
    }

    @Override
    public CartesianCoordinates cartesToCartes(final long lID,
            final CartesianCoordinates cartes, final boolean[][] constraints,
            final boolean isConstricted, final BondInfo bonds) throws Exception{

        final String sCP2KFolder = "cp2k" + lID;
        final String sCP2KInput = sCP2KFolder + System.getProperty("file.separator") +"cp2k" + lID + ".inp";

        final float[] faCharges = cartes.getAllCharges();
        final short[] iaSpins = cartes.getAllSpins();

        /*
         * Determine the cell size we need
         */

        cartes.moveCoordsToCOM();

        final double[] daMaxDims = new double[3];
        final float[] faMaxDims = new float[3];
        final double[][] daXYZ = cartes.getAllXYZCoord();
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < daXYZ[0].length; j++){
                daMaxDims[i] = Math.max(Math.abs(daXYZ[i][j]), daMaxDims[i]);
            }
            daMaxDims[i] *= 2.0;
            daMaxDims[i] += iCellIncrement;
            // make a float out of it
            faMaxDims[i] = (float) daMaxDims[i];
        }

        /*
         * write the output aka input for the to be called program
         */
        try{
            writeCP2KInput(sCP2KFolder, sCP2KInput, whichFunctional.name(), iMaxIter, cartes.getAllAtomTypes(),
                    daXYZ, faMaxDims, saAuxInfos, cartes.getTotalCharge(), cartes.getTotalSpin());
        } catch(Exception e){
            throw new ConvergenceException("Problem in writing geometry for cp2k input (local optimization).",e);
        }

        /*
         * copy the basis set and potential file
         */
        final String sBasisFileOrig = "BASIS_SET";
        final String sPotentFileOrig = "POTENTIAL";
        final String sBasisFileCopy =  sCP2KFolder + System.getProperty("file.separator") +sBasisFileOrig;
        final String sPotentFileCopy = sCP2KFolder + System.getProperty("file.separator") + sPotentFileOrig;

        OutputPrimitives.createLink(sBasisFileOrig,sBasisFileCopy);
        OutputPrimitives.createLink(sPotentFileOrig,sPotentFileCopy);

        /*
         * call cp2k
         */
        //String[] saCP2KStream;
        Process proc = null;
        try {
            final Runtime rt = Runtime.getRuntime();
            final String sLocalInp = sCP2KFolder + ".inp";
            final String sCommand = "cp2k.sopt " + sLocalInp;
            final File dir = new File(sCP2KFolder);

            proc = rt.exec(sCommand,null,dir);

            // any error message?
            final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

            // any output?
            final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

            // kick them off
            errorGobbler.start();
            outputGobbler.start();

            // any error???
            final int iExitValue = proc.waitFor();
            if (iExitValue != 0) {
                throw new ConvergenceException("CP2K returns non-zero return value.");
            } else {
                // cp2k should(!) have completed normally...
                //saCP2KStream = outputGobbler.getData();
            }

        } catch (Exception e) {
            throw new ConvergenceException("CP2K has a problem (local optimization).", e);
        } finally{
            if(proc != null){
                proc.destroy();
            }
        }

        /*if(saCP2KStream.length == 0){
            System.err.println("CP2K stream does not contain any info.");
            throw new ConvergenceException();
        }*/

        /*
         * read the output
         */

        final int iNoOfAts = cartes.getNoOfAtoms();
        final String sCartesFile = sCP2KFolder + System.getProperty("file.separator") + "OGOLEM-pos-1.xyz";

        final String[] saCartCont = Input.ReadFile(sCartesFile);

        int iStartLine = -1;
        for(int i = saCartCont.length-1; i > -1; i--){
            if(saCartCont[i].trim().equalsIgnoreCase(""  + iNoOfAts)){
                // we found the start of our proper xyz coordinates
                iStartLine = i;
            }
        }

        /*
         * now read the cartesian coordinates in
         */
        final CartesianCoordinates endCartes = new CartesianCoordinates(cartes);

        if(iStartLine == -1){
            throw new ConvergenceException("Can't handle CP2K output");
        }

        final double[][] daNewXYZ = new double[3][iNoOfAts];
        double dEnergy = FixedValues.NONCONVERGEDENERGY;

        // read the energy
        final int iEnergyIndex = saCartCont[iStartLine+1].lastIndexOf("=");
        final String sEnergyLine = saCartCont[iStartLine+1].substring(iEnergyIndex+1).trim();
        try{
            // CP2K is so nice to provide the energy in hartree :-)
            dEnergy = Double.parseDouble(sEnergyLine);
        } catch(Exception e){
            throw new ConvergenceException("Couldn't parse energy in output.",e);
        }
        endCartes.setEnergy(dEnergy);

        for(int i = iStartLine+2; i < iStartLine+2+iNoOfAts; i++){
            // the cartesian coordinates
            String sCurr = saCartCont[i].trim();
            sCurr = sCurr.substring(sCurr.indexOf(" ")).trim();

            // now we get to the first real coordinate
            int iIndex = sCurr.indexOf(" ");
            String sFirst = sCurr.substring(0, iIndex).trim();
            sCurr = sCurr.substring(iIndex).trim();
            
            iIndex = sCurr.indexOf(" ");
            String sSecond = sCurr.substring(0, iIndex).trim();
            String sThird = sCurr.substring(iIndex);

            // now we try to parse it
            try{
                daNewXYZ[0][i-iStartLine-2] = Double.parseDouble(sFirst)*Constants.ANGTOBOHR;
                daNewXYZ[1][i-iStartLine-2] = Double.parseDouble(sSecond)*Constants.ANGTOBOHR;
                daNewXYZ[2][i-iStartLine-2] = Double.parseDouble(sThird)*Constants.ANGTOBOHR;
            } catch(Exception e){
                throw new ConvergenceException("Couldn't parse CP2Ks coordinates",e);
            }
        }

        // put it into the cartesian
        endCartes.setAllXYZ(daNewXYZ);

        endCartes.setAllCharges(faCharges);
        endCartes.setAllSpins(iaSpins);

        /*
         * delete the folder
         */
        try{
            ManipulationPrimitives.remove(sCP2KFolder);
        }catch(Exception e){
            System.err.println("Problem cleaning cp2k files up." + e.toString());
        }

        
        /*
         * return the results
         */

        return endCartes;
    }
    
    private static void writeCP2KInput(final String sFolderName, final String sInputPath, final String sFunctional,
            final int iMaxIters, final String[] saAtomTypes, final double[][] daXYZ, final float[] faCellSize,
            final String[][] saBasisAndPots, final int iTotalCharge, final int iTotalSpin)
            throws IOException {

        // check whether the folder already exists
        final boolean bFolderExists = InquiryPrimitives.doesFileExist(sFolderName);

        if(!bFolderExists){
            // then create the folder
            try {
                OutputPrimitives.createAFolder(sFolderName);
            } catch (Exception e) {
                System.err.println("Couldn't create folder. " + e.toString());
            }
        }

        final List<String> llInput = new LinkedList<>();

        // first some general things
        llInput.add("&FORCE_EVAL");
        llInput.add("METHOD Quickstep");
        llInput.add("&DFT");
        llInput.add("CHARGE " + iTotalCharge);
        final int iMultiplicity = 2 * iTotalSpin + 1;
        llInput.add("MULTIPLICITY " + iMultiplicity);
        // we do not need to adjust the grid ATM and if, then we would just
        // do that using the aux file, that is where such thing belongs
        //llInput.add("&MGRID");
        //llInput.add("CUTOFF 300);
        //llInput.add("&END MGRID");
        //Usage hint : 'NGRIDS 5' (section &MGRID) will deal more efficiently (2x speedup) with the diffuse nature of the basis.
        // TODO perhaps we should enable UKS for this using LSD keyword ? -> apparently not for dublett.

        llInput.add("&QS");
        llInput.add("&END QS");
        llInput.add("&SCF");
        llInput.add("SCF_GUESS ATOMIC");
        llInput.add("&END SCF");
        llInput.add("&XC");
        llInput.add("&XC_FUNCTIONAL " + sFunctional);
        llInput.add("&END XC_FUNCTIONAL");
        llInput.add("&END XC");
        llInput.add("&END DFT");
        llInput.add("&SUBSYS");
        llInput.add("&CELL");
        llInput.add("ABC " + faCellSize[0] * (float) BOHRTOANG + " " + faCellSize[1] * (float) BOHRTOANG + " " + faCellSize[2] * (float) BOHRTOANG);
        llInput.add("&END CELL");
        llInput.add("&COORD");

        // put the coordinates in
        final DecimalFormat decForm = new DecimalFormat("0.000000");
        for (int i = 0; i < saAtomTypes.length; i++) {
            llInput.add("" + saAtomTypes[i] + "   " + decForm.format(daXYZ[0][i] * BOHRTOANG) 
                    + "   " + decForm.format(daXYZ[1][i] * BOHRTOANG)
                    + "   " + decForm.format(daXYZ[2][i] * BOHRTOANG));
        }
        llInput.add("&END COORD");

        // now the bloody basis sets and potentials
        for (final String[] saBasisAndPot : saBasisAndPots) {
            llInput.add("&KIND " + saBasisAndPot[0]);
            llInput.add("BASIS_SET " + saBasisAndPot[1]);
            llInput.add("POTENTIAL " + saBasisAndPot[2]);
            llInput.add("&END KIND");
        }

        llInput.add("&END SUBSYS");
        llInput.add("&END FORCE_EVAL");

        llInput.add("&GLOBAL");
        llInput.add("");
        llInput.add("PROJECT OGOLEM");
        llInput.add("RUN_TYPE GEO_OPT");
        llInput.add("PRINT_LEVEL LOW");
        llInput.add("&END GLOBAL");
        llInput.add("&MOTION");
        llInput.add("&GEO_OPT");
        llInput.add("MAX_ITER " + iMaxIters);
        llInput.add("OPTIMIZER BFGS");
        llInput.add("&END GEO_OPT");
        llInput.add("&END MOTION");

        // transform the linked list
        final String[] saInput = new String[llInput.size()];
        for (int i = 0; i < llInput.size(); i++) {
            saInput[i] = llInput.get(i);
        }

        // write the input
        final String sInputFile = sFolderName + System.getProperty("file.separator") + sFolderName + ".inp";

        // we do not catch this exception on purpose but instead throw it upwards
        OutputPrimitives.writeOut(sInputFile, saInput, false);
    }
}
