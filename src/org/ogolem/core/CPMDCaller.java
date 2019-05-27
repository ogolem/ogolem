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
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Calls CPMD for local optimizations.
 * @author Johannes Dieterich
 * @version 2015-04-21
 */
final class CPMDCaller extends AbstractLocOpt {

    // the ID
    private static final long serialVersionUID = (long) 20110225;
    
    static enum FUNCTIONAL{SONLY,LDA,BONLY,BP,BLYP,XLYP,GGA,PBE,PBES,REVPBE,HCTH,OPTX,OLYP,TPSS,PBE0,B1LYP,B3LYP,X3LYP,CUSTOM};

    // which funtional to be used out of the broad variety provided
    private final FUNCTIONAL whichFunctional;

    // we use this as the convergence threshhold for the locopt
    private final double dThreshGeom;

    // theshhold for the orbital convergence
    private final double dThreshOrbitals;

    private final double dCutOff;

    /*
     * The auxilliary info: first entry in each row is the atom type, then the
     * pseudopotential file and then the maximum angular momentum.
     */
    private final String[][] saAuxInfos;

    private final String[] customBlocks;

    private double[] cellSize;


    /**
     * Constructor for use as a local optimization backend.
     * @param globconf
     */
    CPMDCaller(GlobalConfig globconf, FUNCTIONAL whichFunctional){
        super(globconf);
        this.whichFunctional = whichFunctional;
        this.dThreshGeom = globconf.threshLocOptCoord;
        this.dThreshOrbitals = globconf.threshLocOptGradient;

        final String sWhichAuxFile = globconf.outputFolder + "-cpmd.aux";
        String[] saAuxInput;

        try {
            saAuxInput = Input.ReadFile(sWhichAuxFile);
        } catch (Exception e) {
            System.err.println("Error in reading the mandatory aux-file!" + e.toString());
            saAuxInput = null;
        }

        if (saAuxInput == null || !saAuxInput[0].equalsIgnoreCase("###OGOLEMAUX###")) {
            saAuxInfos = null;
            dCutOff = 70;
        } else {
            // read the cutoff in
            double dCut;
            try {
                dCut = Double.parseDouble(saAuxInput[1]);
            } catch (Exception e) {
                System.err.println("ERROR: Couldn't parse the cutoff in the aux info. Using 70 Rydberg now. " + e.toString());
                dCut = 70;
            }
            dCutOff = dCut;

            // try to parse the cell size
            try {
                cellSize = new double[3];
                final String[] sa = saAuxInput[2].trim().split("\\;");
                cellSize[0] = Double.parseDouble(sa[0].trim());
                cellSize[1] = Double.parseDouble(sa[1].trim());
                cellSize[2] = Double.parseDouble(sa[2].trim());
            } catch (Exception e) {
                System.err.println("WARNING: Cell size could not be read, estimating always. " + e.toString());
                cellSize = null;
            }

            // read in, which basis sets and potentials are to be taken
            String[][] saTempInfo = new String[saAuxInput.length - 3][3];
            try {
                for (int i = 3; i < saAuxInput.length; i++) {

                    if(saAuxInput[i].trim().isEmpty()) continue;

                    // always three strings: for which element which basisset and which lmax
                    final String[] sa = saAuxInput[i].trim().split("\\s+");
                    if(sa.length != 3){
                        throw new Exception("ERROR: Aux information not valid! Words in line: " + sa.length);
                    }
                    saTempInfo[i - 3] = sa;
                }
            } catch (Exception e) {
                // we just catch anything that comes along
                System.err.println("ERROR: Can't read in aux file properly. " + e.toString());
                saTempInfo = null;
            }
            saAuxInfos = saTempInfo;
        }
        if(whichFunctional == FUNCTIONAL.CUSTOM){
            
            final String inpFile = globconf.outputFolder + "-cpmdinp.aux";

            String[] custInput;
            try {
                custInput = Input.ReadFile(inpFile);
            } catch (Exception e) {
                System.err.println("Error in reading the custom input file!" + e.toString());
                throw new RuntimeException("Couldn't read custom input file for CPMD.",e);
            }

            // custom method: split the input
            String s = "";
            for(final String t : custInput){
                s+= t + System.getProperty("line.separator");
            }
            customBlocks = s.split("OGOCOORDS");
            // cutoff doesn't matter
        } else{
            customBlocks = null;
        }
    }

    private CPMDCaller(final CPMDCaller orig){
        super(orig);
        this.dThreshGeom = orig.dThreshGeom;
        this.dThreshOrbitals = orig.dThreshOrbitals;
        this.dCutOff = orig.dCutOff;
        this.whichFunctional = orig.whichFunctional;
        if(orig.cellSize != null) this.cellSize = orig.cellSize.clone();
        else this.cellSize = null;
        if(orig.customBlocks != null) this.customBlocks = orig.customBlocks.clone();
        else this.customBlocks = null;
        // again, a shallow copy should be enough
        this.saAuxInfos = orig.saAuxInfos.clone();
    }

    @Override
    public CPMDCaller clone(){
        return new CPMDCaller(this);
    }

    @Override
    public String myIDandMethod(){
        return "CPMD: " + whichFunctional;
    }
    
    @Override
    public String myID(){
        return "CPMD: " + whichFunctional;
    }

    @Override
    public CartesianCoordinates cartesToCartes(final long lID,
            final CartesianCoordinates cartes, final boolean[][] constraints,
            final boolean isConstricted, final BondInfo bonds) throws InitIOException,
            IOException, ConvergenceException, CastException{

        // the foldername
        final String sCPMDFolder = "cpmd" + lID;

        final float[] faCharges = cartes.getAllCharges();
        final short[] iaSpins = cartes.getAllSpins();

        // create the folder, potentially throws an exception
        OutputPrimitives.createAFolder(sCPMDFolder);

        // copy the .psp files
        for(int i = 0; i < saAuxInfos.length; i++){

            final String sPSPFile = saAuxInfos[i][1];
            final String sNewPSP = sCPMDFolder + System.getProperty("file.separator") + sPSPFile;

            // potentially throws an exception
            OutputPrimitives.createLink(sPSPFile, sNewPSP);
        }

        // move the COM to zero
        cartes.moveCoordsToCOM();

        double[] daCell;
        if(cellSize == null || cartes.getNoOfMolecules() == 1){
            // determine the cell size
            daCell = determineCellSize(cartes.getAllXYZCoord(), cartes.getTotalCharge());
        } else{
            daCell = cellSize;
        }

        // write the input, potentially throws an exception
        final List<List<Integer>> llKey = writeCPMDInput(sCPMDFolder,
                dThreshGeom, dThreshOrbitals, whichFunctional.name(), cartes.getTotalCharge(),
                cartes.getTotalSpin(), cartes.getAllAtomTypes(),
                cartes.getAllXYZCoord(), saAuxInfos, daCell, dCutOff, customBlocks);

        // try to run cpmd
        Process proc = null;
        String[] saOutput;
        try{
            Runtime rt = Runtime.getRuntime();
            String sCommand = "cpmd.x run.in";
            String[] saEnvp = new String[0];
            File dir = new File(sCPMDFolder);
            proc = rt.exec(sCommand,saEnvp,dir);

            // any error message?
            StreamGobbler errorGobbler = new
                StreamGobbler(proc.getErrorStream(), "ERROR");

            // any output?
            StreamGobbler outputGobbler = new
                StreamGobbler(proc.getInputStream(), "OUTPUT");

            // kick them off
            errorGobbler.start();
            outputGobbler.start();


            // any error???
            int iExitValue = proc.waitFor();
            if(iExitValue != 0){
                throw new ConvergenceException("CPMD returns non-zero return value (local optimization): " + iExitValue + ".");
            } else{
                // cpmd should(!) have completed normally...
                saOutput = outputGobbler.getData();
            }

        } catch(Exception e){
            System.err.println(e.toString());
            throw new ConvergenceException("CPMD has a problem (local optimization).",e);
        } finally{
            if(proc != null){
                proc.destroy();
            }
        }

        // get the optimized geometry
        CartesianCoordinates optCartes = new CartesianCoordinates(cartes);
        try{
            optCartes = transformCPMDOutput(saOutput, cartes.getNoOfAtoms(),
                    cartes.getNoOfMolecules(), cartes.getAllAtomsPerMol(),
                    sCPMDFolder, iaSpins, faCharges, llKey, cartes.getAllAtomTypes());
        } catch(Exception e){
            throw new ConvergenceException("Couldn't read CPMDs output.",e);
        }


        // delete the folder
        try{
            ManipulationPrimitives.remove(sCPMDFolder);
        } catch(Exception e){
            System.err.println("WARNING: Couldn't remove the cpmd folder. " + e.toString());
        }
        
        // we do not need to set spins or charges anymore, that has been fixed for us :-)

        return optCartes;
    }

    private static CartesianCoordinates transformCPMDOutput(final String[] saOutput,
            final int iNoOfAtoms, final int iNoOfMolecules, final int[] iaAtsPerMol,
            final String sFolderName, final short[] iaSpins, final float[] faCharges,
            final List<List<Integer>> llKey, final String[] saRefAtoms)
            throws Exception{

        // first the coordinates
        final String sWhichCoordFile = sFolderName + System.getProperty("file.separator")
                + "GEOMETRY.xyz";

        // potentially might throw exceptions
        CartesianCoordinates cartes = Input.readCPMDCartesFromFile(sWhichCoordFile, iNoOfMolecules, iaAtsPerMol, iaSpins, faCharges);

        // then the energy
        int iEnergyLine = -1;
        for(int i = saOutput.length -1; i >= 0; i--){
            final String sCurr = saOutput[i];
            if(sCurr.trim().startsWith("=              END OF GEOMETRY OPTIMIZATION")){
                iEnergyLine = i - 5;
                break;
            }
        }

        if(iEnergyLine == -1){
            // we didn't find the energy line
            throw new ConvergenceException("No end of optimization found.");
        }

        // now actually read the energy line in
        final String sEnergyLine = saOutput[iEnergyLine];
        String sEnergy = sEnergyLine.substring(46).trim();
        sEnergy = sEnergy.substring(0,sEnergy.length()-4);

        // parse the energy, potentially throws an exception
        final double dEnergy = Double.parseDouble(sEnergy);

        // set it in
        cartes.setEnergy(dEnergy);

        // actually, due to the way in which CPMD expects the input, we need to
        // reorder the cartesians
        final CartesianCoordinates orderedCartes = new CartesianCoordinates(cartes);
        orderedCartes.setAllAtomTypes(saRefAtoms);

        final double[][] daXYZ = cartes.getAllXYZCoord();
        final double[][] daOrderedXYZ = new double[3][iNoOfAtoms];

        // turn the key into a long int[]
        final int[] iaKey = new int[orderedCartes.getNoOfAtoms()];
        
        int iCounter = 0;
        for(final List<Integer> outer : llKey){
            for(final int inner : outer){
                iaKey[iCounter] = inner;
                iCounter++;
            }
        }

        // now we "just" move coordinates
        for(int i = 0; i < iaKey.length; i++){
            daOrderedXYZ[0][iaKey[i]] = daXYZ[0][i];
            daOrderedXYZ[1][iaKey[i]] = daXYZ[1][i];
            daOrderedXYZ[2][iaKey[i]] = daXYZ[2][i];
        }

        orderedCartes.setAllXYZ(daOrderedXYZ);

        return orderedCartes;
    }

    private static double[] determineCellSize(final double[][] daXYZ, final int iTotalCharge){

        final double[] daCell = new double[3];

        double dMinX = Double.POSITIVE_INFINITY;
        double dMinY = Double.POSITIVE_INFINITY;
        double dMinZ = Double.POSITIVE_INFINITY;

        double dMaxX = Double.NEGATIVE_INFINITY;
        double dMaxY = Double.NEGATIVE_INFINITY;
        double dMaxZ = Double.NEGATIVE_INFINITY;

        // figure the MIN/MAX for x out
        for(double dValue : daXYZ[0]){
            dMinX = Math.min(dValue, dMinX);
            dMaxX = Math.max(dValue, dMaxX);
        }

        // figure the MIN/MAX for y out
        for(double dValue : daXYZ[1]){
            dMinY = Math.min(dValue, dMinY);
            dMaxY = Math.max(dValue, dMaxY);
        }

        // figure the MIN/MAX for z out
        for(double dValue : daXYZ[2]){
            dMinZ = Math.min(dValue, dMinZ);
            dMaxZ = Math.max(dValue, dMaxZ);
        }

        // now compute the cell size and increment
        if(iTotalCharge == 0){
            daCell[0] = Math.abs(dMaxX - dMinX) + 15.0;
            daCell[1] = Math.abs(dMaxY - dMinY) + 15.0;
            daCell[2] = Math.abs(dMaxZ - dMinZ) + 15.0;
        } else{
            // charged species, make it bigger
            daCell[0] = Math.abs(dMaxX - dMinX) + 25.0;
            daCell[1] = Math.abs(dMaxY - dMinY) + 25.0;
            daCell[2] = Math.abs(dMaxZ - dMinZ) + 25.0;
        }

        return daCell;
    }
    
    private static List<List<Integer>> writeCPMDInput(final String sCPMDFolder,
            final double dThreshGeom, final double dThreshOrbitals,
            final String sFunctional, final int iTotCharge, final int iTotSpin,
            final String[] saAtoms, final double[][] daXYZ, final String[][] saAux,
            final double[] daCell, final double dCutOff, final String[] customBlocks)
            throws IOException{

        final List<String> llInput = new LinkedList<>();
        final List<List<Integer>> llKey = new LinkedList<>();

        final boolean custom = sFunctional.equalsIgnoreCase("CUSTOM");

        if(!custom){
            // ask for a geometry optimization
            llInput.add("&CPMD");
            llInput.add("OPTIMIZE GEOMETRY XYZ");
            llInput.add("CONVERGENCE ORBITALS");
            llInput.add("" + dThreshOrbitals);
            llInput.add("CONVERGENCE GEOMETRY");
            llInput.add("" + dThreshGeom);

            if (iTotSpin != 0) {
                // CPMD manual says that not all functionals are implemented for this!
                llInput.add("LOCAL SPIN DENSITY");
            }

            llInput.add("INITIALIZE WAVEFUNCTION RANDOM");
            llInput.add("CENTER MOLECULE ON");
            llInput.add("HESSIAN UNIT");
            llInput.add("&END");
            llInput.add("");

            // specify the system
            llInput.add("&SYSTEM");
            llInput.add("SYMMETRY");
            llInput.add("ISOLATED SYSTEM");

            // the cell size
            llInput.add("CELL");

            // compute the cell size
            final double dCellX = daCell[0];
            final double dCellY = daCell[1] / daCell[0];
            final double dCellZ = daCell[2] / daCell[0];

            final DecimalFormat dec = new DecimalFormat("0.0");

            // we ALWAYS use a rectangular cell
            llInput.add(" " + dec.format(dCellX) + " " + dec.format(dCellY) + " " + dec.format(dCellZ)
                    + " 0 0 0");

            llInput.add("CHARGE");
            llInput.add(" " + iTotCharge);

            if (iTotSpin != 0) {
                llInput.add("MULTIPLICITY");
                llInput.add(" " + (iTotSpin * 2 + 1));
            }

            llInput.add("CUTOFF");
            llInput.add(" " + dCutOff);
            llInput.add("&END");
            llInput.add("");

            // specify the functional to be used
            llInput.add("&DFT");
            llInput.add("FUNCTIONAL " + sFunctional);
            llInput.add("&END");

            // specify the system
            llInput.add("&ATOMS");

            for (int i = 0; i < saAux.length; i++) {
                llInput.add("*" + saAux[i][1]);
                llInput.add("LMAX=" + saAux[i][2]);

                // now the atoms
                final List<Integer> llCurrentAtoms = new LinkedList<>();

                for (int j = 0; j < saAtoms.length; j++) {
                    final String sAtom = saAtoms[j];
                    if (sAtom.equalsIgnoreCase(saAux[i][0])) {
                        llCurrentAtoms.add(j);
                    }
                }

                llInput.add("" + llCurrentAtoms.size());

                // finally the coordinates
                for (int j = 0; j < llCurrentAtoms.size(); j++) {
                    llInput.add("  " + (daXYZ[0][llCurrentAtoms.get(j)]) + "  "
                            + (daXYZ[1][llCurrentAtoms.get(j)]) + "  "
                            + (daXYZ[2][llCurrentAtoms.get(j)]));
                }
                llKey.add(llCurrentAtoms);
            }
            llInput.add("&END");
        } else{
            // custom input
            llInput.add(customBlocks[0]);
            for (int i = 0; i < saAux.length; i++) {
                llInput.add("*" + saAux[i][1]);
                llInput.add("LMAX=" + saAux[i][2]);

                // now the atoms
                final List<Integer> llCurrentAtoms = new LinkedList<>();

                for (int j = 0; j < saAtoms.length; j++) {
                    final String sAtom = saAtoms[j];
                    if (sAtom.equalsIgnoreCase(saAux[i][0])) {
                        llCurrentAtoms.add(j);
                    }
                }

                llInput.add("" + llCurrentAtoms.size());

                // finally the coordinates
                for (int j = 0; j < llCurrentAtoms.size(); j++) {
                    llInput.add("  " + (daXYZ[0][llCurrentAtoms.get(j)]) + "  "
                            + (daXYZ[1][llCurrentAtoms.get(j)]) + "  "
                            + (daXYZ[2][llCurrentAtoms.get(j)]));
                }
                llKey.add(llCurrentAtoms);
            }
            llInput.add(customBlocks[1]);
        }
    
        // write it to run.in in the specified folder
        final String sInputFile = sCPMDFolder + System.getProperty("file.separator")
                + "run.in";
        final String[] saInputData = new String[llInput.size()];
        for(int i = 0; i < llInput.size(); i++){
            saInputData[i] = llInput.get(i);
        }

        // this might throw an exception
        OutputPrimitives.writeOut(sInputFile, saInputData, false);

        // return the key
        return llKey;
    }
}
