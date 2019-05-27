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
package org.ogolem.core;

import java.io.File;
import java.io.IOException;
import java.util.Locale;
import static org.ogolem.core.Constants.ANGTOBOHR;
import static org.ogolem.core.Constants.BOHRTOANG;
import static org.ogolem.core.Constants.KCALTOHARTREE;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.InquiryPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * This calls the amber program suite (namely sander) to locally optimize a geometry.
 * @author Johannes Dieterich
 * @version 2014-09-06
 */
final class AmberCaller extends AbstractLocOpt {

    // the ID
    private static final long serialVersionUID = (long) 20110113;
    private final int noOfIterations;
    private final boolean implicitWater;
    private final String name;

    AmberCaller(final GlobalConfig globconf, final boolean implicitWater){
        super(globconf);
        this.noOfIterations = globconf.maxIterLocOpt;
        this.name = globconf.outputFolder;
        this.implicitWater = implicitWater;
    }

    private AmberCaller(final AmberCaller orig){
        super(orig);
        this.noOfIterations = orig.noOfIterations;
        this.name = orig.name;
        this.implicitWater = orig.implicitWater;
    }

    @Override
    public AmberCaller clone(){
        return new AmberCaller(this);
    }
    
    @Override
    public String myID(){
        if(implicitWater) {return "amber: Implicit water";}
        else {return "amber";}
    }

    @Override
    public String myIDandMethod(){
        if(implicitWater) {return "amber: Implicit water";}
        else {return "amber";}
    }

    @Override
    public Molecule localOptimization(Molecule mStartMolecule) {
        
        // this caller is not capable of molecular optimization due to the amber binary ff format
        System.out.println("WARNING: AMBER caller can't locally optimize a molecule. Returning non-optimized molecule.");
        mStartMolecule.setEnergy(FixedValues.NONCONVERGEDENERGY);

        return mStartMolecule;
    }

    @Override
    public CartesianCoordinates cartesToCartes(long lID, CartesianCoordinates cartes,
           boolean[][] constraints, boolean isConstricted, final BondInfo bonds) throws ConvergenceException{
        
        final String amberFolder = "amber" + lID;
        final float[] charges = cartes.getAllCharges();
        final short[] spins = cartes.getAllSpins();

        /*
         * write the output aka input for the to be called program
         */
        try{
            writeAmberInput(amberFolder, cartes.getAllXYZCoord(), noOfIterations, implicitWater);
        } catch( InitIOException | IOException e){
            throw new ConvergenceException("Problem in writing geometry for amber input (local optimization).",e);
        }

        final String paramFile = System.getProperty("user.dir") + System.getProperty("file.separator") + name +".prmtop";
        final String copiedParamFile = System.getProperty("user.dir") + System.getProperty("file.separator") +
                amberFolder + System.getProperty("file.separator") + name +".prmtop";
        try{
            OutputPrimitives.createLink(paramFile, copiedParamFile);
        } catch(IOException e){
            throw new ConvergenceException("Problem in copying prmtop for amber input (local optimization).",e);
        }


        /*
         * call amber
         */
        boolean success = false;
        try{
            success = callAmber(amberFolder);
        } catch(ConvergenceException e){
            throw e;
        }

        if(!success){
            // try again, it might be a "sander bomb" in amber
            try{
                success = callAmber(amberFolder);
            } catch(ConvergenceException e){
                throw e;
            }
        }

        if(!success){
            // this time we really throw an exception
            throw new ConvergenceException("Amber has a problem in the local optimization.");
        }


        /*
         * read ambers output
         */
        final CartesianCoordinates optCartes;
        try{
            optCartes = readXYZAmberOutput(amberFolder, cartes.getAllAtomTypes(), cartes.getNoOfAtoms(), cartes.getNoOfMolecules(), cartes.getAllAtomsPerMol());
        } catch(Exception e){
            throw new ConvergenceException("Problem in reading the output of amber.",e);
        }

        optCartes.setAllCharges(charges);
        optCartes.setAllSpins(spins);

        /*
         * clean up
         */
        try{
            ManipulationPrimitives.remove(amberFolder);
        }catch(Exception e){
            System.err.println("Problem cleaning amber files up." + e.toString());
        }

        /*
         * return it
         */
        
        return optCartes;
    }

    private boolean callAmber(final String amberFolder) throws ConvergenceException{
        
	boolean success = false;
        
        /*
         * call amber
         */
        Process proc = null;
        try{
            final Runtime rt = Runtime.getRuntime();
            final String localPRMTOP = name +".prmtop";
            final String command = "sander -O -i optimization.in -o optimization.out -c coordinates.inpcrd -p " + localPRMTOP + " -r optimization.rst";
            final File dir = new File(amberFolder);
            proc = rt.exec(command,null,dir);

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
            success = (proc.waitFor() == 0);
        } catch(Exception e){
            e.printStackTrace(System.err);
            throw new ConvergenceException("Amber has a problem (local optimization).",e);
        } finally{
            if(proc != null){
                proc.destroy();
            }
        }

        return success;
    }

    private static CartesianCoordinates readXYZAmberOutput(final String amberFolder,
            final String[] atomNames, final int noOfAtoms, final int noOfMolecules,
            final int[] atomsPerMol)
            throws ConvergenceException{
        
        final CartesianCoordinates cartes = new CartesianCoordinates(noOfAtoms, noOfMolecules, atomsPerMol);
        cartes.setAllAtomTypes(atomNames);

        // read the cartesian coordinates in
        String[] coordData;
        try{
            coordData = InputPrimitives.readFileIn(amberFolder + System.getProperty("file.separator") + "optimization.rst");
        }catch(Exception e){
            throw new ConvergenceException("Couldn't read in ambers coordinate output.",e);
        }

        int atom = 0;
        final double[][] xyz = cartes.getAllXYZCoord();
        for(int i = 2; i < coordData.length; i++){
            // split line up
            final String[] saTemp = coordData[i].trim().split("\\s+");
         
            // the first atom
            for(int j = 0;j < 3; j++){
                try{
                    xyz[j][atom] = Double.parseDouble(saTemp[j]) * ANGTOBOHR;
                } catch(Exception e){
                    throw new ConvergenceException("Couldn't cast the coordinate",e);
                }
            }
            atom++;

            if(atom == noOfAtoms){
                break;
            }

            // the second atom
            for(int j = 0;j < 3; j++){
                try{
                    xyz[j][atom] = Double.parseDouble(saTemp[j+3]) * ANGTOBOHR;
                } catch(Exception e){
                    throw new ConvergenceException("Couldn't cast the coordinate",e);
                }
            }
            atom++;

            if(atom == noOfAtoms){
                break;
            }
        }


        // read the energy in
        String[] output;
        try{
            output = InputPrimitives.readFileIn(amberFolder+System.getProperty("file.separator") + "optimization.out");
        }catch(Exception e){
            throw new ConvergenceException("Couldn't read ambers output",e);
        }

        int energyLine = -1;
        for(int i = 0; i < output.length; i++){
            if(output[i].contains("FINAL RESULTS")){
                energyLine = i+5;
                break;
            }
        }

        final String tmp = output[energyLine].trim().split("\\s+")[1];
        try{
            final double energy = Double.parseDouble(tmp) * KCALTOHARTREE;
            cartes.setEnergy(energy);
        } catch(Exception e){
            throw new ConvergenceException("Couldn't cast energy.",e);
        }

        return cartes;
    }
    
    private static void writeAmberInput(final String folder, final double[][] xyz,
            final int noLocIter, final boolean implicitWater) throws InitIOException,IOException {

        // create the folder, if it doesn't exist
        if (!InquiryPrimitives.doesFileExist(folder)) {
            try {
                OutputPrimitives.createAFolder(folder, false);
            } catch (IOException e) {
                throw e;
            }
        }

        // write the .in file
        final String[] input = new String[9];
        input[0] = "Automatically generated by OGOLEM.";
        input[1] = "&cntrl";
        input[2] = "imin = 1,";
        input[3] = "maxcyc = " + noLocIter + ",";
        input[4] = "ncyc = " + noLocIter / 2 + ",";
        input[5] = "ntb = 0,";
        if (!implicitWater) {
            input[6] = "igb = 0,";
        } else {
            input[6] = "igb = 2,";
        }
        input[7] = "cut = 9999";
        input[8] = "/";


        final String inputFile = folder + System.getProperty("file.separator") + "optimization.in";
        try {
            OutputPrimitives.writeOut(inputFile, input, false);
        } catch (IOException e) {
            throw new InitIOException(e);
        }

        writeAmberCoordinates(xyz, folder);
    }
    
    static void writeAmberCoordinates(final double[][] xyz, final String folder)
            throws IOException {

        // create and write the coordinates file
        final double dNoOfLines = xyz[0].length / 2.;
        final int iNoOfLines = (int) Math.round(dNoOfLines) + 2;
        final String[] coords = new String[iNoOfLines];
        coords[0] = "";
        coords[1] = String.format(Locale.US, "%6d", xyz[0].length);
        
        boolean first = true;
        int line = 2;
        for (int atom = 0; atom < xyz[0].length; atom++) {
            if(first){
                coords[line]  = String.format(Locale.US, "%12.7f",xyz[0][atom]*BOHRTOANG);
                coords[line] += String.format(Locale.US, "%12.7f",xyz[1][atom]*BOHRTOANG);
                coords[line] += String.format(Locale.US, "%12.7f",xyz[2][atom]*BOHRTOANG);
                first = false;
            } else{
                coords[line] += String.format(Locale.US, "%12.7f",xyz[0][atom]*BOHRTOANG);
                coords[line] += String.format(Locale.US, "%12.7f",xyz[1][atom]*BOHRTOANG);
                coords[line] += String.format(Locale.US, "%12.7f",xyz[2][atom]*BOHRTOANG);
                line++;
                first = true;
            }
        }
        
        final String coordFile = folder + System.getProperty("file.separator") + "coordinates.inpcrd";
        OutputPrimitives.writeOut(coordFile, coords, false);
    }
}
