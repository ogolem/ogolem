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

import java.io.IOException;
import java.util.Locale;
import static org.ogolem.core.Constants.ANGTOBOHR;
import static org.ogolem.core.Constants.BOHRTOANG;
import static org.ogolem.core.Constants.EVTOHARTREE;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * This calls the program suite Mopac for local geometry optimizations.
 * @author Johannes Dieterich
 * @version 2015-04-21
 */
final class MopacCaller extends AbstractLocOpt {

    // the ID
    private static final long serialVersionUID = (long) 20150421;

    static enum METHOD{CUSTOM,MNDO,AM1,PM3,PM5,PM6,RM1,PM7};
    
    private final METHOD whichBackend;
    private final int iNoOfCycles;
    private final boolean bUseCOSMO;
    private final float fSolventRadius;
    private final float fDielectricConst;
    private final String mopLine;

    /**
     * The constructor for usage as a local optimizing engine.
     * @param globconf A complete configuration set.
     */
    MopacCaller(final GlobalConfig globconf, final METHOD whichMethod) throws Exception{
        super(globconf);
        if(whichMethod == METHOD.CUSTOM){
                mopLine = Input.ReadFile(globconf.outputFolder + "-mopac.aux")[0];
        } else {
            mopLine = null;
        }
        this.whichBackend = whichMethod;
        this.iNoOfCycles = globconf.maxIterLocOpt;
        this.bUseCOSMO = false;
        this.fSolventRadius = Float.NaN;
        this.fDielectricConst = Float.NaN;
    }

    /**
     * The constructor for usage as a local optimizing engine with COSMO solvation.
     * @param globconf A complete configuration set.
     * @param whichMethod Defines which semiempirical ansatz to choose.
     * @param solventRadius Radius of the solvent molecule for COSMO.
     * @param dielectricConstant Dielectric constant of the solvent.
     */
    MopacCaller(final GlobalConfig globconf, final METHOD whichMethod,
            final float solventRadius, final float dielectricConstant){
        super(globconf);
        this.mopLine = null;
        this.whichBackend = whichMethod;
        this.iNoOfCycles = globconf.maxIterLocOpt;
        this.bUseCOSMO = true;
        this.fSolventRadius = solventRadius;
        this.fDielectricConst = dielectricConstant;
    }

    private MopacCaller(final MopacCaller orig){
        super(orig);
        this.mopLine = orig.mopLine;
        this.whichBackend = orig.whichBackend;
        this.iNoOfCycles = orig.iNoOfCycles;
        this.bUseCOSMO = orig.bUseCOSMO;
        this.fSolventRadius = orig.fSolventRadius;
        this.fDielectricConst = orig.fDielectricConst;
    }

    @Override
    public MopacCaller clone(){
        return new MopacCaller(this);
    }
    
    @Override
    protected String myID(){
        return "Mopac caller";
    }

    @Override
    public String myIDandMethod(){

        String method = whichBackend.name();

        if(!Float.isNaN(fSolventRadius) && !Float.isNaN(fDielectricConst)){
            method += "COSMO solvation, solvent radius: " + fSolventRadius
                    + " dielectric constant: " + fDielectricConst;
        }

        return "MOPAC" + method;
    }
    
    @Override
    public CartesianCoordinates cartesToCartes(final long lID, CartesianCoordinates cartes,
            final boolean[][] baConstraints, boolean isConstricted, final BondInfo bonds) throws ConvergenceException{
        String sMopacInput = "mopac" + lID + ".dat";
        String sMopacBasis = "mopac" + lID;

        float[] faCharges = cartes.getAllCharges();
        short[] iaSpins = cartes.getAllSpins();
        
        /*
         * write the output aka input for the to be called program
         */
        try{
            writeMopacInput(sMopacInput, whichBackend.name(), cartes.getAllXYZCoord(), cartes.getAllAtomTypes(),
                            cartes.getTotalCharge(), cartes.getTotalSpin(), iNoOfCycles, bUseCOSMO,
                            fSolventRadius, fDielectricConst, baConstraints,
                            mopLine);
        } catch(InitIOException e){
            throw new ConvergenceException("Problem in writing geometry for mopac input (local optimization).",e);
        }

        /*
         * call mopac
         */
        try{
            Runtime rt = Runtime.getRuntime();
            String sMOPACCmd = System.getenv("OGO_MOPACCMD");
            if(sMOPACCmd == null){
                // default
                sMOPACCmd = "mopac";
            }
            Process proc = rt.exec(sMOPACCmd + " " + sMopacBasis);

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
                if(DEBUG){
                    final String[] err = errorGobbler.getData();
                    final String[] out = outputGobbler.getData();
                    for(final String s : err){
                        System.err.print(s);
                    }
                    for(final String s : out){
                        System.err.print(s);
                    }
                }
                throw new ConvergenceException("Mopac returns non-zero return value (local optimization).");
            } // mopac should(!) have completed normally...

        } catch(Exception e){
            throw new ConvergenceException("Mopac has a problem (local optimization).",e);
        }


        /*
         * read mopacs output
         */
        String sMopacOutput = sMopacInput.substring(0,sMopacInput.indexOf(".")) + ".out";
        try{
            cartes = ReadXYZMopacOutput(sMopacOutput, cartes.getNoOfAtoms(), cartes.getNoOfMolecules(), cartes.getAllAtomsPerMol(),
                    cartes.getAllAtomTypes());
        } catch(Exception e){
            e.printStackTrace(System.err);
            throw new ConvergenceException("Problem in reading the output of mopac.",e);
        }
        
        /*
         * clean up
         */
        try{
            Input.RemoveMopacFiles(sMopacInput);
        }catch(Exception e){
            System.err.println("Problem cleaning mopac files up." + e.toString());
        }

        /*
         * return it
         */
        cartes.setAllCharges(faCharges);
        cartes.setAllSpins(iaSpins);

        return cartes;
    }
    
    private static CartesianCoordinates ReadXYZMopacOutput(final String sMopacOutput,
            final int iNoOfAtoms, final int iNoOfMolecules, final int[] iaNoAtsPerMol,
            final String[] atoms)
            throws InitIOException, CastException{
        
        /*
         * get the string array
         */
        String[] saData;
        try{
            saData = InputPrimitives.readFileIn(sMopacOutput);
        } catch(Exception e){
            throw new InitIOException("Error in reading mopac's output file.",e);
        }

        /*
         * get the starting point for the geometry
         */
        CartesianCoordinates cartesians = new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaNoAtsPerMol);
        cartesians.setAllAtomTypes(atoms);
        int iGeometryStart = 0;
        for(int i = 0; i < saData.length; i++){
            // search for the last occurence of "CARTESIAN COORDINATES"
            if(saData[i].contains("CARTESIAN COORDINATES")){
                iGeometryStart = i+2;
            }
        }

        // get the actual information...
        final double[][] xyz = cartesians.getAllXYZCoord();
        for(int i= iGeometryStart; i < iNoOfAtoms+iGeometryStart; i++){
            
            final String[] tmp = saData[i].trim().split("\\s+");
            final String thisAt = atoms[i-iGeometryStart];
            if(!thisAt.equalsIgnoreCase(tmp[1])){throw new CastException("Atoms mismatch "
                    + (i-iGeometryStart) + " " + tmp[0] + " vs " + thisAt);}

            // now the actual set of coordinates
            for(int coord = 0; coord < 3; coord++){
                try {
                    xyz[coord][ i-iGeometryStart] = Double.parseDouble(tmp[coord+2]) * ANGTOBOHR;
                } catch (Exception e) {
                    throw new CastException("Failure in mopac coordinate casting.", e);
                }
            }
        }

        // get the line of the energy
        int iEnergyLine = 0;
        for(int i = 0; i < saData.length; i++){
            if(saData[i].contains("TOTAL ENERGY")){
                iEnergyLine = i;
            }
        }

        // actually get to the energy, the 26 is really mopac output specific.
        final String tmp = saData[iEnergyLine].trim().substring(26).trim();
        // get rid of the EV in the end
        final String[] tmp2 = tmp.split("\\s+");
        double energy = 0.0;
        try{
            energy = Double.parseDouble(tmp2[0]) * EVTOHARTREE;
        }catch(Exception e){
            System.err.println("Problems casting the energy of the mopac output.");
            throw new CastException(e);
        }

        // set the energy
        cartesians.setEnergy(energy);

        return cartesians;
    }

    private static void writeMopacInput(final String sMopacInput, final String sMethod, final double[][] daXYZ,
            final String[] saAtoms, final int iTotalCharge, final int iTotalSpin, final int iNoOfCycles,
            final boolean bUseCOSMO, final float fSolventRadius, final float fDielectricConst,
            final boolean[][] baConstraints, final String mopacAux)
            throws InitIOException {

        // case switching for the spins
        String sSpin;

        // the multiplicity (since we operate with full spin units, just 1s+1)
        final int iMultiplicity = 1 * iTotalSpin + 1;

        switch (iMultiplicity) {
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
                sSpin = "";
                break;
        }

        final int iNoOfAtoms = saAtoms.length;
        final String[] saOutput = new String[iNoOfAtoms + 3];
        // first line: the method and all that
        if(mopacAux == null){
            saOutput[0] = "XYZ NOLOG PRECISE GEO-OK T=100H CYCLES=" + iNoOfCycles + " " + sMethod + " charge=" + iTotalCharge + " " + sSpin;
            if(bUseCOSMO){
                saOutput[0] += " EPS=" + fDielectricConst + " RSOLV=" + fSolventRadius;
            }
        } else {
            saOutput[0] = mopacAux;
        }

        // the next line: a tag
        saOutput[1] = "CREATED BY OGOLEM";

        // third line: empty
        saOutput[2] = "";

        // the geometry in the wonderful MOPAC xyz format...
        for (int i = 3; i < iNoOfAtoms + 3; i++) {
            //System.out.println("DEBUG: Atom coords constraints " + (i-3) + "\t" + baConstraints[0][i-3] + "\t" + baConstraints[1][i-3] + "\t" + baConstraints[2][i-3]);
            final int i1 = (baConstraints[0][i-3]) ? 0:1;
            final int i2 = (baConstraints[1][i-3]) ? 0:1;
            final int i3 = (baConstraints[2][i-3]) ? 0:1;

            saOutput[i] = saAtoms[i - 3] + " " + String.format(Locale.US,"%20.6f",(daXYZ[0][i - 3] * BOHRTOANG)) + " " +  i1
                    + " " + String.format(Locale.US,"%20.7f",(daXYZ[1][i - 3] * BOHRTOANG)) + " " + i2
                    + " " + String.format(Locale.US,"%20.7f",(daXYZ[2][i - 3] * BOHRTOANG)) + " " + i3;
        }

        //write it actually
        try {
            OutputPrimitives.writeOut(sMopacInput, saOutput, false);
        } catch (IOException e) {
            throw new InitIOException(e);
        }
    }
}
