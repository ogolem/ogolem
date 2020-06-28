/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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
package org.ogolem.interfaces;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.AbstractLocOpt;
import org.ogolem.core.CartesianFullBackend;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.Constants;
import static org.ogolem.core.Constants.ANGTOBOHR;
import static org.ogolem.core.Constants.BOHRTOANG;
import org.ogolem.core.ConvergenceException;
import org.ogolem.core.FixedValues;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Gradient;
import org.ogolem.core.Input;
import org.ogolem.core.NumericalGradients;
import org.ogolem.core.StreamGobbler;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * This calls Orca as a backend for local optimization and dipole calculation.
 * @author Johannes Dieterich
 * @version 2017-03-03
 */
public final class OrcaCaller extends AbstractLocOpt implements CartesianFullBackend {

    private static final long serialVersionUID = (long) 20121101;
    private final int iWhichMethod;
    private final String sCustomInp;
    private final String orcaCmd;
    private final boolean fileOut;

    /**
     * Constructor for use as a energy/ backend.
     * @param globconf
     */
/*    public OrcaCaller(final GlobalConfig globconf, final int whichMethod, final boolean fileOut,
            final String auxFilePath){
        super(globconf, globconf.bonds);
        if(whichMethod < 0){
            this.sCustomInp = readCustomInp(auxFilePath);
        } else{
            this.sCustomInp = null;
        }

        this.iWhichMethod = whichMethod;
        
        final String tmp = System.getenv("OGO_ORCACMD");
        if(tmp == null){
            orcaCmd = "orca";
        } else {
            orcaCmd = tmp;
        }
        this.fileOut = fileOut;
    }*/
    
    /**
     * Constructor for use as a local optimization backend.
     * @param globconf
     */
    public OrcaCaller(final GlobalConfig globconf, final int whichMethod, final boolean fileOut,
            final String auxFilePath){
        super(globconf);
        if(whichMethod < 0){
            this.sCustomInp = readCustomInp(auxFilePath);
        } else{
            this.sCustomInp = null;
        }

        this.iWhichMethod = whichMethod;
        
        final String tmp = System.getenv("OGO_ORCACMD");
        if(tmp == null){
            orcaCmd = "orca";
        } else {
            orcaCmd = tmp;
        }
        this.fileOut = fileOut;
    }

    private OrcaCaller(final OrcaCaller orig){
        super(orig);
        this.iWhichMethod = orig.iWhichMethod;
        this.orcaCmd = orig.orcaCmd;
        if(orig.sCustomInp != null){
            this.sCustomInp = orig.sCustomInp;
        } else{
            this.sCustomInp = null;
        }
        this.fileOut = orig.fileOut;
    }

    @Override
    public OrcaCaller clone(){
        return new OrcaCaller(this);
    }
    
    @Override
    public String getMethodID(){
        return "Orca interface: " + iWhichMethod;
    }
    
    @Override
    public String myID(){
        return "Orca";
    }

    @Override
    public String myIDandMethod(){

        String method = "";
        switch(iWhichMethod){
            case -2:
                method = "custom, no locopt"; break;
            case -1:
                method = "custom"; break;
            case 0:
                method = "MNDO"; break;
            case 1:
                method = "AM1"; break;
            case 2:
                method = "PM3"; break;
            case 3:
                method = "B3LYP/vdz"; break;
            case 4:
                method = "RI-BP86/svp"; break;
            case 5:
                method = "RI-BP86/tzvp"; break;
            case 6:
                method = "RI-UBP86/tzvp"; break;
            case 7:
                method = "B3LYP/6-31+G**"; break;
            case 8:
                method = "B3LYP/6-31+G** w/ vdW correction"; break;
            case 40:
                method = "RI-MP2/aug-cc-pVDZ"; break;
            case 50:
                method = "CCSD(T)/tzvpp"; break;
            default:
                break;
        }

        return "Orca: " + method;
    }

    @Override
    public CartesianCoordinates cartesToCartes(long lID, CartesianCoordinates cartes,
            final boolean[][] baConstraints, final boolean isConstricted, final BondInfo bonds)
            throws IOException, ConvergenceException, CastException{

        final String orcaInput = "orca" + lID + ".inp";
        final String orcaBasis = "orca" + lID;

        /*
         * write the output aka input for the to be called program
         */
        try{
            writeInput(orcaInput,cartes.getAllXYZCoord(),cartes.getAllAtomTypes(),
                    cartes.getTotalCharge(), cartes.getTotalSpin(), baConstraints, isConstricted, true);
        } catch(IOException e){
            throw new ConvergenceException("Problem in writing geometry for orca input (local optimization).",e);
        }

        /*
         * call orca
         */
        String[] output;
        try {
            output = runOrca(orcaCmd, orcaInput);
        } catch(Exception e){
            e.printStackTrace(System.err);
            throw new ConvergenceException("Orca has a problem (local optimization).");
        }

        /*
         * read orcas output
         */
        
        if(fileOut){
            // read in the output file
            try{
                output = InputPrimitives.readFileIn(orcaBasis + ".out");
            } catch(Exception e){
                System.err.println("ERROR: Couldn't read in mandatory output file " + orcaBasis + ".out");
                throw new ConvergenceException("Couldn't read in mandatory output file " + orcaBasis + ".out",e);
            }
        }

        if(this.iWhichMethod != -2){
            cartes = createCartesFromOutput(output, cartes);
        } else{
            // just parse the energy
            double energy = FixedValues.NONCONVERGEDENERGY;
            try{
                energy = parseEnergy(output);
            } catch(Exception e){
                System.err.println("WARNING: Couldn't parse the energy from orca output. " + e.toString());
            }
            cartes.setEnergy(energy);
        }

        /*
         * clean up
         */
        try{
            Input.RemoveOrcaFiles(orcaBasis);
        }catch(Exception e){
            System.err.println("Problem cleaning orca files up." + e.toString());
        }

        /*
         * return it
         */
        
        return cartes;
    }
    

    
    @Override
    public void gradientCalculation(long lID, int iIteration, double[] xyz1D, String[] saAtomTypes,
            short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges,
            short[] spins, BondInfo bonds, Gradient grad, final boolean hasRigidEnv) {
        
        final Gradient numGrad = NumericalGradients.numericalGradient(lID, iIteration, xyz1D, saAtomTypes, atomNos,
                atsPerMol, energyparts, iNoOfAtoms, faCharges, spins, bonds, this, hasRigidEnv);
        
        grad.copyDataIn(numGrad);
    }

    @Override
    public double energyCalculation(long lID, int iIteration, double[] xyz1D, String[] saAtoms,
            short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges,
            short[] iaSpins, final BondInfo bonds, final boolean hasRigidEnv){
        
        final String orcaInput = "orca" + lID + ".inp";
        final String orcaBasis = "orca" + lID;

        /*
         * write the output aka input for the to be called program
         */
        final double[][] xyz = new double[3][iNoOfAtoms];
        System.arraycopy(xyz1D, 0, xyz[0], 0, iNoOfAtoms);
        System.arraycopy(xyz1D, iNoOfAtoms, xyz[0], 0, iNoOfAtoms);
        System.arraycopy(xyz1D, 2*iNoOfAtoms, xyz[0], 0, iNoOfAtoms);
        int totalSpin = 0;
        float totalCharge = 0.0f;
        for (int i = 0; i < iNoOfAtoms; i++){
            totalSpin += iaSpins[i];
            totalCharge += faCharges[i];
        }
        
        try{
            writeInput(orcaInput,xyz,saAtoms,(int)totalCharge, totalSpin, null, false, true);
        } catch(IOException e){
            System.err.println("WARNING: Problem in writing geometry for orca input (local optimization)." + e.toString());
            return FixedValues.NONCONVERGEDENERGY;
        }

        /*
         * call orca
         */
        String[] output;
        try {
            output = runOrca(orcaCmd, orcaInput);
        } catch(Exception e){
            e.printStackTrace(System.err);
            System.err.println("Orca has a problem (local optimization).");
            return FixedValues.NONCONVERGEDENERGY;
        }

        /*
         * read orcas output
         */
        
        if(fileOut){
            // read in the output file
            try{
                output = InputPrimitives.readFileIn(orcaBasis + ".out");
            } catch(Exception e){
                System.err.println("ERROR: Couldn't read in mandatory output file " + orcaBasis + ".out " + e.toString());
                return FixedValues.NONCONVERGEDENERGY;
            }
        }

        // just parse the energy
        double energy = FixedValues.NONCONVERGEDENERGY;
        try{
            energy = parseEnergy(output);
        } catch(Exception e){
            System.err.println("WARNING: Couldn't parse the energy from orca output. " + e.toString());
        }

        /*
         * clean up
         */
        try{
            Input.RemoveOrcaFiles(orcaBasis);
        }catch(Exception e){
            System.err.println("Problem cleaning orca files up." + e.toString());
        }
        
        return energy;
    }
    
    private void writeInput(String sOrcaInput, double[][] daXYZ, String[] saAtoms,
            int iTotalCharge, int iTotalSpin, boolean[][] baConstraints, boolean isConstricted,
            final boolean doesLocopt)
            throws IOException{

        switch(iWhichMethod){
            case -2:
                // custom, same case
            case -1:
                // custom, glue input together
                final String[] geom = new String[saAtoms.length+2];
                geom[0] = "*xyz " + iTotalCharge + "  " + (iTotalSpin+1);
                for(int i = 0; i < saAtoms.length; i++){
                    geom[i+1] = saAtoms[i] + "    " + daXYZ[0][i]*Constants.BOHRTOANG
                            + "    " + daXYZ[1][i]*Constants.BOHRTOANG
                            + "    " + daXYZ[2][i]*Constants.BOHRTOANG;
                }
                geom[geom.length-1] = "*";
                final String[] out = glue(sCustomInp, geom);
                try{
                    OutputPrimitives.writeOut(sOrcaInput, out, false);
                } catch(Exception e){
                    throw new IOException(e);
                }
                break;
            case 0:
                String sBasis = "";
                try{
                    WriteOrcaInput(sOrcaInput, "mndo", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                }catch(IOException e){
                    throw e;
                };
                break;
            case 1:
                sBasis = "";
                try{
                    WriteOrcaInput(sOrcaInput, "am1", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                }catch(IOException e){
                    throw e;
                };
                break;
            case 2:
                sBasis = "";
                try{
                    WriteOrcaInput(sOrcaInput, "pm3", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                }catch(IOException e){
                    throw e;
                };
                break;
            case 3:
                sBasis = "vdz";
                try{
                    WriteOrcaInput(sOrcaInput, "b3lyp", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                }catch(IOException e){
                    throw e;
                };
                break;
            case 4:
                sBasis = "svp";
                try{
                    WriteOrcaInput(sOrcaInput, "bp86 ri sv/j", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                } catch(IOException e){
                    throw e;
                };
                break;
            case 5:
                sBasis = "tzvp";
                try{
                    WriteOrcaInput(sOrcaInput, "bp86 ri tzv/j", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                } catch(IOException e){
                    throw e;
                };
                break;
            case 6:
                sBasis = "tzvp";
                try{
                    WriteOrcaInput(sOrcaInput, "uks bp86 ri tzv/j", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                } catch(IOException e){
                    throw e;
                };
                break;
            case 7:
                sBasis = "6-31+G**";
                try{
                    WriteOrcaInput(sOrcaInput, "b3lyp", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                } catch(IOException e){
                    throw e;
                };
                break;
            case 8:
                sBasis = "6-31+G** VDW";
                try{
                    WriteOrcaInput(sOrcaInput, "b3lyp", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                } catch(IOException e){
                    throw e;
                };
                break;
            case 40:
                sBasis = "aug-cc-pVDZ";
                try{
                    WriteOrcaInput(sOrcaInput, "ri-mp2", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                } catch(IOException e){
                    throw e;
                };
                break;
            case 50:
                sBasis = "tzvpp";
                try{
                    WriteOrcaInput(sOrcaInput, "ccsd(t)", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                } catch(IOException e){
                    throw e;
                };
                break;
            default:
                sBasis = "";
                try{
                    WriteOrcaInput(sOrcaInput, "am1", sBasis, daXYZ, saAtoms, iTotalCharge, iTotalSpin, isConstricted, baConstraints, doesLocopt);
                }catch(IOException e){
                    throw e;
                };
                break;
        }
    }

    private double parseEnergy(String[] out) throws Exception{
        
        // loop an search for last occurence of FINAL SINGLE POINT
        int line = -1;
        for(int i = out.length-1; i >= 0; i++){
            if(out[i].contains("FINAL SINGLE POINT ENERGY")){
                line = i;
                break;
            }
        }

        // parse energy
        final String temp = out[line].substring(26).trim();
        double energy = FixedValues.NONCONVERGEDENERGY;
        try{
            energy = Double.parseDouble(temp);
        } catch(Exception e){
            throw new CastException(e);
        }

        return energy;
    }

    private CartesianCoordinates createCartesFromOutput(final String[] output, final CartesianCoordinates ref) throws CastException{
        
        final CartesianCoordinates cartes = new CartesianCoordinates(ref);
        /*
         * check where the coordinates start from:
         * Last occurence of "CARTESIAN COORDINATES (ANGSTROEM)"
         * check where the energy is:
         * Last occurence of "FINAL SINGLE POINT ENERGY"
         */
        int cartesStart = 0;
        int energyIndex = 0;
        for(int i = 0; i < output.length; i++){
            if(output[i].contains("CARTESIAN COORDINATES (ANGSTROEM)")){
                cartesStart = i+2;
            } else if(output[i].contains("FINAL SINGLE POINT ENERGY")){
                energyIndex = i;
            }
        }
        
        // get all coordinates
        final double[][] xyz = cartes.getAllXYZCoord();
        for(int i = 0; i < ref.getNoOfAtoms(); i++){
            // we can ignore the atom types, as these should be set already
            final String[] sa = output[i+cartesStart].trim().split("\\s+");
            try{
                xyz[0][i] = Double.parseDouble(sa[1])*ANGTOBOHR;
                xyz[1][i] = Double.parseDouble(sa[2])*ANGTOBOHR;
                xyz[2][i] = Double.parseDouble(sa[3])*ANGTOBOHR;
            } catch(Exception e){
                throw new CastException(e);
            }            
        }

        // we know where the energy is from before
        final String temp = output[energyIndex].substring(26).trim();

        double energy = FixedValues.NONCONVERGEDENERGY;
        try{
            energy = Double.parseDouble(temp);
        } catch(Exception e){
            throw new CastException(e);
        }

        cartes.setEnergy(energy);
        return cartes;
    }

    private static String readCustomInp(final String file){

        String[] sa;
        try{
            sa = InputPrimitives.readFileIn(file);
        } catch(Exception e){
            System.err.println("ERROR: Couldn't read in reference input from file " + file + ". This will fail! " + e.getMessage());
            return "";
        }

        String stot = "";
        for(final String s : sa){
            stot += s + System.getProperty("line.separator");
        }

        return stot;
    }

    private static String[] glue(final String ref, final String[] geom){

        final String[] res = new String[geom.length+1];
        res[0] = ref;
        System.arraycopy(geom, 0, res, 1, geom.length);

        return res;
    }
    
    private static void WriteOrcaInput(final String sOrcaInput, final String sMethod, final String sBasis,
            final double[][] xyz, final String[] saAtoms, final int totCharge,
            final int totSpin, final boolean isConstricted,
            final boolean[][]constr, final boolean doesLocopt) throws IOException {

        final int noOfAtoms = saAtoms.length;
        final List<String> inpData = new ArrayList<>(noOfAtoms+5);

        if(sMethod.equalsIgnoreCase("CCSD(T)") && doesLocopt){
            inpData.add("!" + sMethod + " " + sBasis + " TightSCF NumGrad OPT");
        } else if(!sMethod.equalsIgnoreCase("CCSD(T)") && doesLocopt) {
            inpData.add("!" + sMethod + " " + sBasis + " TightSCF OPT");
        } else{
            inpData.add("!" + sMethod + " " + sBasis + " TIghtSCF");
        }

        // we need more than the standard 256 MB memory. For almost everything.
        inpData.add("%MaxCore 1000");

        // now constraints
        if(isConstricted){
            inpData.add("%geom Constraints");
            for(int i = 0; i < noOfAtoms; i++){
                if(constr[0][i] && constr[1][i] && constr[2][i]){
                    inpData.add("{ C " + i + "C }");
                } else if(!constr[0][i] && !constr[1][i] && !constr[2][i]){
                    // no constraint
                    continue;
                } else {
                    // mixed, unknown case
                    System.err.println("WARNING: Mixed constraint case detected for Orca. We are not capable of doing this." +
                            " Continuing w/o constraint.");
                    continue;
                }
            }
            inpData.add("end");
            inpData.add("end");
        }

        // orca wants the multiplicity to be there (since we operate with full spin units, just 1s+1)
        final int multiplicity = totSpin + 1;

        inpData.add("*xyz " + totCharge + " " + multiplicity);
        for (int i = 3; i < noOfAtoms + 3; i++) {
            inpData.add("\t" + saAtoms[i - 3] + "\t" + xyz[0][i - 3] * BOHRTOANG
                    + "\t" + xyz[1][i - 3] * BOHRTOANG + "\t" + xyz[2][i - 3] * BOHRTOANG);
        }
        inpData.add("*");

        // copy to String[]
        final String[] orcaInputData = new String[inpData.size()];
        for(int i = 0; i < inpData.size(); i++ ){
            orcaInputData[i] = inpData.get(i);
        }

        try {
            OutputPrimitives.writeOut(sOrcaInput, orcaInputData, false);
        } catch (Exception e) {
            throw new IOException(e);
        }
    }
    
    
    private static String[] runOrca(final String orcaCmd, final String orcaInput) throws Exception {
        
        final Runtime rt = Runtime.getRuntime();
        final Process proc = rt.exec(orcaCmd + " " + orcaInput);

        // any error message?
        final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

        // any output?
        final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

        // kick them off
        errorGobbler.start();
        outputGobbler.start();

        // any error???
        if (proc.waitFor() != 0) {
            throw new Exception("Orca returns non-zero return value (local optimization).");
        } else {
            // orca should(!) have completed normally...
            return outputGobbler.getData();
        }
    }
    
    @Override
    public CartesianFullBackend getBackend(){
        return new OrcaCaller(this);
    }
}
