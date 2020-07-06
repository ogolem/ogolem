/**
Copyright (c) 2013-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
              2017-2020, J. M. Dieterich and B. Hartke
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

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.ogolem.core.AbstractLocOpt;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CartesianFullBackend;
import org.ogolem.core.Constants;
import static org.ogolem.core.Constants.BOHRTOANG;
import org.ogolem.core.FixedValues;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Gradient;
import org.ogolem.core.StreamGobbler;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Calls Stefan Grimme's xtb program for local geometry optimization. Has partial
 * support as a CartesianFullBackend (energy evaluation) which is why it is not
 * being exposed.
 * @author Bernd Hartke
 * @author Johannes Dieterich
 * @version 2020-07-04
 */
public class XTBCaller extends AbstractLocOpt implements CartesianFullBackend {
    
    private static final long serialVersionUID = (long) 20200704;
    
    public static enum METHOD {GFNFF, GFN0XTB, GFN1XTB, GFN2XTB};
    public static enum OPTLEVEL {CRUDE, SLOPPY, LOOSE, NORMAL, TIGHT, VERYTIGHT};
    
    private final String xtbCmd;
    private final METHOD xtbMeth;
    private final String[] xtbMethod;
    private final String[] xtbOpt;
    private final String[] xtbOtherOptions;
    private final boolean setEnvironment;
    private final String xControlFileOrig;
        
    public XTBCaller(final GlobalConfig globconf, final METHOD method, final OPTLEVEL level,
    	    final boolean setEnvironment, final String xControlFileOrig) throws Exception {
        super(globconf);
        
        final String cmd = System.getenv("OGO_XTBCMD");
        this.xtbCmd = (cmd == null) ? "xtb" : cmd;
        
        String[] meth = new String[2];
        switch(method) {
            case GFNFF: meth[0] = "--gfnff"; meth[1] = ""; break;
            case GFN0XTB: meth[0] = "--gfn"; meth[1] = "0"; break;
            case GFN1XTB: meth[0] = "--gfn"; meth[1] = "1"; break;
            case GFN2XTB: meth[0] = "--gfn"; meth[1] = "2"; break;
        }
        this.xtbMethod = meth;
        this.xtbMeth = method;
        
        String[] optLevel = new String[2];
        switch(level) {
            case CRUDE: optLevel[0] = "--opt"; optLevel[1] = "crude"; break;
            case SLOPPY: optLevel[0] = "--opt"; optLevel[1] = "sloppy"; break;
            case LOOSE: optLevel[0] = "--opt"; optLevel[1] = "loose"; break;
            case NORMAL: optLevel[0] = "--opt"; optLevel[1] = "normal"; break;
            case TIGHT: optLevel[0] = "--opt"; optLevel[1] = "tight"; break;
            case VERYTIGHT: optLevel[0] = "--opt"; optLevel[1] = "verytight"; break;
        }
        this.xtbOpt = optLevel;
        
        final String otherOpts = System.getenv("OGO_XTBOTHEROPTS");
        this.xtbOtherOptions = (otherOpts == null) ? null : otherOpts.split("\\s+");
        this.setEnvironment = setEnvironment;
        this.xControlFileOrig = xControlFileOrig;
    }
    
    private XTBCaller(final XTBCaller orig){
        super(orig);
        this.xtbCmd = orig.xtbCmd;
        this.xtbMethod = orig.xtbMethod;
        this.xtbMeth = orig.xtbMeth;
        this.xtbOpt = orig.xtbOpt;
        this.xtbOtherOptions = orig.xtbOtherOptions;
        this.setEnvironment = orig.setEnvironment;
        this.xControlFileOrig = orig.xControlFileOrig;
    }
    
    @Override
    public XTBCaller clone() {
        return new XTBCaller(this);
    }

    @Override
    public String myID() {
        return "Grimme xtb";
    }

    @Override
    public String myIDandMethod() {
        return "Grimme xtb " + xtbMeth.name();
    }

    @Override
    public String getMethodID() {
        return "Grimme xtb " + xtbMeth.name();
    }

    @Override
    public CartesianFullBackend getBackend(){
        return new XTBCaller(this);
    }

    @Override
    public CartesianCoordinates cartesToCartes(final long lID, final CartesianCoordinates cartes, 
        final boolean[][] constraints, final boolean isConstricted, final BondInfo bonds) throws Exception {
        
        // create a directory
        final String dirName = "xtblocopt-" + lID + "_" + System.currentTimeMillis();
        final String xyzFile = dirName + File.separator + "input.xyz";
        final String xyzRes = dirName + File.separator + "xtbopt.xyz";
        OutputPrimitives.createAFolder(dirName);
        
        // put the xyz info there
        final String[] printXYZ = cartes.createPrintableCartesians();
        OutputPrimitives.writeOut(xyzFile, printXYZ, false);
        
        // execute the command
        final int charge = cartes.getTotalCharge();
        final int spin = cartes.getTotalSpin();
        
        final List<String> cmdList = new ArrayList<>();
        cmdList.add(xtbCmd);
        cmdList.add(xtbMethod[0]);
        if(!xtbMethod[1].isEmpty()) cmdList.add(xtbMethod[1]);
        cmdList.add(xtbOpt[0]);
        cmdList.add(xtbOpt[1]);
        if(charge != 0){
            cmdList.add("--chrg");
            cmdList.add("" + charge);
        }
        if(spin != 0){
            cmdList.add("--uhf");
            cmdList.add("" + spin);
        }

        if(xControlFileOrig != null) {
            // copy the xcontrol file
            OutputPrimitives.copyFile(xControlFileOrig, dirName + File.separator + "xtb.inp");
            cmdList.add("--input");
            cmdList.add("xtb.inp");
        }

        if(xtbOtherOptions != null){
            for(final String s : xtbOtherOptions) cmdList.add(s);
        }

        cmdList.add("input.xyz");

        final ProcessBuilder pb = new ProcessBuilder(cmdList);
        pb.directory(new File(dirName));
        if(setEnvironment) {
            // ensure that xtb is only running with one thread
            final Map<String, String> envMap = pb.environment();
            envMap.put("OMP_NUM_THREADS", "1");
            envMap.put("OMP_MAX_ACTIVE_LEVELS", "1");
            envMap.put("MKL_NUM_THREADS", "1");
        }

        final Process proc = pb.start();

        // any error message?
        final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

        // any output?
        final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

        // kick them off
        errorGobbler.start();
        outputGobbler.start();

        // any error???
        final int errCode = proc.waitFor();
        if (errCode != 0) {
            throw new Exception("xtb returns non-zero return value (local optimization). Error code " + errCode);
        } 
        if (! new File(dirName + File.separator + ".xtboptok").isFile()){
            throw new Exception("xtb locopt had problems.");
        }
        // cmd should(!) have completed normally...
        
        // read the output back in
        final CartesianCoordinates res = cartes.clone();
        try{
            final String[] resDat = InputPrimitives.readFileIn(xyzRes);
            final String[] atoms = res.getAllAtomTypes();
            final double[][] xyz = res.getAllXYZCoord();
            final int noAtoms = Integer.parseInt(resDat[0].trim());
            if(noAtoms != res.getNoOfAtoms()){
                throw new Exception("Wrong number of atoms in result xyz.");
            }
            for(int i = 2; i < noAtoms+2; i++){
                final String[] line = resDat[i].trim().split("\\s+");
                // check the atom type
                if(!line[0].equalsIgnoreCase(atoms[i-2])){
                    throw new Exception("Wrong atom in output. " + line[0] + " should be " + atoms[i-2]);
                }
                // parse the coordinates
                xyz[0][i-2] = Double.parseDouble(line[1])*Constants.ANGTOBOHR;
                xyz[1][i-2] = Double.parseDouble(line[2])*Constants.ANGTOBOHR;
                xyz[2][i-2] = Double.parseDouble(line[3])*Constants.ANGTOBOHR;
            }
            
            
            // find the total energy in the output
            double totalE = FixedValues.NONCONVERGEDENERGY;
            final String[] out = outputGobbler.getData();
            for(int i = out.length - 1; i >= 0; i--) {
            	if(out[i].contains("TOTAL ENERGY")) {
            		final String[] line = out[i].trim().split("\\s+");
            		totalE = Double.parseDouble(line[3]);
            		break;
            	}
            }
            res.setEnergy(totalE);
            
        } catch(Exception e){
            throw new Exception("Failure to read xtb output",e);
        }
        
        // cleanup
        ManipulationPrimitives.remove(dirName);
        
        return res;
    }

    @Override
    public void gradientCalculation(long lID, int iIteration, double[] xyz1D, String[] saAtomTypes,
            short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges,
            short[] spins, final BondInfo bonds, final Gradient grad, final boolean hasRigidEnvironment) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double energyCalculation(long lID, int iIteration, double[] xyz1D, String[] saAtomTypes,
            short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges,
            short[] spins, final BondInfo bonds, final boolean hasRigidEnvironment){

        // create a directory
        final String dirName = "xtbenergy-" + lID + "_" + System.currentTimeMillis();
        final String xyzFile = dirName + File.separator + "input.xyz";
        try{
            OutputPrimitives.createAFolder(dirName);
        } catch(Exception e){
            System.err.println("Failure to create subdirectory in XTBCaller.energyCalculation");
            e.printStackTrace(System.err);
            return FixedValues.NONCONVERGEDENERGY;
        }

        // put the xyz info there
        final String[] printXYZ = new String[iNoOfAtoms + 2];

        printXYZ[0] = Integer.toString(iNoOfAtoms);
        printXYZ[1] = " created by OGOLEM";

        for (int i = 2; i < iNoOfAtoms + 2; i++) {
            final double dXValue = xyz1D[(i-2)] * BOHRTOANG;
            final double dYValue = xyz1D[iNoOfAtoms+(i-2)] * BOHRTOANG;
            final double dZValue = xyz1D[2*iNoOfAtoms+(i-2)] * BOHRTOANG;
            printXYZ[i] = saAtomTypes[i - 2] + "   " + String.format(Locale.US, "%20.7f", dXValue) + "   " 
                          + String.format(Locale.US, "%20.7f", dYValue) + "   " + String.format(Locale.US, "%20.7f", dZValue);
        }

        try {
            OutputPrimitives.writeOut(xyzFile, printXYZ, false);
        } catch (Exception e) {
            System.err.println("Failure to write xyz file in XTBCaller.energyCalculation");
            e.printStackTrace(System.err);
            return FixedValues.NONCONVERGEDENERGY;

        }

        // get charge/spin
        float totCharge = 0.0f;
        for(int i = 0; i < faCharges.length; i++){
            totCharge += faCharges[i];
        }
        final int charge = Math.round(totCharge);
        int spin = 0;
        for(int i = 0; i < spins.length; i++){
            spin += spins[i];
        }

        final List<String> cmdList = new ArrayList<>();
        cmdList.add(xtbCmd);
        cmdList.add(xtbMethod[0]);
        if(!xtbMethod[1].isEmpty()) cmdList.add(xtbMethod[1]);
        if(charge != 0){
            cmdList.add("--chrg");
            cmdList.add("" + charge);
        }
        if(spin != 0){
            cmdList.add("--uhf");
            cmdList.add("" + spin);
        }

        if(xControlFileOrig != null) {
            // copy the xcontrol file
            try {
                OutputPrimitives.copyFile(xControlFileOrig, dirName + File.separator + "xtb.inp");
            } catch (Exception e){
                System.err.println("Failure to copy xtb input file in XTBCaller.energyCalculation");
                e.printStackTrace(System.err);
                return FixedValues.NONCONVERGEDENERGY;
            }
            cmdList.add("--input");
            cmdList.add("xtb.inp");
        }

        if(xtbOtherOptions != null){
            for(final String s : xtbOtherOptions) cmdList.add(s);
        }

        cmdList.add("input.xyz");

        final ProcessBuilder pb = new ProcessBuilder(cmdList);
        pb.directory(new File(dirName));
        if(setEnvironment) {
            // ensure that xtb is only running with one thread
            final Map<String, String> envMap = pb.environment();
            envMap.put("OMP_NUM_THREADS", "1"); 
            envMap.put("OMP_MAX_ACTIVE_LEVELS", "1");
            envMap.put("MKL_NUM_THREADS", "1");
        }

        String[] xtbOut = null;
        try {
            final Process proc = pb.start();
            
            // any error message?
            final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

            // any output?
            final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

            // kick them off
            errorGobbler.start();
            outputGobbler.start();

            // any error???
            if (proc.waitFor() != 0) {
                throw new Exception("xtb returns non-zero return value (energy calculation).");
            }
            xtbOut = outputGobbler.getData(); 
        } catch (Exception e) {
            System.err.println("Failure in cmd execution in XTBCaller.energyCalculation");
            e.printStackTrace(System.err);
            return FixedValues.NONCONVERGEDENERGY;
        }

        double energy = FixedValues.NONCONVERGEDENERGY;
        try{
            // find the total energy in the output
            double totalE = FixedValues.NONCONVERGEDENERGY;
            for(int i = xtbOut.length - 1; i >= 0; i--) {
                if(xtbOut[i].contains("TOTAL ENERGY")) {
                    final String[] line = xtbOut[i].trim().split("\\s+");
                    energy = Double.parseDouble(line[3]);
                    break;
                }
            }
        } catch(Exception e){
            System.err.println("Failure in reading energy in XTBCaller.energyCalculation");
            e.printStackTrace(System.err);
            return FixedValues.NONCONVERGEDENERGY;
        }

        //cleanup
        try {
            ManipulationPrimitives.remove(dirName);
        } catch (Exception e) {
            System.err.println("failure to remove subdir in XTBCaller.energyCalculation");
            e.printStackTrace(System.err);
        }

        return energy;
    }
}
