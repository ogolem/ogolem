/**
Copyright (c) 2010-2013, J. M. Dieterich
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
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Calls the Local Optimization with a TRivial Force Field code of
 * Bernd Hartke.
 * @author Johannes Dieterich
 * @version 2013-09-23
 */
final class LotriffCaller extends AbstractLocOpt {

    // the ID
    private static final long serialVersionUID = (long) 20101108;
    private String[] saAuxData;

    LotriffCaller(final GlobalConfig globConf){
        super(globConf);
        String sWhichAuxFile = globConf.outputFolder + "-lotriff.aux";
        // read auxfile
        try{
            saAuxData = Input.ReadFile(sWhichAuxFile);
            if(!saAuxData[0].equalsIgnoreCase("###OGOLEMAUX###")){
                throw new ConvergenceException("The read in auxiliary file is not valid!");
            }
        }catch(Exception e){
            System.err.println("Error in reading the mandatory aux-file! " + e.toString());
            saAuxData = null;
        }
    }

    private LotriffCaller(final LotriffCaller orig){
        super(orig);
        if(orig.saAuxData == null){
            this.saAuxData = null;
        } else{
            this.saAuxData = orig.saAuxData.clone();
        }
    }

    @Override
    public LotriffCaller clone(){
        return new LotriffCaller(this);
    }

    @Override
    public String myIDandMethod(){
        return "Lotriff";
    }
    
    @Override
    public String myID(){
        return "Lotriff";
    }

    @Override
    public CartesianCoordinates cartesToCartes(long lID, CartesianCoordinates cartes,
            final boolean[][] baConstraints, final boolean isConstricted, final BondInfo bonds)
            throws InitIOException, IOException, ConvergenceException, CastException{

        final String sFolder = "lotriff" + lID;
        final String sInputXYZ = sFolder + System.getProperty("file.separator")
                + "ingeo.xyz";
        final String sOutputXYZ = sFolder + System.getProperty("file.separator") + "outgeo.xyz";
        

        // we want to run these in folders, so first create a folder (might throw exception)
        OutputPrimitives.createAFolder(sFolder);

        // create the actual input file
        String[] saXYZ = cartes.createPrintableCartesians();

        // now add the atom definitions (start at second line in aux file)
        int iAuxLine = 1;
        for(int i = 0; i < cartes.getNoOfAtoms(); i++){
            if(baConstraints[0][i] && baConstraints[1][i] && baConstraints[2][i]){
                // add atom type
                saXYZ[i+2] += "   " + saAuxData[iAuxLine];
                iAuxLine++;
            } else if(!baConstraints[0][i] && !baConstraints[1][i] && !baConstraints[2][i]){
                continue;
            } else{
                // error
                System.err.println("ERROR: Unknown constraints for atom " + i + ". "
                        + "Lotriff only understands full or no constraint per atom.");
            }
        }

        // write the adjusted xyz (exception possible)
        OutputPrimitives.writeOut(sInputXYZ, saXYZ, false);

        
        try{
            Runtime rt = Runtime.getRuntime();
            String sLotriffCmd = System.getenv("OGO_LOTRIFFCMD");
            if(sLotriffCmd == null){
                sLotriffCmd = "lotriff";
            }
            File dir = new File(sFolder);
            Process proc = rt.exec(sLotriffCmd, null,dir);

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
                throw new ConvergenceException("Lotriff returns non-zero return value (local optimization).");
            } else{
                // lotriff should(!) have completed normally...
            }


        } catch(Exception e){
            e.printStackTrace(System.err);
            throw new ConvergenceException("Lotriff has a problem (local optimization).");
        }

        // parse the output (exception possible)
        final String[] saNewXYZ = org.ogolem.io.InputPrimitives.readFileIn(sOutputXYZ);
        // first the energy in the second line
        final String[] saTmpEnergy = saNewXYZ[1].trim().split("\\s+");
        final double dEnergy = Double.parseDouble(saTmpEnergy[0].trim())* Constants.KJTOHARTREE;
        final double[][] daNewXYZ = new double[3][cartes.getNoOfAtoms()];

        for(int i = 0; i < cartes.getNoOfAtoms(); i++){
            // split at whitespace
            final String[] saTemp = saNewXYZ[i+2].trim().split("\\s+");
            for(int j = 0; j < 3; j++){
                daNewXYZ[j][i] = Double.parseDouble(saTemp[j+1].trim())
                        * Constants.ANGTOBOHR;
            }
        }

        // put it into cartes
        cartes.setEnergy(dEnergy);
        cartes.setAllXYZ(daNewXYZ);

        // delete the folder (might throw an exception)
        try{
            ManipulationPrimitives.remove(sFolder);
        } catch(Exception e){
            System.err.println("WARNING: Couldn't remove folder " + sFolder + ". "
                    + "Ignoring for now. " + e.toString());
        }
        
        return cartes;
    }
}
