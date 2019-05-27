/**
Copyright (c) 2013, J. M. Dieterich
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
package org.ogolem.scanner;

import java.io.File;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Geometry;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Input;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Scans a set of explicit degrees of freedom and (if wished) quench the resulting
 * structures.
 * @author Johannes Dieterich
 * @version 2013-09-22
 */
public class MainScanner {
    
    public static void run(String args[]) {
        
        if(args[0].equalsIgnoreCase("help")){
            System.out.println("This is the scanning functionality. We only need a global config as the argument.");
            System.out.println("Inside the global config, we need:");
            System.out.println(" * ONE molecule which must be flexible");
            System.out.println(" * a <SCANCONF> block containing at least one DoF scanning specification.");
            return;
        }
        
        final String globConf = args[0];
        final int dotPart = args[0].lastIndexOf(".ogo");
        final String outFolder = args[0].substring(0,dotPart) + "-scan";
        final String outputPrefix = outFolder + File.separator + "scangeom";
        GlobalConfig config = null;
        try{
            config = Input.ConfigureMe(globConf);
        } catch(Exception e){
            System.err.println("ERROR: Couldn't read global config.");
            e.printStackTrace(System.err);
            System.exit(1);
        }
        
        // create a reference geometry
        final Geometry refGeom = new Geometry(config.geoConfCopy());
        if(!refGeom.isThereAFlexy()){
            System.err.println("ERROR: No flexible molecule in parsed config. Therefore: no scanning.");
            System.exit(2);
        }
        
        if(refGeom.getNumberOfIndieParticles() != 1){
            System.err.println("ERROR: Unfortunately, we do not support geometries with more than one molecule for scanning currently. Inform author(s) if you are stuck here.");
            System.exit(3);
        }
        
        // create our scanner config
        ScannerConfig scanConf = null;
        try{
            final String[] inp = InputPrimitives.readFileIn(globConf);
            scanConf = new ScannerConfig(inp);
        } catch(Exception e){
            System.err.println("ERROR: Failure to configure scanner.");
            e.printStackTrace(System.err);
            System.exit(4);
        }
        
        // get our cartesians and all the other stuff
        final CartesianCoordinates cartes = refGeom.getCartesians();
        final boolean[][] constraints = refGeom.getAllConstraintsXYZ(true);
        final boolean isConstricted = refGeom.isThereAConstraint();
        
        // create the folder for our output
        try {
            final File f = new File(outFolder);
            if(f.exists() && f.isDirectory()){
                System.out.println("INFO: Output folder exists. This will overwrite whatever is in there!");
            } else{
                OutputPrimitives.createAFolder(outFolder);
            }
        } catch(Exception e){
            System.out.println("ERROR: Couldn't create output folder for scanned geometries.");
            e.printStackTrace(System.err);
            System.exit(5);
        }
        
        // dunk the first one (yeah, here is our limitation) into the scanner core
        ScannerCore.scan(scanConf, cartes.giveMolecularCartes(0, true), config.getRefNewton().clone(), outputPrefix, constraints, isConstricted, refGeom.getBondInfo());
    }
}
