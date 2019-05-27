/**
Copyright (c) 2014, J. M. Dieterich
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
package org.ogolem.dimerizer;

import java.io.File;
import java.util.List;
import org.ogolem.core.FixedValues;
import org.ogolem.core.Geometry;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Input;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.InquiryPrimitives;
import org.ogolem.io.OutputPrimitives;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A very basic grid search for a full 6D hypersurface of dimerization.
 * @author Johannes Dieterich
 * @version 2014-03-07
 */
public class MainDimerizer {
    
    private static final Logger log = LoggerFactory.getLogger(MainDimerizer.class);
    
    public static void run(final String[] args){
        
        if(args == null || args[0].equalsIgnoreCase("help")){
            System.out.println("This does a full 6D grid search for a dimerization. Required arguments:");
            System.out.println(" * an ogo file for the basic stuff including a geometry containing TWO molecules, locopt etc pp");
            System.out.println(" * an dimerizing config for the actual stuff (grid/locopt/molecules)");
            System.out.println(" * the output directory");
            System.exit(0);
        }
        
        // setup
        GlobalConfig conf = null;
        DimerizerConfig dimConf = null;
        String outputDir = null;
        try{
            final String ogoFile = args[0];
            final String dimConfFile = args[1];
            outputDir = args[2];
            
            // first get the global config
            conf = Input.ConfigureMe(ogoFile);
            
            // now parse the dimer config
            final String[] dat = InputPrimitives.readFileIn(dimConfFile);
            dimConf = new DimerizerConfig(dat,conf);
            
        } catch(Exception e){
            System.err.println("ERROR: Couldn't setup configuration.");
            e.printStackTrace(System.err);
            System.exit(1);
        }
        
        log.info("Done with the config parsing!");
        
        // run the actual dimerization kernel
        assert(conf != null);
        assert(dimConf != null);
        assert(dimConf.grid != null);
        assert(dimConf.geom != null);
        assert(dimConf.locopt != null);
        DimerizerCore.runAllPoints(dimConf.grid, dimConf.geom, dimConf.locopt,
                conf.getBlowFacBondDetect());
        
        log.info("Completed the dimerizer core!");
        
        // print out the results
        final List<double[]> allGridPoints = dimConf.grid.getAllGridPoints();
        final List<Geometry> allOptGeoms = dimConf.grid.getListOfAllAttachments();
        
        try{
            OutputPrimitives.createAFolder(outputDir);
            final String sep = File.separator;
        
            log.info("Number of points to be plotted " + allGridPoints.size());
            
            for(int i = 0; i < allGridPoints.size(); i++){
                                
                final double[] p = allGridPoints.get(i);
                final Geometry g = allOptGeoms.get(i);
                final String outputPath = outputDir + sep + "geom-relaxed_" + 
                        + p[0] + "_" + p[1] + "_" + p[2] + "_" + p[3] + "_"
                        + p[4] + "_" + p[5] + ".xyz";
                
                final boolean fileExist = InquiryPrimitives.doesFileExist(outputPath);
                if(fileExist){
                    System.err.println("File " + outputPath + " exits.");
                    //System.exit(1);
                }
                
                if(g == null){
                    System.err.println("No data for " + outputPath);
                    continue;
                }
                
                if(g.getFitness() >= FixedValues.NONCONVERGEDENERGY){
                    System.err.println("No useful data for " + outputPath);
                    continue;
                }
                
                final String[] geomDat = g.makePrintableAbsoluteCoord(true);
                OutputPrimitives.writeOut(outputPath, geomDat, true);
            }
        } catch(Exception e){
            System.err.println("ERROR: Something went wrong in printing results... :-(");
            e.printStackTrace(System.err);
        }
    }
}
