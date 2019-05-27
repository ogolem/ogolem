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
package org.ogolem.dimeranalyzer;

import contrib.edu.princeton.eac.Grid;
import java.io.File;
import java.util.LinkedList;
import java.util.List;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Constants;
import org.ogolem.core.CoordTranslation;
import org.ogolem.core.Geometry;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Input;
import org.ogolem.dimerizer.DimerizerConfig;
import org.ogolem.helpers.Tuple;
import org.ogolem.io.FileFilter;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.OutputPrimitives;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Very basic analysis of dimerizer output.
 * @author Johannes Dieterich
 * @version 2014-05-12
 */
public class MainDimerAnalyzer {
    
    private static final double RMSDCUTOFF = 0.05;
    private static final Logger log = LoggerFactory.getLogger(MainDimerAnalyzer.class);
    
    public static void run (final String[] args){
        
        if(args.length == 0 || args[0].equalsIgnoreCase("help")){
            System.out.println("This provides BASIC analysis of dimerizing runs.");
            System.out.println("Needed input:");
            System.out.println(" * the folder containing the output data of the dimerizer");
            System.out.println(" * the energy of subsystem 1 in Hartree");
            System.out.println(" * the energy of subsystem 2 in Hartree");
            System.out.println(" * an ogo file for the basic stuff including a geometry containing TWO molecules, locopt etc pp");
            System.out.println(" * an dimerizing config used for the actual stuff (grid/locopt/molecules)");
            System.exit(0);
        }
        
        final String folder = args[0];
        final double eSys1 = Double.parseDouble(args[1]);
        final double eSys2 = Double.parseDouble(args[2]);
        
        GlobalConfig conf = null;
        DimerizerConfig dimConf = null;
        try{
            final String ogoFile = args[3];
            final String dimConfFile = args[4];
            
            // first get the global config
            conf = Input.ConfigureMe(ogoFile);
            
            // now parse the dimer config
            final String[] dat = InputPrimitives.readFileIn(dimConfFile);
            dimConf = new DimerizerConfig(dat,conf,2); // ONLY x and y   
        } catch(Exception e){
            System.err.println("ERROR: Couldn't setup configuration.");
            e.printStackTrace(System.err);
            System.exit(1);
        }
        
        final File fol = new File(folder);
        if(!fol.exists() || !fol.isDirectory()){
            System.err.println("Output folder does not exist or is not a directory.");
            System.exit(2);
        }
        
        /*
         * let's go over our grid
         */
        final Grid<Geometry> grid = dimConf.getGrid();
        final List<String> avgFormation = new LinkedList<>();
        final List<String> maxFormation = new LinkedList<>();
        double lastXVal = Double.NaN;
        for(final double[] p : grid){
            
            log.info("Working on grid point " + p[0] + "\t" + p[1]);
            
            /*
             * assemble prefix and
             * find a list of all the files with that
             */
            final String prefix = "geom-relaxed_" + p[0] + "_" + p[1];
            final String[] allFiles = fol.list(new FileFilter(prefix,true));
            
            final List<CartesianCoordinates> cartes = new LinkedList<>();
            for(final String file : allFiles){
                try{
                    final CartesianCoordinates c = Input.readCartesFromFile(folder + File.separator + file);
                    
                    // read in the energy as well (yes, somewhat suboptimal)
                    final String[] sa = InputPrimitives.readFileIn(folder + File.separator + file);
                    final String[] lineCont = sa[1].trim().split("\\s+");
                    final double e = Double.parseDouble(lineCont[1])*Constants.KJTOHARTREE;
                    
                    // at least, directly correct to a "heat of formation"
                    c.setEnergy(e-eSys1-eSys2);
                    
                    tryToAppend(c,cartes, RMSDCUTOFF);
                } catch(Exception e){
                    System.err.println("Something seriously cocked up in your output folder. File was " + (folder + File.separator + file));
                    e.printStackTrace(System.err);
                }
            }
            
            // now we have a list of all the (non-redundant) cartesian coordinates
            if(p[0] != lastXVal){
                avgFormation.add("");
                maxFormation.add("");
                lastXVal = p[0];
            }
            
            double max = 0.0;
            double avg = 0.0;
            final int no = cartes.size();
            for(final CartesianCoordinates c : cartes){
                final double e = c.getEnergy();
                avg += e;
                max = Math.min(max,e);
            }
            avg /= no;
            
            avgFormation.add(" " + p[0] + "\t" + p[1] + "\t" + avg);
            maxFormation.add(" " + p[0] + "\t" + p[1] + "\t" + max);
        }
        
        // well, just output left
        final String avgFile = folder + "-avg.dat";
        final String maxFile = folder + "-max.dat";
        
        try{
            OutputPrimitives.writeOut(maxFile, maxFormation, false);
            OutputPrimitives.writeOut(avgFile, avgFormation, false);
        } catch(Exception e){
            System.err.println("Error writing analysis out.");
            e.printStackTrace(System.err);
        }
    }
    
    private static void tryToAppend(final CartesianCoordinates c, final List<CartesianCoordinates> cartes,
            final double rmsdCut) throws Exception {
        
        // loop over and check for doubles
        for(final CartesianCoordinates c2 : cartes){
            // since we can assume that the two cartesian coordinates are in the same atom order, we should be able to use
            // a standard Kearlsey overlap
            final Tuple<CartesianCoordinates, Double> tup = CoordTranslation.alignTwoCartesians(c, c2);
            
            // check the rmsd
            if(tup.getObject2() < rmsdCut){
                return;
            }
        }

        // if we made it here, we are cool
        cartes.add(c);
    }
}
