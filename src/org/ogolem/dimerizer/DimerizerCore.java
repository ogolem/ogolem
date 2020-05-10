/**
Copyright (c) 2014, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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

import contrib.edu.princeton.eac.Grid;
import java.util.Arrays;
import java.util.List;
import org.ogolem.core.CollisionDetection;
import org.ogolem.core.CollisionDetectionEngine;
import org.ogolem.core.Geometry;
import org.ogolem.core.Newton;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * The actual working core of the dimerizer.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
class DimerizerCore {
    
    private static final Logger log = LoggerFactory.getLogger(DimerizerCore.class);
    
    static void runAllPoints(final Grid<Geometry> grid, final Geometry geom, final Newton newton,
            final double blowFac){
    
        final CollisionDetectionEngine cd = new CollisionDetection(CollisionDetection.CDTYPE.SIMPLEPAIRWISE);
        
        // setup XYZ constraints on the first atom of each of the molecules.
        final Geometry work = geom.copy();
        final boolean[][] constr1 = work.getMoleculeAtPosition(0).getConstraints();
        constr1[0][0] = true;
        constr1[1][0] = true;
        constr1[2][0] = true;
        final boolean[][] constr2 = work.getMoleculeAtPosition(1).getConstraints();
        constr2[0][0] = true;
        constr2[1][0] = true;
        constr2[2][0] = true;
    
        final List<double[]> allGridPoints = grid.getAllGridPoints();
        
        log.info("Total number of grid points to be done " + allGridPoints.size());
        
        long counter = 0l;
        for(final double[] point : allGridPoints){
            
            log.info("Running " + counter + " " + Arrays.toString(point));
            
            final Geometry thisG = runForOne(point, work, newton, cd, blowFac, counter);
            // put in
            try{
                log.info("Trying to add " + counter);
                grid.attachToPoint(point, thisG);
            } catch(Exception e){
                System.err.println("Something failed for grid point " + point[0] + "\t" + point[1]
                    + "\t" + point[2] + "\t" + point[3] + "\t" + point[4] + "\t" + point[5]);
                e.printStackTrace(System.err);
            }
            counter++;
        }
    }
    
    private static Geometry runForOne(final double[] point, final Geometry geom, final Newton newton,
            final CollisionDetectionEngine cd, final double blowFac, final long counter){
            
        // create one just for us
        final Geometry ggg = geom.copy();
        
        // rotate, translate
        final double[] com2 = ggg.getCOM(1);
        final double[] extCoords = point.clone();
        extCoords[0] += com2[0];
        extCoords[1] += com2[1];
        extCoords[2] += com2[2];
        
        ggg.setExtCoordMolecule(extCoords, 1);
        ggg.setID(counter);
        
        // check CD (makes sense, really!)
        final boolean hasColl = cd.checkOnlyForCollision(ggg.getCartesians(), blowFac, ggg.getBondInfo());
        if(hasColl){
            return null;
        }
    
        // relax
        final Geometry opt = newton.localOptimization(ggg);
    
        // return
        return opt;
    }
}
