/**
Copyright (c) 2012-2014, J. M. Dieterich
              2015, J. M. Dieterich and B.Hartke
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
package org.ogolem.tests;

import java.util.List;
import org.ogolem.core.*;
import org.ogolem.core.CollisionInfo.Collision;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericGlobalOptimization;

/**
 * Enables some simple tests mainly for development.
 *
 * @author Johannes Dieterich
 * @version 2015-07-25
 */
public class MainTests {

    private static final boolean VERBOSE = true;
    
    public static void run(final String[] args) {

        /*
         * PARSE THE ARGUMENTS
         * 1) the ogo file
         * 2) the geom file
         * 3) globopt?
         * 4) secondary file (if applicable)
         */
        final String ogoFile = args[0];
        final String geomFile1 = args[1];
        final boolean globopt = Boolean.parseBoolean(args[2]);
        String geomFile2 = null;
        if(globopt) {geomFile2 = args[3];}

        // get the global config
        GlobalConfig conf = null;
        try{
            conf = Input.ConfigureMe(ogoFile);
        } catch(Exception e){
            System.err.println("ERROR: Reading configuration failed.");
            e.printStackTrace(System.err);
            System.exit(234);
        }
        
        assert(conf != null);
        System.out.println("INFO: Configuration successfully parsed and loaded.");
        
        // read the geometry in
        Geometry geom1 = null;
        Geometry geom2 = null;
        try{
            geom1 = SimpleSeedingInit.readSeedIn(geomFile1, conf.geoConfCopy(), 0);
            if(globopt) {geom2 = SimpleSeedingInit.readSeedIn(geomFile2, conf.geoConfCopy(), 0);}
        } catch(Exception e){
            System.err.println("ERROR: Reading geometry failed.");
            e.printStackTrace(System.err);
            System.exit(235);
        }
        
        assert(geom1 != null);
                
        if(VERBOSE){
            System.out.println("INFO: READ IN GEOMETRY COMING");
            final String[] sa = geom1.makePrintableAbsoluteCoord(true);
            for(final String s : sa) {System.out.println(s);}
            if(globopt){
                assert(geom2 != null);
                System.out.println("INFO: READ IN GEOMETRY2 COMING");
                final String[] sa2 = geom2.makePrintableAbsoluteCoord(true);
                for(final String s : sa2) {System.out.println(s);}
            }
        }
        
        // CD and DD
        final CollisionDetectionEngine cd = new CollisionDetection(conf.getWhichCollisionEngine());
        
        final CartesianCoordinates c1 = geom1.getCartesians();
        final CollisionInfo info1 = cd.checkForCollision(c1, conf.getBlowFacBondDetect(), conf.geoConfCopy().bonds);
        System.out.println("INFO: Collision found for geom 1? " + info1.hasCollision());
        if(info1.hasCollision()){
            final List<Collision> collisions = info1.getCollisions();
            String s = "";
            for(final Collision coll : collisions){
                s += " between " + coll.getAtomOne() + " and " + coll.getAtomTwo();
            }
            System.out.println("INFO: Collision between (at least) " + s);
        }
        
        final boolean diss1 = DissociationDetection.checkForDissociation(info1.getPairWiseDistances(), c1.getAllAtomTypes(),
                c1.getAllAtomNumbers(), conf.getBlowFacDissocDetect(), conf.getWhichDissociationEngine());
        System.out.println("INFO: Dissociation found for geom 1? " + diss1);
        
        if(globopt){
            assert(geom2 != null);
            final CartesianCoordinates c2 = geom2.getCartesians();
            final CollisionInfo info2 = cd.checkForCollision(c2, conf.getBlowFacBondDetect(), conf.geoConfCopy().bonds);
            System.out.println("INFO: Collision found for geom 2? " + info2.hasCollision());
            if(info2.hasCollision()){
                final List<Collision> collisions = info1.getCollisions();
                String s = "";
                for(final Collision coll : collisions){
                    s += " between " + coll.getAtomOne() + " and " + coll.getAtomTwo();
                }
                System.out.println("INFO: Collision between (at least) " + s);
            }
            final boolean diss2 = DissociationDetection.checkForDissociation(info2.getPairWiseDistances(), c2.getAllAtomTypes(),
                c2.getAllAtomNumbers(), conf.getBlowFacDissocDetect(), conf.getWhichDissociationEngine());
            System.out.println("INFO: Dissociation found for geom 2? " + diss2);
        }
        
        // locopt
        GenericFitnessFunction<Molecule,Geometry> opt = null;
        try{
            opt = conf.getFitnessFunction();
        } catch(Exception e){
            System.err.println("Unable to obtain the fitness function for this geometry.");
            e.printStackTrace(System.err);
            System.exit(42);
        }
        geom1.setFitness(FixedValues.NONCONVERGEDENERGY);
        final Geometry geomL1 = opt.fitness(geom1, false);
        if(globopt){
            assert(geom2 != null);
            geom2.setFitness(FixedValues.NONCONVERGEDENERGY);
        }
        final Geometry geomL2 = (globopt) ? opt.fitness(geom2, false) : null;
        final double finalFitness1 = geomL1.getFitness();
        System.out.println("full fitness: " + finalFitness1 + " hartree");
        System.out.println("full fitness: " + (finalFitness1 * Constants.HARTREETOKJ) + " kJ/mol");
        System.out.println("INFO: LOCALLY OPTIMIZED GEOMETRY COMING");
        final String[] sa = geomL1.makePrintableAbsoluteCoord(true);
        for(final String s : sa) {System.out.println(s);}
        if(globopt){
            System.out.println("INFO: LOCALLY OPTIMIZED GEOMETRY2 COMING");
            assert(geomL2 != null);
            final String[] sa2 = geomL2.makePrintableAbsoluteCoord(true);
            for(final String s : sa2) {System.out.println(s);}
        }

        // globopt?
        if(globopt){
            GenericGlobalOptimization<Molecule,Geometry> glob = null;
            try{
                glob = conf.getGlobalOptimization();
            } catch(Exception e){
                System.err.println("Couldn't instantiate global optimization.");
                e.printStackTrace(System.err);
                System.exit(42);
            }
            final String globoptName = glob.getMyID();
            System.out.println("The globopt algo is " + globoptName);
            final Geometry geom = glob.globalOptimization(2, geomL1, geomL2);
            System.out.println("INFO: GLOBALLY OPTIMIZED GEOMETRY COMING");
            final String[] sa2 = geom.makePrintableAbsoluteCoord(true);
            for(final String s : sa2) {System.out.println(s);}
        }
    }
}
