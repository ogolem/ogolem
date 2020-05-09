/**
Copyright (c) 2014-2020, J. M. Dieterich and B. Hartke
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

import contrib.bobyqa.AbstractBOBYQAMethod;
import contrib.bobyqa.BOBYQAOptimizer;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.generic.Copyable;
import org.ogolem.generic.GenericFitnessBackend;
import org.ogolem.generic.GenericLocOpt;
import org.ogolem.helpers.Tuple;

/**
 * A collection of optimization strategies for GDM.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class DirMutOptStrategies implements Serializable {
    
    private static final boolean DEBUG = false;
    private static final long serialVersionUID = (long) 20200429;
    
    static interface PointOptStrategy extends Copyable, Serializable {
        
    	@Override
        PointOptStrategy copy();
        
        String getMyID();
        
        /**
         * Do (optionally) an optimization of the geometry at a certain point.
         * @param geom the to be optimized geometry.
         * @param locopt the local optimization to be used. Can just be a fitness function.
         * @param moveMol the molecule to move
         * @param movePartner the partner of the movable molecule
         * @return a tuple of the fitness and the geometry. The geometry may be identical to geom!
         */
        Tuple<Double,Geometry> pointOpt(final Geometry geom, final GenericLocOpt<Molecule,Geometry> locopt, final int moveMol, final int movePartner);
    }
    
    static PointOptStrategy parseStrategy (final String configString,
            final CollisionDetectionEngine collDetect, final double blowColl, final boolean fullyRelaxed)
            throws Exception {
        
        PointOptStrategy strategy = null;
        if(configString.equalsIgnoreCase("energyonly")){
            strategy = new EnergyOnlyEvaluator();
        } else if(configString.startsWith("euleropt")){
            // XXX make configurable
            final RigidBodyOptimization.GDMBOBYQAConfig config = new RigidBodyOptimization.GDMBOBYQAConfig();
            strategy = new EulerOnlyOptimization(config,collDetect,blowColl,fullyRelaxed);
        } else if(configString.startsWith("rigidopt")){
            
            final String[] tokens = configString.substring(8).trim().split("/");
            for(int i = 0; i < tokens.length; i++){tokens[i] = tokens[i].trim();}
            
            double maxCOMBorders = 10.0*Constants.ANGTOBOHR;
            
            for(final String token : tokens){
                if(token.isEmpty()){continue;}
                if(token.startsWith("comborders=")){
                    maxCOMBorders = Double.parseDouble(token.substring(11).trim());
                } else {
                    throw new RuntimeException("Unknown option " + token + " in GDM rigidopt.");
                }
            }
            
            // XXX make configurable
            final RigidBodyOptimization.GDMBOBYQAConfig config = new RigidBodyOptimization.GDMBOBYQAConfig();
            strategy = new RigidBodyOptimization(maxCOMBorders,config,collDetect,blowColl,fullyRelaxed);
        } else if(configString.startsWith("fullopt")){
            strategy = new FullOptimization();
        } else {
            throw new RuntimeException("Unknown GDM strategy: " + configString);
        }
        
        return strategy;
    }
    
    static class FullOptimization implements PointOptStrategy {

        private static final long serialVersionUID = (long) 20140526;
        
        @Override
        public PointOptStrategy copy() {
            return new FullOptimization();
        }

        @Override
        public String getMyID() {
            return "fully flexible strategy";
        }

        @Override
        public Tuple<Double, Geometry> pointOpt(final Geometry geom, final GenericLocOpt<Molecule,Geometry> locopt, final int moveMol, final int movePartner) {
            
            final Geometry gOpt = locopt.fitness(geom, false);
            
            return new Tuple<>(gOpt.getFitness(),gOpt);
        }
    }
    
    static class EulerOnlyOptimization implements PointOptStrategy {
        
        private static final long serialVersionUID = (long) 20140526;

        private final CollisionDetectionEngine collDetect;
        private final RigidBodyOptimization.GDMBOBYQAConfig bobConf;
        private final double blowColl;
        private final boolean fullyRelaxed;
        
        EulerOnlyOptimization(final RigidBodyOptimization.GDMBOBYQAConfig config,
                final CollisionDetectionEngine collDetect, final double blowColl,
                final boolean fullyRelaxed){
            this.bobConf = config;
            this.blowColl = blowColl;
            this.fullyRelaxed = fullyRelaxed;
            this.collDetect = collDetect;
        }
        
        @Override
        public PointOptStrategy copy() {
            return new EulerOnlyOptimization(this.bobConf, this.collDetect.clone(), this.blowColl, this.fullyRelaxed);
        }

        @Override
        public String getMyID() {
            return "Euler only optimization";
        }

        @Override
        public Tuple<Double, Geometry> pointOpt(final Geometry geom, final GenericLocOpt<Molecule,Geometry> locopt, final int moveMol, final int movePartner) {
            
            final double eBefore = geom.getFitness();
            if(DEBUG){System.out.println("DEBUG: Starting with energy " + eBefore);}
            
            // use bobyqa for a quick optimization of that one moved molecule
            
            final double[] eulerFirst = geom.getEulers(moveMol); // this we will tinker with
            final double[] guess = new double[3];
            
            final double[][] boundaries = new double[2][3];
            // phi
            boundaries[0][0] = -Math.PI;
            boundaries[1][0] = Math.PI;
            // omega
            boundaries[0][1] = -0.5*Math.PI;
            boundaries[1][1] = 0.5*Math.PI;
            // psi
            boundaries[0][2] = -Math.PI;
            boundaries[1][2] = Math.PI;
            
            final BOBYQAOptimizer bobyqa = new BOBYQAOptimizer(boundaries, 7,
                    bobConf.initialTrust, bobConf.stoppingTrust, bobConf.maxIter,
                    true); // we do normalize
            
            // do something...
            guess[0] = eulerFirst[0];
            guess[1] = eulerFirst[1];
            guess[2] = eulerFirst[2];
            
            final GDMAdaptor adap = new GDMAdaptor(geom, moveMol,
                locopt, collDetect, guess, boundaries, blowColl, fullyRelaxed,
                true);
            
            final double[] init = new double[3];
            adap.normalizeFromBounds(boundaries, guess, init);
            
            final Tuple<Double, double[]> opt = bobyqa.doOptimize(init, adap);

            final double currBestE = opt.getObject1();
            
            if (currBestE >= FixedValues.NONCONVERGEDENERGY) {
                System.out.println("INFO: Could not find any place that was collision free and had a reasonable energy. Complaining.");
                return new Tuple<>(FixedValues.NONCONVERGEDENERGY, geom.copy());
            }

            // work is already setup properly, except for the energy :-)
            final Geometry res = adap.getCurrentGeom().copy();
            res.setFitness(currBestE);
            
            return new Tuple<>(res.getFitness(),res);
        }
        
    }
    
    static class RigidBodyOptimization implements PointOptStrategy {
        
        private static final long serialVersionUID = (long) 20140526;
        
        private final double halfDiff;
        private final CollisionDetectionEngine collDetect;
        private final GDMBOBYQAConfig bobConf;
        private final double blowColl;
        private final boolean fullyRelaxed;
        
        RigidBodyOptimization(final double maxDiff, final GDMBOBYQAConfig config,
                final CollisionDetectionEngine collDetect, final double blowColl,
                final boolean fullyRelaxed){
            this.halfDiff = maxDiff;
            this.bobConf = config;
            this.blowColl = blowColl;
            this.fullyRelaxed = fullyRelaxed;
            this.collDetect = collDetect;
        }

        @Override
        public PointOptStrategy copy() {
            return new RigidBodyOptimization(this.halfDiff, this.bobConf,
                this.collDetect.clone(), this.blowColl, this.fullyRelaxed);
        }

        @Override
        public String getMyID() {
            return "rigid body optimization of all external coordinates";
        }

        @Override
        public Tuple<Double, Geometry> pointOpt(final Geometry geom, final GenericLocOpt<Molecule,Geometry> locopt, final int moveMol, final int movePartner) {
            
            final double eBefore = geom.getFitness();
            if(DEBUG){System.out.println("DEBUG: Starting with energy " + eBefore);}
            
            // use bobyqa for a quick optimization of that one least connected molecule
            
            final double[] comMe = geom.getCOM(moveMol);
            final double[] guess = new double[6];
            
            final double[][] boundaries = new double[2][6];
            boundaries[0][0] = -halfDiff + comMe[0];
            boundaries[1][0] = halfDiff + comMe[0];
            
            boundaries[0][1] = -halfDiff + comMe[1];
            boundaries[1][1] = halfDiff + comMe[1];
            
            boundaries[0][2] = -halfDiff + comMe[2];
            boundaries[1][2] = halfDiff + comMe[2];
            // phi
            boundaries[0][3] = -Math.PI;
            boundaries[1][3] = Math.PI;
            // omega
            boundaries[0][4] = -0.5*Math.PI;
            boundaries[1][4] = 0.5*Math.PI;
            // psi
            boundaries[0][5] = -Math.PI;
            boundaries[1][5] = Math.PI;
            
            final BOBYQAOptimizer bobyqa = new BOBYQAOptimizer(boundaries, bobConf.interPoints,
                    bobConf.initialTrust, bobConf.stoppingTrust, bobConf.maxIter,
                    true); // we do normalize
            
            // do something...
            final double[] com = geom.getCOM(moveMol);
            final double[] eulers =  geom.getEulers(moveMol);
            
            guess[0] = com[0];
            guess[1] = com[1];
            guess[2] = com[2];
            guess[3] = eulers[0];
            guess[4] = eulers[1];
            guess[5] = eulers[2];

            final GDMAdaptor adap = new GDMAdaptor(geom, moveMol,
                locopt, collDetect, guess, boundaries, blowColl, fullyRelaxed,
                false);
            
            final double[] init = new double[6];
            adap.normalizeFromBounds(boundaries, guess, init);
            if(DEBUG){
                for(int i = 0; i < 6; i++){
                    assert(init[i] >= 0.0 && init[i] <= 1.0);
                }
            }
            
            final Tuple<Double, double[]> opt = bobyqa.doOptimize(init, adap);

            final double currBestE = opt.getObject1();
            
            if (currBestE >= FixedValues.NONCONVERGEDENERGY) {
                System.out.println("INFO: Could not find any place that was collision free and had a reasonable energy. Complaining.");
                return new Tuple<>(FixedValues.NONCONVERGEDENERGY, geom.copy());
            }

            // work is already setup properly, except for the energy :-)
            final Geometry res = adap.getCurrentGeom().copy();
            res.setFitness(currBestE);
            
            return new Tuple<>(res.getFitness(),res);
        }
        
        static class GDMBOBYQAConfig implements Serializable {
        
            private static final long serialVersionUID = (long) 20140525;
        
            // XXX this should ideally be configurable
            final double initialTrust = 1E-1;
            final double stoppingTrust = 1E-5;
            final int maxIter = 500;
            final int interPoints = 13; // 2N+1 w/ 6 degrees of freedom
        }
    }
    
    private static class GDMAdaptor extends AbstractBOBYQAMethod {
        
        private static final long serialVersionUID = (long) 20140526;
        private Geometry work;
        private final int leastConn;
        private final GenericLocOpt<Molecule,Geometry> locopt;
        private final GenericFitnessBackend<Molecule,Geometry> back;
        private final CollisionDetectionEngine coll;
        private final double[] point;
        private final double[][] bounds;
        private final double blowCD;
        private final boolean fullyRelaxed;
        private final boolean eulerOnly;
        private int cIter = 0;
        
        GDMAdaptor(final Geometry work, final int leastConn,
                final GenericLocOpt<Molecule,Geometry> locopt, final CollisionDetectionEngine coll,
                final double[] point, final double[][] bounds, final double blowCD,
                final boolean fullyRelaxed, final boolean eulerOnly){
            super();
            this.coll = coll;
            this.locopt = locopt;
            this.back = locopt.getFitnessBackend();
            assert(back != null);
            this.leastConn = leastConn;
            this.work = work;
            this.bounds = bounds;
            this.point = point;
            this.blowCD = blowCD;
            this.eulerOnly = eulerOnly;
            this.fullyRelaxed = fullyRelaxed;
        }

        @Override
        public double computeObjectiveValue(final double[] normalized) {
            
            cIter++;
            
            // denormalize
            denormalizeToBounds(bounds, normalized, point);
            
            // put in
            if(eulerOnly){
                assert(point.length == 3);
                work.getMoleculeAtPosition(leastConn).setOrientation(point);
            } else {
                assert(point.length == 6);
                work.setExtCoordMolecule(point, leastConn);
            }
            
            // run a quick CD
            final boolean hasColl = CollisionDetection.checkOnlyForCollision(work, blowCD);
            if(hasColl){
                if(DEBUG) {System.out.println("DEBUG: Collision found.");}
                return FixedValues.NONCONVERGEDENERGY;
            }
            
            // now try the backend
            if(fullyRelaxed){
                
                work = locopt.fitness(work, false);
                final double e = work.getFitness();
                if(DEBUG){System.out.println("DEBUG: new relaxed energy is " + e);}
                
                return e;
            } else {
                final double[] coords = back.getActiveCoordinates(work.copy());
                final double e = back.fitness(coords, cIter);
                if(DEBUG){System.out.println("DEBUG: new energy is " + e);}
                
                return e;
            }
        }

        @Override
        public boolean doesNormalize() {
            return true;
        }
        
        Geometry getCurrentGeom (){
            return work;
        }
    }
    
    static class EnergyOnlyEvaluator implements PointOptStrategy {
        
        private static final long serialVersionUID = (long) 20140528;
        
        private int evalCounter = 0;
        private double[] coordsCache;
        private List<double[][]> rotCache;
        
        EnergyOnlyEvaluator(){
        }

        @Override
        public PointOptStrategy copy() {
            return new EnergyOnlyEvaluator();
        }

        @Override
        public String getMyID() {
            return "energy only evaluator";
        }

        @Override
        public Tuple<Double, Geometry> pointOpt(final Geometry geom, final GenericLocOpt<Molecule,Geometry> locopt, final int moveMol, final int movePartner) {
            
            final GenericFitnessBackend<Molecule,Geometry> back = locopt.getFitnessBackend();
            if(back == null){
                throw new RuntimeException("Backend is null.");
            }
            
            // special clause to reduce memory pressure
            if(back instanceof FullyCartesianCoordinates){
                if(coordsCache == null){
                    
                    this.rotCache = new ArrayList<>();
                    final int noMols = geom.getNumberOfIndieParticles();
                    for(int part = 0; part < noMols; part++){
                        
                        final Molecule mol = geom.getMoleculeAtPosition(part);
                        final double[][] cache = new double[3][mol.getNumberOfAtoms()];
                        rotCache.add(cache);
                    }
                    
                    // clone() is necessary as object is potentially mutable
                    final double[] currCoords = back.getActiveCoordinates(geom.copy());
                    coordsCache = currCoords;
                    evalCounter++;
                    final double e = back.fitness(currCoords, evalCounter);
                    geom.setFitness(e);
                    
                    return new Tuple<>(e,geom);
                }
                
                // check if this is the right length
                final int noAtoms = geom.getNumberOfAtoms();
                if(coordsCache.length != 3*noAtoms){
                    // do the full dance
                    System.err.println("Something definitely changed for the Cartesian coordinates. Doing full evaluation.");
                    final Geometry evald = locopt.fitness(geom, true);
                    evalCounter++;
                    
                    return new Tuple<>(evald.getFitness(),evald);
                }
                
                // copy new coordinates in
                final int noMols = geom.getNumberOfIndieParticles();
                int offset = 0;
                for(int part = 0; part < noMols; part++){
                    
                    final Molecule mol = geom.getMoleculeAtPosition(part);
                    final int atsInMol = mol.getNumberOfAtoms();
                    final double[][] rotXYZ = rotCache.get(part);
                    if(rotXYZ[0].length != atsInMol){
                        System.err.println("Something wrong with the rotation matrix cache. Doing full evaluation.");
                        final Geometry evald = locopt.fitness(geom, true);
                        evalCounter++;
                    
                        return new Tuple<>(evald.getFitness(),evald);
                    }
                    
                    mol.giveRotTransCartesians(rotXYZ);
                    
                    // layout of the 1D array (yes, this IS ugly and makes use of knowledge on internal state): first all x, then all y, then all z
                    // XXX ugly solution
                    System.arraycopy(rotXYZ[0], 0, coordsCache,           offset, atsInMol);
                    System.arraycopy(rotXYZ[1], 0, coordsCache,   noAtoms+offset, atsInMol);
                    System.arraycopy(rotXYZ[2], 0, coordsCache, 2*noAtoms+offset, atsInMol);
                    
                    offset += atsInMol;
                }
                
                // do the evaluation
                evalCounter++;
                final double e = back.fitness(coordsCache, evalCounter);
                geom.setFitness(e);
                    
                return new Tuple<>(e,geom);
            }
            
            final Geometry evald = locopt.fitness(geom, true);
            
            evalCounter++;
            
            // clone() is NOT necessary as consuming code knows this may be the same object by API definition
            return new Tuple<>(evald.getFitness(),evald);
        }
    }    
}
