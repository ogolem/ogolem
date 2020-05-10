/**
Copyright (c) 2012-2013, J. M. Dieterich
              2016-2017, J. M. Dieterich and B. Hartke
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
package org.ogolem.heat;

import java.io.Serializable;
import java.util.List;
import org.ogolem.core.*;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericInitializer;
import org.ogolem.heat.LocalHeatPulses.Configuration.CHOOSEMODE;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * Implements local heat pulses (potentially slightly changed)
 * as invented by A. Moebius et al.
 * @author Johannes Dieterich
 * @version 2017-06-27
 */
public class LocalHeatPulses implements Serializable {
    
    private static final long serialVersionUID = (long) 20170625;
    private static final boolean DEBUG = false;

    public static Geometry cycle(final Geometry initial, final GenericFitnessFunction<Molecule,Geometry> fitness,
            final CollisionDetectionEngine cd, final double blowCD, final double blowDD,
            final BondInfo bonds, LocalHeatPulses.Configuration config, final double acceptableFitness,
            final GenericInitializer<Molecule,Geometry> initer){
        
        final Lottery r = Lottery.getInstance();        
        Geometry work = new Geometry(initial);
        
        double radCluster = 0.0;
        if(config.chooseMode == CHOOSEMODE.INSPHERE || config.chooseMode == CHOOSEMODE.ONSPHERE){
            // figure out the radius of the work geometry
            final CartesianCoordinates c = work.getCartesians();
            c.moveCoordsToCOM();
            final double[][] xyz = c.getAllXYZCoord();
            for(int i = 0; i < c.getNoOfAtoms(); i++){
                final double dist = Math.sqrt(xyz[0][i]*xyz[0][i] + xyz[1][i]*xyz[1][i] + xyz[2][i]*xyz[2][i]);
                radCluster = Math.max(radCluster,dist);
            }
        }
        
        final double[] mvVec = new double[3];
        Geometry prev;
        Geometry best = null;
        double bestE = Double.MAX_VALUE;
        double currAmpl = config.startAmplitude;
        double currTemp = config.temperature;
        double currEuler = config.startEulerStrength;
        int noProgress = 0;
        for(int iter = 0; iter < config.iters; iter++){
            
            if(bestE <= acceptableFitness){
                System.out.println("INFO: We found the acceptable fitness " + acceptableFitness);
                System.out.println("INFO: Number of heat pulses needed to get here: " + iter);
                break;
            }
            
            prev = new Geometry(work);
            
	    if(iter > 0){
	      final CartesianCoordinates c = work.getCartesians();
              
              if(config.chooseMode == CHOOSEMODE.INSPHERE || config.chooseMode == CHOOSEMODE.ONSPHERE || config.chooseMode == CHOOSEMODE.INCENTER){
                  double[] pulsePoint = new double[3];
                  if(null != config.chooseMode)switch (config.chooseMode) {
                      case INSPHERE:
                          pulsePoint = new double[3];
                          RandomUtils.randomPointInSphere(radCluster, pulsePoint);
                          if(DEBUG){System.out.println("DEBUG: Pulsepoint is (in sphere) " + pulsePoint[0] + " " + pulsePoint[1] + " " + pulsePoint[2]);}
                          break;
                      case ONSPHERE:
                          pulsePoint = LocalHeatPulses.findPulsePointOnSurface(c);
                          if(DEBUG){System.out.println("DEBUG: Pulsepoint is (on sphere) " + pulsePoint[0] + " " + pulsePoint[1] + " " + pulsePoint[2]);}
                          break;
                      case INCENTER:
                          pulsePoint = new double[3];
                          break;
                      default:
                          System.err.println("ERROR: In local heat pulses. How did you end up here?!");
                          pulsePoint = new double[3];
                          break;
                  }
                  
                  // now figure out where and how much to move
                  if(config.moveMode == Configuration.MOVEMODE.MOVECOMS){
                      int countMoved = 0;
                      for(int com = 0; com < work.getNumberOfIndieParticles(); com++){
                          final double[] comP = work.getCOM(com);
                          final double dX = comP[0]-pulsePoint[0];
                          final double dY = comP[1]-pulsePoint[1];
                          final double dZ = comP[2]-pulsePoint[2];
                          
                          final double gauss = currAmpl*Math.exp(-(dX*dX/config.sigmaX + dY*dY/config.sigmaY + dZ*dZ/config.sigmaZ));
                          if(DEBUG){System.out.println("DEBUG: Gauss is evaluated to " + gauss + " for COM " + com);}
                          if(gauss > 1E-2){
                              // only then we actually go through the pain of heating this
                              RandomUtils.randomVector(mvVec, gauss);
                              final double rotStrength = Math.min(currEuler*gauss, 1.0); // need to make sure that this is <= 1.0, since this is assert-enforced in RandomUtils used by move!
                              move(work, com, mvVec, rotStrength); // scale not just COM-movement but also rotation by Gauss-centrality!
                              countMoved++;
                          }
                      }
                      if(config.printVerbose){
                          System.out.println("DEBUG: Moved " + countMoved + " COMs from spheric heat pulse.");
                      }
                  } else if(config.moveMode == Configuration.MOVEMODE.MOVECARTESIAN){
                      final double[][] xyz = c.getAllXYZCoord();
                      int countMoved = 0;
                      for(int atom = 0; atom < c.getNoOfAtoms(); atom++){
                          final double dX = xyz[0][atom]-pulsePoint[0];
                          final double dY = xyz[1][atom]-pulsePoint[1];
                          final double dZ = xyz[2][atom]-pulsePoint[2];
                          
                          final double gauss = currAmpl*Math.exp(-(dX*dX/config.sigmaX + dY*dY/config.sigmaY + dZ*dZ/config.sigmaZ));
                          if(DEBUG){System.out.println("DEBUG: Gauss is evaluated to " + gauss + " for atom " + atom);}
                          if(gauss > 1E-2){
                              // only then we actually go through the pain of heating this
                              RandomUtils.randomVector(mvVec, gauss);
                              move(c, atom, mvVec);
                              countMoved++;
                          }
                      }
                      if(config.printVerbose){
                          System.out.println("DEBUG: Moved " + countMoved + " atoms from spheric heat pulse.");
                      }
                  }
              } else{
                    final int noMove = disturbHowManyAtoms(config.chooseMode, (config.moveMode == Configuration.MOVEMODE.MOVECARTESIAN) ? c.getNoOfAtoms() : c.getNoOfMolecules(), r);
                    final List<Integer> mvs = RandomUtils.listOfPoints(noMove, (config.moveMode == Configuration.MOVEMODE.MOVECARTESIAN) ? c.getNoOfAtoms() : c.getNoOfMolecules());
            
                    for (final int mv : mvs) {

                        RandomUtils.randomVector(mvVec, currAmpl);
                        // move COM or atom
                        if (config.moveMode == Configuration.MOVEMODE.MOVECOMS) {
                            move(work, mv, mvVec, currEuler);
                        } else if(config.moveMode == Configuration.MOVEMODE.MOVECARTESIAN){
                            move(c, mv, mvVec);
                        }
                    }
                  
                    if(config.printVerbose){
                        System.out.println("DEBUG: Moved " + noMove + " from distributed heat pulse.");
                    }
              }

                if (config.moveMode == Configuration.MOVEMODE.MOVECARTESIAN) {
                    // translate c back to work as this is needed because of moved atoms
                    final boolean hasEnv = work.containsEnvironment();
                    final long id = work.getID();
                    work = CoordTranslation.cartesianToGeometry(c, c.getNoOfMolecules(), c.getAllAtomsPerMol(),
                            work.getAllFlexies(), work.getExplicitDoFs(), work.getAllConstraints(hasEnv), work.getAllConstraintsXYZ(hasEnv),
                            c.getAllAtomTypes(), work.getBondInfo());
                    work.setID(id); // the translation looses the ID
                }

	      // check CD/DD
	      if(config.doCDDD){
                if(checkCDandDD(work, cd, blowCD, blowDD, bonds)){
                    if(DEBUG){System.out.println("DEBUG: Checked CD/DD, insane geometry, resetting to previous.");}
                    work = prev;
                    continue;
                }
	      }
	    }
            
            if(DEBUG){
                System.err.println("DEBUG: Geometry " + work.getID() + " at iter " + iter + " BEFORE local optimization");
                for(final String s : work.makePrintableAbsoluteCoord(true)){
                    System.out.println(s);
                }
            }
            
            // locally optimize
            work = fitness.fitness(work, false); //XXX should this be tunable?
            
            if(DEBUG){
                System.err.println("DEBUG: Geometry " + work.getID() + " at iter " + iter + " AFTER local optimization");
                for(final String s : work.makePrintableAbsoluteCoord(true)){
                    System.out.println(s);
                }
            }
            
            if(config.useMetropolis){
                if(work.getFitness() > prev.getFitness()){
                    final double beta = 1.0 / (8.31447215 / 1000.0 * currTemp); // for sake of simplicity in kj/mol and K
                    final double boltz = Math.exp(-beta * (work.getFitness() - prev.getFitness()) * Constants.HARTREETOKJ);
                    final double trial = r.nextDouble();

                    if (boltz < trial) {
                        // reject
                        if(DEBUG) {System.out.println("DEBUG: Rejecting geometry after Metropolis. This fitness: " + work.getFitness() + " vs " + prev.getFitness());}                        
                        work = prev;
                        continue;
                    } else {
                        if(DEBUG) {System.out.println("DEBUG: Accepting geometry after Metropolis. This fitness: " + work.getFitness() + " vs " + prev.getFitness());}
                    }
                } else {
                    if(DEBUG) {System.out.println("DEBUG: Accepting geometry. This fitness: " + work.getFitness() + " vs " + prev.getFitness());}
                }
            } else{
                if(work.getFitness() > prev.getFitness()){
                    if(DEBUG) {System.out.println("DEBUG: Rejecting geometry w/o Metropolis. This fitness: " + work.getFitness() + " vs " + prev.getFitness());}
                    work = prev;
                    continue;
                } else {
                    if(DEBUG) {System.out.println("DEBUG: Accepting geometry w/o Metropolis. This fitness: " + work.getFitness() + " vs " + prev.getFitness());}
                }
            }
            
            // are we better?
            if(work.getFitness() < bestE){
		if(DEBUG) {
                    System.out.println("DEBUG: At iter " + iter + " new best " + work.getFitness());
                }
                bestE = work.getFitness();
                best = new Geometry(work);
                noProgress = 0;
            } else {
                noProgress++;
            }
            
            if(config.resetAfterNoProgress && noProgress >= config.resetNoProgressIters){
                if(!config.resetToRandom){
                    work = new Geometry(initial);
                } else {
                    work = initer.initializeOnly(work, initial.getID());
                    work.setFather(initial.getFatherID());
                    work.setMother(initial.getMotherID());
                    work.setFitness(FixedValues.NONCONVERGEDENERGY);
                }
                noProgress = 0;
            }
            
            if(iter > config.eqIter) {
                if(DEBUG) {System.out.println("DEBUG: Rescaling amplitudes, temperatures and Euler moves: " + currAmpl + " " + currTemp + " " + currEuler);}
                currAmpl *= config.scaleFac;
                currTemp *= config.scaleFac;
                currEuler *= config.scaleFac;
            }
        }
        
        if(DEBUG) {System.out.println("DEBUG: Returning " + bestE);}
                
        return (best == null) ? work : best;
    }
    
    private static void move(final Geometry g, final int which, final double[] move,
            final double eulerStrength){
        
        final double[] com = g.getCOM(which);
        com[0] += move[0];
        com[1] += move[1];
        com[2] += move[2];
    
        RandomUtils.randomEulerIncrements(g.getEulers(which), eulerStrength);
    }
    
    private static void move(final CartesianCoordinates c, final int which, final double[] move){
        
        final double[][] xyz = c.getAllXYZCoord();
        xyz[0][which] += move[0];
        xyz[1][which] += move[1];
        xyz[2][which] += move[2];
    }
    
    private static int disturbHowManyAtoms(final Configuration.CHOOSEMODE mode, final int noAtoms,
            final Lottery r){
        
        switch(mode){
            case PICK5: return 5;
            case PERCENT10: return (int) Math.round(0.1*noAtoms);
            case UPTO5: return r.nextInt(5);
            case UPTO10PERCENT: return (int) Math.round(r.nextDouble()*0.1*noAtoms);
            default: return 5;
        }
    }
        
    private static boolean checkCDandDD(final Geometry g, final CollisionDetectionEngine cd,
            final double blowCD, final double blowDD, final BondInfo bonds){
        
        final CartesianCoordinates c = g.getCartesians();
        
        final CollisionInfo ci = cd.checkForCollision(c, blowCD, bonds);
        if (ci.hasCollision()) {
            return true;
        }
        
        return DissociationDetection.checkForDissociation(ci.getPairWiseDistances(),c.getAllAtomTypes(), c.getAllAtomNumbers(),
                blowDD, DissociationDetection.DEFAULTDD);
    }
    
    private static double[] findPulsePointOnSurface(final CartesianCoordinates c){
        
        final double angThresh = 0.1745; // approx 10 deg
        
        final double[] full = new double[3];
        RandomUtils.randomSphericals(full, 1.0);
        
        final double[][] xyz = c.getAllXYZCoord();
        final double[][] sphers = CoordTranslation.cartesianToSphericalCoord(xyz);
        
        double maxR = 0.0;
        double bestPhi = Double.MAX_VALUE;
        double bestOmega = Double.MAX_VALUE;
        for(int i = 0; i < c.getNoOfAtoms(); i++){
            
            if(Math.abs(sphers[1][i]-full[1]) < Math.abs(bestPhi-full[1])
                    && Math.abs(sphers[1][i]-bestPhi) > angThresh
                    && Math.abs(sphers[2][i]-full[2]) < Math.abs(bestOmega-full[2])
                    && Math.abs(sphers[2][i]-bestOmega) > angThresh
                    && sphers[0][i] > maxR){
                maxR = sphers[0][i];
                bestPhi = sphers[1][i];
                bestOmega = sphers[2][i];
            }
        }
        
        full[0] = maxR;
        
        final double[] pulse = new double[3];
        CoordTranslation.sphericalToCartesianCoord(full, pulse);
        
        return pulse;
    }
        
    public static class Configuration implements Serializable {
        
        private static final long serialVersionUID = (long) 20170313;
        
        public static enum CHOOSEMODE{PICK5, PERCENT10, UPTO5, UPTO10PERCENT, INCENTER, ONSPHERE, INSPHERE};
        
        public static enum MOVEMODE{MOVECARTESIAN, MOVECOMS};
        
        public boolean doCDDD = false;
        public CHOOSEMODE chooseMode = CHOOSEMODE.PICK5;
        public MOVEMODE moveMode = MOVEMODE.MOVECOMS;
        public double startAmplitude = 2.0*org.ogolem.core.Constants.ANGTOBOHR;
        public int iters = 10;
        public int eqIter = 3;
        public double scaleFac = 0.9;
        public boolean useMetropolis = true;
        public double temperature = 50;
        public double sigmaX = 0.1;
        public double sigmaY = 0.1;
        public double sigmaZ = 0.1;
        public double startEulerStrength = 0.5;
        public boolean printVerbose = false;
        public boolean resetAfterNoProgress = false;
        public int resetNoProgressIters = 1000;
        public boolean resetToRandom = false;
        
        public Configuration(){}
        
        public Configuration(final Configuration orig){
            this.chooseMode = orig.chooseMode;
            this.doCDDD = orig.doCDDD;
            this.eqIter = orig.eqIter;
            this.iters = orig.iters;
            this.moveMode = orig.moveMode;
            this.scaleFac = orig.scaleFac;
            this.startAmplitude = orig.startAmplitude;
            this.useMetropolis = orig.useMetropolis;
            this.temperature = orig.temperature;
            this.sigmaX = orig.sigmaX;
            this.sigmaY = orig.sigmaY;
            this.sigmaZ = orig.sigmaZ;
            this.startEulerStrength = orig.startEulerStrength;
            this.printVerbose = orig.printVerbose;
            this.resetAfterNoProgress = orig.resetAfterNoProgress;
            this.resetNoProgressIters = orig.resetNoProgressIters;
            this.resetToRandom = orig.resetToRandom;
        }
        
        public String printConfig(){
            
            String config = "";
            config += "choose mode:       " + chooseMode.name();
            config += "\ndo CD/DD:        " + doCDDD;
            config += "\neq iters:        " + eqIter;
            config += "\niters:           " + iters;
            config += "\nmove mode:       " + moveMode.name();
            config += "\nscale factor:    " + scaleFac;
            config += "\nstart amplitude: " + startAmplitude;
            config += "\nuse Metropolis:  " + useMetropolis;
            config += "\ntemperature:     " + temperature;
            config += "\n sigma X:        " + sigmaX;
            config += "\n sigma Y:        " + sigmaY;
            config += "\n sigma Z:        " + sigmaZ;
            config += "\n print verbose:  " + printVerbose;
            config += "\n reset for no progress: " + resetAfterNoProgress;
            config += "\n reset after X iterations (if applicable): " + resetNoProgressIters;
            config += "\n reset from random init and not initial structure (if applicable): " + resetToRandom;
            
            return config;
        }
    }
}
