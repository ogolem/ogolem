/**
Copyright (c) 2014-2015, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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

import static org.ogolem.core.DirMutOptStrategies.PointOptStrategy;
import static org.ogolem.core.DirMutPointProviders.PointProvider;
import org.ogolem.generic.GenericFitnessBackend;
import org.ogolem.generic.GenericLocOpt;
import org.ogolem.generic.GenericMutation;
import org.ogolem.helpers.Tuple;

/**
 * A directed mutation that is somewhat user-modifiable through plug'n'play.
 * Essentially a stripped down version of the AdvGDM so far. Eventually, this will
 * though implement AdvGDM.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class PluggableDirMut implements GenericMutation<Molecule,Geometry>{
    
    private static final long serialVersionUID = (long) 20200429;

    private final boolean DEBUG;
    private final double blowBonds;
    private final CollisionDetectionEngine collDetect;
    private final double blowColl;
    private final GenericLocOpt<Molecule,Geometry> locopt;
    private final GenericFitnessBackend<Molecule,Geometry> back;
    private final boolean doLocOpt;
    private final boolean doCDFirst;
    private final double realCD;
    private final double realDD;
    private final PointProvider pointProv;
    private final PointOptStrategy pointOpt;
    private final boolean fullyRelaxed;
    private final boolean markMovedMolsUnmovable;
    private final boolean doCDForEveryTrial;
    private final int noMoved;
    
    /**
     * A pluggable dirmut constructor.
     * @param config the GDM configuration
     * @param collDetect the collision detection engine to be used
     * @param locopt the local optimization engine (including a backend!) to be used
     * @param realCDBlow the CD blow factor used in the rest of the program, for debugging purposes
     * @param realDDBlow the DD blow factor used in the rest of the program, for debugging purposes
     * @param doLocOptFirst do a local optimization first, prior to the GDM
     * @param doCDFirst do a CD first, AFTER the local optimization (if defined) and the GDM
     * @param provider the point provider (e.g. a grid)
     * @param strategy the strategy what to do on each grid point (e.g. an optimization or only a rigid body opt or...)
     */
    public PluggableDirMut(final AdvancedGraphBasedDirMut.GDMConfiguration config, final CollisionDetectionEngine collDetect,
            final GenericLocOpt<Molecule,Geometry> locopt, final double realCDBlow, final double realDDBlow,
            final boolean doLocOptFirst, final boolean doCDFirst, final PointProvider provider,
            final PointOptStrategy strategy){
        this.pointProv = provider;
        this.pointOpt = strategy;
        this.DEBUG = (GlobalConfig.DEBUGLEVEL > 1);
        this.blowBonds = config.blowBonds;
        this.collDetect = collDetect;
        this.blowColl = config.blowCD;
        this.doLocOpt = doLocOptFirst;
        this.doCDFirst = doCDFirst;
        if(locopt.getBackend() == null){
            throw new RuntimeException("GDM must have a non-null Backend from Newton!");
        }
        this.locopt = locopt;
        this.back = locopt.getBackend();
        this.realCD = realCDBlow;
        this.realDD = realDDBlow;
        this.fullyRelaxed = config.fullOptEachMove;
        this.markMovedMolsUnmovable = config.markMovedMolsUnmoveable;
        this.noMoved = config.noOfMolsToMove;
        this.doCDForEveryTrial = config.doCDEveryTrial;
        
        assert(blowBonds > 0.0);
        assert(blowColl > 0.0);
    }
    
    private PluggableDirMut(final PluggableDirMut orig){
        this.DEBUG = (GlobalConfig.DEBUGLEVEL > 1);
        this.blowBonds = orig.blowBonds;
        this.collDetect = orig.collDetect.clone();
        this.blowColl = orig.blowColl;
        this.back = orig.back.copy();
        this.locopt = orig.locopt.copy();
        this.realCD = orig.realCD;
        this.realDD = orig.realDD;
        this.doLocOpt = orig.doLocOpt;
        this.fullyRelaxed = orig.fullyRelaxed;
        this.doCDFirst = orig.doCDFirst;
        this.pointProv = orig.pointProv.clone();
        this.pointOpt = orig.pointOpt.copy();
        this.markMovedMolsUnmovable = orig.markMovedMolsUnmovable;
        this.noMoved = orig.noMoved;
        this.doCDForEveryTrial = orig.doCDForEveryTrial;
    }
    
    @Override
    public GenericMutation<Molecule, Geometry> clone() {
        return new PluggableDirMut(this);
    }

    @Override
    public String getMyID() {
        
        return "pluggable directed mutation: \n\t point provider: " + pointProv.getMyID()
                + "\n\t opt strategy? " + pointOpt.getMyID()
                + "\n\t blow factor bonds " + blowBonds + "\n\t blow factor collision " + blowColl 
                + "\n\t CD first? " + doCDFirst + "\n\t opt first? " + doLocOpt
                + "\n\t number of moved mols: " + noMoved + "\n\t mark moved mols unmoveable? " + markMovedMolsUnmovable
                + "\n\t fully relaxed? " + fullyRelaxed + "\n\t do CD for every trial? " + doCDForEveryTrial;
    }

    @Override
    public Geometry mutate(final Geometry geom) {
        
        if(DEBUG){System.out.println("DEBUG: Running our pluggable directed mutation.");}
                
        Geometry work = new Geometry(geom);
        
        if(doLocOpt){
            if(DEBUG){System.out.println("DEBUG: doing a local optimization first.");}
            work = locopt.fitness(work, false);
        }
        
        final BondInfo bonds = work.getBondInfo();
        
        // get the energy before we manipulate
        final double[] currCoords = back.getActiveCoordinates(work.copy());
        final double eBefore = back.fitness(currCoords, 42);
        final double[] eparts = new double[work.getNumberOfIndieParticles()];
        
        for(int moveCounter = 0; moveCounter < noMoved; moveCounter++){
            
            if(fullyRelaxed){
                if(DEBUG){System.out.println("DEBUG: Optimizing at step " + moveCounter);}
                work = locopt.fitness(work, false);
                if(DEBUG){System.out.println("DEBUG: Post optimization energy: " + work.getFitness());}
            }
                        
            // hand this info over to the point provider in exchange for a list of trial points
            Geometry bestFoundGeom = null;
            double bestFoundFitness = Double.MAX_VALUE;
            while(pointProv.hasNextPoint()){
                if(DEBUG){System.out.println("DEBUG: Working on new trial point.");}
                
                final int moved = pointProv.putNextPointIn(work, null, -1, -1);
                if(doCDForEveryTrial){
                    final boolean hasColl = collDetect.checkOnlyForCollision(work.getCartesians(), blowColl, bonds);
                    if(hasColl){
                        if(DEBUG){System.out.println("DEBUG: Found a collision for this trial point.");}
                        continue;
                    }
                }
                
                // now do "something" whatever the optimization strategy is, really
                final Tuple<Double,Geometry> res = pointOpt.pointOpt(work, locopt, 
                        moved, -1);
                if(bestFoundGeom == null){
                    if(DEBUG){System.out.println("DEBUG: Unconditionally making this best found.");}
                    bestFoundGeom = res.getObject2().copy();
                    bestFoundFitness = res.getObject1();
                    continue;
                }
                if(res.getObject1() < bestFoundFitness){
                    if(DEBUG){System.out.println("DEBUG: Conditionally making this best found.");}
                    bestFoundGeom = res.getObject2().copy();
                    bestFoundFitness = res.getObject1();
                }
            }
            
            pointProv.reset(); // reset so that we can get a new set of trial points :-)
            
            if(bestFoundGeom != null){
            
                // well, just continue with the best we found :-)
                final Geometry mutationResult = bestFoundGeom.copy();
                mutationResult.setFitness(bestFoundFitness);
            
                work = mutationResult;
            }
            // continue the loop over molecules :-)
        }
        
        if (work.getFitness() >= (eBefore-1E-8)) {
            System.out.println("INFO: Pluggable DM tried hard but could not find a substantially better configuration than "
                    + eBefore + " vs best found " + work.getFitness() + " continuing but complaining.");
            /*
             * the usual BXH disclaimer applies that this makes sense anyways because it adds new genome etc pp.
             */
        } else {
            final double percImprove = (eBefore-work.getFitness())*100/eBefore;
            System.out.println("INFO: Advanced GDM did find a better configuration: " + String.format("%6.2f", percImprove) + " % better.");
        }
        
        return work;
    }
}
