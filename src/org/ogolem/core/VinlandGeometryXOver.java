/**
Copyright (c) 2016, J. M. Dieterich and B. Hartke
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

import java.util.List;
import org.ogolem.core.CollisionInfo.Collision;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.Lottery;

/**
 * Xovers with an environment: Vinland.
 * @author Johannes Dieterich
 * @version 2016-12-18
 */
class VinlandGeometryXOver implements GenericCrossover<Molecule,Geometry> {
    
    static enum MOVEMODE {FULLRANDOM, SURFACESTYLE};
    static final double DEFAULTINCRBOHR = 0.2;
    static final int DEFAULTMAXMOVETRIES = 20; // max: 4 bohr move!
    static final int DEFAULTMAXTRIES = 200;
    
    private static final long serialVersionUID = (long) 20160403;
    private static final boolean DEBUG = false;
    
    private final GenericCrossover<Molecule,Geometry> xover;
    private final CollisionDetectionEngine cd;
    private final double blowFac;
    private final Lottery rnd;
    private final MOVEMODE mode;
    private final double movePerStep;
    private final int maxMoveTries;
    private final int maxTries;
    
    VinlandGeometryXOver(final GenericCrossover<Molecule,Geometry> xover, final double blowFac,
            final MOVEMODE mode, final double movePerStep, final int maxMoveTries, final int maxTries){
        assert(maxMoveTries > 0);
        assert(maxTries > 0);
        this.xover = xover;
        this.cd = (mode == MOVEMODE.SURFACESTYLE) ? new CollisionDetection(CollisionDetection.CDTYPE.SIMPLEPAIRWISE) : new CollisionDetection(CollisionDetection.CDTYPE.ADVANCEDPAIRWISE);
        this.blowFac = blowFac;
        this.rnd = Lottery.getInstance();
        this.mode = mode;
        this.movePerStep = movePerStep;
        this.maxMoveTries = maxMoveTries;
        this.maxTries = maxTries;
    }
    
    VinlandGeometryXOver(final VinlandGeometryXOver orig){
        this.cd = orig.cd.clone();
        this.blowFac = orig.blowFac;
        this.xover = orig.xover;
        this.rnd = Lottery.getInstance();
        this.mode = orig.mode;
        this.movePerStep = orig.movePerStep;
        this.maxMoveTries = orig.maxMoveTries;
        this.maxTries = orig.maxTries;
    }
    
    @Override
    public VinlandGeometryXOver clone() {
        return new VinlandGeometryXOver(this);
    }

    @Override
    public String getMyID() {
        return "VINLAND GEOMETRY XOVER:\n" + xover.getMyID();
    }

    @Override
    public Tuple<Geometry, Geometry> crossover(final Geometry mother, final Geometry father, final long futureID) {
        
        boolean geomHasCollision = true;
        int tries = 0;
        
        Geometry g1 = null;
        Geometry g2 = null;
        
        do {
            
            final Tuple<Geometry, Geometry> children = xover.crossover(mother, father, futureID);
            final Geometry gTmp1 = children.getObject1();
            final Geometry gTmp2 = children.getObject2();
            
            if(gTmp1 != null){
                final CartesianCoordinates c1 = gTmp1.getCartesians();
                final BondInfo bonds = gTmp1.getBondInfo();
                final CollisionInfo cI1 = cd.checkForCollision(c1, blowFac, bonds);
                final boolean hasColl1 = cI1.hasCollision();
                if(!hasColl1){
                    if(g1 == null){
                        g1 = gTmp1;
                    } else {
                        g1 = gTmp1;
                        // we still could have dissociation, but we do not care...
                        geomHasCollision = false; // got two now
                    }
                }
            }
            
            if(gTmp2 != null){
                final CartesianCoordinates c2 = gTmp2.getCartesians();
                final BondInfo bonds = gTmp2.getBondInfo();
                final CollisionInfo cI2 = cd.checkForCollision(c2, blowFac, bonds);
                final boolean hasColl2 = cI2.hasCollision();
                if(!hasColl2){
                    if(g2 == null){
                        g2 = gTmp2;
                    } else {
                        g2 = gTmp2;
                        // we still could have dissociation, but we do not care...
                        geomHasCollision = false; // got two now
                    }
                }
            }
            
            tries++;
        } while(geomHasCollision && tries < maxTries);
        
        if(tries >= maxTries){
            return new Tuple<>(g1,g2); // they may be partially null
        }
        
        if(g1 != null){
        
            boolean hasEnvColl = true;
            int lastColls = 0;
            int dirMove = -1;
            boolean positiveMove = true;
            int myTries = 0;
            do {
            
                final Tuple<CartesianCoordinates, BondInfo> tup = g1.getCartesiansAndBondsWithEnvironment();
                final CartesianCoordinates c = tup.getObject1();
                final BondInfo bonds = tup.getObject2();
                final CollisionInfo cI = cd.checkForCollision(c, blowFac, bonds);
            
                final boolean hasColl = cI.hasCollision();
            
                if(!hasColl){
                    // we still could have dissociation, but we do not care...
                    break;
                }
            
                final List<Collision> colls = cI.getCollisions();
                final int thisColls = colls.size();
               
                switch(mode){
                    case FULLRANDOM:
                        if(thisColls > lastColls){
                            // we do need a new direction
                            dirMove = rnd.nextInt(3);
                            positiveMove = rnd.nextBoolean();
                        } // if not: we are doing fine -> continue in that direction
                
                        // let's move the whole cluster a bit
                        final double move = (positiveMove) ? movePerStep : -movePerStep;
                        g1.getEnvironment().moveCluster(dirMove, move);
                        break;
                    case SURFACESTYLE:
                        g1.getEnvironment().moveCluster(2, movePerStep); // always in z, always positive
                        break;
                }
            
                lastColls = thisColls;
                tries++;
                myTries++;
            } while(hasEnvColl && tries < maxTries && myTries < maxMoveTries);
        }
        
        if(g2 != null){
        
            boolean hasEnvColl = true;
            int lastColls = 0;
            int dirMove = -1;
            boolean positiveMove = true;
            int myTries = 0;
            do {
            
                final Tuple<CartesianCoordinates, BondInfo> tup = g2.getCartesiansAndBondsWithEnvironment();
                final CartesianCoordinates c = tup.getObject1();
                final BondInfo bonds = tup.getObject2();
                final CollisionInfo cI = cd.checkForCollision(c, blowFac, bonds);
            
                final boolean hasColl = cI.hasCollision();
            
                if(!hasColl){
                    // we still could have dissociation, but we do not care...
                    break;
                }
            
                final List<Collision> colls = cI.getCollisions();
                final int thisColls = colls.size();
               
                switch(mode){
                    case FULLRANDOM:
                        if(thisColls > lastColls){
                            // we do need a new direction
                            dirMove = rnd.nextInt(3);
                            positiveMove = rnd.nextBoolean();
                        } // if not: we are doing fine -> continue in that direction
                
                        // let's move the whole cluster a bit
                        final double move = (positiveMove) ? movePerStep : -movePerStep;
                        g2.getEnvironment().moveCluster(dirMove, move);
                        break;
                    case SURFACESTYLE:
                        g2.getEnvironment().moveCluster(2, movePerStep); // always in z, always positive
                        break;
                }
            
                lastColls = thisColls;
                tries++;
                myTries++;
            } while(hasEnvColl && tries < maxTries && myTries < maxMoveTries);
        }
        
        return new Tuple<>(g1,g2);
    }
    
    @Override
    public short hasPriority() {
        return -1;
    }
}
