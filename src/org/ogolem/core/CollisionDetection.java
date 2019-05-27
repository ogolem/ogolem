/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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

import java.util.ArrayList;
import java.util.List;

/**
 * The decorator for collision detections in and between molecules.
 * @author Johannes Dieterich
 * @version 2015-07-23
 */
public class CollisionDetection implements CollisionDetectionEngine {

    private static final long serialVersionUID = (long) 20111007;
    
    public static enum CDTYPE{SIMPLEPAIRWISE, ADVANCEDPAIRWISE, SIMPLEGRID, ADVANCEDGRID};
    public static final CDTYPE DEFAULTCD = CDTYPE.SIMPLEPAIRWISE;
    
    private final CDTYPE which;
    protected CollisionDetectionEngine detection;

    public CollisionDetection(final CDTYPE whichCollisionEngine) {

        this.which = whichCollisionEngine;
        
        // XXX: dummy
        final DummyCollisionStrengthComputer dummy = new DummyCollisionStrengthComputer();
        
        switch(whichCollisionEngine){
            case SIMPLEPAIRWISE:
                detection = new AdvancedPairWise(true, dummy);
                break;
            case ADVANCEDPAIRWISE:
                detection = new AdvancedPairWise(false, dummy);
                break;
            case SIMPLEGRID:
                detection = new GridCollisionDetection(true, dummy);
                break;
            case ADVANCEDGRID:
                detection = new GridCollisionDetection(false, dummy);
                break;
            default:
                System.err.println("There is no such collision detection engine representing " + whichCollisionEngine + ".");
                System.err.println("Please fix this ASAP! Using fallback alternative AdvancedPairWise checking.");
                System.err.println("This is computationally expensive!");
                detection = new AdvancedPairWise(false,dummy);
        }
    }

    @Override
    public CollisionDetection clone(){
        return new CollisionDetection(this.which);
    }
    
    @Override
    public CollisionInfo checkForCollision(final CartesianCoordinates cartes,
        final double blowFactor, final BondInfo bonds) {
        return detection.checkForCollision(cartes, blowFactor, bonds);
    }

    @Override
    public boolean checkOnlyForCollision(final CartesianCoordinates cartesians,
            final double blowFactor, final BondInfo bonds) {
        return detection.checkOnlyForCollision(cartesians, blowFactor, bonds);
    }
    
    @Override
    public void checkForCollision(final CartesianCoordinates cartesians, final double blowFactor, final BondInfo bonds, final CollisionInfo info) {
        detection.checkForCollision(cartesians, blowFactor, bonds, info);
    }
    
    static boolean checkOnlyForCollision(final Geometry geom, final double blowFactor){
        
        final BondInfo bonds = geom.getBondInfo();        
        final int noMols = geom.getNumberOfIndieParticles();
        
        final List<double[][]> rotTransXYZ = new ArrayList<>();
        for(int part = 0; part < noMols; part++){
            final Molecule mol = geom.getMoleculeAtPosition(part);
            final int noAts = mol.getNumberOfAtoms();
            if(noAts == 1){
                rotTransXYZ.add(null);
            } else {
                rotTransXYZ.add(mol.giveRotTransCartesians());
            }
        }
        
        int off1 = 0;
        for(int part1 = 0; part1 < noMols; part1++){
            
            final Molecule mol1 = geom.getMoleculeAtPosition(part1);
            final double[][] xyz1 = rotTransXYZ.get(part1);
            final short[] atomNos1 = mol1.getAtomNumbers();
            final int noAts1 = mol1.getNumberOfAtoms();
            
            if(noAts1 > 1){
                // inside each molecule:
                for (int i = 0; i < noAts1 - 1; i++) {

                    final double rad1 = AtomicProperties.giveRadius(atomNos1[i]);

                    for (int j = i + 1; j < noAts1; j++) {
                        final double radii = blowFactor * (rad1
                            + AtomicProperties.giveRadius(atomNos1[j]));
                        final double radiiSq = radii*radii;

                        final double dx = xyz1[0][i] - xyz1[0][j];
                        final double dy = xyz1[1][i] - xyz1[1][j];
                        final double dz = xyz1[2][i] - xyz1[2][j];
                
                        final double distanceSq = dx*dx+dy*dy+dz*dz;
                        if(distanceSq <= radiiSq && !bonds.hasBond(i+off1, j+off1)){
                            return true;
                        }
                    }
                }
            }
            
            // in between molecules
            for(int part2 = part1+1; part2 < noMols; part2++){
                
                final Molecule mol2 = geom.getMoleculeAtPosition(part2);
                final double[][] xyz2 = rotTransXYZ.get(part2);
                final short[] atomNos2 = mol2.getAtomNumbers();
                final int noAts2 = mol2.getNumberOfAtoms();
                
                if(noAts1 > 1 && noAts2 > 1){
                    for (int i = 0; i < noAts1; i++) {
                        
                        final double rad1 = AtomicProperties.giveRadius(atomNos1[i]);

                        for (int j = 0; j < noAts2; j++) {
                            final double radii = blowFactor * (rad1
                                + AtomicProperties.giveRadius(atomNos2[j]));
                            final double radiiSq = radii*radii;

                            final double dx = xyz1[0][i] - xyz2[0][j];
                            final double dy = xyz1[1][i] - xyz2[1][j];
                            final double dz = xyz1[2][i] - xyz2[2][j];
                        
                            final double distanceSq = dx*dx+dy*dy+dz*dz;
                        
                            if(distanceSq <= radiiSq){
                                // by definition no bonds in between molecules!
                                return true;
                            }
                        }
                    }
                } else if(noAts1 == 1 && noAts2 == 1){
                    final double rad1 = AtomicProperties.giveRadius(atomNos1[0]);
                    final double radii = blowFactor * (rad1
                                + AtomicProperties.giveRadius(atomNos2[0]));
                    final double radiiSq = radii*radii;
                    
                    final double[] com1 = mol1.getExternalCenterOfMass();
                    final double[] com2 = mol2.getExternalCenterOfMass();
                    
                    final double dx = com1[0] - com2[0];
                    final double dy = com1[1] - com2[1];
                    final double dz = com1[2] - com2[2];
                        
                    final double distanceSq = dx*dx+dy*dy+dz*dz;
                        
                    if(distanceSq <= radiiSq){
                        // by definition no bonds in between molecules!
                        return true;
                    }
                } else if(noAts1 == 1 && noAts2 > 1){
                    final double rad1 = AtomicProperties.giveRadius(atomNos1[0]);
                    final double[] com1 = mol1.getExternalCenterOfMass();
                    
                    for (int j = 0; j < noAts2; j++) {
                        final double radii = blowFactor * (rad1
                            + AtomicProperties.giveRadius(atomNos2[j]));
                        final double radiiSq = radii*radii;

                        final double dx = com1[0] - xyz2[0][j];
                        final double dy = com1[1] - xyz2[1][j];
                        final double dz = com1[2] - xyz2[2][j];
                        
                        final double distanceSq = dx*dx+dy*dy+dz*dz;
                        
                        if(distanceSq <= radiiSq){
                            // by definition no bonds in between molecules!
                            return true;
                        }
                    }
                } else if(noAts1 > 1 && noAts2 == 1){
                    final double rad2 = AtomicProperties.giveRadius(atomNos2[0]);
                    final double[] com2 = mol2.getExternalCenterOfMass();
                    
                    for (int j = 0; j < noAts1; j++) {
                        final double radii = blowFactor * (rad2
                            + AtomicProperties.giveRadius(atomNos1[j]));
                        final double radiiSq = radii*radii;

                        final double dx = com2[0] - xyz1[0][j];
                        final double dy = com2[1] - xyz1[1][j];
                        final double dz = com2[2] - xyz1[2][j];
                        
                        final double distanceSq = dx*dx+dy*dy+dz*dz;
                        
                        if(distanceSq <= radiiSq){
                            // by definition no bonds in between molecules!
                            return true;
                        }
                    }
                }
            }
            
            off1 += noAts1;
        }
        
        return false;
    }
    
    public static CDTYPE parseType(final String type) throws Exception {
        
        if(type.equalsIgnoreCase("simplepairwise")){
            return CDTYPE.SIMPLEPAIRWISE;
        } else if(type.equalsIgnoreCase("advancedpairwise")){
            return CDTYPE.ADVANCEDPAIRWISE;
        } else if(type.equalsIgnoreCase("simplegrid")){
            System.err.println("WARNING: GRID COLLISION DETECTION ENGINE UNTESTED. NOT RECOMMENDED!");
            return CDTYPE.SIMPLEGRID;
        } else if(type.equalsIgnoreCase("advancedgrid")){
            System.err.println("WARNING: GRID COLLISION DETECTION ENGINE UNTESTED. NOT RECOMMENDED!");
            return CDTYPE.ADVANCEDGRID;
        }
        
        throw new Exception("Illegal collision detection " + type + ".");
    }
}
