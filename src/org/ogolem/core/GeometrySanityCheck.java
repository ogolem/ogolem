/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

import org.ogolem.generic.GenericSanityCheck;

/**
 * Checks a given geometry for post-local-optimization-sanity.
 * @author Johannes Dieterich
 * @version 2014-09-06
 */
public final class GeometrySanityCheck implements GenericSanityCheck<Molecule,Geometry>{

    private static final long serialVersionUID = (long) 20140401;
    private static final boolean DEBUG = false;
    private final double blowBonds;
    private final boolean checkCollisions;
    private final boolean doDD;
    private final double blowDiss;
    
    public GeometrySanityCheck(final double blowBonds, final double blowDiss, final boolean checkColl,
            final boolean doDD){
        this.blowBonds = blowBonds;
        this.blowDiss = blowDiss;
        this.checkCollisions = checkColl;
        this.doDD = doDD;
    }
    
    GeometrySanityCheck(final GeometrySanityCheck orig){
        this.blowBonds = orig.blowBonds;
        this.blowDiss = orig.blowDiss;
        this.doDD = orig.doDD;
        this.checkCollisions = orig.checkCollisions;
    }
    
    @Override
    public GeometrySanityCheck clone() {
        return new GeometrySanityCheck(this);
    }

    @Override
    public boolean isSane(final Geometry individual) {
        
        final CartesianCoordinates cartes = individual.getCartesians();
        final BondInfo bonds = individual.getBondInfo();
        
        return checkSanity(cartes, bonds, blowBonds, checkCollisions, doDD, blowDiss);
    }

    
    /**
     * Checks whether the provided geometry still has all intramolecular bonds as before.
     * Scales O(N^2) with N as the number of atoms. The arrays may be of different length,
     * then the bond array decides which interactions are checked. This is e.g. of
     * interest when the cartesian still contains an environment but the bond array doesn't.
     * @param cartes A cartesian.
     * @param bonds The bonds in this cartesian.
     * @param blowBonds A blow factor.
     * @return true if all bonds still exist, false otherwise.
     */
    public static boolean checkSanity(final CartesianCoordinates cartes, final BondInfo bonds,
            final double blowBonds){
        
        return checkSanity(cartes, bonds, blowBonds, false, false, 4.0);
    }

    /**
     * Checks whether the provided geometry still has all intramolecular bonds as before.
     * Scales O(N^2) with N as the number of atoms. The arrays may be of different length,
     * then the bond array decides which interactions are checked. This is e.g. of
     * interest when the cartesian still contains an environment but the bond array doesn't.
     * @param cartes A cartesian.
     * @param bonds The bonds in this cartesian
     * @param doDD do dissociation detection?
     * @param blowBonds A blow factor.
     * @param checkCollisions If also collisions shall be checked.
     * @param blowDiss A blow factor for the dissociation check.
     * @return true if all bonds still exist, false otherwise.
     */
    public static boolean checkSanity(final CartesianCoordinates cartes, final BondInfo bonds,
            final double blowBonds, final boolean checkCollisions, final boolean doDD,
            final double blowDiss){
        
        if(DEBUG){
            System.out.println("DEBUG: Sanity check working on the following cartesian...");
            final String[] sa = cartes.createPrintableCartesians();
            for(final String s : sa){
                System.out.println(s);
            }
        }

        // see above why we use the bond lengths
        final int noOfAtoms = cartes.getNoOfAtoms();
        final short[] nos = cartes.getAllAtomNumbers();
        final double[][] xyz = cartes.getAllXYZCoord();
        final boolean[][] adjacency = (!doDD) ? null : new boolean[noOfAtoms][noOfAtoms];
        
        for(int i = 0; i < noOfAtoms-1; i++){
            final double rad1 = AtomicProperties.giveRadius(nos[i]);
            for(int j = i+1; j < noOfAtoms; j++){
                final double rad2 = AtomicProperties.giveRadius(nos[j]);
                final double radii = blowBonds*(rad1+rad2);
                final double radiiSq = radii*radii;
                
                final double distX = xyz[0][i]-xyz[0][j];
                final double distY = xyz[1][i]-xyz[1][j];
                final double distZ = xyz[2][i]-xyz[2][j];
                
                final double distSq = distX*distX + distY*distY + distZ*distZ;
                final boolean bond = bonds.hasBond(i, j);
                if(bond && distSq > radiiSq){
                    // there is no bond where there should be one
                    if(DEBUG) {
                        System.out.println("DEBUG: Backing out because of no bond. " + i + "\t" + j);
                        System.out.println("DEBUG: dist " + Math.sqrt(distSq) + " bigger than " + radii);
                    }
                    return false;
                } else if(checkCollisions && !bond && distSq < radiiSq){
                    // there is a bond where there should be none
                    if(DEBUG) {
                        System.out.println("DEBUG: Backing out because of a bond. " + i + "\t" + j);
                        System.out.println("DEBUG: dist " + Math.sqrt(distSq) + " smaller than " + radii);
                    }
                    return false;
                }
                
                if(doDD){
                    final double summedRadii = (rad1 + rad2) * blowDiss;
                    final double summedRadiiSq = summedRadii*summedRadii;
                    if (summedRadiiSq >= distSq) {
                        adjacency[i][j] = true;
                        adjacency[j][i] = true;
                    } else{
                        // XXX not really needed...
                        adjacency[i][j] = false;
                        adjacency[j][i] = false;
                    }
                }
            }
        }
        
        if (!doDD) {
            return true;
        }
        
        // check for DD using DFS
        final boolean isConnected = DistanceCalc.dfsReachability(adjacency);
        if(DEBUG){
            if(!isConnected){System.out.println("DEBUG: Cluster is not connected.");}
        }
        
        return isConnected;
    }
}
