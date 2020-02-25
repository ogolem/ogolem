/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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
package org.ogolem.core;

/**
 * This, in difference to the SimplePairWise collision detection engine does NOT 
 * exit on the first collision detected but goes on and gathers more information together.
 * The algorithm is in both cases a pairwise checking.
 * @author Johannes Dieterich
 * @version 2020-02-01
 */
public class AdvancedPairWise implements CollisionDetectionEngine {

    private static final long serialVersionUID = (long) 20200201;
    
    private final boolean exit;
    private final CollisionStrengthComputer comp;
    
    public AdvancedPairWise(final boolean exitOnFirst, final CollisionStrengthComputer comp){
        assert(comp != null);
        this.exit = exitOnFirst;
        this.comp = comp;
    }
    
    private AdvancedPairWise(final AdvancedPairWise orig){
        this.exit = orig.exit;
        this.comp = orig.comp.clone();
    }
    
    @Override
    public AdvancedPairWise clone(){
        return new AdvancedPairWise(this);
    }
        
    /**
     * 
     * @param cartesians A complete set of cartesian coordinates which should
     * be checked for collisions.
     * @param blowFactor Bond detection works with this blow factor.
     * @param bonds Information on existing (and therefore wanted) bonds in
     * the molecule.
     * @return info The returned {@code CollisionInfo}, an object keeping
     * information not only on collisions but also on the calculated pairwise
     * distances.
     */
    @Override
    public CollisionInfo checkForCollision(final CartesianCoordinates cartesians, final double blowFactor, final BondInfo bonds) {
        
        final CollisionInfo info = (exit) ? new SingleCollisionInfo() : new MultiCollisionInfo();
        checkForCollision(cartesians, blowFactor, bonds, info);
        
        return info;
    }
    
    @Override
    public void checkForCollision(final CartesianCoordinates cartesians, final double blowFactor, final BondInfo bonds, final CollisionInfo info) {

        if (bonds == null) {
            // apparently no one has cared to initialize this
            System.err.println("WARNING: No bond information found in AdvancedPairWise.");
            System.err.println("WARNING: Perhaps it is not initialized yet?");
            System.err.println("WARNING: Returning null'd collision information object.");
            return;
        }
        
        if (cartesians == null) {
            System.err.println("WARNING: Cartesians were null'd in AdvancedPairWise!");
            return;
        }
        
        info.resizeDistsAndClearState(cartesians.getNoOfAtoms());
        
        final int noOfAtoms = cartesians.getNoOfAtoms();
        final double[][] dists = info.getPairWiseDistances();

        final double[][] xyz = cartesians.getAllXYZCoord();
        final short[] numbers = cartesians.getAllAtomNumbers();
        
        assert(xyz != null);
        assert(numbers != null);
        assert(numbers.length >= noOfAtoms);
        assert(xyz.length == 3);
        assert(xyz[0].length >= noOfAtoms);
        assert(xyz[1].length >= noOfAtoms);
        assert(xyz[2].length >= noOfAtoms);
        
        if(!bonds.bondMatrixFast()){
            System.err.println("WARNING: The retrival of the full bond matrix is slow, reconsider to NOT use AdvancedPairWise or change implementation!");
        }
        final boolean[][] bondMat = bonds.getFullBondMatrix();
        
        // prefetch the radii in O(N) - the allocation here is not ideal
        // but caching it would make this whole object thread-unsafe
        final double[] radii = new double[noOfAtoms];
        for(int i = 0; i < noOfAtoms; i++){
            radii[i] = AtomicProperties.giveRadius(numbers[i]);
        }
        
        boolean distCompl = true;
        Outer: for (int i = 0; i < noOfAtoms - 1; i++) {
            final double rad1 = radii[i];
            final double x = xyz[0][i]; final double y = xyz[1][i]; final double z = xyz[2][i];
            for (int j = i + 1; j < noOfAtoms; j++) {
                
                final double rad2 = radii[j];
                final double radiiAdd = blowFactor * (rad1 + rad2);
                final double dX = x - xyz[0][j];
                final double dY = y - xyz[1][j];
                final double dZ = z - xyz[2][j];
                final double dist = Math.sqrt(dX * dX + dY * dY + dZ * dZ);
                dists[i][j] = dist;
                dists[j][i] = dist;
                if (dist < radiiAdd && !bondMat[i][j]) {
                    // collision
                    final double strength = comp.calculateCollisionStrength(i, j, dist, radiiAdd);
                    final boolean succ = info.reportCollision(i, j, strength);
                    if(!succ){
                        System.err.println("No success setting collision!");
                    }
                    
                    // increment the noOfCollisions AFTER setting the collision info since arrays start from 0.
                    if(exit){
                        distCompl = false;
                        break Outer;
                    }
                }
            }
        }

        // set the pairwise distances to our collision info object
        info.setPairWiseDistances(dists, distCompl);
    }

    @Override
    public boolean checkOnlyForCollision(final CartesianCoordinates cartesians, final double blowFactor,
            final BondInfo bonds) {
        
        return checkOnlyForCollision(cartesians, blowFactor, bonds, 0, cartesians.getNoOfAtoms());
    }
    
    @Override
    public boolean checkOnlyForCollision(final CartesianCoordinates cartesians,
                final double blowFactor, final BondInfo bonds, final int offset,
                final int endset){
        
        if (bonds == null) {
            // apparently no one has cared to initialize this
            System.err.println("WARNING: No bond information found in AdvancedPairWise.");
            System.err.println("WARNING: Perhaps it is not initialized yet?");
            System.err.println("WARNING: Returning null'd collision information object.");
            throw new RuntimeException("No bonds object given to AdvancedPairWise.");
        }
        
        if (cartesians == null) {
            System.err.println("WARNING: Cartesians were null'd in AdvancedPairWise!");
            throw new RuntimeException("No Cartesian coordinates given to AdvancedPairWise.");
        }
        
        final int noOfAtoms = cartesians.getNoOfAtoms();
        final double[][] xyz = cartesians.getAllXYZCoord();
        final short[] numbers = cartesians.getAllAtomNumbers();
        
        assert(xyz != null);
        assert(numbers != null);
        assert(numbers.length >= noOfAtoms);
        assert(xyz.length == 3);
        assert(xyz[0].length >= noOfAtoms);
        assert(xyz[1].length >= noOfAtoms);
        assert(xyz[2].length >= noOfAtoms);
        assert(offset >= 0);
        assert(endset <= noOfAtoms);
        
        if(!bonds.bondMatrixFast()){
            System.err.println("WARNING: The retrival of the full bond matrix is slow, reconsider to NOT use AdvancedPairWise or change implementation!");
        }
        final boolean[][] bondMat = bonds.getFullBondMatrix();
        
        // prefetch the radii in O(N) - the allocation here is not ideal
        // but caching it would make this whole object thread-unsafe
        final double[] radii = new double[noOfAtoms];
        for(int i = 0; i < noOfAtoms; i++){
            radii[i] = AtomicProperties.giveRadius(numbers[i]);
        }
        
        final int end = Math.min(noOfAtoms, endset);
        for (int i = 0; i < end; i++) { // note: this is by spec
            final double rad1 = radii[i];
            final double x = xyz[0][i]; final double y = xyz[1][i]; final double z = xyz[2][i];
            final int start = Math.max(offset, i+1);
            for (int j = start; j < noOfAtoms; j++) {                
                final double rad2 = radii[j];
                final double radiiAdd = blowFactor * (rad1 + rad2);
                final double dX = x - xyz[0][j];
                final double dY = y - xyz[1][j];
                final double dZ = z - xyz[2][j];
                final double distSq = dX * dX + dY * dY + dZ * dZ;
                if (distSq < radiiAdd*radiiAdd && !bondMat[i][j]) {
                    // collision
                    return true;
                }
            }
        }

        return false;
    }
}
