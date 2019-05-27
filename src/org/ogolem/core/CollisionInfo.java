/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2012, J. M. Dieterich
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

import java.io.Serializable;
import java.util.List;

/**
 * Interface for a value object holding all information on detected collision(s). It
 * keeps the pairwise distances.
 * @author Johannes Dieterich
 * @version 2015-07-23
 */
public interface CollisionInfo extends Serializable {
    
    public static final double DEFAULTSTRENGTH = 42;

    /**
     * Report a collision into this info object
     * @param atom1 the first atom in the collision
     * @param atom2 the second atom in the collision
     * @param strength the collision strength, can be
     * @return true if the report was successful, false otherwise
     */
    boolean reportCollision(final int atom1, final int atom2, final double strength);
    
    /**
     * Whether a collision was reported.
     * @return True if a collision was reported, false otherwise.
     */
    boolean hasCollision();
    
    /**
     * Get the number of reported collisions.
     * @return the number of reported collision.
     */
    int getNumberOfStoredCollisions();
    
    /**
     * Get all collisions that were reported.
     * @return list of reported collisions. If there are no collisions, an empty list will be returned.
     */
    List<Collision> getCollisions();
    
    /**
     * Whether or not the pairwise distances stored are complete.
     * @return True if they are all complete and usable, false otherwis.
     */
    boolean pairWiseDistsComplete();
    
    /**
     * Set the pairwise distances.
     * @param distances the pairwise distances as a full matrix
     * @param areComplete whether the pairwise distances are complete.
     * @return true if the set was successful, false otherwise
     */
    boolean setPairWiseDistances(final double[][] distances, final boolean areComplete);
    
    /**
     * Get the pairwise distances.
     * @return the pairwise distances. May be not complete. May be null if no storage of distances happened.
     */
    double[][] getPairWiseDistances();
    
    /**
     * Resize the pairwise matrix and clean state for reuse. Will NOT clean the pairwise distance matrix but only the collision status..
     * @param size the size of the pairwise matrix. Will reallocate and clean internal datastructures if different than before, otherwise will only clean.
     */
    void resizeDistsAndClearState(final int size);
    
    public static class Collision implements Serializable {
        
        private static final long serialVersionUID = (long) 20150719;
        
        private final int at1;
        private final int at2;
        private final double strength;
        
        Collision(final int at1, final int at2, final double strength){
            this.at1 = at1;
            this.at2 = at2;
            this.strength = strength;
        }
        
        public int getAtomOne(){
            return at1;
        }
        
        public int getAtomTwo(){
            return at2;
        }
        
        public double getCollisionStrength(){
            return strength;
        }
    }
}
