/**
Copyright (c) 2014, J. M. Dieterich
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
 * An interface for rigid body backends.
 * @author Johannes Dieterich
 * @version 2014-12-19
 */
public interface RigidBodyBackend extends Serializable, Cloneable {
    
    /**
     * Clones this backend.
     * @return a clone, i.e., a deep copy.
     */
    RigidBodyBackend clone();
    
    /**
     * Name myself.
     * @return this method's ID as a nice String.
     */
    String getMethodID();
    
    /**
     * Checks if this implementation is suitable for the system under study. May
     * set some internal state unique to this geometry!
     * @param ref the reference geometry, must not be changed!
     * @return true if this backend is happy with ALL aspects of this geometry, false otherwise
     */
    boolean suitableForThisBackend(final Geometry ref);
    
    /**
     * Adjusts a molecule to put it in line with the assumptions of the force field behind it. May NOT deliberately change position or orientation beyond mere sanitization (or any other
     * fields in the structure).
     * @param mol The molecule which we want to rigidify.
     */
    void rigidify(final Molecule mol);
    
    /**
     * Adjusts a CartesianCoordinates set for this particular backend. May add or delete atoms, may alter coordinates, symbols, ...
     * @param ref the reference geometry
     * @param unadjusted the unadjusted coordinates
     * @param molID the molecule ID with respect to the Geometry object
     * @return a CartesianCoordinates object (may be the same unadjusted object) with a fitting coordinate representation for this backend
     */
    CartesianCoordinates adjustCartesians(final Geometry ref, final CartesianCoordinates unadjusted, final int molID);
    
    /**
     * Calculates the energy of this geometry (as a rigid body based off the 
     * molecule units and their orientations and translations). Must not change
     * the geometry!
     * @param ref the reference geometry
     * @param cartes a list of Cartesian coordinate objects that were pre-rotated, i.e., the respective Euler angle rotations where applied
     * @param coms the COM translations of all the Cartesian coordinates
     * @param counter a unique counter per evaluation w.r.t. this ref individual
     * @return the energy based on the method in use
     */
    double energy(final Geometry ref, final List<CartesianCoordinates> cartes,
            final double[][] coms, final int counter);
    
    /**
     * Calculates the energy and gradient of this geometry (as a rigid body based off the 
     * molecule units and their orientations and translations).
     * @param ref the reference geometry
     * @param cartes a list of Cartesian coordinate objects that were pre-rotated, i.e., the respective Euler angle rotations where applied
     * @param coms the COM translations of all the Cartesian coordinates
     * @param gradient the gradient for the external coordinates of the molecules in this geometry. Will be changed on return, does not need to be set to zero prior to calling. Contains: per molecule a double array with 3 Cartesians per atom as defined in the adjusted Cartesian set.
     * @param counter a unique counter per evaluation w.r.t. this ref individual
     * @return the energy based on the method in use
     */
    double gradient(final Geometry ref, final List<CartesianCoordinates> cartes,
            final double[][] coms, final List<double[][]> gradient, final int counter);
}
