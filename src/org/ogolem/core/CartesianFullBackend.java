/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
              2017, J. M. Dieterich and B. Hartke
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

/**
 * This interface defines what all the various energy and gradient delivering
 * backbones need to know.
 * ATTENTION: All input coordinates are actually 1D! Therefore chopping them in
 * the appropriate spots is needed to make them 3D.
 * @author Johannes Dieterich
 * @version 2017-03-03
 */
public interface CartesianFullBackend extends Cloneable, Serializable{

    CartesianFullBackend clone();

    String getMethodID();

    /**
     * A energy and gradient evaluation.
     * @param lID the ID of the individual used for this calculation, used for making e.g. file names unique
     * @param iIteration the iteration ID for this calculation, again used for uniqueness
     * @param xyz1D the Cartesian coordinates as a 1D array. First all x, then all y, then all z coordinates.
     * @param saAtomTypes the atom types as strings
     * @param atomNos the atom numbers associated with the atom types
     * @param atsPerMol array of length molecules, each entry defining how many atoms are in this molecule
     * @param energyparts an array for energyparts of length molecules, may be used to store an energy decomposition in
     * @param iNoOfAtoms the total number of atoms
     * @param faCharges the charges as floats, of length number of atoms
     * @param spins the spins as integers, of length number of atoms
     * @param bonds the bonds in this "atom assembly" (i.e., geometry)
     * @param grad the gradient, changed on exit with energy and gradient. Must be pre-allocated and of right size. May be filled with rubbish.
     * @param hasRigidEnvironment whether the coordinate set contains a rigid environment (per definition the last "molecule" whose internal energy contribution can be ignored (if the method allows for that).
     */
    void gradientCalculation(long lID, int iIteration, double[] xyz1D, String[] saAtomTypes,
            short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges,
            short[] spins, final BondInfo bonds, final Gradient grad, final boolean hasRigidEnvironment);

    /**
     * An energy evaluation
     * @param lID the ID of the individual used for this calculation, used for making e.g. file names unique
     * @param iIteration the iteration ID for this calculation, again used for uniqueness
     * @param xyz1D the Cartesian coordinates as a 1D array. First all x, then all y, then all z coordinates.
     * @param saAtomTypes the atom types as strings
     * @param atomNos the atom numbers associated with the atom types
     * @param atsPerMol array of length molecules, each entry defining how many atoms are in this molecule
     * @param energyparts an array for energyparts of length molecules, may be used to store an energy decomposition in
     * @param iNoOfAtoms the total number of atoms
     * @param faCharges the charges as floats, of length number of atoms
     * @param spins the spins as integers, of length number of atoms
     * @param bonds the bonds in this "atom assembly" (i.e., geometry)
     * @return the energy of this geometry
     * @param hasRigidEnvironment whether the coordinate set contains a rigid environment (per definition the last "molecule" whose internal energy contribution can be ignored (if the method allows for that).
     */
    double energyCalculation(long lID, int iIteration, double[] xyz1D, String[] saAtomTypes,
            short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges,
            short[] spins, final BondInfo bonds, final boolean hasRigidEnvironment);
}
