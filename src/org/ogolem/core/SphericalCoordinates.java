/**
Copyright (c) 2009, J. M. Dieterich and B. Hartke
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
 * A value object storing spherical coordinates.
 * @author Johannes Dieterich
 * @version 2009-06-08
 */
class SphericalCoordinates implements Serializable{

    // serialVersion so that serialization works flawlessly
	private static final long serialVersionUID = (long) 20090608;

    private final int iNoOfMolecules;

	private final int iNoOfAtoms;

	private final int[] iaAtomsPerMolecule;

	private final double[][] daCoordinates;

	private String [] saAtoms;

    /**
     * A specific constructor to initialize all sizes correct.
     * @param iNoOfAtoms The total number of atoms in the geometry.
     * @param iNoOfMolecules The number of molecules in the geometry.
     * @param iaAtPerMol The number of atoms per molecule.
     */
    SphericalCoordinates(final int iNoOfAtoms, final int iNoOfMolecules,final int[] iaAtPerMol){
        this.iaAtomsPerMolecule = iaAtPerMol;
		this.iNoOfAtoms = iNoOfAtoms;
		this.iNoOfMolecules = iNoOfMolecules;
		this.daCoordinates = new double[3][iNoOfAtoms];
		this.saAtoms = new String[iNoOfAtoms];
    }

    /**
     * Constructing the spherical coordinates by using cartesian ones.
     * @param cartes The cartesian coordinates object.
     */
    SphericalCoordinates(final CartesianCoordinates cartes){
        this.iaAtomsPerMolecule = cartes.getAllAtomsPerMol();
        this.iNoOfAtoms = cartes.getNoOfAtoms();
        this.iNoOfMolecules = cartes.getNoOfMolecules();
        this.daCoordinates = CoordTranslation.cartesianToSphericalCoord(cartes.getAllXYZCoord());
    }

    /**
     * Translated the sperical coordinates into cartesian ones and returns them as an object. It does NOT
     * though set ANY reference z-matrices, energies, methods etc. . These things need to be set afterwards
     * if needed.
     * @return CartesianCoordinates object.
     */
    CartesianCoordinates translateToCartes(){
        CartesianCoordinates cartes = new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaAtomsPerMolecule);
        cartes.setAllAtomTypes(saAtoms);
        double[][] daXYZ = CoordTranslation.sphericalToCartesianCoord(daCoordinates);
        cartes.setAllXYZ(daXYZ);
        return cartes;
    }
}
