/**
Copyright (c) 2014, J. M. Dieterich
              2016-2020, J. M. Dieterich and B. Hartke
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

import static org.ogolem.core.Constants.ANGTOBOHR;
import org.ogolem.generic.IndividualReader;
import org.ogolem.io.InputPrimitives;

/**
 * Individual reader implementation for geometries.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class GeometryReader implements IndividualReader<Geometry>{

    private static final long serialVersionUID = (long) 20200429;
    private static final boolean DEBUG = false;
    
    @Override
    public GeometryReader copy() {
        return new GeometryReader();
    }

    @Override
    public void populateIndividualFromFile(final Geometry individual, final String file) throws Exception {
        
        final String[] data = InputPrimitives.readFileIn(file);
        final CartesianCoordinates cartes = individual.getCartesiansWithEnvironment();
        final String[] atoms = cartes.getAllAtomTypes();
        final double[][] xyz = cartes.getAllXYZCoord();
        final int noAtoms = cartes.getNoOfAtoms();
        
        if(noAtoms > data.length-2){throw new RuntimeException("Not enough data in xyz file " + file);}
        
        for(int i = 2; i < noAtoms+2; i++){
            final String[] sa = data[i].trim().split("\\s+");
            if(!sa[0].equalsIgnoreCase(atoms[i-2])){throw new RuntimeException("Wrong atoms. Should be " + atoms[i-2] + " is " + sa[0]);}
            
            xyz[0][i-2] = Double.parseDouble(sa[1])*ANGTOBOHR;
            xyz[1][i-2] = Double.parseDouble(sa[2])*ANGTOBOHR;
            xyz[2][i-2] = Double.parseDouble(sa[3])*ANGTOBOHR;
            
            // just to make sure...
            assert(!Double.isNaN(xyz[0][i-2]) && !Double.isInfinite(xyz[0][i-2]));
            assert(!Double.isNaN(xyz[1][i-2]) && !Double.isInfinite(xyz[1][i-2]));
            assert(!Double.isNaN(xyz[2][i-2]) && !Double.isInfinite(xyz[2][i-2]));
        }
        
        // translate back...
        final GeometryConfig geoConf = individual.returnMyConfig();
        final Geometry geom = CoordTranslation.cartesianToGeometry(cartes, cartes.getNoOfMolecules(),
                cartes.getAllAtomsPerMol(), individual.getAllFlexies(), individual.getExplicitDoFs(),
                individual.getAllConstraints(true), individual.getAllConstraintsXYZ(true),
                individual.getSIDs(), geoConf.bonds.clone());
        
        individual.setAllCOMs(geom.getAllCOMs());
        individual.setAllExtCoords(geom.getAllExtCoords());
        for(int i = 0; i < cartes.getNoOfMolecules(); i++){
            final Molecule mol = geom.getMoleculeAtPosition(i);
            individual.setMoleculeAtPosition(i, new Molecule(mol));
        }
        
        if(individual.containsEnvironment()){
            individual.setEnvironment(geom.getEnvironment().clone());
        }
        
        if(DEBUG){
            System.out.println("DEBUG: Individual from file " + file);
            for(final String s : individual.makePrintableAbsoluteCoord(true)){
                System.out.println(s);
            }
        }
    }
}
