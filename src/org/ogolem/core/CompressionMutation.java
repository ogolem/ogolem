/**
Copyright (c) 2012, J. M. Dieterich
              2016, J. M. Dieterich and B. Hartke
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
import java.util.Random;

/**
 * Compresses the cluster in order to mutate it a bit. Cannot be used with constraints.
 * @author Johannes Dieterich
 * @version 2016-09-03
 */
public class CompressionMutation implements Serializable {
    
    private static final long serialVersionUID = (long) 20120612;
    private static final Random random = new Random();
    
    public static Molecule mutate(final Molecule mol, final double maxCompression){
        
        if(mol.isConstricted()){
            System.err.println("WARNING: Compression mutation not suitable for constraint molecules. Returning.");
            return mol;
        }
        
        CartesianCoordinates cartes = mol.getCartesians();
        final ZMatrix[] zmatVec = cartes.getZMatrices();

        cartes = cartesToCartes(cartes, maxCompression);

        // we anyway set the old zmatrix first in, so that we have something in there
        cartes.setZMatrices(zmatVec);

       
        if (zmatVec[0] != null) {
            final CartesianCoordinates cartesTemp = cartes.giveMolecularCartes(0, true);
            zmatVec[0] = cartesTemp.calculateZMatrix();
        } // else: nothing needs to happen, unflexible molecule

        cartes.setZMatrices(zmatVec);
        
        final Molecule mEnd = new Molecule (cartes, mol.getMolPosition(), mol.getSID(),
                mol.getFlexy(), mol.getDegreesOfFreedom(), mol.isConstricted(), mol.getConstraints());
        mEnd.setID(mol.getMolPosition());
        mEnd.setSID(mol.getSID());
        
        return mEnd;
    }
    
    public static Geometry mutate(final Geometry g, final double maxCompression){
        
        if(g.isThereAConstraint()){
            System.err.println("WARNING: Compression mutation not suitable for constraint geoemtries. Returning.");
            return g;
        }
        
        CartesianCoordinates cartes = (g.containsEnvironment()) ? g.getCartesiansWithEnvironment(): g.getCartesians();

        final ZMatrix[] zmatVec = cartes.getZMatrices();
        final Environment refEnv = cartes.getReferenceEnvironmentCopy();

        cartes = cartesToCartes(cartes, maxCompression);

        // we anyway set the old zmatrices first in, so that we have something in there
        cartes.setZMatrices(zmatVec);

        for(int i = 0; i < zmatVec.length; i++){
            if(zmatVec[i] != null){
                final CartesianCoordinates cartesTemp = cartes.giveMolecularCartes(i, true);
                zmatVec[i] = cartesTemp.calculateZMatrix();
            } // else: nothing needs to happen, unflexible molecule
        }

        cartes.setZMatrices(zmatVec);
        cartes.setRefEnvironment(refEnv);
        cartes.recalcAtomNumbers();

        final Geometry gEnd = new Geometry(cartes, g.getID(), g.getNumberOfIndieParticles(),
                cartes.getAllAtomsPerMol(), g.getAllFlexies(), g.getExplicitDoFs(), g.getAllConstraints(false),
                g.getAllConstraintsXYZ(false), g.getSIDs(),g.getBondInfo().clone());
        
        gEnd.setFitness(cartes.getEnergy());
        gEnd.setFather(g.getFatherID());
        gEnd.setMother(g.getMotherID());
        gEnd.setLocalOptimized(true);
        
        return gEnd;
    }
    
    private static CartesianCoordinates cartesToCartes(final CartesianCoordinates start, final double maxCompression){
        
        final CartesianCoordinates work = new CartesianCoordinates(start);
        final double[][] xyz = work.getAllXYZCoord();
        
        final double compr = 1-(maxCompression*random.nextDouble());
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < xyz[0].length; j++) xyz[i][j] *= compr;
        }
        
        return work;
    }
    // TODO hook up and test
    // TODO add partial compression (as in: part of a sphere), requires transformation of the coordinates
}
