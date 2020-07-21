/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

import java.util.ArrayList;
import java.util.List;
/**
 * Use seeding geometries to create an initial geometry.
 * @author Johannes Dieterich
 * @version 2016-09-03
 */
public final class SimpleSeedingInit {

    private static final boolean DEBUG = true;
    
    public static Geometry readSeedIn(final String sSeedFile, final GeometryConfig refConf,
            final long lID) throws InitIOException, CastException, Exception {

        // set the needed things up
        final boolean hasEnv = (refConf.env != null);        
        
        final int iNoOfMolecules = (hasEnv) ? refConf.noOfParticles + 1 : refConf.noOfParticles;
        final int[] iaAtsPerMol = new int[iNoOfMolecules];
        final boolean[] baFlexies = new boolean[iNoOfMolecules];
        final boolean[] baConstraints = new boolean[iNoOfMolecules];
        final String[] sids = new String[refConf.noOfParticles];

        int iTotNoOfAts = 0;

        boolean anyFlexy = false;
        for(int i = 0; i < refConf.noOfParticles; i++){
            iaAtsPerMol[i] = refConf.geomMCs.get(i).noOfAtoms;
            baFlexies[i] = refConf.geomMCs.get(i).flexy;
            baConstraints[i] = refConf.geomMCs.get(i).constricted;
            sids[i] = refConf.geomMCs.get(i).sID;
            if(refConf.geomMCs.get(i).flexy) anyFlexy = true;
            iTotNoOfAts += iaAtsPerMol[i];
        }
        
        if(hasEnv){
            final int noEnvAtoms = refConf.env.atomsInEnv();
            iTotNoOfAts += noEnvAtoms;
            iaAtsPerMol[refConf.noOfParticles] = noEnvAtoms;
        }

        final short[] iaSpins = new short[iTotNoOfAts];
        final float[] faCharges = new float[iTotNoOfAts];
        final boolean[][] baConstraintsXYZ = new boolean[3][iTotNoOfAts];
        final List<boolean[][]> allFlexies = new ArrayList<>();
        
        int iCounter = 0;
        for(int i = 0; i < refConf.noOfParticles; i++){
            // first get the relevant arrays
            final float[] faCurrCharges = refConf.geomMCs.get(i).charges;
            final short[] iaCurrSpins =  refConf.geomMCs.get(i).spins;
            final boolean[][] baCurrConstr = refConf.geomMCs.get(i).constraints;
            final boolean[][] myFlexies = refConf.geomMCs.get(i).degreesOfFreedom;
            
            if(refConf.geomMCs.get(i).flexy){
                
                final boolean[][] flexies = new boolean[myFlexies.length][];
                for(int x = 0; x < myFlexies.length; x++){
                    flexies[i] = myFlexies[i].clone();
                }
                
                allFlexies.add(flexies);
            } else {
                allFlexies.add(null);
            }

            // copy the stuff
            for(int j = 0; j < faCurrCharges.length; j++){
                faCharges[iCounter] = faCurrCharges[j];
                iaSpins[iCounter] = iaCurrSpins[j];
                if(baCurrConstr != null){
                    baConstraintsXYZ[0][iCounter] = baCurrConstr[0][j];
                    baConstraintsXYZ[1][iCounter] = baCurrConstr[1][j];
                    baConstraintsXYZ[2][iCounter] = baCurrConstr[2][j];
                } else{
                    // which is the case if there is no constraint
                    baConstraintsXYZ[0][iCounter] = false;
                    baConstraintsXYZ[1][iCounter] = false;
                    baConstraintsXYZ[2][iCounter] = false;
                }
                iCounter++;
            }
        }

        // read it in, this might throw exceptions
        CartesianCoordinates cartes;
        boolean wasXYZ = false;
        if(sSeedFile.endsWith(".xyz")){
            cartes = Input.readCartesFromFile(sSeedFile,
                iNoOfMolecules, iaAtsPerMol, iaSpins, faCharges);
            wasXYZ = true;
        } else if(sSeedFile.endsWith(".zmat")){
            final ZMatrix zmat = Input.readZmatFromFile(sSeedFile);
            cartes = CoordTranslation.zMatToCartesians(zmat);
            cartes.setAllCharges(faCharges);
            cartes.setAllSpins(iaSpins);
        } else{
            System.err.println("WARNING: Seeding file " + sSeedFile + " has an unknown extension. Trying as xyz...");
            cartes = Input.readCartesFromFile(sSeedFile,
                iNoOfMolecules, iaAtsPerMol, iaSpins, faCharges);
            wasXYZ = true;
        }
        
        // if it is flexible and xyz, we MUST adapt the cartesian object
        if(wasXYZ && anyFlexy){
            final ZMatrix[] allZMat = new ZMatrix[iNoOfMolecules];
            int count = -1;
            for(final MoleculeConfig mc : refConf.geomMCs){
                count++;
                if(!mc.flexy){
                    allZMat[count] = null;
                    continue;
                }
                
                final ZMatrix zmat = CoordTranslation.cartesToZMat(cartes.giveMolecularCartes(count, true), mc.zmat);
                allZMat[count] = zmat;
            }
            
            cartes.setZMatrices(allZMat);
        }
        
        // consider the environment
        if(refConf.env != null){
            Environment newEnv = refConf.env.clone();
            cartes.setRefEnvironment(newEnv);
        }
        
        int[] myAtsPerMol = iaAtsPerMol;
        if(hasEnv){
            myAtsPerMol = new int[refConf.noOfParticles];
            System.arraycopy(iaAtsPerMol, 0, myAtsPerMol, 0, refConf.noOfParticles);
        }

        // translate it, this might also fail with an exception
        final Geometry geom = CoordTranslation.cartesianToGeometry(cartes, refConf.noOfParticles,
                myAtsPerMol, baFlexies, allFlexies, baConstraints, baConstraintsXYZ, sids, refConf.bonds.clone());

        // finally set the ID in
        geom.setID(lID);
        
        if(DEBUG) System.out.println("DEBUG: Successfully read in " + sSeedFile + " as geometry " + lID);

        return geom;
    }
}
