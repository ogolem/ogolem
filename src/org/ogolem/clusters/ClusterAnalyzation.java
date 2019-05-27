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
package org.ogolem.clusters;

import contrib.jama.*;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.CollisionDetection;
import org.ogolem.core.CollisionInfo;
import org.ogolem.core.DissociationDetection;

import static java.lang.Math.*;
import org.ogolem.core.BondInfo;
import org.ogolem.core.SimpleBondInfo;

/**
 * Provides atomics for the actual analyzation of geometries.
 * @author Johannes Dieterich
 * @version 2015-07-20
 */
public class ClusterAnalyzation {

    public enum Shape {
        SPHERICAL, OBLATE, PROLATE, SCALENE
    }

    /**
     * Calculates the moments of inertia of a provided set of cartesian coordinates, orders them by
     * increasing size and scales them to the biggest one.
     * @param cartes
     * @return moments of inertia
     */
    public static double[] calcMomentsOfInertia(final CartesianCoordinates cartes){

        // the number of atoms
        int iNoOfAtoms = cartes.getNoOfAtoms();
        
        // all atom types
        String[] saAtoms = cartes.getAllAtomTypes();

        // all cartesian coordinates
        cartes.moveCoordsToCOM();
        final double[][] daCoords = cartes.getAllXYZCoordsCopy();

        /*
         * now calculate the actual moments of inertia
         */
        double[] daMasses = new double[iNoOfAtoms];
        double[] daDistances1 = new double[iNoOfAtoms];
        double[] daDistances2 = new double[iNoOfAtoms];

        // get all masses and all distances to the COM
        for(int i = 0; i < iNoOfAtoms; i++){
            daMasses[i] = AtomicProperties.giveWeight(saAtoms[i]);
            daDistances2[i] = Math.sqrt(daCoords[0][i]*daCoords[0][i]
                    + daCoords[1][i]*daCoords[1][i]
                    + daCoords[2][i]*daCoords[2][i]);
        }

        // calculate the moment of inertia tensor for both sets
        final Matrix tensor = new Matrix(3,3);
        final double[][] daTensor = tensor.getArray();

        for(int i = 0; i < 3; i++){
            for(int j = i; j < 3; j++){
                if(i == 0 && j == 0){
                    for(int k = 0; k < iNoOfAtoms; k++){
                        daTensor[i][j] += daMasses[k] * (daCoords[1][k]*daCoords[1][k] + daCoords[2][k]*daCoords[2][k]);
                    }
                } else if(i == 1 && j == 1){
                    for(int k = 0; k < iNoOfAtoms; k++){
                        daTensor[i][j] += daMasses[k] * (daCoords[0][k]*daCoords[0][k] + daCoords[2][k]*daCoords[2][k]);
                    }
                } else if(i == 2 && j == 2){
                    for(int k = 0; k < iNoOfAtoms; k++){
                        daTensor[i][j] += daMasses[k] * (daCoords[0][k]*daCoords[0][k] + daCoords[1][k]*daCoords[1][k]);
                    }
                } else{
                    for(int k = 0; k < iNoOfAtoms; k++){
                        daTensor[i][j] -= daMasses[k] * daCoords[i][k] * daCoords[j][k];
                    }
                    // the tensor is actually symmetric
                    daTensor[j][i] = daTensor[i][j];
                }
            }
        }

        /*for(int i = 0; i < 3; i++){
            System.err.println(daTensor[i][i]/daTensor[2][2]);
        }*/

        // calculate the eigenvectors of the tensors
        final EigenvalueDecomposition eigen = new EigenvalueDecomposition(tensor,true);
        //calculate the eigenvalues
        final double[] daEigenVals = eigen.getRealEigenvalues();

        final double[] daMoments = new double[3];
        daMoments[0] = daEigenVals[0] / daEigenVals[2];
        daMoments[1] = daEigenVals[1] / daEigenVals[2];
        daMoments[2] = daEigenVals[2] / daEigenVals[2];

        //for(double d : daMoments){
        //    System.err.println(d);
        //}

        return daMoments;
    }

    private static double[][] substractValues(double[][] daCoords, final double[] daValues){
        // substract the COM values from all XYZ coordinates
        for(int i = 0; i < daValues.length; i++){
            for(int j = 0; j < daCoords[0].length; j++){
                daCoords[i][j] = daCoords[i][j] - daValues[i];
            }
        }
        return daCoords;
    }

    /**
     * Estimates the shape of the cluster from its moments of inertia.
     * @param daMomentsOfInertia
     * @return the shape
     */
    public static Shape shapeOfCluster(final double[] daMomentsOfInertia, final double dThreshSame,
            final double dThreshDiff){
        // i'm feeling lazy.. ;-)
        double[] da = daMomentsOfInertia;
        double dTS = dThreshSame;
        double dTD = dThreshDiff;

        if((da[0]-da[1]) < dTS && (da[1]-da[2]) < dTS){
            // spherical
            return Shape.SPHERICAL;
        } else if(((da[0]-da[1]) < dTS && (da[2] - da[0]) > dTD)
                || ((da[1]-da[2]) < dTS && (da[0] - da[1]) > dTD)
                || ((da[0]-da[2]) < dTS && (da[1] - da[0]) > dTD)
                ){
            // oblate
            return Shape.OBLATE;
        } else if (((da[0]-da[1]) < dTS && (da[2] - da[0]) < dTD)
                || ((da[1]-da[2]) < dTS && (da[0] - da[1]) < dTD)
                || ((da[0]-da[2]) < dTS && (da[1] - da[0]) < dTD)
                ){
            // prolate
            return Shape.PROLATE;
        } else {
            // non-defined
            return Shape.SCALENE;
        }
    }

    /**
     * Estimates the shape of the cluster from its moments of inertia. Works with
     * a single threshold.
     * @param daMomentsOfInertia
     * @return the shape
     */
    public static Shape shapeOfCluster(final double[] daMomentsOfInertia, final double dThreshSame){
        // i'm feeling lazy.. ;-)
        double[] da = daMomentsOfInertia;
        double ts = dThreshSame;

        if(abs(da[2]-da[1]) < ts && abs(da[2]-da[0]) < ts){
            // spherical
            return Shape.SPHERICAL;
        } else if((abs(da[2]-da[1]) < ts && abs(da[2]-da[0]) > ts)){
            // oblate
            return Shape.OBLATE;
        } else if((abs(da[2]-da[0]) > ts && abs(da[0]-da[1]) < ts)){
            // prolate
            return Shape.PROLATE;
        } else if(abs(da[2]-da[1]) > ts && abs(da[1]-da[0]) > ts){
            // scalene
            return Shape.SCALENE;
        } else{
            System.err.println("Error in figuring ellipsoid out.");
            return Shape.SCALENE;
        }
    }

    /**
     * Detects dissociations in a set of cartesian coordinates. Quite expensive with bad scaling,
     * two times O(N^2).
     * @param cartes
     * @param dBlowDissoc
     * @return true if a dissociation has been detected
     */
    static boolean detectDissociations(CartesianCoordinates cartes, double dBlowDissoc){
        
        // we need the advanced one, so full detection.
        final CollisionDetection collDetect = new CollisionDetection(CollisionDetection.DEFAULTCD);

        // the  bonds is a fake, of course...
        final BondInfo bonds = new SimpleBondInfo(cartes.getNoOfAtoms());

        final CollisionInfo info = collDetect.checkForCollision(cartes, 1.0, bonds);

        // now do the dissociation detection (use DFS)
        final boolean bDissoc = DissociationDetection.checkForDissociation(info.getPairWiseDistances(), cartes.getAllAtomTypes(), cartes.getAllAtomNumbers(), dBlowDissoc, DissociationDetection.DEFAULTDD);

        return bDissoc;
    }

    /**
     * Calculates all COMs of all molecules relative to the overall COM.
     * @param cartes
     * @return an array of the distances of all molecules to the COM
     */
    static double[][] calculateDistanceToCOM(CartesianCoordinates cartes){

        // first, again, the overall COM
        double[] daOverallCOM = cartes.calculateTheCOM();

        double[][] daCoordsCopy = cartes.getAllXYZCoordsCopy();

        daCoordsCopy = substractValues(daCoordsCopy, daOverallCOM);

        CartesianCoordinates cartesCopy = new CartesianCoordinates(cartes);

        // put the new xyz values in
        cartesCopy.setAllXYZ(daCoordsCopy);

        // now create all the relative COMs of all indies
        int iNoOfIndies = cartesCopy.getNoOfMolecules();

        double[][] daAllCOMs = new double[iNoOfIndies][];

        for(int i = 0; i < iNoOfIndies; i++){
            daAllCOMs[i] = cartesCopy.giveMolecularCartes(i, false).calculateTheCOM();
        }

        return daAllCOMs;
    }

    static double ljNearestNeighbors(final CartesianCoordinates cartes){

        // the atomic types
        final String[] saAtomTypes = cartes.getAllAtomTypes();
        // all cartesian coordinates
        final double[][] daXYZ = cartes.getAllXYZCoordsCopy();

        final int iNoOfAtoms = cartes.getNoOfAtoms();

        final double dFactor = Math.pow(2.0, (1.0/6.0));

        double dEnergy = 0.0;

        for(int i = 0; i < iNoOfAtoms - 1; i++){
            // get the first lj parameter
            final double dSigma1 = AtomicProperties.giveLennardJonesSigma(saAtomTypes[i]);
            final double dEpsilon1 = AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[i]);

            for(int j = i + 1; j < iNoOfAtoms; j++){
                // second lj parameter
                final double dSigma2 = AtomicProperties.giveLennardJonesSigma(saAtomTypes[j]);
                final double dEpsilon2 = AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[j]);

                // mix them using Lorentz-Berthelot mixing rules
                final double dEpsilon = Math.sqrt(dEpsilon1 * dEpsilon2);
                final double dSigma = 0.5 * (dSigma1 + dSigma2);

                // calculate the interatomic distance
                final double dDist = Math.sqrt(Math.pow((daXYZ[0][i] - daXYZ[0][j]), 2.0)
                        + Math.pow((daXYZ[1][i] - daXYZ[1][j]), 2.0)
                        + Math.pow((daXYZ[2][i] - daXYZ[2][j]), 2.0));
                
                // check whether it is below the cutoff
                if(dDist < dSigma*dFactor){
                    dEnergy += dEpsilon;
                }
            }
        }

        return dEnergy;
    }

    static double ljStrainEnergy(final CartesianCoordinates cartes){
        
        // the atomic types
        final String[] saAtomTypes = cartes.getAllAtomTypes();
        // all cartesian coordinates
        final double[][] daXYZ = cartes.getAllXYZCoordsCopy();

        final int iNoOfAtoms = cartes.getNoOfAtoms();

        final double dFactor = Math.pow(2.0, (1.0/6.0));

        double dEnergy = 0.0;

        for(int i = 0; i < iNoOfAtoms - 1; i++){
            // get the first set of lj parameters
            final double dEpsilon1 = AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[i]);
            final double dSigma1 = AtomicProperties.giveLennardJonesSigma(saAtomTypes[i]);

            for(int j = i + 1; j < iNoOfAtoms; j++){
                // second set of lj parameters
                final double dEpsilon2 = AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[j]);
                final double dSigma2 = AtomicProperties.giveLennardJonesSigma(saAtomTypes[j]);

                // mix them using Lorentz-Berthelot mixing rules
                final double dEpsilon = Math.sqrt(dEpsilon1 * dEpsilon2);
                final double dSigma = 0.5 * (dSigma1 + dSigma2);

                // calculate the interatomic distance
                final double dDist = Math.sqrt(Math.pow((daXYZ[0][i] - daXYZ[0][j]), 2)
                        + Math.pow((daXYZ[1][i] - daXYZ[1][j]), 2)
                        + Math.pow((daXYZ[2][i] - daXYZ[2][j]), 2));

                // check whether it is below the cutoff
                if(dDist < dSigma*dFactor){
                    // calculate the LJ contribution
                    final double dDistSquared = Math.pow(dDist, 2);
                    final double dInvRPow2 = dSigma * dSigma / dDistSquared;
                    final double dInvRPow6 = dInvRPow2 * dInvRPow2 * dInvRPow2;
                    final double dInvRPow12 = dInvRPow6 * dInvRPow6;
                    final double dTempEnergy = 4.0 * dEpsilon * (dInvRPow12 - dInvRPow6);

                    // calculate the complete thing
                    dEnergy += dTempEnergy - dEpsilon;
                }
            }
        }

        return dEnergy;
    }
    
    static double ljEnergyNonNearest(final CartesianCoordinates cartes){

        // the atomic types
        final String[] saAtomTypes = cartes.getAllAtomTypes();
        // all cartesian coordinates
        final double[][] daXYZ = cartes.getAllXYZCoordsCopy();

        final int iNoOfAtoms = cartes.getNoOfAtoms();

        final double dFactor = Math.pow(2.0, (1.0/6.0));

        double dEnergy = 0.0;

        for(int i = 0; i < iNoOfAtoms - 1; i++){
            // get the first set of lj parameters
            final double dEpsilon1 = AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[i]);
            final double dSigma1 = AtomicProperties.giveLennardJonesSigma(saAtomTypes[i]);

            for(int j = i + 1; j < iNoOfAtoms; j++){
                // second set of lj parameters
                final double dEpsilon2 = AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[j]);
                final double dSigma2 = AtomicProperties.giveLennardJonesSigma(saAtomTypes[j]);

                // mix them using Lorentz-Berthelot mixing rules
                final double dEpsilon = Math.sqrt(dEpsilon1 * dEpsilon2);
                final double dSigma = 0.5 * (dSigma1 + dSigma2);

                // calculate the interatomic distance
                final double dDist = Math.sqrt(Math.pow((daXYZ[0][i] - daXYZ[0][j]), 2)
                        + Math.pow((daXYZ[1][i] - daXYZ[1][j]), 2)
                        + Math.pow((daXYZ[2][i] - daXYZ[2][j]), 2));

                // check whether it is below the cutoff
                if(dDist >= dSigma*dFactor){
                    // calculate the LJ contribution
                    final double dDistSquared = Math.pow(dDist, 2);
                    final double dInvRPow2 = dSigma * dSigma / dDistSquared;
                    final double dInvRPow6 = dInvRPow2 * dInvRPow2 * dInvRPow2;
                    final double dInvRPow12 = dInvRPow6 * dInvRPow6;
                    final double dTempEnergy = 4.0 * dEpsilon * (dInvRPow12 - dInvRPow6);

                    // add it up
                    dEnergy += dTempEnergy;
                }
            }
        }

        return dEnergy;
    }
}
