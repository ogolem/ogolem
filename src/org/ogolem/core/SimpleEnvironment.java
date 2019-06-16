/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Random;
import contrib.jama.*;

/**
 * Holds as the father of all environments for the cluster.
 * @author Johannes Dieterich
 * @version 2015-07-20
 */
public final class SimpleEnvironment implements Environment{

    // serialVersion so that serialization works flawlessly
    private static final long serialVersionUID = (long) 20101107;

    private final boolean[][] envBonds;

    private final boolean flexySurface;

    //TODO howto use this? and what is with a more fine-grained rigid-part?
    private final List<Integer> listSecAtoms;

    private final CollisionDetection.CDTYPE whichCollDetect;

    private final double blowColl;

    private final Atom[] referencePoints;

    private final CartesianCoordinates cartEnv;

    private double[] distanceCOMs;

    private final double[] eulersToCluster;

    private final AllowedSpace space;
    
    /**
     * The constructor for the initial creation.
     * @param environment
     * @param flexy
     * @param bonds
     */
    SimpleEnvironment(final CartesianCoordinates environment, final boolean flexy,
            final boolean[][] bonds, final List<Integer> secondaryAtoms,
            final CollisionDetection.CDTYPE whichColl, final double collBlow, final Atom[] referenceAtoms,
            final AllowedSpace allowedSpace){
        this.cartEnv = environment;
        this.envBonds = bonds;
        this.flexySurface = flexy;
        this.distanceCOMs = new double[3];
        this.eulersToCluster = new double[3];
        this.listSecAtoms = secondaryAtoms;
        this.blowColl = collBlow;
        this.whichCollDetect = whichColl;
        this.referencePoints = referenceAtoms;
        this.space = allowedSpace;
    }

   /**
    * A standard copy constructor.
    * @param original
    */
   private SimpleEnvironment(final SimpleEnvironment original){
        this.cartEnv = new CartesianCoordinates(original.cartEnv);

        final boolean[][] baOrigBonds = original.envBonds;
        this.envBonds = new boolean[baOrigBonds.length][];

        // make a deep copy. yes, this is expensive.
        for(int i = 0; i < baOrigBonds.length; i++){
            envBonds[i] = baOrigBonds[i].clone();
        }

        this.distanceCOMs = original.distanceCOMs;
        this.eulersToCluster = original.eulersToCluster.clone();
        this.flexySurface = original.flexySurface;
        this.whichCollDetect = original.whichCollDetect;
        this.blowColl = original.blowColl;

        // now the List of secondary atoms
        this.listSecAtoms = new ArrayList<>(original.listSecAtoms.size());
        final Iterator<Integer> itSecOrig = original.listSecAtoms.iterator();
        while(itSecOrig.hasNext()){
            int iTemp = itSecOrig.next();
            listSecAtoms.add(iTemp);
        }

        this.referencePoints = original.referencePoints;
        this.space = original.space.clone();
   }

    @Override
    public SimpleEnvironment clone(){
        return new SimpleEnvironment(this);
    }

    @Override
    public boolean doesItFit(final CartesianCoordinates clusterCartes){

        if(!space.isPointInSpace(distanceCOMs)){
            // even if it might potentially fit, we just do not want it!
            return false;
        }

        // marry them
        final CartesianCoordinates married = marryThem(clusterCartes);

        // check for CD (no DD required)
        final CollisionDetection detect = new CollisionDetection(whichCollDetect);

        final BondInfo marriedBonds = this.marryBonds(clusterCartes.getNoOfAtoms(),
                (clusterCartes.getNoOfAtoms()+cartEnv.getNoOfAtoms()));
        final CollisionInfo info = detect.checkForCollision(married, blowColl, marriedBonds);
        /*
         * important side-note: due to us actually *wanting* interactions between
         * the environment and the cluster, the blow factor above needs to be very
         * small, below 1.0, to actually just alert *real* clashes.
         */


        return info.hasCollision();
    }

    @Override
    public CartesianCoordinates marryThem(final CartesianCoordinates clusterCartes){

        if(clusterCartes.containsEnvironment()){
            System.err.println("ERROR: Cluster cartesian already contains an environment. This is a bug. Contact author(s).");
            return clusterCartes;
        }

        /*
         * we initially use the saved eulers and distance of COMs
         */
        clusterCartes.moveCoordsToCOM();
        final double[][] xyzCluster = clusterCartes.getAllXYZCoord();

        // rotation
        CoordTranslation.rotateXYZ(xyzCluster, eulersToCluster);

        // translation
        final double[] trans = new double[3];
        for(int i = 0; i < 3; i++){
            trans[i] = distanceCOMs[i] + referencePoints[0].getPosition()[i];
        }
        clusterCartes.moveCoordsToPoint(trans);

        // put the surface and the cluster together
        final int iNoOfAtoms = clusterCartes.getNoOfAtoms() + cartEnv.getNoOfAtoms();
        final int iNoOfMols = clusterCartes.getNoOfMolecules() + 1;
        final int[] iaAtsPerMol = new int[iNoOfMols];
        final int[] iaAtsPerMolCluster = clusterCartes.getAllAtomsPerMol();
        System.arraycopy(iaAtsPerMolCluster, 0, iaAtsPerMol, 0, iaAtsPerMolCluster.length);
        iaAtsPerMol[iaAtsPerMol.length-1] = cartEnv.getNoOfAtoms();

        CartesianCoordinates completeCartes = new CartesianCoordinates(iNoOfAtoms, iNoOfMols, iaAtsPerMol);

        final double[][] xyzEnv = cartEnv.getAllXYZCoord();
        final double[][] xyzComplete = new double[3][iNoOfAtoms];

        // copy around a little
        System.arraycopy(xyzCluster[0], 0, xyzComplete[0], 0, xyzCluster[0].length);
        System.arraycopy(xyzCluster[1], 0, xyzComplete[1], 0, xyzCluster[0].length);
        System.arraycopy(xyzCluster[2], 0, xyzComplete[2], 0, xyzCluster[0].length);
        System.arraycopy(xyzEnv[0], 0, xyzComplete[0], xyzCluster[0].length, xyzEnv[0].length);
        System.arraycopy(xyzEnv[1], 0, xyzComplete[1], xyzCluster[0].length, xyzEnv[0].length);
        System.arraycopy(xyzEnv[2], 0, xyzComplete[2], xyzCluster[0].length, xyzEnv[0].length);
        
        completeCartes.setAllXYZ(xyzComplete);

        // the atom types
        final String[] saAtomsCompl = new String[iNoOfAtoms];
        final String[] saAtomsCluster = clusterCartes.getAllAtomTypes();
        final String[] saAtomsEnv = cartEnv.getAllAtomTypes();

        System.arraycopy(saAtomsCluster,0,saAtomsCompl,0,saAtomsCluster.length);
        System.arraycopy(saAtomsEnv,0,saAtomsCompl,saAtomsCluster.length,saAtomsEnv.length);

        completeCartes.setAllAtomTypes(saAtomsCompl);
        
        // spins and charges
        final short[] iaSpinsCompl = new short[completeCartes.getNoOfAtoms()];
        final float[] faChargesCompl = new float[completeCartes.getNoOfAtoms()];

        final short[] iaSpinsCluster = clusterCartes.getAllSpins();
        final float[] faChargesCluster = clusterCartes.getAllCharges();

        System.arraycopy(iaSpinsCluster, 0, iaSpinsCompl, 0, iaSpinsCluster.length);
        System.arraycopy(faChargesCluster, 0, faChargesCompl, 0, faChargesCluster.length);

        completeCartes.setAllSpins(iaSpinsCompl);
        completeCartes.setAllCharges(faChargesCompl);

        // zmats
        final ZMatrix[] allZmatsCluster = clusterCartes.getZMatrices();
        final ZMatrix[] allZmatsCompl = new ZMatrix[allZmatsCluster.length + 1];

        for(int i = 0; i < allZmatsCluster.length; i++){
            final ZMatrix currZMat = allZmatsCluster[i];
            if(currZMat == null) allZmatsCompl[i] = null;
            else allZmatsCompl[i] = new ZMatrix(currZMat);
        }
        
        // add one for this "molecule"
        allZmatsCompl[allZmatsCompl.length-1] = null;

        completeCartes.setZMatrices(allZmatsCompl);

        // and of course we should also add a clone of this
        completeCartes.setRefEnvironment(clone());

        
        return completeCartes;
    }

    @Override
    public CartesianCoordinates divorceThem(final CartesianCoordinates completeCartes){

        if(!completeCartes.containsEnvironment()){
            System.err.println("ERROR: This is fatal and a bug. Contact the author(s).");
            return null;
        }

        final int[] iaComplAtsPerMol = completeCartes.getAllAtomsPerMol();
        final int[] iaClusterAtsPerMol = new int[iaComplAtsPerMol.length-1];
        int iNoOfAtomsCluster = 0;
        for(int i = 0; i < iaClusterAtsPerMol.length; i++){
            iaClusterAtsPerMol[i] = iaComplAtsPerMol[i];
            iNoOfAtomsCluster += iaComplAtsPerMol[i];
        }

        /*
         * who get's the car and the house?
         * affine transformation
         */

        // first find out where our references are now
        double[] daRef1New = completeCartes.getXYZCoordinatesOfAtom(referencePoints[0].getID()+iNoOfAtomsCluster);
        double[] daRef2New = completeCartes.getXYZCoordinatesOfAtom(referencePoints[1].getID()+iNoOfAtomsCluster);
        double[] daRef3New = completeCartes.getXYZCoordinatesOfAtom(referencePoints[2].getID()+iNoOfAtomsCluster);

        double[] daRef1Old = referencePoints[0].getPosition();
        double[] daRef2Old = referencePoints[1].getPosition();
        double[] daRef3Old = referencePoints[2].getPosition();

        // now build up two matrices
        double[][] daMatT1 = new double[3][3];
        daMatT1[0][0] = daRef1Old[0];
        daMatT1[1][0] = daRef1Old[1];
        daMatT1[2][0] = daRef1Old[2];
        daMatT1[0][1] = daRef2Old[0];
        daMatT1[1][1] = daRef2Old[1];
        daMatT1[2][1] = daRef2Old[2];
        daMatT1[0][2] = daRef3Old[0];
        daMatT1[1][2] = daRef3Old[1];
        daMatT1[2][2] = daRef3Old[2];


        double[][] daMatT2 = new double[3][3];
        daMatT2[0][0] = daRef1New[0];
        daMatT2[1][0] = daRef1New[1];
        daMatT2[2][0] = daRef1New[2];
        daMatT2[0][1] = daRef2New[0];
        daMatT2[1][1] = daRef2New[1];
        daMatT2[2][1] = daRef2New[2];
        daMatT2[0][2] = daRef3New[0];
        daMatT2[1][2] = daRef3New[1];
        daMatT2[2][2] = daRef3New[2];
        
        Matrix matT1 = new Matrix(daMatT1);
        Matrix matT2 = new Matrix(daMatT2);

        // calculate the transformation matrix
        Matrix matTrans = matT2.times(matT1.inverse());

        // transform the cartesians
        Matrix matCartes = new Matrix(completeCartes.getAllXYZCoord());
        Matrix matNewCartes = matTrans.times(matCartes);

        completeCartes.setAllXYZ(matNewCartes.getArrayCopy());


        /*
         * now we can divorce
         */
        final CartesianCoordinates clusterCartes = new CartesianCoordinates(iNoOfAtomsCluster,
                iaClusterAtsPerMol.length, iaClusterAtsPerMol);

        // the spins and charges
        final short[] iaSpinsCluster = new short[iNoOfAtomsCluster];
        final float[] faChargesCluster = new float[iNoOfAtomsCluster];
        final short[] iaSpinsCompl = completeCartes.getAllSpins();
        final float[] faChargesCompl = completeCartes.getAllCharges();

        System.arraycopy(iaSpinsCompl,0,iaSpinsCluster,0,iNoOfAtomsCluster);
        System.arraycopy(faChargesCompl,0,faChargesCluster,0,iNoOfAtomsCluster);

        clusterCartes.setAllSpins(iaSpinsCompl);
        clusterCartes.setAllCharges(faChargesCompl);

        // the zmats
        final ZMatrix[] allComplZmats = completeCartes.getZMatrices();
        final ZMatrix[] allClusterZmats = new ZMatrix[allComplZmats.length-1];
        for(int i = 0; i < allComplZmats.length-1; i++){
            final ZMatrix zmatTemp = allComplZmats[i];
            if(zmatTemp == null) allClusterZmats[i] = null;
            else allClusterZmats[i] = new ZMatrix(zmatTemp);
        }
        clusterCartes.setZMatrices(allClusterZmats);

        // coordinates and atom types
        final String[] saAtomsCluster = new String[iNoOfAtomsCluster];
        final String[] saAtomsComplete = completeCartes.getAllAtomTypes();

        System.arraycopy(saAtomsComplete,0,saAtomsCluster,0,iNoOfAtomsCluster);
        clusterCartes.setAllAtomTypes(saAtomsCluster);

        final double[][] daClusterCoords = new double[3][iNoOfAtomsCluster];
        final double[][] daComplCoords = completeCartes.getAllXYZCoord();

        System.arraycopy(daComplCoords[0], 0, daClusterCoords[0], 0, iNoOfAtomsCluster);
        System.arraycopy(daComplCoords[1], 0, daClusterCoords[1], 0, iNoOfAtomsCluster);
        System.arraycopy(daComplCoords[2], 0, daClusterCoords[2], 0, iNoOfAtomsCluster);
        clusterCartes.setAllXYZ(daClusterCoords);

        // the euler angles can stay zero, if (and ONLY if) nothing gets rotated afterwards
        eulersToCluster[0] = 0.0;
        eulersToCluster[1] = 0.0;
        eulersToCluster[2] = 0.0;

        // measure the translation from reference point [0] (daRef1Old) to the COM of the cluster
        // therefore, we first need to find out where the COM now is
        double[] daNewCOM = clusterCartes.calculateTheCOM();
        for(int i = 0; i < 3; i++){
            distanceCOMs[i] = daNewCOM[i] - daRef1Old[i];
        }

        if(flexySurface){
            // we need to update the cartesian in here, since something might have changed
            int iAtSurf = cartEnv.getNoOfAtoms();
            double[][] daSurfXYZ = cartEnv.getAllXYZCoord();
            System.arraycopy(daComplCoords[0], iNoOfAtomsCluster, daSurfXYZ[0], 0, iAtSurf);
            System.arraycopy(daComplCoords[1], iNoOfAtomsCluster, daSurfXYZ[1], 0, iAtSurf);
            System.arraycopy(daComplCoords[2], iNoOfAtomsCluster, daSurfXYZ[2], 0, iAtSurf);
        }

        // some last minor modifications
        clusterCartes.setEnergy(completeCartes.getEnergy());

        return clusterCartes;
    }

    @Override
    public void initializeConnections(final CartesianCoordinates clusterCartes){

        final Random random = new Random();

        final CollisionDetection collDetect = new CollisionDetection(whichCollDetect);

        boolean hasColl;
        do {
            // first the eulers
            double dr;
            boolean br = random.nextBoolean();

            // first acting on phi
            while (true) {
                dr = random.nextDouble();
                if (dr == 1.0 && br == true) {
                    // need a new one
                } else {
                    break;
                }
            }
            if (br) {
                dr *= -Math.PI;
            } else {
                dr *= Math.PI;
            }
            eulersToCluster[0] = dr;

            // then omega
            br = random.nextBoolean();
            while (true) {
                dr = random.nextDouble();
                if (dr == 1.0 && br == true) {
                    // need a new one
                } else {
                    break;
                }
            }
            if (br) {
                dr *= -Math.PI / 2;
            } else {
                dr *= Math.PI / 2;
            }
            eulersToCluster[1] = dr;

            // then psi
            br = random.nextBoolean();
            while (true) {
                dr = random.nextDouble();
                if (dr == 1.0 && br == true) {
                    // need a new one
                } else {
                    break;
                }
            }
            if (br) {
                dr *= -Math.PI;
            } else {
                dr *= Math.PI;
            }
            eulersToCluster[2] = dr;

            // now the distances
            distanceCOMs = space.getPointInSpace();

            // check for collisions
            hasColl = collDetect.checkForCollision(marryThem(clusterCartes), blowColl,
                    marryBonds(clusterCartes.getNoOfAtoms(), (clusterCartes.getNoOfAtoms()+cartEnv.getNoOfAtoms()))).hasCollision();
        } while (hasColl);
    }

    @Override
    public void mutateConnections(final CartesianCoordinates clusterCartes){
        
        // first we want to figure out where to mutate
        final Random random = new Random();
        final int mutPos = random.nextInt(6);

        double dr;
        boolean br = random.nextBoolean();
        while (true) {
            dr = random.nextDouble();
            if (dr == 1.0 && br == true) {
                // we need to take a new double
            } else {
                break;
            }
        }

        switch(mutPos){
            case 0:
                // COM-COM first
                distanceCOMs[0] *= 1.5*dr;
                break;
            case 1:
                // COM-COM second
                distanceCOMs[1] *= 1.5*dr;
                break;
            case 2:
                // COM-COM third
                distanceCOMs[2] *= 1.5*dr;
                break;
            case 3:
                // phi
                if (br) {
                    dr = -dr;
                    dr = dr * Math.PI;
                } else {
                    dr = dr * Math.PI;
                }
                eulersToCluster[0] = dr;
                break;
            case 4:
                // omega
                if (br) {
                    dr = -dr;
                    dr = dr * Math.PI / 2;
                } else {
                    dr = dr * Math.PI / 2;
                }
                eulersToCluster[1] = dr;
                break;
            case 5:
                // psi
                if (br) {
                    dr = -dr;
                    dr = dr * Math.PI;
                } else {
                    dr = dr * Math.PI;
                }
                eulersToCluster[2] = dr;
                break;
            default:
                System.err.println("ERROR: You should NEVER end up here in SimpleEnvironment. Please contact the author(s).");
                break;
        }
    }

    @Override
    public List<Environment> createOffspring(final Environment father){

        final ArrayList<Environment> alChildren = new ArrayList<>(2);

        final Environment surfChild1 = father.clone();
        final Environment surfChild2 = clone();
        
        final double[] daFatherGenom = father.returnMyGenom();
        final double[] daMotherGenom = returnMyGenom();

        final Random random = new Random();
        // we make a linearily distributed cut
        final int iCut = random.nextInt(6);

        // copy over
        final double[] daChild1 = new double[6];
        final double[] daChild2 = new double[6];

        for(int i = 0; i < iCut; i++){
            daChild1[i] = daFatherGenom[i];
            daChild2[i] = daMotherGenom[i];
        }

        for(int i = iCut; i < daChild2.length; i++){
            daChild1[i] = daMotherGenom[i];
            daChild2[i] = daFatherGenom[i];
        }

        // put the new genes into the children
        surfChild1.putGenomIn(daChild1);
        surfChild2.putGenomIn(daChild2);

        // return the children
        alChildren.add(surfChild1);
        alChildren.add(surfChild2);

        return alChildren;
    }

    @Override
    public List<Integer> whichSecondaryAtoms(){
        
        final List<Integer> secAtoms = new ArrayList<>();

        final Iterator<Integer> itSecs = listSecAtoms.iterator();
        while(itSecs.hasNext()){
            int iSec = itSecs.next();
            secAtoms.add(iSec);
        }

        return secAtoms;
    }

    @Override
    public boolean isEnvironmentRigid(){
        return !flexySurface;
    }

    @Override
    public double[] returnMyGenom(){
        
        final double[] genom = new double[6];
        genom[0] = distanceCOMs[0];
        genom[1] = distanceCOMs[1];
        genom[2] = distanceCOMs[2];
        genom[3] = eulersToCluster[0];
        genom[4] = eulersToCluster[1];
        genom[5] = eulersToCluster[2];

        return genom;
    }

    @Override
    public void putGenomIn(final double[] genom){
        
        distanceCOMs[0] = genom[0];
        distanceCOMs[1] = genom[1];
        distanceCOMs[2] = genom[2];
        eulersToCluster[0] = genom[3];
        eulersToCluster[1] = genom[4];
        eulersToCluster[2] = genom[5];
    }

    @Override
    public int atomsInEnv(){
        return cartEnv.getNoOfAtoms();
    }

    private BondInfo marryBonds(final int noClusterAts, final int noTotalAts){

        final BondInfo bondsCompl = new SimpleBondInfo(noTotalAts);

        for(int i = 0; i < noClusterAts; i++){
            for(int j = 0; j < noClusterAts; j++){
                // it doesn't matter, what we fill in here since it's not about these distances anyways (but there must be bonds)
                bondsCompl.setBond(i, j, BondInfo.UNCERTAIN);
            }
        }

        for(int i = noClusterAts; i < noTotalAts; i++){
            for(int j = noClusterAts; j < noTotalAts; j++){
                if(envBonds[i-noClusterAts][j-noClusterAts]){
                    bondsCompl.setBond(i, j, BondInfo.UNCERTAIN);
                } // default is no bond, so this is fine
            }
        }

        return bondsCompl;
    }
}
