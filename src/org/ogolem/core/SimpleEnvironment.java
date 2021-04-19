/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

import contrib.jama.*;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.CollisionInfo.Collision;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * Holds as the father of all environments for the cluster.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class SimpleEnvironment implements Environment {

  private static final long serialVersionUID = (long) 20200719;
  protected static final boolean DEBUG = false;

  public static enum FITMODE {
    FULLEXTERNAL,
    LAYERONLY
  };

  public static final double THRESHENPOSSAME = 1e-6;

  protected final Lottery random;

  protected final BondInfo envBonds;

  protected final boolean flexySurface;

  // TODO howto use this? and what is with a more fine-grained rigid-part?
  protected final List<Integer> listSecAtoms;

  protected final CollisionDetection.CDTYPE whichCollDetect;

  protected final double blowColl;

  protected final Atom[] referencePoints;

  protected final CartesianCoordinates cartEnv;

  protected double[] distanceCOMs;

  protected final double[] eulersToCluster;

  protected final FITMODE mode;

  protected final ENVIRONMENTTYPE type;

  protected final AllowedSpace space;

  /**
   * The constructor for the initial creation.
   *
   * @param environment
   * @param flexy
   * @param bonds
   */
  SimpleEnvironment(
      final CartesianCoordinates environment,
      final boolean flexy,
      final BondInfo bonds,
      final List<Integer> secondaryAtoms,
      final CollisionDetection.CDTYPE whichColl,
      final double collBlow,
      final Atom[] referenceAtoms,
      final AllowedSpace allowedSpace,
      final FITMODE fit,
      final ENVIRONMENTTYPE type) {
    this.cartEnv = environment;
    this.envBonds = bonds;
    this.flexySurface = flexy;
    this.distanceCOMs = new double[3];
    this.eulersToCluster = (fit == FITMODE.LAYERONLY) ? null : new double[3];
    this.listSecAtoms = secondaryAtoms;
    this.blowColl = collBlow;
    this.whichCollDetect = whichColl;
    this.referencePoints = referenceAtoms;
    this.space = allowedSpace;
    this.mode = fit;
    this.random = Lottery.getInstance();
    this.type = type;
  }

  /**
   * A standard copy constructor.
   *
   * @param original
   */
  protected SimpleEnvironment(final SimpleEnvironment original) {
    this.cartEnv =
        (original.flexySurface)
            ? new CartesianCoordinates(original.cartEnv)
            : original
                .cartEnv; // only need to copy for flexible surfaces, otherwise the Cartesian is
    // frozen
    this.envBonds =
        original.envBonds
            .shallowCopy(); // the bonds are effectively immutable, so this is allowed to safe on GC
    this.distanceCOMs = original.distanceCOMs.clone();
    this.eulersToCluster =
        (original.mode == FITMODE.LAYERONLY) ? null : original.eulersToCluster.clone();
    this.flexySurface = original.flexySurface;
    this.whichCollDetect = original.whichCollDetect;
    this.blowColl = original.blowColl;

    // now the List of secondary atoms
    this.listSecAtoms = new ArrayList<>(original.listSecAtoms.size());
    for (final int i : listSecAtoms) {
      listSecAtoms.add(i);
    }

    this.referencePoints = new Atom[3];
    this.referencePoints[0] = original.referencePoints[0].copy();
    this.referencePoints[1] = original.referencePoints[1].copy();
    this.referencePoints[2] = original.referencePoints[2].copy();

    this.space = original.space.copy();
    this.mode = original.mode;
    this.random = Lottery.getInstance();
    this.type = original.type;
  }

  @Override
  public SimpleEnvironment copy() {
    return new SimpleEnvironment(this);
  }

  @Override
  public boolean doesItFit(final CartesianCoordinates clusterCartes, final BondInfo clusterBonds) {

    if (!space.isPointInSpace(distanceCOMs)) {
      // even if it might potentially fit, we just do not want it!
      return false;
    }

    // marry them
    final CartesianCoordinates married = marryThem(clusterCartes);
    if (DEBUG) {
      System.out.println("DEBUG: Married trial cartesian:");
      for (final String s : married.createPrintableCartesians()) {
        System.out.println(s);
      }
    }

    // check for CD (no DD required)
    final CollisionDetection detect = new CollisionDetection(whichCollDetect);

    final BondInfo marriedBonds =
        this.marryBonds(
            clusterCartes.getNoOfAtoms(), (clusterCartes.getNoOfAtoms() + cartEnv.getNoOfAtoms()));
    final boolean hasCollision = detect.checkOnlyForCollision(married, blowColl, marriedBonds);
    /*
     * important side-note: due to us actually *wanting* interactions between
     * the environment and the cluster, the blow factor above needs to be very
     * small, below 1.0, to actually just alert *real* clashes.
     */

    return !hasCollision;
  }

  @Override
  public Tuple<CartesianCoordinates, BondInfo> marryThem(
      final CartesianCoordinates clusterCartes, final BondInfo clusterBonds) {

    if (clusterCartes.containsEnvironment()) {
      System.err.println(
          "ERROR: Cluster cartesian already contains an environment. This is a bug. Contact author(s).");
      return new Tuple<>(clusterCartes, clusterBonds);
    }

    final CartesianCoordinates completeCartes = marryThem(clusterCartes);

    final int iNoOfAtoms = clusterCartes.getNoOfAtoms() + cartEnv.getNoOfAtoms();
    final BondInfo completeBonds =
        marryBonds(clusterBonds, clusterCartes.getNoOfAtoms(), iNoOfAtoms);

    return new Tuple<>(completeCartes, completeBonds);
  }

  @Override
  public CartesianCoordinates marryThem(final CartesianCoordinates clusterCartes) {

    if (clusterCartes.containsEnvironment()) {
      System.err.println(
          "ERROR: Cluster cartesian already contains an environment. This is a bug. Contact author(s).");
      return clusterCartes;
    }

    /*
     * we initially use the saved eulers and distance of COMs
     */
    clusterCartes.moveCoordsToCOM();
    final double[][] xyzCluster = clusterCartes.getAllXYZCoord();

    // rotation
    if (mode != FITMODE.LAYERONLY) {
      CoordTranslation.rotateXYZ(xyzCluster, eulersToCluster);
    }

    // translation
    final double[] trans = new double[3];
    for (int i = 0; i < 3; i++) {
      trans[i] = distanceCOMs[i] + referencePoints[0].getPosition()[i];
    }
    clusterCartes.moveCoordsToPoint(trans);

    // put the surface and the cluster together
    final int iNoOfAtoms = clusterCartes.getNoOfAtoms() + cartEnv.getNoOfAtoms();
    final int iNoOfMols = clusterCartes.getNoOfMolecules() + 1;
    final int[] iaAtsPerMol = new int[iNoOfMols];
    final int[] iaAtsPerMolCluster = clusterCartes.getAllAtomsPerMol();
    System.arraycopy(iaAtsPerMolCluster, 0, iaAtsPerMol, 0, iaAtsPerMolCluster.length);
    iaAtsPerMol[iaAtsPerMol.length - 1] = cartEnv.getNoOfAtoms();

    CartesianCoordinates completeCartes =
        new CartesianCoordinates(iNoOfAtoms, iNoOfMols, iaAtsPerMol);

    final double[][] xyzEnv = cartEnv.getAllXYZCoord();
    final double[][] xyzComplete = completeCartes.getAllXYZCoord();

    // copy around a little
    System.arraycopy(xyzCluster[0], 0, xyzComplete[0], 0, xyzCluster[0].length);
    System.arraycopy(xyzCluster[1], 0, xyzComplete[1], 0, xyzCluster[0].length);
    System.arraycopy(xyzCluster[2], 0, xyzComplete[2], 0, xyzCluster[0].length);
    System.arraycopy(xyzEnv[0], 0, xyzComplete[0], xyzCluster[0].length, xyzEnv[0].length);
    System.arraycopy(xyzEnv[1], 0, xyzComplete[1], xyzCluster[0].length, xyzEnv[0].length);
    System.arraycopy(xyzEnv[2], 0, xyzComplete[2], xyzCluster[0].length, xyzEnv[0].length);

    // the atom types
    final String[] saAtomsCompl = completeCartes.getAllAtomTypes();
    final String[] saAtomsCluster = clusterCartes.getAllAtomTypes();
    final String[] saAtomsEnv = cartEnv.getAllAtomTypes();

    System.arraycopy(saAtomsCluster, 0, saAtomsCompl, 0, saAtomsCluster.length);
    System.arraycopy(saAtomsEnv, 0, saAtomsCompl, saAtomsCluster.length, saAtomsEnv.length);

    // spins and charges
    final short[] iaSpinsCompl = completeCartes.getAllSpins();
    final float[] faChargesCompl = completeCartes.getAllCharges();

    final short[] iaSpinsCluster = clusterCartes.getAllSpins();
    final float[] faChargesCluster = clusterCartes.getAllCharges();

    System.arraycopy(iaSpinsCluster, 0, iaSpinsCompl, 0, iaSpinsCluster.length);
    System.arraycopy(faChargesCluster, 0, faChargesCompl, 0, faChargesCluster.length);

    // zmats
    final ZMatrix[] allZmatsCluster = clusterCartes.getZMatrices();
    final ZMatrix[] allZmatsCompl = new ZMatrix[allZmatsCluster.length + 1];

    for (int i = 0; i < allZmatsCluster.length; i++) {
      final ZMatrix currZMat = allZmatsCluster[i];
      if (currZMat == null) allZmatsCompl[i] = null;
      else allZmatsCompl[i] = new ZMatrix(currZMat);
    }

    // add one for this "molecule"
    allZmatsCompl[allZmatsCompl.length - 1] = null;

    completeCartes.setZMatrices(allZmatsCompl);

    // and of course we should also add a clone of this
    completeCartes.setRefEnvironment(copy());

    // recalc atom numbers
    completeCartes.recalcAtomNumbers();

    return completeCartes;
  }

  @Override
  public CartesianCoordinates divorceThem(final CartesianCoordinates completeCartes) {

    if (DEBUG) {
      final String[] printCartes = completeCartes.createPrintableCartesians();
      System.out.println("ENVIRONMENT AND CLUSTER PRIOR TO DIVORCE:");
      for (final String s : printCartes) {
        System.out.println(s);
      }
    }

    if (!completeCartes.containsEnvironment()) {
      System.err.println("ERROR: This is fatal and a bug. Contact the author(s).");
      return null;
    }

    final int[] iaComplAtsPerMol = completeCartes.getAllAtomsPerMol();
    final int[] iaClusterAtsPerMol = new int[iaComplAtsPerMol.length - 1];
    int iNoOfAtomsCluster = 0;
    for (int i = 0; i < iaClusterAtsPerMol.length; i++) {
      iaClusterAtsPerMol[i] = iaComplAtsPerMol[i];
      iNoOfAtomsCluster += iaComplAtsPerMol[i];
    }

    /*
     * who get's the car and the house?
     * affine transformation
     */

    // first find out where our references are now
    final double[] daRef1New =
        completeCartes.getXYZCoordinatesOfAtom(referencePoints[0].getID() + iNoOfAtomsCluster);
    final double[] daRef2New =
        completeCartes.getXYZCoordinatesOfAtom(referencePoints[1].getID() + iNoOfAtomsCluster);
    final double[] daRef3New =
        completeCartes.getXYZCoordinatesOfAtom(referencePoints[2].getID() + iNoOfAtomsCluster);

    final double[] daRef1Old = referencePoints[0].getPosition();
    final double[] daRef2Old = referencePoints[1].getPosition();
    final double[] daRef3Old = referencePoints[2].getPosition();

    final double diff1X = daRef1New[0] - daRef1Old[0];
    final double diff1Y = daRef1New[1] - daRef1Old[1];
    final double diff1Z = daRef1New[2] - daRef1Old[2];
    final double diff1 = Math.abs(Math.sqrt(diff1X * diff1X + diff1Y * diff1Y + diff1Z + diff1Z));

    final double diff2X = daRef2New[0] - daRef2Old[0];
    final double diff2Y = daRef2New[1] - daRef2Old[1];
    final double diff2Z = daRef2New[2] - daRef2Old[2];
    final double diff2 = Math.abs(Math.sqrt(diff2X * diff2X + diff2Y * diff2Y + diff2Z + diff2Z));

    final double diff3X = daRef3New[0] - daRef3Old[0];
    final double diff3Y = daRef3New[1] - daRef3Old[1];
    final double diff3Z = daRef3New[2] - daRef3Old[2];
    final double diff3 = Math.abs(Math.sqrt(diff3X * diff3X + diff3Y * diff3Y + diff3Z + diff3Z));

    // check if the environment position has changed significantly. should typically not be the
    // case, I think (i.e., when we run the environment constraint)
    boolean adjusted = false;
    if (diff1 > THRESHENPOSSAME || diff2 > THRESHENPOSSAME || diff3 > THRESHENPOSSAME) {

      // now build up two matrices
      final double[][] daMatT1 = new double[3][3];
      daMatT1[0][0] = daRef1Old[0];
      daMatT1[1][0] = daRef1Old[1];
      daMatT1[2][0] = daRef1Old[2];
      daMatT1[0][1] = daRef2Old[0];
      daMatT1[1][1] = daRef2Old[1];
      daMatT1[2][1] = daRef2Old[2];
      daMatT1[0][2] = daRef3Old[0];
      daMatT1[1][2] = daRef3Old[1];
      daMatT1[2][2] = daRef3Old[2];

      final double[][] daMatT2 = new double[3][3];
      daMatT2[0][0] = daRef1New[0];
      daMatT2[1][0] = daRef1New[1];
      daMatT2[2][0] = daRef1New[2];
      daMatT2[0][1] = daRef2New[0];
      daMatT2[1][1] = daRef2New[1];
      daMatT2[2][1] = daRef2New[2];
      daMatT2[0][2] = daRef3New[0];
      daMatT2[1][2] = daRef3New[1];
      daMatT2[2][2] = daRef3New[2];

      final Matrix matT1 = new Matrix(daMatT1);
      final Matrix matT2 = new Matrix(daMatT2);

      // calculate the transformation matrix
      final Matrix matTrans = matT2.times(matT1.inverse());

      // transform the cartesians
      final Matrix matCartes = new Matrix(completeCartes.getAllXYZCoord());
      final Matrix matNewCartes = matTrans.times(matCartes);

      completeCartes.setAllXYZ(matNewCartes.getArrayCopy());

      // adjust the reference atoms
      final double[] xyzRef1 =
          completeCartes.getXYZCoordinatesOfAtom(referencePoints[0].getID() + iNoOfAtomsCluster);
      daRef1Old[0] = xyzRef1[0];
      daRef1Old[1] = xyzRef1[1];
      daRef1Old[2] = xyzRef1[2];

      final double[] xyzRef2 =
          completeCartes.getXYZCoordinatesOfAtom(referencePoints[1].getID() + iNoOfAtomsCluster);
      daRef2Old[0] = xyzRef2[0];
      daRef2Old[1] = xyzRef2[1];
      daRef2Old[2] = xyzRef2[2];

      final double[] xyzRef3 =
          completeCartes.getXYZCoordinatesOfAtom(referencePoints[2].getID() + iNoOfAtomsCluster);
      daRef3Old[0] = xyzRef3[0];
      daRef3Old[1] = xyzRef3[1];
      daRef3Old[2] = xyzRef3[2];

      adjusted = true;
    }

    /*
     * now we can divorce
     */
    final CartesianCoordinates clusterCartes =
        new CartesianCoordinates(iNoOfAtomsCluster, iaClusterAtsPerMol.length, iaClusterAtsPerMol);

    // the spins and charges
    final short[] iaSpinsCluster = clusterCartes.getAllSpins();
    final float[] faChargesCluster = clusterCartes.getAllCharges();
    final short[] iaSpinsCompl = completeCartes.getAllSpins();
    final float[] faChargesCompl = completeCartes.getAllCharges();

    System.arraycopy(iaSpinsCompl, 0, iaSpinsCluster, 0, iNoOfAtomsCluster);
    System.arraycopy(faChargesCompl, 0, faChargesCluster, 0, iNoOfAtomsCluster);

    // the zmats
    final ZMatrix[] allComplZmats = completeCartes.getZMatrices();
    final ZMatrix[] allClusterZmats = new ZMatrix[allComplZmats.length - 1];
    for (int i = 0; i < allComplZmats.length - 1; i++) {
      final ZMatrix zmatTemp = allComplZmats[i];
      if (zmatTemp == null) allClusterZmats[i] = null;
      else allClusterZmats[i] = new ZMatrix(zmatTemp);
    }
    clusterCartes.setZMatrices(allClusterZmats);

    // coordinates and atom types
    final String[] saAtomsCluster = clusterCartes.getAllAtomTypes();
    final String[] saAtomsComplete = completeCartes.getAllAtomTypes();

    System.arraycopy(saAtomsComplete, 0, saAtomsCluster, 0, iNoOfAtomsCluster);

    final double[][] clusterCoords = clusterCartes.getAllXYZCoord();
    final double[][] completeCoords = completeCartes.getAllXYZCoord();

    System.arraycopy(completeCoords[0], 0, clusterCoords[0], 0, iNoOfAtomsCluster);
    System.arraycopy(completeCoords[1], 0, clusterCoords[1], 0, iNoOfAtomsCluster);
    System.arraycopy(completeCoords[2], 0, clusterCoords[2], 0, iNoOfAtomsCluster);

    if (mode != FITMODE.LAYERONLY) {
      // the euler angles can stay zero, if (and ONLY if) nothing gets rotated afterwards
      eulersToCluster[0] = 0.0;
      eulersToCluster[1] = 0.0;
      eulersToCluster[2] = 0.0;
    }

    // measure the translation from reference point [0] (daRef1Old) to the COM of the cluster
    // therefore, we first need to find out where the COM now is
    final double[] newCOM = clusterCartes.calculateTheCOM();
    for (int i = 0; i < 3; i++) {
      if (adjusted) {
        distanceCOMs[i] = newCOM[i] - daRef1Old[i];
      } else {
        distanceCOMs[i] = newCOM[i] - daRef1New[i];
      }
    }

    if (flexySurface) {
      // we need to update the cartesian in here, since something might have changed
      final int noSurfAts = cartEnv.getNoOfAtoms();
      final double[][] surfaceXYZ = cartEnv.getAllXYZCoord();
      System.arraycopy(completeCoords[0], iNoOfAtomsCluster, surfaceXYZ[0], 0, noSurfAts);
      System.arraycopy(completeCoords[1], iNoOfAtomsCluster, surfaceXYZ[1], 0, noSurfAts);
      System.arraycopy(completeCoords[2], iNoOfAtomsCluster, surfaceXYZ[2], 0, noSurfAts);
    }

    // some last minor modifications
    clusterCartes.setEnergy(completeCartes.getEnergy());
    clusterCartes.moveCoordsToCOM();

    return clusterCartes;
  }

  @Override
  public void initializeConnections(
      final CartesianCoordinates clusterCartes, final BondInfo clusterBonds) {

    final CollisionDetection collDetect = new CollisionDetection(whichCollDetect);
    final BondInfo marriedBonds =
        marryBonds(
            clusterCartes.getNoOfAtoms(), (clusterCartes.getNoOfAtoms() + cartEnv.getNoOfAtoms()));

    long emergCounter = 0;
    boolean hasColl;
    do {

      if (DEBUG) {
        System.out.println(
            "DEBUG: Trying to initialize connections for " + (emergCounter + 1) + " time");
      }

      if (mode != FITMODE.LAYERONLY) {
        // first the eulers
        RandomUtils.randomEulers(eulersToCluster);
      }

      // now the distances
      distanceCOMs = space.getPointInSpace();

      // check for collisions. DO NOT check for collisions within cluster!
      final CartesianCoordinates c = marryThem(clusterCartes);
      hasColl = collDetect.checkOnlyForCollision(c, blowColl, marriedBonds);
      emergCounter++;

      if (hasColl && DEBUG) {
        // get more info
        final CollisionInfo collInfo = collDetect.checkForCollision(c, blowColl, marriedBonds);
        final List<Collision> colls = collInfo.getCollisions();
        for (final Collision coll : colls) {
          System.out.println(
              "DEBUG: Collision between " + coll.getAtomOne() + " and " + coll.getAtomTwo());
        }
      }

    } while (hasColl && emergCounter < FixedValues.MAXTOEMERGENCY);

    assert (space.isPointInSpace(distanceCOMs));

    if (DEBUG) {
      System.out.println(
          "Took " + emergCounter + " attempts to initialize connections of cluster with env.");
    }
  }

  @Override
  public void mutateConnections(final CartesianCoordinates clusterCartes) {

    if (mode == FITMODE.LAYERONLY) {

      final double[] rndPoint = space.getPointInSpace();
      distanceCOMs[0] = rndPoint[0];
      distanceCOMs[1] = rndPoint[1];
      distanceCOMs[2] = rndPoint[2];

      return;
    }

    // first we want to figure out where to mutate
    final int mutPos = random.nextInt(6);

    double dr;
    final boolean br = random.nextBoolean();
    while (true) {
      dr = random.nextDouble();
      if (!(dr == 1.0 && br)) {
        break;
      } // else: we need to take a new double
    }

    switch (mutPos) {
      case 0:
      case 1:
      case 2:
        // COM-COM coordinates
        final double[] rndPoint = space.getPointInSpace();
        distanceCOMs[0] = rndPoint[0];
        distanceCOMs[1] = rndPoint[1];
        distanceCOMs[2] = rndPoint[2];
        break;
      case 3:
        // phi
        eulersToCluster[0] = (br) ? -dr * Math.PI : dr * Math.PI;
        break;
      case 4:
        // omega
        eulersToCluster[1] = (br) ? -dr * Math.PI / 2 : dr * Math.PI / 2;
        break;
      case 5:
        // psi
        eulersToCluster[2] = (br) ? -dr * Math.PI : dr * Math.PI;
        break;
      default:
        System.err.println(
            "ERROR: You should NEVER end up here in SimpleEnvironment. Please contact the author(s).");
        break;
    }
  }

  @Override
  public List<Environment> createOffspring(final Environment father) {

    final List<Environment> children = new ArrayList<>(2);

    final Environment envChild1 = father.copy();
    final Environment envChild2 = copy();

    final double[] fatherGenom = father.returnMyGenom();
    final double[] motherGenom = returnMyGenom();

    // we make a linearily distributed cut
    final int cut = random.nextInt(6);

    // copy over
    final double[] child1 = new double[6];
    final double[] child2 = new double[6];

    for (int i = 0; i < cut; i++) {
      child1[i] = fatherGenom[i];
      child2[i] = motherGenom[i];
    }

    for (int i = cut; i < child2.length; i++) {
      child1[i] = motherGenom[i];
      child2[i] = fatherGenom[i];
    }

    // put the new genes into the children
    envChild1.putGenomIn(child1);
    envChild2.putGenomIn(child2);

    // return the children
    children.add(envChild1);
    children.add(envChild2);

    return children;
  }

  @Override
  public List<Integer> whichSecondaryAtoms() {

    final List<Integer> secAtoms = new ArrayList<>();
    for (final int sec : listSecAtoms) {
      secAtoms.add(sec);
    }

    return secAtoms;
  }

  @Override
  public boolean isEnvironmentRigid() {
    return !flexySurface;
  }

  @Override
  public double[] returnMyGenom() {

    if (mode == FITMODE.LAYERONLY) {

      final double[] genom = new double[6];
      genom[0] = distanceCOMs[0];
      genom[1] = distanceCOMs[1];
      genom[2] = distanceCOMs[2];
      genom[3] = 0.0;
      genom[4] = 0.0;
      genom[5] = 0.0;

      return genom;
    } else {

      final double[] genom = new double[6];
      genom[0] = distanceCOMs[0];
      genom[1] = distanceCOMs[1];
      genom[2] = distanceCOMs[2];
      genom[3] = eulersToCluster[0];
      genom[4] = eulersToCluster[1];
      genom[5] = eulersToCluster[2];

      return genom;
    }
  }

  @Override
  public void putGenomIn(final double[] genom) {

    distanceCOMs[0] = genom[0];
    distanceCOMs[1] = genom[1];
    distanceCOMs[2] = genom[2];
    if (mode != FITMODE.LAYERONLY) {
      eulersToCluster[0] = genom[3];
      eulersToCluster[1] = genom[4];
      eulersToCluster[2] = genom[5];
    }
  }

  @Override
  public int atomsInEnv() {
    return cartEnv.getNoOfAtoms();
  }

  @Override
  public void moveCluster(final int direction, final double move) {
    assert (direction >= 0);
    assert (direction < 3);
    this.distanceCOMs[direction] += move;
  }

  @Override
  public CartesianCoordinates getEnvironmentCartes() {
    return this.cartEnv;
  }

  @Override
  public ENVIRONMENTTYPE getEnvironmentType() {
    return type;
  }

  /**
   * This is for internal use for the CD, it does NOT set the cluster bonds up properly.
   *
   * @param noClusterAts the number of cluster atoms
   * @param noTotalAts the total number of atoms in the system
   * @return a bond info object useful for CD applications.
   */
  protected BondInfo marryBonds(final int noClusterAts, final int noTotalAts) {

    final BondInfo bondsCompl = new SimpleBondInfo(noTotalAts);

    for (int i = 0; i < noClusterAts; i++) {
      for (int j = 0; j < noClusterAts; j++) {
        // XXX could be done with fill
        // it doesn't matter, what we fill in here since it's not about these distances anyways (but
        // there must be bonds)
        bondsCompl.setBond(i, j, BondInfo.UNCERTAIN);
      }
    }

    for (int i = noClusterAts; i < noTotalAts; i++) {
      for (int j = i; j < noTotalAts; j++) {
        if (envBonds.hasBond(i - noClusterAts, j - noClusterAts)) {
          final short bond = envBonds.bondType(i - noClusterAts, j - noClusterAts);
          bondsCompl.setBond(i, j, bond);
          bondsCompl.setBond(j, i, bond);
        } // default is no bond, so this is fine
      }
    }

    return bondsCompl;
  }

  protected BondInfo marryBonds(
      final BondInfo clusterBonds, final int noClusterAts, final int noTotalAts) {

    final BondInfo bondsCompl = new SimpleBondInfo(noTotalAts);

    for (int i = 0; i < noClusterAts; i++) {
      for (int j = i; j < noClusterAts; j++) {
        if (clusterBonds.hasBond(i, j)) {
          bondsCompl.setBond(i, j, clusterBonds.bondType(i, j));
          bondsCompl.setBond(j, i, clusterBonds.bondType(i, j));
        }
      }
    }

    // WE DO NOT SET ANY BONDS IN BETWEEN ENV AND CLUSTER. IF THIS IS WANTED, YOU'LL NEED TO ADJUST!

    for (int i = noClusterAts; i < noTotalAts; i++) {
      for (int j = i; j < noTotalAts; j++) {
        if (envBonds.hasBond(i - noClusterAts, j - noClusterAts)) {
          final short bond = envBonds.bondType(i - noClusterAts, j - noClusterAts);
          bondsCompl.setBond(i, j, bond);
          bondsCompl.setBond(j, i, bond);
        } // default is no bond, so this is fine
      }
    }

    return bondsCompl;
  }

  protected BondInfo divorceBonds(final BondInfo complete, final int noClusterAts) {

    final BondInfo clusterBonds = new SimpleBondInfo(noClusterAts);

    for (int i = 0; i < noClusterAts; i++) {
      for (int j = i; j < noClusterAts; j++) {
        if (complete.hasBond(i, j)) {
          final short bond = complete.bondType(i, j);
          clusterBonds.setBond(i, j, bond);
          clusterBonds.setBond(j, i, bond);
        }
      }
    }

    return clusterBonds;
  }
}
