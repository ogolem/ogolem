/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

import org.ogolem.generic.GenericSanityCheck;

/**
 * Checks a given geometry for post-local-optimization-sanity.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public final class GeometrySanityCheck implements GenericSanityCheck<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20160121;
  private static final boolean DEBUG = false;
  private final double blowBonds;
  private final double blowBondsEnv;
  private final boolean checkCollisions;
  private final boolean doDD;
  private final double blowDiss;

  public GeometrySanityCheck(
      final double blowBonds,
      final double blowBondsEnv,
      final double blowDiss,
      final boolean checkColl,
      final boolean doDD) {
    this.blowBonds = blowBonds;
    this.blowBondsEnv = blowBondsEnv;
    this.blowDiss = blowDiss;
    this.checkCollisions = checkColl;
    this.doDD = doDD;
  }

  GeometrySanityCheck(final GeometrySanityCheck orig) {
    this.blowBonds = orig.blowBonds;
    this.blowBondsEnv = orig.blowBondsEnv;
    this.blowDiss = orig.blowDiss;
    this.doDD = orig.doDD;
    this.checkCollisions = orig.checkCollisions;
  }

  @Override
  public GeometrySanityCheck copy() {
    return new GeometrySanityCheck(this);
  }

  @Override
  public boolean isSane(final Geometry individual) {

    final CartesianCoordinates cartes =
        (individual.containsEnvironment())
            ? individual.getCartesiansWithEnvironment()
            : individual.getCartesians();
    final BondInfo bonds = individual.getBondInfo();

    if (individual.containsEnvironment()) {
      if (!individual.doesFitWithEnvironment()) {
        return false;
      }
    }

    return checkSanity(cartes, bonds, blowBonds, blowBondsEnv, checkCollisions, doDD, blowDiss);
  }

  /**
   * Checks whether the provided geometry still has all intramolecular bonds as before. Scales
   * O(N^2) with N as the number of atoms. The arrays may be of different length, then the bond
   * array decides which interactions are checked. This is e.g. of interest when the cartesian still
   * contains an environment but the bond array doesn't.
   *
   * @param cartes A cartesian.
   * @param bonds The bonds in this cartesian.
   * @param blowBonds A blow factor.
   * @param blowBondsEnv A blow factor w.r.t. cluster/environment clashes (if there is one)
   * @return true if all bonds still exist, false otherwise.
   */
  public static boolean checkSanity(
      final CartesianCoordinates cartes,
      final BondInfo bonds,
      final double blowBonds,
      final double blowBondsEnv) {

    return checkSanity(cartes, bonds, blowBonds, blowBondsEnv, false, false, 4.0);
  }

  /**
   * Checks whether the provided geometry still has all intramolecular bonds as before. Scales
   * O(N^2) with N as the number of atoms. The arrays may be of different length, then the bond
   * array decides which interactions are checked. This is e.g. of interest when the cartesian still
   * contains an environment but the bond array doesn't.
   *
   * @param cartes A cartesian.
   * @param bonds The bonds in this cartesian
   * @param doDD do dissociation detection?
   * @param blowBonds A blow factor.
   * @param blowBondsEnv A blow factor w.r.t. cluster/environment clashes (if there is one)
   * @param checkCollisions If also collisions shall be checked.
   * @param blowDiss A blow factor for the dissociation check.
   * @return true if all bonds still exist, false otherwise.
   */
  public static boolean checkSanity(
      final CartesianCoordinates cartes,
      final BondInfo bonds,
      final double blowBonds,
      final double blowBondsEnv,
      final boolean checkCollisions,
      final boolean doDD,
      final double blowDiss) {

    if (!checkCollisions && !doDD) {
      // we are doing neither: do not waste time in here.
      return true;
    }

    if (DEBUG) {
      System.out.println("DEBUG: Sanity check working on the following cartesian...");
      final String[] sa = cartes.createPrintableCartesians();
      for (final String s : sa) {
        System.out.println(s);
      }
    }

    // see above why we use the bond lengths
    final int noEnvAtoms =
        (cartes.containsEnvironment()) ? cartes.getReferenceEnvironmentCopy().atomsInEnv() : 0;
    final int noOfAtomsCluster = cartes.getNoOfAtoms() - noEnvAtoms;
    final short[] nos = cartes.getAllAtomNumbers();
    final double[][] xyz = cartes.getAllXYZCoord();
    final boolean[][] adjacency = (!doDD) ? null : new boolean[noOfAtomsCluster][noOfAtomsCluster];

    final double[] radiiCluster = new double[noOfAtomsCluster];
    for (int i = 0; i < noOfAtomsCluster; i++) {
      radiiCluster[i] = AtomicProperties.giveRadius(nos[i]);
    }

    // first check CD within the cluster
    for (int i = 0; i < noOfAtomsCluster - 1; i++) {
      final double rad1 = radiiCluster[i];
      final double x0 = xyz[0][i];
      final double y0 = xyz[1][i];
      final double z0 = xyz[2][i];
      for (int j = i + 1; j < noOfAtomsCluster; j++) {
        final double rad2 = radiiCluster[j];
        final double radii = blowBonds * (rad1 + rad2);
        final double radiiSq = radii * radii;

        final double distX = x0 - xyz[0][j];
        final double distY = y0 - xyz[1][j];
        final double distZ = z0 - xyz[2][j];

        final double distSq = distX * distX + distY * distY + distZ * distZ;
        final boolean bond = bonds.hasBond(i, j);
        if (bond && distSq > radiiSq) {
          // there is no bond where there should be one
          if (DEBUG) {
            System.out.println("DEBUG: Backing out because of no bond. " + i + "\t" + j);
            System.out.println("DEBUG: dist " + Math.sqrt(distSq) + " bigger than " + radii);
          }
          return false;
        } else if (checkCollisions && !bond && distSq < radiiSq) {
          // there is a bond where there should be none
          if (DEBUG) {
            System.out.println("DEBUG: Backing out because of a bond. " + i + "\t" + j);
            System.out.println("DEBUG: dist " + Math.sqrt(distSq) + " smaller than " + radii);
          }
          return false;
        }

        if (doDD) {
          final double summedRadii = (rad1 + rad2) * blowDiss;
          final double summedRadiiSq = summedRadii * summedRadii;
          if (summedRadiiSq >= distSq) {
            adjacency[i][j] = true;
            adjacency[j][i] = true;
          } else {
            // XXX not really needed...
            adjacency[i][j] = false;
            adjacency[j][i] = false;
          }
        }
      }
    }

    // then the same for cluster -> env (NOT within env!)
    if (checkCollisions && noEnvAtoms > 0) {

      final double[] radiiEnv = new double[noEnvAtoms];
      for (int j = 0; j < noEnvAtoms; j++) {
        radiiEnv[j] = AtomicProperties.giveRadius(nos[noOfAtomsCluster + j]);
      }

      for (int i = 0; i < noOfAtomsCluster; i++) {
        final double rad1 = radiiCluster[i];
        final double x0 = xyz[0][i];
        final double y0 = xyz[1][i];
        final double z0 = xyz[2][i];
        for (int j = 0; j < noEnvAtoms; j++) {
          final double rad2 = radiiEnv[j];
          final double radii = blowBondsEnv * (rad1 + rad2);
          final double radiiSq = radii * radii;

          final double distX = x0 - xyz[0][noOfAtomsCluster + j];
          final double distY = y0 - xyz[1][noOfAtomsCluster + j];
          final double distZ = z0 - xyz[2][noOfAtomsCluster + j];

          final double distSq = distX * distX + distY * distY + distZ * distZ;

          // by definition, there are no *actual* bonds defined between cluster and env
          if (distSq < radiiSq) {
            // there is a clash
            if (DEBUG) {
              System.out.println("DEBUG: Backing out because of a clash. " + i + "\t" + j);
              System.out.println("DEBUG: dist " + Math.sqrt(distSq) + " smaller than " + radii);
            }
            return false;
          }
        }
      }
    }

    if (!doDD) {
      return true;
    }

    // check for DD using DFS (only within the cluster, not w.r.t. any potential environment)
    final boolean isConnected = DistanceCalc.dfsReachability(adjacency);
    if (DEBUG) {
      if (!isConnected) {
        System.out.println("DEBUG: Cluster is not connected.");
      }
    }

    return isConnected;
  }
}
