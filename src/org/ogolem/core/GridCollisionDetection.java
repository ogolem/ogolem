/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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

import java.util.Arrays;
import java.util.LinkedList;
import org.ogolem.math.SymmetricMatrixNoDiag;

/**
 * This is a O(N) scaling method working with a grid to find collisions in the geometry.
 *
 * @author Johannes Dieterich
 * @version 2020-12-21
 */
class GridCollisionDetection implements CollisionDetectionEngine {

  private static final long serialVersionUID = (long) 20101209;

  // TODO test and debug
  private final boolean bExitOnFirstClash;

  // TODO needs to be made tunable, even if this requires a major code rewrite
  private final double dCellDimension = 7;

  private final CollisionStrengthComputer comp;

  GridCollisionDetection(
      final boolean bExitOnFirstCollision, final CollisionStrengthComputer comp) {
    this.bExitOnFirstClash = bExitOnFirstCollision;
    this.comp = comp;
  }

  private GridCollisionDetection(final GridCollisionDetection orig) {
    this.comp = orig.comp.copy();
    this.bExitOnFirstClash = orig.bExitOnFirstClash;
  }

  @Override
  public GridCollisionDetection copy() {
    return new GridCollisionDetection(this);
  }

  @Override
  public CollisionInfo checkForCollision(
      final CartesianCoordinates cartes, final double dBlowFactor, final BondInfo bonds) {

    final CollisionInfo info =
        (bExitOnFirstClash) ? new SingleCollisionInfo() : new MultiCollisionInfo();
    checkForCollision(cartes, dBlowFactor, bonds, info);

    return info;
  }

  @Override
  public void checkForCollision(
      final CartesianCoordinates cartes,
      final double dBlowFactor,
      final BondInfo bonds,
      final CollisionInfo info) {

    final double[][] daAllCartes = cartes.getAllXYZCoord();
    /*
     * first we need to set the grid up, for that we need the maximal and
     * minimal cartesian coordinates.
     */
    double dMaxX = Double.NEGATIVE_INFINITY;
    double dMinX = Double.POSITIVE_INFINITY;
    double dMaxY = Double.NEGATIVE_INFINITY;
    double dMinY = Double.POSITIVE_INFINITY;
    double dMaxZ = Double.NEGATIVE_INFINITY;
    double dMinZ = Double.POSITIVE_INFINITY;
    for (int i = 0; i < daAllCartes[0].length; i++) {
      // x
      if (daAllCartes[0][i] > dMaxX) {
        dMaxX = daAllCartes[0][i];
      }
      if (daAllCartes[0][i] < dMinX) {
        dMinX = daAllCartes[0][i];
      }

      // y
      if (daAllCartes[1][i] > dMaxY) {
        dMaxY = daAllCartes[1][i];
      }
      if (daAllCartes[1][i] < dMinY) {
        dMinY = daAllCartes[1][i];
      }

      // z
      if (daAllCartes[2][i] > dMaxZ) {
        dMaxZ = daAllCartes[2][i];
      }
      if (daAllCartes[2][i] < dMinZ) {
        dMinZ = daAllCartes[2][i];
      }
    }

    /*
     * now we set the small cells up
     */
    final double dLengthX = Math.abs(dMaxX - dMinX);
    final double dLengthY = Math.abs(dMaxY - dMinY);
    final double dLengthZ = Math.abs(dMaxZ - dMinZ);

    final int iNoOfCellsX = (int) Math.ceil(dLengthX / dCellDimension);
    final int iNoOfCellsY = (int) Math.ceil(dLengthY / dCellDimension);
    final int iNoOfCellsZ = (int) Math.ceil(dLengthZ / dCellDimension);

    // set the container up
    final Object[][][] oaCells = new Object[iNoOfCellsX][iNoOfCellsY][iNoOfCellsZ];

    // populate the cells
    for (int x = 0; x < iNoOfCellsX; x++) {
      for (int y = 0; y < iNoOfCellsY; y++) {
        for (int z = 0; z < iNoOfCellsZ; z++) {
          oaCells[x][y][z] = new LinkedList<>();
        }
      }
    }

    /*
     * assign the atoms to the cells
     */
    for (int i = 0; i < cartes.getNoOfAtoms(); i++) {
      // TODO the assigning might be broken and might need bettering.
      final int x = (int) Math.round(Math.abs(daAllCartes[0][i] - dMinX) / dCellDimension);
      final int y = (int) Math.round(Math.abs(daAllCartes[1][i] - dMinY) / dCellDimension);
      final int z = (int) Math.round(Math.abs(daAllCartes[2][i] - dMinZ) / dCellDimension);
      // put it in
      final Atom atom = cartes.getAtomAtPos(i);
      // we can safely suppress the unchecked warning
      @SuppressWarnings(value = "unchecked")
      final LinkedList<Atom> llTemp = (LinkedList<Atom>) oaCells[x][y][z];
      llTemp.add(atom);
      oaCells[x][y][z] = llTemp;
    }

    final int iNoOfAtoms = cartes.getNoOfAtoms();
    info.resizeDistsAndClearState(iNoOfAtoms);
    final SymmetricMatrixNoDiag daDistances = info.getPairWiseDistances();

    // first "initialize" all distances with big values
    final double[] buffer = daDistances.underlyingStorageBuffer();
    Arrays.fill(buffer, 1E6);

    /*
     * now compute the actual distances of atoms in the same box or in neighboring ones
     * first in the same box
     */
    int iNoOfCollisions = 0;

    for (int i = 0; i < iNoOfCellsX; i++) {
      for (int j = 0; j < iNoOfCellsY; j++) {
        for (int k = 0; k < iNoOfCellsZ; k++) {
          // get all atoms as the linked list
          // we can safely suppress the unchecked warning
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTemp = (LinkedList<Atom>) oaCells[i][j][k];

          final int iNoOfAtsInCell = llTemp.size();
          for (int at = 0; at < iNoOfAtsInCell; at++) {
            // the running index is on purpose one too small
            Atom atOne = llTemp.get(at);
            Atom atTwo = llTemp.get(at + 1);

            final double[] daPosOne = atOne.getPosition();
            final double[] daPosTwo = atTwo.getPosition();

            final double dDistTemp =
                Math.sqrt(
                    Math.pow((daPosOne[0] - daPosTwo[0]), 2)
                        + Math.pow((daPosOne[1] - daPosTwo[1]), 2)
                        + Math.pow((daPosOne[2] - daPosTwo[2]), 2));
            final double dRadii = (atOne.getRadius() + atTwo.getRadius()) * dBlowFactor;

            // put the distance in
            final int iIDOne = atOne.getID();
            final int iIDTwo = atTwo.getID();
            daDistances.setElement(iIDOne, iIDTwo, dDistTemp);

            // check for collision
            if (dDistTemp <= dRadii && !bonds.hasBond(iIDOne, iIDTwo)) {
              // we have a clash
              if (bExitOnFirstClash) {
                // put the collision in and return the info
                info.setPairWiseDistances(daDistances, false);
                final double strength =
                    comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                info.reportCollision(iIDOne, iIDTwo, strength);

                return;
              } else {
                // put a collision in
                final double strength =
                    comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                info.reportCollision(iIDOne, iIDTwo, strength);
                iNoOfCollisions++;
              }
            }
          }
        }
      }
    }

    // now all to the right in x dimension (but the last one)
    for (int i = 0; i < iNoOfCellsX - 1; i++) {
      for (int j = 0; j < iNoOfCellsY; j++) {
        for (int k = 0; k < iNoOfCellsZ; k++) {
          // get all atoms in box one as the linked list
          // we can safely suppress the unchecked warning
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempOne = (LinkedList<Atom>) oaCells[i][j][k];
          // same procedure for the box to the right
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempTwo = (LinkedList<Atom>) oaCells[i + 1][j][k];

          final int iNoOfAtsCell1 = llTempOne.size();
          final int iNoOfAtsCell2 = llTempTwo.size();
          // calculate all pairwise distances
          for (int at1 = 0; at1 < iNoOfAtsCell1 + 1; at1++) {
            for (int at2 = 0; at2 < iNoOfAtsCell2 + 1; at2++) {
              Atom atOne = llTempOne.get(at1);
              Atom atTwo = llTempTwo.get(at2);
              final double[] daPosOne = atOne.getPosition();
              final double[] daPosTwo = atTwo.getPosition();

              final double dDistTemp =
                  Math.sqrt(
                      Math.pow((daPosOne[0] - daPosTwo[0]), 2)
                          + Math.pow((daPosOne[1] - daPosTwo[1]), 2)
                          + Math.pow((daPosOne[2] - daPosTwo[2]), 2));
              final double dRadii = (atOne.getRadius() + atTwo.getRadius()) * dBlowFactor;

              // put the distance in
              final int iIDOne = atOne.getID();
              final int iIDTwo = atTwo.getID();
              daDistances.setElement(iIDOne, iIDTwo, dDistTemp);

              // check for collision
              if (dDistTemp <= dRadii && !bonds.hasBond(iIDOne, iIDTwo)) {
                // we have a clash
                if (bExitOnFirstClash) {
                  // put the collision in and return the info
                  info.setPairWiseDistances(daDistances, false);
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);

                  return;
                } else {
                  // put a collision in
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);
                  iNoOfCollisions++;
                }
              }
            }
          }
        }
      }
    }

    // now all to the right in y dimension
    for (int i = 0; i < iNoOfCellsX; i++) {
      for (int j = 0; j < iNoOfCellsY - 1; j++) {
        for (int k = 0; k < iNoOfCellsZ; k++) {
          // get all atoms in box one as the linked list
          // we can safely suppress the unchecked warning
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempOne = (LinkedList<Atom>) oaCells[i][j][k];
          // same procedure for the box to the right
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempTwo = (LinkedList<Atom>) oaCells[i][j + 1][k];

          final int iNoOfAtsCell1 = llTempOne.size();
          final int iNoOfAtsCell2 = llTempTwo.size();
          // calculate all pairwise distances
          for (int at1 = 0; at1 < iNoOfAtsCell1 + 1; at1++) {
            for (int at2 = 0; at2 < iNoOfAtsCell2 + 1; at2++) {
              Atom atOne = llTempOne.get(at1);
              Atom atTwo = llTempTwo.get(at2);
              final double[] daPosOne = atOne.getPosition();
              final double[] daPosTwo = atTwo.getPosition();

              final double dDistTemp =
                  Math.sqrt(
                      Math.pow((daPosOne[0] - daPosTwo[0]), 2)
                          + Math.pow((daPosOne[1] - daPosTwo[1]), 2)
                          + Math.pow((daPosOne[2] - daPosTwo[2]), 2));
              final double dRadii = (atOne.getRadius() + atTwo.getRadius()) * dBlowFactor;

              // put the distance in
              final int iIDOne = atOne.getID();
              final int iIDTwo = atTwo.getID();
              daDistances.setElement(iIDOne, iIDTwo, dDistTemp);

              // check for collision
              if (dDistTemp <= dRadii && !bonds.hasBond(iIDOne, iIDTwo)) {
                // we have a clash
                if (bExitOnFirstClash) {
                  // put the collision in and return the info
                  info.setPairWiseDistances(daDistances, false);
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);
                  return;
                } else {
                  // put a collision in
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);
                  iNoOfCollisions++;
                }
              }
            }
          }
        }
      }
    }

    // now all to the right in z dimension
    for (int i = 0; i < iNoOfCellsX; i++) {
      for (int j = 0; j < iNoOfCellsY; j++) {
        for (int k = 0; k < iNoOfCellsZ - 1; k++) {
          // get all atoms in box one as the linked list
          // we can safely suppress the unchecked warning
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempOne = (LinkedList<Atom>) oaCells[i][j][k];
          // same procedure for the box to the right
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempTwo = (LinkedList<Atom>) oaCells[i][j][k + 1];

          final int iNoOfAtsCell1 = llTempOne.size();
          final int iNoOfAtsCell2 = llTempTwo.size();
          // calculate all pairwise distances
          for (int at1 = 0; at1 < iNoOfAtsCell1 + 1; at1++) {
            for (int at2 = 0; at2 < iNoOfAtsCell2 + 1; at2++) {
              Atom atOne = llTempOne.get(at1);
              Atom atTwo = llTempTwo.get(at2);
              final double[] daPosOne = atOne.getPosition();
              final double[] daPosTwo = atTwo.getPosition();

              final double dDistTemp =
                  Math.sqrt(
                      Math.pow((daPosOne[0] - daPosTwo[0]), 2)
                          + Math.pow((daPosOne[1] - daPosTwo[1]), 2)
                          + Math.pow((daPosOne[2] - daPosTwo[2]), 2));
              final double dRadii = (atOne.getRadius() + atTwo.getRadius()) * dBlowFactor;

              // put the distance in
              final int iIDOne = atOne.getID();
              final int iIDTwo = atTwo.getID();
              daDistances.setElement(iIDOne, iIDTwo, dDistTemp);

              // check for collision
              if (dDistTemp <= dRadii && !bonds.hasBond(iIDOne, iIDTwo)) {
                // we have a clash
                if (bExitOnFirstClash) {
                  // put the collision in and return the info
                  info.setPairWiseDistances(daDistances, false);
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);
                  return;
                } else {
                  // put a collision in
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);
                  iNoOfCollisions++;
                }
              }
            }
          }
        }
      }
    }

    // slightly more difficult: all diogonal shifted (so [x+1][y+1][z])
    for (int i = 0; i < iNoOfCellsX - 1; i++) {
      for (int j = 0; j < iNoOfCellsY - 1; j++) {
        for (int k = 0; k < iNoOfCellsZ; k++) {
          // get all atoms in box one as the linked list
          // we can safely suppress the unchecked warning
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempOne = (LinkedList<Atom>) oaCells[i][j][k];
          // same procedure for the box to the right
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempTwo = (LinkedList<Atom>) oaCells[i + 1][j + 1][k];

          final int iNoOfAtsCell1 = llTempOne.size();
          final int iNoOfAtsCell2 = llTempTwo.size();
          // calculate all pairwise distances
          for (int at1 = 0; at1 < iNoOfAtsCell1 + 1; at1++) {
            for (int at2 = 0; at2 < iNoOfAtsCell2 + 1; at2++) {
              Atom atOne = llTempOne.get(at1);
              Atom atTwo = llTempTwo.get(at2);
              final double[] daPosOne = atOne.getPosition();
              final double[] daPosTwo = atTwo.getPosition();

              final double dDistTemp =
                  Math.sqrt(
                      Math.pow((daPosOne[0] - daPosTwo[0]), 2)
                          + Math.pow((daPosOne[1] - daPosTwo[1]), 2)
                          + Math.pow((daPosOne[2] - daPosTwo[2]), 2));
              final double dRadii = (atOne.getRadius() + atTwo.getRadius()) * dBlowFactor;

              // put the distance in
              final int iIDOne = atOne.getID();
              final int iIDTwo = atTwo.getID();
              daDistances.setElement(iIDOne, iIDTwo, dDistTemp);

              // check for collision
              if (dDistTemp <= dRadii && !bonds.hasBond(iIDOne, iIDTwo)) {
                // we have a clash
                if (bExitOnFirstClash) {
                  // put the collision in and return the info
                  info.setPairWiseDistances(daDistances, false);
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);

                  return;
                } else {
                  // put a collision in
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);
                  iNoOfCollisions++;
                }
              }
            }
          }
        }
      }
    }

    // also slightly more difficult: all diogonal up (so [x][y+1][z+1])
    for (int i = 0; i < iNoOfCellsX; i++) {
      for (int j = 0; j < iNoOfCellsY - 1; j++) {
        for (int k = 0; k < iNoOfCellsZ - 1; k++) {
          // get all atoms in box one as the linked list
          // we can safely suppress the unchecked warning
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempOne = (LinkedList<Atom>) oaCells[i][j][k];
          // same procedure for the box to the right
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempTwo = (LinkedList<Atom>) oaCells[i][j + 1][k + 1];

          final int iNoOfAtsCell1 = llTempOne.size();
          final int iNoOfAtsCell2 = llTempTwo.size();
          // calculate all pairwise distances
          for (int at1 = 0; at1 < iNoOfAtsCell1 + 1; at1++) {
            for (int at2 = 0; at2 < iNoOfAtsCell2 + 1; at2++) {
              Atom atOne = llTempOne.get(at1);
              Atom atTwo = llTempTwo.get(at2);
              final double[] daPosOne = atOne.getPosition();
              final double[] daPosTwo = atTwo.getPosition();

              final double dDistTemp =
                  Math.sqrt(
                      Math.pow((daPosOne[0] - daPosTwo[0]), 2)
                          + Math.pow((daPosOne[1] - daPosTwo[1]), 2)
                          + Math.pow((daPosOne[2] - daPosTwo[2]), 2));
              final double dRadii = (atOne.getRadius() + atTwo.getRadius()) * dBlowFactor;

              // put the distance in
              final int iIDOne = atOne.getID();
              final int iIDTwo = atTwo.getID();
              daDistances.setElement(iIDOne, iIDTwo, dDistTemp);

              // check for collision
              if (dDistTemp <= dRadii && !bonds.hasBond(iIDOne, iIDTwo)) {
                // we have a clash
                if (bExitOnFirstClash) {
                  // put the collision in and return the info
                  info.setPairWiseDistances(daDistances, false);
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);

                  return;
                } else {
                  // put a collision in
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);
                  iNoOfCollisions++;
                }
              }
            }
          }
        }
      }
    }

    // also slightly more difficult: all diogonal up (so [x+1][y][z+1])
    for (int i = 0; i < iNoOfCellsX - 1; i++) {
      for (int j = 0; j < iNoOfCellsY; j++) {
        for (int k = 0; k < iNoOfCellsZ - 1; k++) {
          // get all atoms in box one as the linked list
          // we can safely suppress the unchecked warning
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempOne = (LinkedList<Atom>) oaCells[i][j][k];
          // same procedure for the box to the right
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempTwo = (LinkedList<Atom>) oaCells[i + 1][j][k + 1];

          final int iNoOfAtsCell1 = llTempOne.size();
          final int iNoOfAtsCell2 = llTempTwo.size();
          // calculate all pairwise distances
          for (int at1 = 0; at1 < iNoOfAtsCell1 + 1; at1++) {
            for (int at2 = 0; at2 < iNoOfAtsCell2 + 1; at2++) {
              Atom atOne = llTempOne.get(at1);
              Atom atTwo = llTempTwo.get(at2);
              final double[] daPosOne = atOne.getPosition();
              final double[] daPosTwo = atTwo.getPosition();

              final double dDistTemp =
                  Math.sqrt(
                      Math.pow((daPosOne[0] - daPosTwo[0]), 2)
                          + Math.pow((daPosOne[1] - daPosTwo[1]), 2)
                          + Math.pow((daPosOne[2] - daPosTwo[2]), 2));
              final double dRadii = (atOne.getRadius() + atTwo.getRadius()) * dBlowFactor;

              // put the distance in
              final int iIDOne = atOne.getID();
              final int iIDTwo = atTwo.getID();
              daDistances.setElement(iIDOne, iIDTwo, dDistTemp);

              // check for collision
              if (dDistTemp <= dRadii && !bonds.hasBond(iIDOne, iIDTwo)) {
                // we have a clash
                if (bExitOnFirstClash) {
                  // put the collision in and return the info
                  info.setPairWiseDistances(daDistances, false);
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);

                  return;
                } else {
                  // put a collision in
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);
                  iNoOfCollisions++;
                }
              }
            }
          }
        }
      }
    }

    // even more difficult: all diagonal up shifted (so [x+1][y+1][z+1])
    for (int i = 0; i < iNoOfCellsX - 1; i++) {
      for (int j = 0; j < iNoOfCellsY - 1; j++) {
        for (int k = 0; k < iNoOfCellsZ - 1; k++) {
          // get all atoms in box one as the linked list
          // we can safely suppress the unchecked warning
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempOne = (LinkedList<Atom>) oaCells[i][j][k];
          // same procedure for the box to the right
          @SuppressWarnings(value = "unchecked")
          LinkedList<Atom> llTempTwo = (LinkedList<Atom>) oaCells[i + 1][j + 1][k + 1];

          final int iNoOfAtsCell1 = llTempOne.size();
          final int iNoOfAtsCell2 = llTempTwo.size();
          // calculate all pairwise distances
          for (int at1 = 0; at1 < iNoOfAtsCell1 + 1; at1++) {
            for (int at2 = 0; at2 < iNoOfAtsCell2 + 1; at2++) {
              Atom atOne = llTempOne.get(at1);
              Atom atTwo = llTempTwo.get(at2);
              final double[] daPosOne = atOne.getPosition();
              final double[] daPosTwo = atTwo.getPosition();

              final double dDistTemp =
                  Math.sqrt(
                      Math.pow((daPosOne[0] - daPosTwo[0]), 2)
                          + Math.pow((daPosOne[1] - daPosTwo[1]), 2)
                          + Math.pow((daPosOne[2] - daPosTwo[2]), 2));
              final double dRadii = (atOne.getRadius() + atTwo.getRadius()) * dBlowFactor;

              // put the distance in
              final int iIDOne = atOne.getID();
              final int iIDTwo = atTwo.getID();
              daDistances.setElement(iIDOne, iIDTwo, dDistTemp);

              // check for collision
              if (dDistTemp <= dRadii && !bonds.hasBond(iIDOne, iIDTwo)) {
                // we have a clash
                if (bExitOnFirstClash) {
                  // put the collision in and return the info
                  info.setPairWiseDistances(daDistances, false);
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);

                  return;
                } else {
                  // put a collision in
                  final double strength =
                      comp.calculateCollisionStrength(iIDOne, iIDTwo, dRadii, dDistTemp);
                  info.reportCollision(iIDOne, iIDTwo, strength);
                  iNoOfCollisions++;
                }
              }
            }
          }
        }
      }
    }

    // ok, if we reached here we either found no collision, or, we didn't exit on the first one
    info.setPairWiseDistances(daDistances, true);
    System.out.println("DEBUG: DONE WITH GRID COLLDETECT!");
  }

  @Override
  public boolean checkOnlyForCollision(
      CartesianCoordinates cartesians, double blowFactor, BondInfo bonds) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public boolean checkOnlyForCollision(
      final CartesianCoordinates cartesians,
      final double blowFactor,
      final BondInfo bonds,
      final int offset,
      final int endset) {
    throw new UnsupportedOperationException("Not supported yet.");
  }
}
