/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
              2015-2022, J. M. Dieterich and B. Hartke
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

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorSpecies;
import org.ogolem.math.BoolSymmetricMatrixNoDiag;
import org.ogolem.math.SymmetricMatrixNoDiag;

/**
 * This, in difference to the SimplePairWise collision detection engine does NOT exit on the first
 * collision detected but goes on and gathers more information together. The algorithm is in both
 * cases a pairwise checking.
 *
 * @author Johannes Dieterich
 * @version 2022-01-08
 */
public class AdvancedPairWise implements CollisionDetectionEngine {

  private static final long serialVersionUID = (long) 20210422;

  private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

  private final boolean exit;
  private final CollisionStrengthComputer comp;

  public AdvancedPairWise(final boolean exitOnFirst, final CollisionStrengthComputer comp) {
    assert (comp != null);
    this.exit = exitOnFirst;
    this.comp = comp;
  }

  private AdvancedPairWise(final AdvancedPairWise orig) {
    this.exit = orig.exit;
    this.comp = orig.comp.copy();
  }

  @Override
  public AdvancedPairWise copy() {
    return new AdvancedPairWise(this);
  }

  /**
   * @param cartesians A complete set of cartesian coordinates which should be checked for
   *     collisions.
   * @param blowFactor Bond detection works with this blow factor.
   * @param bonds Information on existing (and therefore wanted) bonds in the molecule.
   * @return info The returned {@code CollisionInfo}, an object keeping information not only on
   *     collisions but also on the calculated pairwise distances.
   */
  @Override
  public CollisionInfo checkForCollision(
      final CartesianCoordinates cartesians, final double blowFactor, final BondInfo bonds) {

    final CollisionInfo info = (exit) ? new SingleCollisionInfo() : new MultiCollisionInfo();
    checkForCollision(cartesians, blowFactor, bonds, info);

    return info;
  }

  @Override
  public void checkForCollision(
      final CartesianCoordinates cartesians,
      final double blowFactor,
      final BondInfo bonds,
      final CollisionInfo info) {

    if (bonds == null) {
      // apparently no one has cared to initialize this
      System.err.println("WARNING: No bond information found in AdvancedPairWise.");
      System.err.println("WARNING: Perhaps it is not initialized yet?");
      System.err.println("WARNING: Returning null'd collision information object.");
      return;
    }

    if (cartesians == null) {
      System.err.println("WARNING: Cartesians were null'd in AdvancedPairWise!");
      return;
    }

    info.resizeDistsAndClearState(cartesians.getNoOfAtoms());

    final int noOfAtoms = cartesians.getNoOfAtoms();
    final SymmetricMatrixNoDiag dists = info.getPairWiseDistances();

    final int noOfAtomsBonds = bonds.getNoOfAtoms();
    assert (noOfAtomsBonds >= noOfAtoms);

    final double[][] xyz = cartesians.getAllXYZCoord();
    final short[] numbers = cartesians.getAllAtomNumbers();

    assert (xyz != null);
    assert (numbers != null);
    assert (numbers.length >= noOfAtoms);
    assert (xyz.length == 3);
    assert (xyz[0].length >= noOfAtoms);
    assert (xyz[1].length >= noOfAtoms);
    assert (xyz[2].length >= noOfAtoms);

    if (!bonds.bondMatrixFast()) {
      System.err.println(
          "WARNING: The retrival of the full bond matrix is slow, reconsider to NOT use AdvancedPairWise or change implementation!");
    }
    final BoolSymmetricMatrixNoDiag bondMat = bonds.getFullBondMatrix();
    final boolean[] bondsBuffer = bondMat.underlyingStorageBuffer();

    // prefetch the radii in O(N) - the allocation here is not ideal
    // but caching it would make this whole object thread-unsafe
    final double[] radii = new double[noOfAtoms];
    for (int i = 0; i < noOfAtoms; i++) {
      radii[i] = blowFactor * AtomicProperties.giveRadius(numbers[i]);
    }

    final double[] distsBuffer = dists.underlyingStorageBuffer();
    int distsIdx = 0;

    final int idxPref = (noOfAtomsBonds * (noOfAtomsBonds - 1) / 2);

    boolean distCompl = true;
    Outer:
    for (int i = 0; i < noOfAtoms - 1; i++) {
      final double rad1 = radii[i];
      final var vRad1 = DoubleVector.broadcast(SPECIES, rad1);
      final double x = xyz[0][i];
      final double y = xyz[1][i];
      final double z = xyz[2][i];
      final var vX = DoubleVector.broadcast(SPECIES, x);
      final var vY = DoubleVector.broadcast(SPECIES, y);
      final var vZ = DoubleVector.broadcast(SPECIES, z);

      final int loopBound = SPECIES.loopBound(noOfAtoms - i - 1);
      int j = i + 1;

      int bondsIdx = idxPref - (noOfAtomsBonds - i) * ((noOfAtomsBonds - i) - 1) / 2;

      // vector loop
      for (; j < loopBound + i + 1; j += SPECIES.length()) {

        final var vRad2 = DoubleVector.fromArray(SPECIES, radii, j);
        final var vRadiiAdd = vRad2.add(vRad1);
        final var vDX = DoubleVector.fromArray(SPECIES, xyz[0], j).sub(vX);
        final var vDY = DoubleVector.fromArray(SPECIES, xyz[1], j).sub(vY);
        final var vDZ = DoubleVector.fromArray(SPECIES, xyz[2], j).sub(vZ);
        final var vFirst = vDX.fma(vDX, vDY.mul(vDY));
        final var vDistSq = vDZ.fma(vDZ, vFirst);
        final var vDist = vDistSq.sqrt();
        final var vComp = vDist.lt(vRadiiAdd);
        final var vNotBond = VectorMask.fromArray(SPECIES, bondsBuffer, bondsIdx).not();
        final var vMask = vComp.and(vNotBond);

        vDist.intoArray(distsBuffer, distsIdx);
        distsIdx += SPECIES.length();
        bondsIdx += SPECIES.length();

        final boolean anyColl = vMask.anyTrue();

        if (anyColl) {
          // collision - loop over mask and add collisions
          for (int v = 0; v < SPECIES.length(); v++) {
            if (!vMask.laneIsSet(v)) continue;
            final double strength =
                comp.calculateCollisionStrength(i, j + v, vDist.lane(v), vRadiiAdd.lane(v));
            final boolean succ = info.reportCollision(i, j + v, strength);
            if (!succ) {
              System.err.println("No success setting collision!");
            }
            if (exit) break; // just in case there is more than one collision in this vector lane
          }

          // increment the noOfCollisions AFTER setting the collision info since arrays start from
          // 0.
          if (exit) {
            distCompl = false;
            break Outer;
          }
        }
      }

      // cleanup loop
      for (; j < noOfAtoms; j++) {

        final double rad2 = radii[j];
        final double radiiAdd = rad1 + rad2;
        final double dX = x - xyz[0][j];
        final double dY = y - xyz[1][j];
        final double dZ = z - xyz[2][j];
        final double dist = Math.sqrt(dX * dX + dY * dY + dZ * dZ);
        distsBuffer[distsIdx] = dist;

        if (dist < radiiAdd && !bondsBuffer[bondsIdx]) {
          // collision
          final double strength = comp.calculateCollisionStrength(i, j, dist, radiiAdd);
          final boolean succ = info.reportCollision(i, j, strength);
          if (!succ) {
            System.err.println("No success setting collision!");
          }

          // increment the noOfCollisions AFTER setting the collision info since arrays start from
          // 0.
          if (exit) {
            distCompl = false;
            break Outer;
          }
        }

        distsIdx++;
        bondsIdx++;
      }
    }

    // set the pairwise distances to our collision info object
    info.setPairWiseDistances(dists, distCompl);
  }

  @Override
  public boolean checkOnlyForCollision(
      final CartesianCoordinates cartesians, final double blowFactor, final BondInfo bonds) {

    return checkOnlyForCollision(cartesians, blowFactor, bonds, 0, cartesians.getNoOfAtoms());
  }

  @Override
  public boolean checkOnlyForCollision(
      final CartesianCoordinates cartesians,
      final double blowFactor,
      final BondInfo bonds,
      final int offset,
      final int endset) {

    if (bonds == null) {
      // apparently no one has cared to initialize this
      System.err.println("WARNING: No bond information found in AdvancedPairWise.");
      System.err.println("WARNING: Perhaps it is not initialized yet?");
      System.err.println("WARNING: Returning null'd collision information object.");
      throw new RuntimeException("No bonds object given to AdvancedPairWise.");
    }

    if (cartesians == null) {
      System.err.println("WARNING: Cartesians were null'd in AdvancedPairWise!");
      throw new RuntimeException("No Cartesian coordinates given to AdvancedPairWise.");
    }

    final int noOfAtoms = cartesians.getNoOfAtoms();
    final double[][] xyz = cartesians.getAllXYZCoord();
    final short[] numbers = cartesians.getAllAtomNumbers();

    final int noOfAtomsBonds = bonds.getNoOfAtoms();
    assert (noOfAtomsBonds >= noOfAtoms);

    assert (xyz != null);
    assert (numbers != null);
    assert (numbers.length >= noOfAtoms);
    assert (xyz.length == 3);
    assert (xyz[0].length >= noOfAtoms);
    assert (xyz[1].length >= noOfAtoms);
    assert (xyz[2].length >= noOfAtoms);
    assert (offset >= 0);
    assert (endset <= noOfAtoms);

    if (!bonds.bondMatrixFast()) {
      System.err.println(
          "WARNING: The retrival of the full bond matrix is slow, reconsider to NOT use AdvancedPairWise or change implementation!");
    }
    final BoolSymmetricMatrixNoDiag bondMat = bonds.getFullBondMatrix();
    final boolean[] bondsBuffer = bondMat.underlyingStorageBuffer();

    // prefetch the radii in O(N) - the allocation here is not ideal
    // but caching it would make this whole object thread-unsafe
    // and fold blowfactor in
    final double[] radii = new double[noOfAtoms];
    for (int i = 0; i < noOfAtoms; i++) {
      radii[i] = AtomicProperties.giveRadius(numbers[i]) * blowFactor;
    }

    final int idxPref = (noOfAtomsBonds * (noOfAtomsBonds - 1) / 2);
    final int end = Math.min(noOfAtoms - 1, endset);
    for (int i = 0; i < end; i++) { // note: this is by spec
      final double rad1 = radii[i];
      final var vRad1 = DoubleVector.broadcast(SPECIES, rad1);
      final double x = xyz[0][i];
      final double y = xyz[1][i];
      final double z = xyz[2][i];
      final var vX = DoubleVector.broadcast(SPECIES, x);
      final var vY = DoubleVector.broadcast(SPECIES, y);
      final var vZ = DoubleVector.broadcast(SPECIES, z);
      final int start = Math.max(offset, i + 1);

      int bondsIdx =
          idxPref - (noOfAtomsBonds - i) * ((noOfAtomsBonds - i) - 1) / 2 + start - i - 1;

      final int loopBound = SPECIES.loopBound(noOfAtoms - start);
      int j = start;

      // vector loop
      for (; j < loopBound + start; j += SPECIES.length()) {
        final var vRad2 = DoubleVector.fromArray(SPECIES, radii, j);
        final var vRadiiAdd = vRad2.add(vRad1);
        final var vDX = DoubleVector.fromArray(SPECIES, xyz[0], j).sub(vX);
        final var vDY = DoubleVector.fromArray(SPECIES, xyz[1], j).sub(vY);
        final var vDZ = DoubleVector.fromArray(SPECIES, xyz[2], j).sub(vZ);
        final var vFirst = vDX.fma(vDX, vDY.mul(vDY));
        final var vDistSq = vDZ.fma(vDZ, vFirst);
        final var vComp = vDistSq.lt(vRadiiAdd.mul(vRadiiAdd));
        final var vNotBond = VectorMask.fromArray(SPECIES, bondsBuffer, bondsIdx).not();
        bondsIdx += SPECIES.length();
        final var vMask = vComp.and(vNotBond);
        // one might think andNot() would be reasonable here - performance craters then (4x in
        // microbench) - this may be a temporary issue, reevaluate later
        // final var vMask = vComp.andNot(VectorMask.fromArray(SPECIES, bondRow, j));
        final boolean anyColl = vMask.anyTrue();
        if (anyColl) return true; // at least one collision
      }

      // cleanup loop
      for (; j < noOfAtoms; j++) {
        final double rad2 = radii[j];
        final double radiiAdd = rad1 + rad2;
        final double dX = x - xyz[0][j];
        final double dY = y - xyz[1][j];
        final double dZ = z - xyz[2][j];
        final double distSq = dX * dX + dY * dY + dZ * dZ;
        if (distSq < radiiAdd * radiiAdd && !bondsBuffer[bondsIdx]) {
          // collision
          return true;
        }
        bondsIdx++;
      }
    }

    return false;
  }
}
