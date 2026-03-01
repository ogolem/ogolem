/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2020-2026, J. M. Dieterich and B. Hartke
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
import java.util.Arrays;
import java.util.List;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This is the toolbox for the global optimization of geometries. It is supposed to help eliminating
 * redundancies in the algorithms. Some redundancy is "normal" but if it should exceed a certain
 * "uncomfortable" level, consider moving parts into functions and (as a next solution) into more
 * classes. Ideally a globoptatomics subpackage.
 *
 * @author Johannes Dieterich
 * @version 2026-02-26
 */
class GlobOptAtomics {

  private static final Logger log = LoggerFactory.getLogger(GlobOptAtomics.class);

  /**
   * Assuming, that the molecular size is almost proportional to the number of atoms in a molecule,
   * this algorithm orders the molecules by their size.
   *
   * @param geom The geometry whos sizes should be assigned.
   * @return Returns the molecular IDs ordered by size.
   */
  static List<Integer> calculateMolecularSizes(final Geometry geom) {

    final int noMols = geom.getNumberOfIndieParticles();

    final List<Integer> sizes = new ArrayList<>(noMols);
    final List<Integer> atomNos = new ArrayList<>(noMols);

    // we put the first molecule in the first spot.
    sizes.add(0);
    atomNos.add(geom.getMoleculeAtPosition(0).getNumberOfAtoms());

    for (int molID = 1; molID < noMols; molID++) {
      // we check for the sizes
      final int noAtoms = geom.getMoleculeAtPosition(molID).getNumberOfAtoms();
      assert (noAtoms > 0);

      final int length = sizes.size();
      for (int j = 0; j < length; j++) {
        if (noAtoms > atomNos.get(j)) {
          // put it before and break the loop
          sizes.add(j, molID);
          atomNos.add(j, noAtoms);
          break;
        }
        if (j == length - 1) {
          // put it in the end
          sizes.add(molID);
          atomNos.add(noAtoms);
        }
      }
    }

    return sizes;
  }

  static int[] assignMolecularTypes(final List<Molecule> molecules) {

    final int noMols = molecules.size();
    final int[] molTypes = new int[noMols];

    final List<String> availTypes = new ArrayList<>();

    int typeCounter = 0;
    for (int i = 0; i < noMols; i++) {

      final String type = molecules.get(i).getSID();
      assert (type != null);
      assert (!type.equalsIgnoreCase("N/A"));

      if (i == 0) {
        // must be a new type
        availTypes.add(type);
        molTypes[i] = typeCounter;
        typeCounter++;
        continue;
      }

      for (int j = 0; j < availTypes.size(); j++) {
        if (availTypes.get(j).equalsIgnoreCase(type)) {
          // no new type
          molTypes[i] = j;
          break;
        }
        if (j == availTypes.size() - 1) {
          // new type
          availTypes.add(type);
          molTypes[i] = typeCounter;
          typeCounter++;
          break;
        }
      }
    }

    return molTypes;
  }

  static enum CUTTINGMODE {
    ZEROZ,
    GAUSSDISTR,
    NORMDISTR
  };

  /**
   * Returns a cutting plane height that is guaranteed to have at least one COM strictly below and
   * one strictly above (no retry loop). Uses sorted COMs and places the plane between two
   * consecutive values. More stable than random sampling when COMs are clustered or when GAUSSDISTR
   * would often sample outside [min,max].
   *
   * @param whichMode cutting mode: ZEROZ uses 0 if it lies between min and max, else deterministic;
   *     GAUSSDISTR/NORMDISTR choose a random split index in [1, n-1] with plane in the gap.
   * @param coms1 COMs[coord][particle] (length 3 for x,y,z).
   * @param random RNG for choosing split index and optional jitter.
   * @param coord axis index (0=x, 1=y, 2=z).
   * @return plane height, or Double.NaN if all COMs have the same coordinate (no valid split).
   */
  static double guaranteedCuttingPlaneHeight(
      final CUTTINGMODE whichMode, final double[][] coms1, final Lottery random, final int coord) {

    assert (coms1 != null);
    assert (coms1.length == 3);
    assert (coord >= 0);
    assert (coord < 3);

    final int n = coms1[coord].length;
    if (n < 2) {
      log.warn(
          "Cutting plane unavailable: fewer than 2 particles (n="
              + n
              + ", axis="
              + axisName(coord)
              + ")");
      return Double.NaN;
    }

    final double[] sorted = coms1[coord].clone();
    Arrays.sort(sorted);
    final double min = sorted[0];
    final double max = sorted[n - 1];
    final double range = max - min;
    if (Math.abs(range) < 1E-8) {
      log.warn(
          "Cutting plane unavailable: all COMs have same coordinate on axis "
              + axisName(coord)
              + " (value="
              + min
              + ", range="
              + range
              + ", tolerance=1E-8)");
      return Double.NaN;
    }

    if (whichMode == CUTTINGMODE.ZEROZ) {
      if (min < 0.0 && 0.0 < max) return 0.0;
    }

    // Choose how many COMs should be below the plane (1 .. n-1).
    int noUnder;
    if (whichMode == CUTTINGMODE.GAUSSDISTR) {
      // Gaussian-like choice around center: map gauss in (-1,1) to index in [1, n-1]
      final double g = RandomUtils.gaussDouble(-1.0, 1.0);
      final double t = (g + 1.0) * 0.5; // in [0,1]
      noUnder = 1 + (int) (t * (n - 1));
      if (noUnder > n - 1) noUnder = n - 1;
      if (noUnder < 1) noUnder = 1;
    } else {
      // NORMDISTR or ZEROZ fallback: uniform split index
      noUnder = 1 + random.nextInt(n - 1);
    }

    final double below = sorted[noUnder - 1];
    final double above = sorted[noUnder];
    final double gap = above - below;
    if (Math.abs(gap) < 1E-8) {
      log.warn(
          "Cutting plane unavailable: degenerate gap at split (axis="
              + axisName(coord)
              + ", noUnder="
              + noUnder
              + ", n="
              + n
              + ", below="
              + below
              + ", above="
              + above
              + ", gap="
              + gap
              + ")");
      return Double.NaN;
    }
    // midpoint
    final double cuttingPlaneHeight = below + gap * 0.5;

    return cuttingPlaneHeight;
  }

  static double guaranteedCuttingZPlane(
      final CUTTINGMODE whichMode, final double[][] coms1, final Lottery random) {
    return guaranteedCuttingPlaneHeight(whichMode, coms1, random, 2);
  }

  static double guaranteedCuttingYPlane(
      final CUTTINGMODE whichMode, final double[][] coms1, final Lottery random) {
    return guaranteedCuttingPlaneHeight(whichMode, coms1, random, 1);
  }

  static double guaranteedCuttingXPlane(
      final CUTTINGMODE whichMode, final double[][] coms1, final Lottery random) {
    return guaranteedCuttingPlaneHeight(whichMode, coms1, random, 0);
  }

  private static String axisName(final int coord) {
    switch (coord) {
      case 0:
        return "x";
      case 1:
        return "y";
      case 2:
        return "z";
      default:
        return "coord" + coord;
    }
  }

  static double findOptimalZPlaneHeight(final double[][] coords, final int noUnder) {
    return findOptimalPlaneHeight(coords, noUnder, 2);
  }

  static double findOptimalYPlaneHeight(final double[][] coords, final int noUnder) {
    return findOptimalPlaneHeight(coords, noUnder, 1);
  }

  static double findOptimalXPlaneHeight(final double[][] coords, final int noUnder) {
    return findOptimalPlaneHeight(coords, noUnder, 0);
  }

  private static double findOptimalPlaneHeight(
      final double[][] coords, final int noUnder, final int coord) {

    assert (noUnder > 0);
    assert (noUnder < coords[2].length);
    assert (coord >= 0);
    assert (coord < 3);

    // order the coordinates by their coord value (so number 2)
    final double[] coordVals = coords[coord].clone();
    Arrays.sort(coordVals);

    // now we figure out where the plane would be
    if (noUnder == 0) {
      // everything is above (should never happen as it is pointless)
      return coordVals[0] - 1.0;
    }

    if (noUnder == coordVals.length) {
      // everything is underneath (should also never happen as it is pointless)
      return coordVals[coordVals.length - 1] + 1.0;
    }

    final double valFirstBelow = coordVals[noUnder - 1];
    final double valFirstAbove = coordVals[noUnder];
    final double diff = valFirstAbove - valFirstBelow;

    if (Math.abs(diff) < 1E-8) {
      // we return NaN since two COMs are too close together (mother / findOptimalPlaneHeight)
      System.err.println(
          "Cutting plane unavailable (mother/optimal): two COMs too close on axis "
              + axisName(coord)
              + " (noUnder="
              + noUnder
              + ", diff="
              + diff
              + ", valFirstBelow="
              + valFirstBelow
              + ", valFirstAbove="
              + valFirstAbove
              + ")");
      System.err.flush();
      return Double.NaN;
    }

    return (valFirstBelow + diff / 2);
  }

  static double findOptimalSphereRadius(final double[] dists, final int noInside) {

    assert (noInside > 0);
    assert (noInside < dists.length);

    // order the distances
    final double[] tmp = dists.clone();
    Arrays.sort(tmp);

    // now we figure out where the plane would be
    if (noInside == 0) {
      // everything is outside (should never happen as it is pointless)
      return tmp[0] - 1.0;
    }

    if (noInside == tmp.length) {
      // everything is underneath (should also never happen as it is pointless)
      return tmp[tmp.length] + 1.0;
    }

    final double valFirstBelow = tmp[noInside - 1];
    final double valFirstAbove = tmp[noInside];
    final double diff = valFirstAbove - valFirstBelow;

    if (Math.abs(diff) < 1E-8) {
      // we return NaN since two COMs are too close together
      log.info("Absolute difference " + Math.abs(diff) + " too small. Tried to find " + noInside);
      if (log.isDebugEnabled()) {
        log.debug(Arrays.toString(tmp));
      }
      return Double.NaN;
    }

    return (valFirstBelow + diff / 2);
  }

  static void randomRotation(final CartesianCoordinates c) {

    final double[][] rotXYZ = RandomUtils.randomRotation(c.getAllXYZCoord());
    c.setAllXYZ(rotXYZ);
  }

  /**
   * Checks if an offspring is still in the relevant properties (number of atoms, number of
   * molecules) the same as the original. Will be extended in the future and is NOT cheap to run.
   * Does NOT claim it checks EVERYTHING.
   *
   * @param original the original
   * @param offspring the offspring
   * @return true if everything it checked is fine, false otherwise
   */
  static boolean areGeometriesStillCorrect(final Geometry original, final Geometry offspring) {

    if (original.getMolecules().size() != offspring.getMolecules().size()) {
      return false;
    }
    if (original.getNumberOfIndieParticles() != offspring.getNumberOfIndieParticles()) {
      return false;
    }

    for (int i = 0; i < original.getNumberOfIndieParticles(); i++) {
      final Molecule mOr = original.getMoleculeAtPosition(i);
      final Molecule mOf = offspring.getMoleculeAtPosition(i);
      if (mOf == null) {
        return false;
      }
      if (mOr.getNumberOfAtoms() != mOf.getNumberOfAtoms()) {
        return false;
      }
      final String[] atomsOr = mOr.getAtomTypes();
      final String[] atomsOf = mOf.getAtomTypes();
      if (atomsOf.length != atomsOr.length) {
        return false;
      }

      for (int at = 0; at < atomsOr.length; at++) {
        if (!atomsOf[at].equalsIgnoreCase(atomsOr[at])) {
          return false;
        }
      }
    }

    return true;
  }
}
