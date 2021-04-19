/*
Copyright (c) 2011-2014, J. M. Dieterich
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
package org.ogolem.adaptive;

import static java.lang.Math.sqrt;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CoordTranslation;
import org.ogolem.core.Gradient;
import org.ogolem.core.Topology;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.io.InputPrimitives;
import org.ogolem.math.TrivialLinearAlgebra;

/**
 * An implementation of AMBERs dihedral term. Please note that to-date it is fixed to two terms per
 * dihedral.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public final class AdaptiveAmberDihedralTerm implements AdaptiveInteractionTerm {
  // XXX do not always compute both sin/cos() depending upon param[0] and param[3]
  private static final long serialVersionUID = (long) 20111021;
  private static final boolean DEBUG = false;
  private static final double LARGE = 1E5;
  private static final double EDGECASE = 0.0001;
  private final boolean useSpecialIDs;
  private final String[] ids;
  private final AdaptiveAmberFF.AmberMath math;
  private final boolean useCaching;
  private int[] paramOffsetCache;

  public AdaptiveAmberDihedralTerm(
      boolean specialIDs, boolean useCache, AdaptiveAmberFF.AmberMath ambermath) throws Exception {
    this.useSpecialIDs = specialIDs;
    this.useCaching = useCache;
    this.math = ambermath;

    if (specialIDs) {
      // read special IDs in
      final String aux = "ids-amberdihedral.aux";
      try {
        final String[] sa = InputPrimitives.readFileIn(aux);
        ids = new String[sa.length];
        for (int i = 0; i < ids.length; i++) {
          final String[] sa2 = sa[i].split("\\s+");
          ids[i] = sa2[1];
        }
      } catch (Exception e) {
        throw new Exception("Couldn't use ids-amberdihedral.aux.", e);
      }
    } else {
      this.ids = null;
    }

    // initialize cache anyways
    this.paramOffsetCache = null;
  }

  private AdaptiveAmberDihedralTerm(AdaptiveAmberDihedralTerm orig) {
    this.useSpecialIDs = orig.useSpecialIDs;
    this.useCaching = orig.useCaching;
    this.math = orig.math.copy();

    if (orig.ids != null) this.ids = orig.ids.clone();
    else this.ids = null;

    if (orig.paramOffsetCache != null) this.paramOffsetCache = orig.paramOffsetCache.clone();
    else this.paramOffsetCache = null;
  }

  @Override
  public AdaptiveAmberDihedralTerm copy() {
    return new AdaptiveAmberDihedralTerm(this);
  }

  @Override
  public double partialInteraction(final Topology topology, final AdaptiveParameters params) {

    final String[] atoms = topology.getAtomNames();
    final double[][] pos = topology.getPositions();
    final int[][] inter14 = topology.get14ContributionsField();

    if (useCaching && paramOffsetCache == null) {
      paramOffsetCache = new int[inter14.length];
      int i = 0;
      for (int[] inter : inter14) {
        // fill cache up a little
        final String s = buildKey(params, atoms, ids, inter, useSpecialIDs);
        paramOffsetCache[i] = params.getStartPointForKey(s);
        i++;
      }
    }

    // a paramobject, per default with all params
    double[] p = params.getAllParamters();
    int offset = 0;
    double energy = 0.0;
    for (int counter = 0; counter < inter14.length; counter++) {

      final int[] inter = inter14[counter];

      if (!useCaching) {
        // get params
        p = getParams(params, atoms, ids, inter, useSpecialIDs);
        if (p == null) {
          System.err.println(
              "WARNING: No parameters for Amber dihedral "
                  + atoms[inter[0]]
                  + inter[0]
                  + atoms[inter[1]]
                  + inter[1]
                  + " "
                  + atoms[inter[2]]
                  + inter[2]
                  + " "
                  + atoms[inter[3]]
                  + inter[3]);
          continue;
        }
      } else offset = paramOffsetCache[counter]; // use the cache

      // compute angle
      final double angle =
          CoordTranslation.calcDihedral(pos, inter[0], inter[1], inter[2], inter[3]);
      if (angle != angle) {
        System.err.println(
            "WARNING: NaN angle in Amber dihedral for "
                + inter[0]
                + " "
                + inter[1]
                + " "
                + inter[2]
                + " "
                + inter[3]
                + ". Continuing.");
        continue;
      }

      // compute energy
      energy +=
          0.5
              * (p[offset] * (1 + math.cos(p[offset + 1] * angle - p[offset + 2]))
                  + p[offset + 3] * (1 + math.cos(p[offset + 4] * angle - p[offset + 5])));
    }

    return energy;
  }

  @Override
  public double partialParamGradient(
      Topology topology, AdaptiveParameters params, final double[] daGrad) {

    final int[] startEnd = params.getStartAndEndForPrefix("amberdihedral:");
    if (startEnd == null) return FixedValues.NONCONVERGEDENERGY;

    final double[][] pos = topology.getPositions();
    final String[] atoms = topology.getAtomNames();
    final int[][] inters = topology.get14ContributionsField();

    if (useCaching && paramOffsetCache == null) {
      // initialize the cache
      paramOffsetCache = new int[inters.length];
      int i = 0;
      for (int[] inter : inters) {
        // fill cache up a little
        final String s = buildKey(params, atoms, ids, inter, useSpecialIDs);
        paramOffsetCache[i] = params.getStartPointForKey(s);
        i++;
      }
    }

    final double[] all = params.getAllParamters();

    // loop over all interactions
    double energy = 0.0;
    for (int counter = 0; counter < inters.length; counter++) {

      final int[] inter = inters[counter];

      // what is the offset for the parameters
      int param;
      if (!useCaching) {

        final String s = buildKey(params, atoms, ids, inter, useSpecialIDs);
        param = params.getStartPointForKey(s);

      } else param = paramOffsetCache[counter];

      // compute angle
      final double angle =
          CoordTranslation.calcDihedral(pos, inter[0], inter[1], inter[2], inter[3]);
      if (angle != angle) {
        System.err.println(
            "WARNING: NaN angle in Amber dihedral for "
                + inter[0]
                + " "
                + inter[1]
                + " "
                + inter[2]
                + " "
                + inter[3]
                + ". Continuing.");
        continue;
      }

      daGrad[param] += 0.5 * (1 + math.cos(all[param + 1] * angle - all[param + 2]));
      daGrad[param + 1] +=
          -0.5 * angle * all[param] * math.sin(all[param + 1] * angle - all[param + 2]);
      daGrad[param + 2] += 0.5 * all[param] * math.sin(all[param + 1] * angle - all[param + 2]);
      daGrad[param + 3] += 0.5 * (1 + math.cos(all[param + 4] * angle - all[param + 5]));
      daGrad[param + 4] +=
          -0.5 * angle * all[param + 3] * math.sin(all[param + 4] * angle - all[param + 5]);
      daGrad[param + 5] += 0.5 * all[param + 3] * math.sin(all[param + 4] * angle - all[param + 5]);
      energy +=
          0.5
              * (all[param] * (1 + math.cos(all[param + 1] * angle - all[param + 2]))
                  + all[param + 3] * (1 + math.cos(all[param + 4] * angle - all[param + 5])));
    }

    return energy;
  }

  @Override
  public Gradient partialCartesianGradient(
      final Topology topology, final AdaptiveParameters params) {

    final String[] atoms = topology.getAtomNames();
    final double[][][] diffs = topology.getPosDiffs();
    final int[][] inter14 = topology.get14ContributionsField();

    final Gradient gradient = new Gradient();
    final double[][] grad = new double[3][atoms.length];
    final double[] c1 = new double[3];
    final double[] c2 = new double[3];
    final double[] c3 = new double[3];

    if (useCaching && paramOffsetCache == null) {
      paramOffsetCache = new int[inter14.length];
      int i = 0;
      for (int[] inter : inter14) {
        // fill cache up a little
        final String s = buildKey(params, atoms, ids, inter, useSpecialIDs);
        paramOffsetCache[i] = params.getStartPointForKey(s);
        i++;
      }
    }

    // a paramobject, per default with all params
    double[] p = params.getAllParamters();
    int offset = 0;
    double energy = 0.0;
    for (int counter = 0; counter < inter14.length; counter++) {

      if (DEBUG) {
        System.out.println("dihedral " + counter);
      }

      final int i = inter14[counter][0];
      final int j = inter14[counter][1];
      final int k = inter14[counter][2];
      final int l = inter14[counter][3];

      if (DEBUG) {
        System.out.println("atoms: " + i + "\t" + j + "\t" + k + "\t" + l);
      }

      if (!useCaching) {
        // get params
        p = getParams(params, atoms, ids, inter14[counter], useSpecialIDs);
        if (p == null) {
          System.err.println(
              "WARNING: No parameters for Amber dihedral "
                  + atoms[i]
                  + i
                  + " "
                  + atoms[j]
                  + j
                  + " "
                  + atoms[k]
                  + k
                  + " "
                  + atoms[l]
                  + l);
          continue;
        }
      } else {
        // use the cache
        offset = paramOffsetCache[counter];
      }

      // direction vectors
      final double[] diffIJ = diffs[i][j];
      final double[] diffIK = diffs[i][k];
      final double[] diffKL = diffs[k][l];
      final double[] diffJK = diffs[j][k];
      final double[] diffLJ = diffs[l][j];

      // compute cross products (then stored in c1-c3)
      TrivialLinearAlgebra.crossProduct(diffIJ, diffJK, c1);
      TrivialLinearAlgebra.crossProduct(diffJK, diffKL, c2);
      TrivialLinearAlgebra.crossProduct(c1, c2, c3);

      final double lengthC1 = sqrt(c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2]);
      final double lengthC2 = sqrt(c2[0] * c2[0] + c2[1] * c2[1] + c2[2] * c2[2]);
      final double lengthC3 = sqrt(c3[0] * c3[0] + c3[1] * c3[1] + c3[2] * c3[2]);
      double val21 = lengthC2 / lengthC1;
      double val12 = lengthC1 / lengthC2;

      // dot product
      double dot = 0.0;
      double dot2 = 0.0;
      boolean isEdgeCase = false;
      double dihedral;
      if (lengthC3 > EDGECASE) {
        for (int z = 0; z < 3; z++) {
          dot += c1[z] * c2[z];
          dot2 += diffJK[z] * c3[z];
        }
        dot /= (lengthC1 * lengthC2);
        dihedral = math.acos(dot);
      } else {
        // take care of numerical inaccuracies
        // XXX might be wrong but should not matter in the context of the symmetric dihedral
        isEdgeCase = true;
        if (DEBUG) {
          System.err.println("WARNING: AMBER FF hits edge case in dihedral.");
        }
        dot = 1.0;
        if (DEBUG) {
          System.out.println("before setting: " + val21 + " " + val12);
        }
        dihedral = 0.0;
      }

      // GRADIENT AFTER THE DIHEDRAL (check minus/plus)
      final double front;
      if (dot2 < 0) { // check for minus sign using dot2
        // plus
        front =
            0.5
                * (p[offset] * p[offset + 1] * math.sin(p[offset + 2] - p[offset + 1] * dihedral)
                    + p[offset + 3]
                        * p[offset + 4]
                        * math.sin(p[offset + 5] - p[offset + 4] * dihedral));
      } else {
        // minus
        front =
            -0.5
                * (p[offset] * p[offset + 1] * math.sin(p[offset + 1] * dihedral + p[offset + 2])
                    + p[offset + 3]
                        * p[offset + 4]
                        * math.sin(p[offset + 4] * dihedral + p[offset + 5]));
        dihedral = -dihedral;
      }

      // GRADIENT AFTER THE CARTESIAN COMPONENTS
      // first the common prefactor
      final double prefac = isEdgeCase ? LARGE : 1 / (sqrt(1 - dot * dot)) / (lengthC1 * lengthC2);
      final double comb = front * prefac;

      // i and l atoms are easier
      final double aX = -diffJK[1] * c2[2] + diffJK[2] * c2[1];
      final double bX = -diffJK[1] * c1[2] + diffJK[2] * c1[1];
      grad[0][i] += comb * (aX - dot * val21 * bX);
      grad[0][l] += comb * (bX - dot * val12 * aX);

      final double aY = -diffJK[2] * c2[0] + diffJK[0] * c2[2];
      final double bY = -diffJK[2] * c1[0] + diffJK[0] * c1[2];
      grad[1][i] += comb * (aY - dot * val21 * bY);
      grad[1][l] += comb * (bY - dot * val12 * aY);

      final double aZ = -diffJK[0] * c2[1] + diffJK[1] * c2[0];
      final double bZ = -diffJK[0] * c1[1] + diffJK[1] * c1[0];
      grad[2][i] += comb * (aZ - dot * val21 * bZ);
      grad[2][l] += comb * (bZ - dot * val12 * aZ);

      // j and k are more difficult
      grad[0][j] +=
          comb
              * (diffKL[2] * c1[1]
                  - diffKL[1] * c1[2]
                  - diffIK[2] * c2[1]
                  + diffIK[1] * c2[2]
                  + dot
                      * (val21 * (diffIK[2] * c1[1] - diffIK[1] * c1[2])
                          + val12 * (-diffKL[2] * c2[1] + diffKL[1] * c2[2])));
      grad[1][j] +=
          comb
              * (diffKL[0] * c1[2]
                  - diffKL[2] * c1[0]
                  - diffIK[0] * c2[2]
                  + diffIK[2] * c2[0]
                  + dot
                      * (val21 * (diffIK[0] * c1[2] - diffIK[2] * c1[0])
                          + val12 * (-diffKL[0] * c2[2] + diffKL[2] * c2[0])));
      grad[2][j] +=
          comb
              * (diffKL[1] * c1[0]
                  - diffKL[0] * c1[1]
                  - diffIK[1] * c2[0]
                  + diffIK[0] * c2[1]
                  + dot
                      * (val21 * (diffIK[1] * c1[0] - diffIK[0] * c1[1])
                          + val12 * (-diffKL[1] * c2[0] + diffKL[0] * c2[1])));
      grad[0][k] +=
          comb
              * (diffLJ[2] * c1[1]
                  - diffLJ[1] * c1[2]
                  + diffIJ[2] * c2[1]
                  - diffIJ[1] * c2[2]
                  + dot
                      * (val21 * (-diffIJ[2] * c1[1] + diffIJ[1] * c1[2])
                          + val12 * (-diffLJ[2] * c2[1] + diffLJ[1] * c2[2])));
      grad[1][k] +=
          comb
              * (diffLJ[0] * c1[2]
                  - diffLJ[2] * c1[0]
                  + diffIJ[0] * c2[2]
                  - diffIJ[2] * c2[0]
                  + dot
                      * (val21 * (-diffIJ[0] * c1[2] + diffIJ[2] * c1[0])
                          + val12 * (-diffLJ[0] * c2[2] + diffLJ[2] * c2[0])));
      grad[2][k] +=
          comb
              * (diffLJ[1] * c1[0]
                  - diffLJ[0] * c1[1]
                  + diffIJ[1] * c2[0]
                  - diffIJ[0] * c2[1]
                  + dot
                      * (val21 * (-diffIJ[1] * c1[0] + diffIJ[0] * c1[1])
                          + val12 * (-diffLJ[1] * c2[0] + diffLJ[0] * c2[1])));

      energy +=
          0.5
              * (p[offset] * (1 + math.cos(p[offset + 1] * dihedral - p[offset + 2]))
                  + p[offset + 3] * (1 + math.cos(p[offset + 4] * dihedral - p[offset + 5])));
    }

    gradient.setGradientTotal(grad);
    gradient.setTotalEnergy(energy);

    return gradient;
  }

  @Override
  public Tuple3D<String[], int[], Integer> requiredParams(
      ArrayList<CartesianCoordinates> cartesians, ArrayList<Topology> topologies, String sMethod) {

    final Topology ref = topologies.get(0);
    final int[][] dihedralContr = ref.get14ContributionsField();
    final String[] atoms = ref.getAtomNames();

    final LinkedList<String> nonRed = new LinkedList<>();
    for (int[] contr : dihedralContr) {
      String s1;
      String s2;
      if (useSpecialIDs) {
        s1 = "amberdihedral:" + ids[contr[0]] + ids[contr[1]] + ids[contr[2]] + ids[contr[3]];
        s2 = "amberdihedral:" + ids[contr[3]] + ids[contr[2]] + ids[contr[1]] + ids[contr[0]];
      } else {
        s1 =
            "amberdihedral:"
                + atoms[contr[0]]
                + atoms[contr[1]]
                + atoms[contr[2]]
                + atoms[contr[3]];
        s2 =
            "amberdihedral:"
                + atoms[contr[3]]
                + atoms[contr[2]]
                + atoms[contr[1]]
                + atoms[contr[0]];
      }
      if (!nonRed.contains(s1) && !nonRed.contains(s2)) {
        nonRed.add(s1);
      }
    }

    final String[] keys = new String[nonRed.size()];
    final int[] paramsPerKey = new int[nonRed.size()];
    int totalSum = 0;
    int counter = 0;
    for (String key : nonRed) {
      keys[counter] = key;
      paramsPerKey[counter] = 6;
      totalSum += 6;
      counter++;
    }

    Tuple3D<String[], int[], Integer> result = new Tuple3D<>(keys, paramsPerKey, totalSum);

    return result;
  }

  @Override
  public double[][] bordersForMyParams(AdaptiveParameters params) {

    // not all of the parameters belong to us
    final List<String> corrKeys = params.getKeysStartingWith("amberdihedral:");

    final int noOfKeys = corrKeys.size();
    final int noOfParams = noOfKeys * 6;
    final double[][] borders = new double[2][noOfParams];

    int counter = 0;
    for (int i = 0; i < noOfKeys; i++) {
      // first is a force constant +- 0.1Eh
      borders[0][counter] = -0.1;
      borders[1][counter] = +0.1;
      counter++;

      // second is a modulation
      borders[0][counter] = 0.0;
      borders[1][counter] = 5.0;
      counter++;

      // second is the dihedral -pi -> pi
      borders[0][counter] = -Math.PI;
      borders[1][counter] = Math.PI;
      counter++;

      // first is a force constant +- 0.1Eh
      borders[0][counter] = -0.1;
      borders[1][counter] = +0.1;
      counter++;

      // second is a modulation
      borders[0][counter] = 0.0;
      borders[1][counter] = 5.0;
      counter++;

      // second is the dihedral -pi -> pi
      borders[0][counter] = -Math.PI;
      borders[1][counter] = Math.PI;
      counter++;
    }

    return borders;
  }

  private static double[] getParams(
      final AdaptiveParameters params,
      final String[] atoms,
      final String[] ids,
      final int[] which,
      final boolean useIds) {

    String s;
    if (useIds) {
      s = "amberdihedral:" + ids[which[0]] + ids[which[1]] + ids[which[2]] + ids[which[3]];
      double[] p = params.getParametersForKey(s);
      if (p != null) return p;

      s = "amberdihedral:" + ids[which[3]] + ids[which[2]] + ids[which[1]] + ids[which[0]];
      // needs to be this one now, in good and in bad...
      return params.getParametersForKey(s);
    } else {
      s = "amberdihedral:" + atoms[which[0]] + atoms[which[1]] + atoms[which[2]] + atoms[which[3]];
      double[] p = params.getParametersForKey(s);
      if (p != null) return p;

      s = "amberdihedral:" + atoms[which[3]] + atoms[which[2]] + atoms[which[1]] + atoms[which[0]];
      // needs to be this one now, in good and in bad...
      return params.getParametersForKey(s);
    }
  }

  private static String buildKey(
      final AdaptiveParameters params,
      final String[] atoms,
      final String[] ids,
      final int[] which,
      final boolean useIds) {

    final StringBuilder sb1 = new StringBuilder(25);
    final StringBuilder sb2 = new StringBuilder(25);
    if (useIds) {
      sb1.append("amberdihedral:");
      sb1.append(ids[which[0]]);
      sb1.append(ids[which[1]]);
      sb1.append(ids[which[2]]);
      sb1.append(ids[which[3]]);
      final double[] p = params.getParametersForKey(sb1.toString());

      if (p != null) return sb1.toString();

      sb2.append("amberdihedral:");
      sb2.append(ids[which[3]]);
      sb2.append(ids[which[2]]);
      sb2.append(ids[which[1]]);
      sb2.append(ids[which[0]]);
      // needs to be this one now. in good and in bad...
      return sb2.toString();
    } else {
      sb1.append("amberdihedral:");
      sb1.append(atoms[which[0]]);
      sb1.append(atoms[which[1]]);
      sb1.append(atoms[which[2]]);
      sb1.append(atoms[which[3]]);
      final double[] p = params.getParametersForKey(sb1.toString());

      if (p != null) return sb1.toString();

      sb2.append("amberdihedral:");
      sb2.append(atoms[which[3]]);
      sb2.append(atoms[which[2]]);
      sb2.append(atoms[which[1]]);
      sb2.append(atoms[which[0]]);
      // needs to be this one now. in good and in bad...
      return sb2.toString();
    }
  }
}
