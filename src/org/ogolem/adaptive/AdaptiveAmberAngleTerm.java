/*
Copyright (c) 2011-2012, J. M. Dieterich
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

/**
 * Implementation of a classical AMBER-style angle term. If caching is used, it is assumed that
 * *all* systems calculated with this object are the same! Please note this!
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public final class AdaptiveAmberAngleTerm implements AdaptiveInteractionTerm {

  private static final long serialVersionUID = (long) 20110113;
  private final boolean useSpecialIDs;
  private final String[] ids;
  private final boolean useCaching;
  private int[] paramOffsetCache;

  private final AdaptiveAmberFF.AmberMath math;

  public AdaptiveAmberAngleTerm(
      boolean specialIDs, boolean useCache, AdaptiveAmberFF.AmberMath ambermath) throws Exception {
    this.useSpecialIDs = specialIDs;
    this.math = ambermath;

    if (specialIDs) {
      // read special IDs in
      final String aux = "ids-amberangle.aux";
      try {
        final String[] sa = InputPrimitives.readFileIn(aux);
        ids = new String[sa.length];
        for (int i = 0; i < ids.length; i++) {
          final String[] sa2 = sa[i].split("\\s+");
          ids[i] = sa2[1];
        }
      } catch (Exception e) {
        throw new Exception("Couldn't use ids-amberangle.aux.", e);
      }
    } else {
      this.ids = null;
    }

    this.useCaching = useCache;
    // initialize a null'd cache anyhow
    paramOffsetCache = null;
  }

  private AdaptiveAmberAngleTerm(AdaptiveAmberAngleTerm orig) {
    this.useSpecialIDs = orig.useSpecialIDs;
    this.useCaching = orig.useCaching;
    this.math = orig.math.copy();

    if (orig.ids != null) this.ids = orig.ids.clone();
    else this.ids = null;

    if (orig.paramOffsetCache != null) this.paramOffsetCache = orig.paramOffsetCache.clone();
    else this.paramOffsetCache = null;
  }

  @Override
  public AdaptiveAmberAngleTerm copy() {
    return new AdaptiveAmberAngleTerm(this);
  }

  @Override
  public double partialInteraction(Topology topology, AdaptiveParameters params) {

    final String[] atoms = topology.getAtomNames();
    final double[][] pos = topology.getPositions();
    final int[][] inter13 = topology.get13ContributionsField();

    if (useCaching && paramOffsetCache == null) {
      paramOffsetCache = new int[inter13.length];
      int i = 0;
      for (int[] inter : inter13) {
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
    for (int counter = 0; counter < inter13.length; counter++) {

      final int[] inter = inter13[counter];

      if (!useCaching) {
        // get params
        p = getParams(params, atoms, ids, inter, useSpecialIDs);
        if (p == null) {
          System.err.println(
              "WARNING: No parameters for Amber angle "
                  + atoms[inter[0]]
                  + inter[0]
                  + atoms[inter[1]]
                  + inter[1]
                  + " "
                  + atoms[inter[2]]
                  + inter[2]);
          continue;
        }
      } else {
        // use the cache
        offset = paramOffsetCache[counter];
      }

      // compute angle
      final double angle = CoordTranslation.calcAngle(pos, inter[0], inter[1], inter[2]);
      if (angle != angle) {
        System.err.println(
            "WARNING: NaN angle in Amber angle for "
                + inter[0]
                + " "
                + inter[1]
                + " "
                + inter[2]
                + ". Continuing.");
        continue;
      }

      // compute energy
      final double d = angle - p[offset + 1];
      energy += 0.5 * p[offset] * d * d;
    }

    return energy;
  }

  @Override
  public double partialParamGradient(
      Topology topology, AdaptiveParameters params, final double[] daGrad) {

    // get the start and endpoint for the parameters of *this* term
    final int[] startEnd = params.getStartAndEndForPrefix("amberangle:");
    if (startEnd == null) return FixedValues.NONCONVERGEDENERGY;

    final double[] daParams = params.getAllParamters();
    final double[][] pos = topology.getPositions();
    final String[] saAtoms = topology.getAtomNames();
    final int[][] inters = topology.get13ContributionsField();

    if (useCaching && paramOffsetCache == null) {
      // initialize the cache
      paramOffsetCache = new int[inters.length];
      int i = 0;
      for (int[] inter : inters) {
        // fill cache up a little
        final String s = buildKey(params, saAtoms, ids, inter, useSpecialIDs);
        paramOffsetCache[i] = params.getStartPointForKey(s);
        i++;
      }
    }

    // loop over all interactions
    double energy = 0.0;
    for (int counter = 0; counter < inters.length; counter++) {

      // what is the offset for the parameters
      int param;
      if (!useCaching) {

        final String s = buildKey(params, saAtoms, ids, inters[counter], useSpecialIDs);
        param = params.getStartPointForKey(s);

      } else {
        param = paramOffsetCache[counter];
      }

      // compute angle
      final double angle =
          CoordTranslation.calcAngle(
              pos, inters[counter][0], inters[counter][1], inters[counter][2]);
      if (angle != angle) {
        System.err.println(
            "WARNING: NaN angle in Amber angle for "
                + inters[counter][0]
                + " "
                + inters[counter][1]
                + " "
                + inters[counter][2]
                + ". Continuing.");
        continue;
      }

      final double d = angle - daParams[param + 1];
      daGrad[param] += 0.5 * d * d;
      daGrad[param + 1] += -1.0 * daParams[param] * d;
      energy += 0.5 * daParams[param] * d * d;
    }

    return energy;
  }

  @Override
  public Gradient partialCartesianGradient(Topology topology, AdaptiveParameters params) {

    final String[] atoms = topology.getAtomNames();
    final double[][][] diffs = topology.getPosDiffs();
    final int[][] inter13 = topology.get13ContributionsField();

    if (useCaching && paramOffsetCache == null) {
      paramOffsetCache = new int[inter13.length];
      int i = 0;
      for (int[] inter : inter13) {
        // fill cache up a little
        final String s = buildKey(params, atoms, ids, inter, useSpecialIDs);
        paramOffsetCache[i] = params.getStartPointForKey(s);
        i++;
      }
    }

    final Gradient gradient = new Gradient();
    final double[][] grad = new double[3][atoms.length];

    double[] p = params.getAllParamters();
    int offset = 0;
    double energy = 0.0;
    for (int counter = 0; counter < inter13.length; counter++) {

      final int i = inter13[counter][0];
      final int j = inter13[counter][1];
      final int k = inter13[counter][2];

      if (!useCaching) {
        // get params
        p = getParams(params, atoms, ids, inter13[counter], useSpecialIDs);
        if (p == null) {
          System.err.println(
              "WARNING: No parameters for Amber angle "
                  + atoms[inter13[counter][0]]
                  + inter13[counter][0]
                  + atoms[inter13[counter][1]]
                  + inter13[counter][1]
                  + " "
                  + atoms[inter13[counter][2]]
                  + inter13[counter][2]);
          continue;
        }
      } else {
        // use the cache
        offset = paramOffsetCache[counter];
      }

      // compute angle

      // vectors between
      final double[] v1 = diffs[inter13[counter][0]][inter13[counter][1]];
      final double[] v2 = diffs[inter13[counter][2]][inter13[counter][1]];

      // the squared norms
      final double v1nSq = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
      final double v2nSq = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];

      // the norms
      final double v1norm = Math.sqrt(v1nSq);
      final double v2norm = Math.sqrt(v2nSq);

      // the scalar product of the non-normalized vectors
      final double v1v2 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
      final double phi1 = -v1v2;

      final double angle = math.acos(v1v2 / (v1norm * v2norm));
      if (angle != angle) {
        System.err.println(
            "WARNING: NaN angle in Amber angle for " + i + " " + j + " " + k + ". Continuing.");
        continue;
      }

      final double angleTerm = p[offset + 1] - angle;
      final double denom = sqrt(1 - phi1 * phi1 / (v1nSq * v2nSq));

      // gradient components
      final double v1nv2n = v1norm * v2norm;
      final double v1nv2nInv = 1.0 / v1nv2n;
      final double phi1m1 = phi1 / (v1nSq * v1nv2n);
      final double phi1m2 = phi1 / (v2nSq * v1nv2n);
      final double t1_0 = v1[0] * phi1m1;
      final double t1_1 = v1[1] * phi1m1;
      final double t1_2 = v1[2] * phi1m1;
      final double t2_0 = v2[0] * phi1m2;
      final double t2_1 = v2[1] * phi1m2;
      final double t2_2 = v2[2] * phi1m2;
      final double mpang = -p[offset] * angleTerm / denom;

      // we do not ptimize the v1[X]*v1nv2nInv and v2[X]*v1nv2nInv away
      // too little gain, too obfuscated code
      grad[0][i] += mpang * (-v2[0] * v1nv2nInv - t1_0);
      grad[1][i] += mpang * (-v2[1] * v1nv2nInv - t1_1);
      grad[2][i] += mpang * (-v2[2] * v1nv2nInv - t1_2);

      grad[0][j] += mpang * ((v1[0] + v2[0]) * v1nv2nInv + t1_0 + t2_0);
      grad[1][j] += mpang * ((v1[1] + v2[1]) * v1nv2nInv + t1_1 + t2_1);
      grad[2][j] += mpang * ((v1[2] + v2[2]) * v1nv2nInv + t1_2 + t2_2);

      grad[0][k] -= mpang * (v1[0] * v1nv2nInv + t2_0);
      grad[1][k] -= mpang * (v1[1] * v1nv2nInv + t2_1);
      grad[2][k] -= mpang * (v1[2] * v1nv2nInv + t2_2);

      // energy (that's easy...)
      energy += 0.5 * p[offset] * angleTerm * angleTerm;
    }

    gradient.setGradientTotal(grad);
    gradient.setTotalEnergy(energy);

    return gradient;
  }

  @Override
  public Tuple3D<String[], int[], Integer> requiredParams(
      ArrayList<CartesianCoordinates> cartesians, ArrayList<Topology> topologies, String sMethod) {

    final Topology ref = topologies.get(0);
    final int[][] angleContr = ref.get13ContributionsField();
    final String[] atoms = ref.getAtomNames();

    final LinkedList<String> nonRed = new LinkedList<>();
    for (int[] contr : angleContr) {
      String s1;
      String s2;
      if (useSpecialIDs) {
        s1 = "amberangle:" + ids[contr[0]] + ids[contr[1]] + ids[contr[2]];
        s2 = "amberangle:" + ids[contr[2]] + ids[contr[1]] + ids[contr[0]];
      } else {
        s1 = "amberangle:" + atoms[contr[0]] + atoms[contr[1]] + atoms[contr[2]];
        s2 = "amberangle:" + atoms[contr[2]] + atoms[contr[1]] + atoms[contr[0]];
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
      paramsPerKey[counter] = 2;
      totalSum += 2;
      counter++;
    }

    Tuple3D<String[], int[], Integer> result = new Tuple3D<>(keys, paramsPerKey, totalSum);

    return result;
  }

  @Override
  public double[][] bordersForMyParams(AdaptiveParameters params) {

    // not all of the parameters belong to us
    final List<String> corrKeys = params.getKeysStartingWith("amberangle:");

    final int noOfKeys = corrKeys.size();
    final int noOfParams = noOfKeys * 2;
    final double[][] borders = new double[2][noOfParams];

    int counter = 0;
    for (int i = 0; i < noOfKeys; i++) {
      // first is a force constant  up to 0.1Eh
      borders[0][counter] = 0.0;
      borders[1][counter] = 0.1;
      counter++;

      // second is an angle 0 -> pi
      borders[0][counter] = 0.0;
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

    final StringBuilder sb1 = new StringBuilder(20);
    final StringBuilder sb2 = new StringBuilder(20);
    if (useIds) {
      sb1.append("amberangle:");
      sb1.append(ids[which[0]]);
      sb1.append(ids[which[1]]);
      sb1.append(ids[which[2]]);
      final double[] p = params.getParametersForKey(sb1.toString());
      if (p != null) return p;

      sb2.append("amberangle:");
      sb2.append(ids[which[2]]);
      sb2.append(ids[which[1]]);
      sb2.append(ids[which[0]]);
      return params.getParametersForKey(sb2.toString());
      // needs to be this one, in good and in bad
    } else {
      sb1.append("amberangle:");
      sb1.append(atoms[which[0]]);
      sb1.append(atoms[which[1]]);
      sb1.append(atoms[which[2]]);
      final double[] p = params.getParametersForKey(sb1.toString());
      if (p != null) return p;

      sb2.append("amberangle:");
      sb2.append(atoms[which[2]]);
      sb2.append(atoms[which[1]]);
      sb2.append(atoms[which[0]]);
      return params.getParametersForKey(sb2.toString());
      // needs to be this one, in good and in bad
    }
  }

  private static String buildKey(
      final AdaptiveParameters params,
      final String[] atoms,
      final String[] ids,
      final int[] which,
      final boolean useIds) {

    final StringBuilder sb1 = new StringBuilder(20);
    final StringBuilder sb2 = new StringBuilder(20);
    if (useIds) {
      sb1.append("amberangle:");
      sb1.append(ids[which[0]]);
      sb1.append(ids[which[1]]);
      sb1.append(ids[which[2]]);
      final double[] p = params.getParametersForKey(sb1.toString());

      if (p != null) return sb1.toString();

      sb2.append("amberangle:");
      sb2.append(ids[which[2]]);
      sb2.append(ids[which[1]]);
      sb2.append(ids[which[0]]);
      // needs to be this one now. in good and in bad...
      return sb2.toString();
    } else {
      sb1.append("amberangle:");
      sb1.append(atoms[which[0]]);
      sb1.append(atoms[which[1]]);
      sb1.append(atoms[which[2]]);
      final double[] p = params.getParametersForKey(sb1.toString());

      if (p != null) return sb1.toString();

      sb2.append("amberangle:");
      sb2.append(atoms[which[2]]);
      sb2.append(atoms[which[1]]);
      sb2.append(atoms[which[0]]);
      // needs to be this one now. in good and in bad...
      return sb2.toString();
    }
  }
}
