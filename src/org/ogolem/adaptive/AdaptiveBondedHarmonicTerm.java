/*
Copyright (c) 2010-2012, J. M. Dieterich
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

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.Topology;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.io.InputPrimitives;

/**
 * Simple harmonic potential for bonded interactions. Using specific terms makes little to no sense
 * in case that the reference cartesians are different systems! Also, if one enables caching, the
 * system is not allowed to change for different invocations of the same object!
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class AdaptiveBondedHarmonicTerm implements AdaptiveInteractionTerm {

  private static final long serialVersionUID = (long) 20101122;

  protected static final double CUTOFF = 100.0;
  protected final boolean specialIDs;
  protected final double distCutoff;
  protected final String[] ids;
  protected final boolean useCaching;
  protected int[] paramOffsetCache;

  public AdaptiveBondedHarmonicTerm(boolean useTermsForEachPair, double distCut, boolean useCache)
      throws Exception {
    this.specialIDs = useTermsForEachPair;
    this.distCutoff = distCut;
    this.useCaching = useCache;
    if (specialIDs) {
      // read special IDs in
      final String aux = "ids-harmonicbond.aux";
      try {
        final String[] sa = InputPrimitives.readFileIn(aux);
        ids = new String[sa.length];
        for (int i = 0; i < ids.length; i++) {
          final String[] sa2 = sa[i].split("\\s+");
          ids[i] = sa2[1];
        }
      } catch (Exception e) {
        throw new Exception("Couldn't use ids-harmonicbond.aux.", e);
      }
    } else {
      this.ids = null;
    }
    // initialize a null'd cache anyhow
    paramOffsetCache = null;
  }

  protected AdaptiveBondedHarmonicTerm(AdaptiveBondedHarmonicTerm orig) {
    this.distCutoff = orig.distCutoff;
    this.specialIDs = orig.specialIDs;
    this.useCaching = orig.useCaching;

    if (orig.ids != null) this.ids = orig.ids.clone();
    else this.ids = null;

    if (orig.paramOffsetCache != null) this.paramOffsetCache = orig.paramOffsetCache.clone();
    else this.paramOffsetCache = null;
  }

  @Override
  public AdaptiveBondedHarmonicTerm copy() {
    return new AdaptiveBondedHarmonicTerm(this);
  }

  @Override
  public double partialInteraction(Topology topology, AdaptiveParameters params) {

    final BondInfo bonds = topology.getBonds();
    final double[][] dists = topology.getDistances();
    final String[] atoms = topology.getAtomNames();
    final int iNoOfAtoms = dists.length;

    if (useCaching && paramOffsetCache == null) {
      int noOfBonds = 0;
      for (int i = 0; i < iNoOfAtoms - 1; i++) {
        for (int j = i + 1; j < iNoOfAtoms; j++) {
          if (bonds.hasBond(i, j)) {
            noOfBonds++;
          }
        }
      }
      paramOffsetCache = new int[noOfBonds];

      // fill cache up on the fly
      int counter = -1;
      for (int i = 0; i < iNoOfAtoms - 1; i++) {
        for (int j = i + 1; j < iNoOfAtoms; j++) {
          if (!bonds.hasBond(i, j)) continue;
          counter++;
          // fill cache up a little
          final String s = buildKey(atoms, ids, i, j, params, specialIDs);
          paramOffsetCache[counter] = params.getStartPointForKey(s);
        }
      }
    }

    double[] p = params.getAllParamters();
    double energy = 0.0;
    int counter = -1;
    int offset = 0;
    for (int i = 0; i < iNoOfAtoms - 1; i++) {
      for (int j = i + 1; j < iNoOfAtoms; j++) {
        if (!bonds.hasBond(i, j)) continue;

        counter++;

        if (!useCaching) {
          // get params
          p = getParamsForPair(atoms, ids, i, j, params, specialIDs);
          if (p == null) continue;
        } else {
          // use the cache
          offset = paramOffsetCache[counter];
        }

        // check wether we are so far apart that we need a cutoff
        if ((dists[i][j] - p[offset + 1]) >= distCutoff) {
          energy += CUTOFF;
          continue;
        }

        final double d = dists[i][j] - p[offset + 1];
        energy += 0.5 * p[offset] * d * d;
      }
    }

    // check for NaN/Infinity cases and catch them
    if (Double.isInfinite(energy) || Double.isNaN(energy)) return FixedValues.NONCONVERGEDENERGY;

    return energy;
  }

  @Override
  public double partialParamGradient(
      Topology topology, AdaptiveParameters params, final double[] daGrad) {

    final BondInfo bonds = topology.getBonds();
    final double[] daParams = params.getAllParamters();
    final double[][] dists = topology.getDistances();
    final String[] saAtoms = topology.getAtomNames();
    final int iNoOfAtoms = saAtoms.length;

    // get the start and endpoint for the parameters of *this* term
    final int[] startEnd = params.getStartAndEndForPrefix("adaptiveharmonicterm:");

    for (int param = startEnd[0]; param <= startEnd[1]; param += 2) {

      // loop over all interactions
      for (int i = 0; i < iNoOfAtoms - 1; i++) {
        for (int j = i + 1; j < iNoOfAtoms; j++) {

          if (!bonds.hasBond(i, j)) continue;

          String sPair1;
          String sPair2;
          if (specialIDs) {
            sPair1 = "adaptiveharmonicterm:" + ids[i] + ids[j];
            sPair2 = "adaptiveharmonicterm:" + ids[j] + ids[i];
          } else {
            sPair1 = "adaptiveharmonicterm:" + saAtoms[i] + saAtoms[j];
            sPair2 = "adaptiveharmonicterm:" + saAtoms[j] + saAtoms[i];
          }

          if (!sPair1.equalsIgnoreCase(params.getForWhichKey(param))
              && !sPair2.equalsIgnoreCase(params.getForWhichKey(param))) continue;

          final double d = dists[i][j] - daParams[param + 1];

          daGrad[param] += 0.5 * d * d;
          daGrad[param + 1] += -1.0 * daParams[param] * (dists[i][j] - daParams[param + 1]);
        }
      }
    }

    // XXX
    return partialInteraction(topology, params);
  }

  @Override
  public Gradient partialCartesianGradient(Topology topology, AdaptiveParameters params) {

    final int noOfAtoms = topology.getNumberOfAtoms();
    final String[] atoms = topology.getAtomNames();
    final BondInfo bonds = topology.getBonds();
    final double[][][] diffs = topology.getPosDiffs();
    final double[][] dists = topology.getDistances();

    final Gradient analyticalGradient = new Gradient();
    final double[][] grad = new double[3][noOfAtoms];

    if (useCaching && paramOffsetCache == null) {
      int noOfBonds = 0;
      for (int i = 0; i < noOfAtoms - 1; i++) {
        for (int j = i + 1; j < noOfAtoms; j++) {
          if (bonds.hasBond(i, j)) {
            noOfBonds++;
          }
        }
      }
      paramOffsetCache = new int[noOfBonds];

      // fill cache up on the fly
      int counter = -1;
      for (int i = 0; i < noOfAtoms - 1; i++) {
        for (int j = i + 1; j < noOfAtoms; j++) {
          if (!bonds.hasBond(i, j)) continue;
          counter++;
          // fill cache up a little
          final String s = buildKey(atoms, ids, i, j, params, specialIDs);
          paramOffsetCache[counter] = params.getStartPointForKey(s);
        }
      }
    }

    double[] p = params.getAllParamters();
    double energy = 0.0;
    int counter = -1;
    int offset = 0;
    for (int i = 0; i < noOfAtoms - 1; i++) {
      for (int j = i + 1; j < noOfAtoms; j++) {

        if (!bonds.hasBond(i, j)) continue;

        counter++;

        if (!useCaching) {
          // get params
          p = getParamsForPair(atoms, ids, i, j, params, specialIDs);
          if (p == null) {
            continue;
          }
        } else {
          // use the cache
          offset = paramOffsetCache[counter];
        }

        // check wether we are so far apart that we need a cutoff
        if ((dists[i][j] - p[offset + 1]) >= distCutoff) {
          energy += CUTOFF;
          continue;
        }

        double tempDeriv = p[offset] * (p[offset + 1] - dists[i][j]) / dists[i][j];

        // normalize NaNs and infinity
        if (Double.isNaN(tempDeriv) || Double.isInfinite(tempDeriv)) {
          tempDeriv = FixedValues.NONCONVERGEDGRADIENT;
        }

        grad[0][i] -= tempDeriv * (diffs[i][j][0]);
        grad[0][j] += tempDeriv * (diffs[i][j][0]);
        grad[1][i] -= tempDeriv * (diffs[i][j][1]);
        grad[1][j] += tempDeriv * (diffs[i][j][1]);
        grad[2][i] -= tempDeriv * (diffs[i][j][2]);
        grad[2][j] += tempDeriv * (diffs[i][j][2]);

        final double d = dists[i][j] - p[offset + 1];
        energy += 0.5 * p[offset] * d * d;
      }
    }

    // normalize NaNs and infinity
    if (Double.isNaN(energy) || Double.isInfinite(energy)) energy = FixedValues.NONCONVERGEDENERGY;

    analyticalGradient.setGradientTotal(grad);
    analyticalGradient.setTotalEnergy(energy);

    return analyticalGradient;
  }

  @Override
  public Tuple3D<String[], int[], Integer> requiredParams(
      ArrayList<CartesianCoordinates> cartesians, ArrayList<Topology> topologies, String sMethod) {

    if (specialIDs) {

      final Topology top = topologies.get(0);
      final BondInfo bonds = top.getBonds();
      final int noOfAtoms = top.getNumberOfAtoms();

      final LinkedList<String> nonRedPairs = new LinkedList<>();
      for (int i = 0; i < noOfAtoms - 1; i++) {
        for (int j = i + 1; j < noOfAtoms; j++) {
          if (bonds.hasBond(i, j)) {
            final String pair1 = ids[i] + ids[j];
            final String pair2 = ids[j] + ids[i];
            if (!nonRedPairs.contains(pair1) && !nonRedPairs.contains(pair2)) {
              nonRedPairs.add(pair1);
            }
          }
        }
      }

      final String[] keys = new String[nonRedPairs.size()];
      final int[] paramsPerKey = new int[nonRedPairs.size()];
      int totalSum = 0;
      int counter = 0;
      for (String pair : nonRedPairs) {
        keys[counter] = "adaptiveharmonicterm:" + pair;
        paramsPerKey[counter] = 2;
        totalSum += 2;
        counter++;
      }

      Tuple3D<String[], int[], Integer> result = new Tuple3D<>(keys, paramsPerKey, totalSum);

      return result;
    } else {
      // loop through all reference geometries and figure a list of non-redundant atom tyes out
      final LinkedList<String> llAtoms = new LinkedList<>();

      for (CartesianCoordinates cartes : cartesians) {

        final String[] saAtoms = cartes.getAllAtomTypes();
        for (final String atom : saAtoms) {
          if (!llAtoms.contains(atom)) {
            // add it to the list
            llAtoms.add(atom);
          }
        }
      }

      int iNoOfPairs = 0;
      for (int i = 0; i < llAtoms.size(); i++) {
        for (int j = 0; j <= i; j++) {
          iNoOfPairs++;
        }
      }

      // each pair gets 2 parameters
      final int iNoOfAtoms = llAtoms.size();
      int iCounter = 0;
      String[] saKeys = new String[iNoOfPairs];

      int[] iaParamsPerKey = new int[iNoOfPairs];
      int iParamSum = 0;

      for (int i = 0; i < iNoOfAtoms; i++) {
        for (int j = i; j < iNoOfAtoms; j++) {

          final String sPair = "adaptiveharmonicterm:" + llAtoms.get(i) + llAtoms.get(j);
          saKeys[iCounter] = sPair;
          iaParamsPerKey[iCounter] = 2;
          iParamSum += iaParamsPerKey[iCounter];
          iCounter++;
        }
      }

      Tuple3D<String[], int[], Integer> result = new Tuple3D<>(saKeys, iaParamsPerKey, iParamSum);

      return result;
    }
  }

  @Override
  public double[][] bordersForMyParams(AdaptiveParameters params) {

    // not all of the parameters belong to us
    final List<String> corrKeys = params.getKeysStartingWith("adaptiveharmonicterm:");

    int iNoOfKeys = corrKeys.size();
    int iNoOfParams = iNoOfKeys * 2;

    double[][] daBorders = new double[2][iNoOfParams];

    int iCounter = 0;
    for (int i = 0; i < iNoOfKeys; i++) {

      /*
       * force constant: 0.0 and 0.5 hartree is a lot
       */
      daBorders[0][iCounter] = 0.0;
      daBorders[1][iCounter] = 0.5;
      iCounter++;

      // optimal distance: between 1.0 and 8.0 a_0
      daBorders[0][iCounter] = 1.0;
      daBorders[1][iCounter] = 8.0;
      iCounter++;
    }

    return daBorders;
  }

  private static double[] getParamsForPair(
      final String[] atoms,
      final String[] ids,
      final int i,
      final int j,
      final AdaptiveParameters parameters,
      final boolean useIds) {

    double[] params;

    if (useIds) {

      final StringBuilder sb = new StringBuilder(25);
      sb.append("adaptiveharmonicterm:");
      sb.append(ids[i]);
      sb.append(ids[j]);
      params = parameters.getParametersForKey(sb.toString());

      if (params != null) return params;

      final StringBuilder sb2 = new StringBuilder(25);
      sb2.append("adaptiveharmonicterm:");
      sb2.append(ids[j]);
      sb2.append(ids[i]);
      params = parameters.getParametersForKey(sb2.toString());

      if (params != null) return params;
      else {
        System.err.println("WARNING: No parameters for key " + sb.toString() + ". Returning null.");
        return null;
      }
    } else {
      final StringBuilder sb = new StringBuilder(25);
      sb.append("adaptiveharmonicterm:");
      sb.append(atoms[i]);
      sb.append(atoms[j]);
      params = parameters.getParametersForKey(sb.toString());

      if (params != null) return params;

      final StringBuilder sb2 = new StringBuilder(25);
      sb2.append("adaptiveharmonicterm:");
      sb2.append(atoms[j]);
      sb2.append(atoms[i]);
      params = parameters.getParametersForKey(sb2.toString());

      if (params != null) return params;
      else {
        System.err.println("WARNING: No parameters for key " + sb.toString() + ". Returning null.");
        return null;
      }
    }
  }

  protected static String buildKey(
      final String[] atoms,
      final String[] ids,
      final int i,
      final int j,
      final AdaptiveParameters parameters,
      final boolean useIds) {

    if (useIds) {

      final StringBuilder sb = new StringBuilder(25);
      sb.append("adaptiveharmonicterm:");
      sb.append(ids[i]);
      sb.append(ids[j]);
      final double[] params = parameters.getParametersForKey(sb.toString());

      if (params != null) return sb.toString();

      final StringBuilder sb2 = new StringBuilder(25);
      sb2.append("adaptiveharmonicterm:");
      sb2.append(ids[j]);
      sb2.append(ids[i]);
      // needs to be this one now. in good and in bad...
      return sb2.toString();
    } else {
      final StringBuilder sb = new StringBuilder(25);
      sb.append("adaptiveharmonicterm:");
      sb.append(atoms[i]);
      sb.append(atoms[j]);
      final double[] params = parameters.getParametersForKey(sb.toString());

      if (params != null) return sb.toString();

      final StringBuilder sb2 = new StringBuilder(25);
      sb2.append("adaptiveharmonicterm:");
      sb2.append(atoms[j]);
      sb2.append(atoms[i]);
      // needs to be this one now. in good and in bad...
      return sb2.toString();
    }
  }
}
