/*
Copyright (c) 2016-2020, J. M. Dieterich and B. Hartke
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
import java.util.List;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.Topology;
import org.ogolem.helpers.Tuple;
import org.ogolem.helpers.Tuple3D;

/**
 * An adaptive SWG FF for, e.g., silicon.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class AdaptiveSWGFF extends AbstractAdaptiveBackend {

  private static final long serialVersionUID = (long) 20200622;

  private static final boolean DEBUG = false;
  private final List<int[]> contr13 = new ArrayList<>();
  private final List<int[]> contr14 = new ArrayList<>();

  private final AdaptiveParameters params;
  private final boolean isAdaptive;
  private final boolean useCaches;
  private final AdaptiveSWG2BodyTerm term2B;
  private final AdaptiveSWG3BodyTerm term3B;

  private ArrayList<Topology> topoCache;

  private double[][] xyz; // used as an object cache when in backend mode.

  public AdaptiveSWGFF(
      final boolean inAdaptive,
      final AdaptiveParameters parameters,
      final boolean useCaches,
      final double blowFacClose) {

    this.isAdaptive = inAdaptive;
    this.useCaches = useCaches;
    this.term2B = new AdaptiveSWG2BodyTerm(useCaches, blowFacClose);
    this.term3B = new AdaptiveSWG3BodyTerm(useCaches);
    if (inAdaptive) {
      this.params = null;
      this.topoCache = new ArrayList<>();
    } else {
      this.params = parameters;
    }
  }

  private AdaptiveSWGFF(final AdaptiveSWGFF orig) {
    this.isAdaptive = orig.isAdaptive;
    this.useCaches = orig.useCaches;
    this.term2B = orig.term2B.copy();
    this.term3B = orig.term3B.copy();
    if (orig.isAdaptive) {
      params = null;
      this.topoCache = new ArrayList<>();
      for (final Topology topo : orig.topoCache) {
        this.topoCache.add(new Topology(topo));
      }
    } else {
      this.params = orig.params.copy();
    }
  }

  @Override
  public AdaptiveSWGFF copy() {
    return new AdaptiveSWGFF(this);
  }

  @Override
  public String getMethodID() {
    return "adaptive SWG force field";
  }

  @Override
  public void gradientCalculation(
      final long lID,
      final int iIteration,
      final double[] xyz1D,
      final String[] saAtomTypes,
      final short[] atomNos,
      final int[] atsPerMol,
      final double[] energyparts,
      final int iNoOfAtoms,
      final float[] faCharges,
      final short[] spins,
      final BondInfo bonds,
      final Gradient grad,
      final boolean hasRigidEnv) {

    if (xyz == null || !useCaches) xyz = new double[3][iNoOfAtoms];
    System.arraycopy(xyz1D, 0, xyz[0], 0, iNoOfAtoms);
    System.arraycopy(xyz1D, iNoOfAtoms, xyz[1], 0, iNoOfAtoms);
    System.arraycopy(xyz1D, iNoOfAtoms * 2, xyz[2], 0, iNoOfAtoms);

    // then the topology
    final Topology topology =
        new Topology(saAtomTypes, xyz, bonds, faCharges, spins, atomNos, contr13, contr14, true);

    final Gradient g2B = term2B.partialCartesianGradient(topology, params);
    final double[][] g2BMat = g2B.getTotalGradient();

    final Gradient g3B = term3B.partialCartesianGradient(topology, params);
    final double[][] g3BMat = g3B.getTotalGradient();

    final double[][] gMat = grad.getTotalGradient();

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < gMat[0].length; j++) {
        final double gelem = g2BMat[i][j] + g3BMat[i][j];
        if (DEBUG)
          System.out.println(
              "Value at "
                  + i
                  + "/"
                  + j
                  + ": two body "
                  + g2BMat[i][j]
                  + ", three body "
                  + g3BMat[i][j]);
        assert (!Double.isNaN(gelem));
        gMat[i][j] = (Double.isInfinite(gelem)) ? FixedValues.NONCONVERGEDGRADIENT : gelem;
        if (DEBUG) System.out.println("Value at " + i + "/" + j + ": total " + gMat[i][j]);
      }
    }

    final double e = g2B.getFunctionValue() + g3B.getFunctionValue();

    assert (!Double.isNaN(e));

    final double swg = (Double.isInfinite(e)) ? FixedValues.NONCONVERGEDENERGY : e;

    grad.setTotalEnergy(swg);
    grad.normalizeGradient();
  }

  @Override
  public double energyCalculation(
      final long lID,
      final int iIteration,
      final double[] xyz1D,
      final String[] saAtomTypes,
      final short[] atomNos,
      final int[] atsPerMol,
      final double[] energyparts,
      final int iNoOfAtoms,
      final float[] faCharges,
      final short[] spins,
      final BondInfo bonds,
      final boolean hasRigidEnv) {

    if (xyz == null || !useCaches) xyz = new double[3][iNoOfAtoms];
    System.arraycopy(xyz1D, 0, xyz[0], 0, iNoOfAtoms);
    System.arraycopy(xyz1D, iNoOfAtoms, xyz[1], 0, iNoOfAtoms);
    System.arraycopy(xyz1D, iNoOfAtoms * 2, xyz[2], 0, iNoOfAtoms);

    final Topology topology =
        new Topology(saAtomTypes, xyz, bonds, faCharges, spins, atomNos, contr13, contr14, false);

    final double e2B = term2B.partialInteraction(topology, params);

    final double e3B = term3B.partialInteraction(topology, params);

    final double swg = e2B + e3B;

    assert (!Double.isNaN(swg));

    return (Double.isInfinite(swg)) ? FixedValues.NONCONVERGEDENERGY : swg;
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {

    Topology topology;
    if (topoCache.size() < geomID + 1) {
      // first the topology
      topology =
          new Topology(
              cartes.getAllAtomTypes(),
              cartes.getAllXYZCoord(),
              bonds,
              cartes.getAllCharges(),
              cartes.getAllSpins(),
              cartes.getAllAtomNumbers(),
              contr13,
              contr14,
              false);
      topoCache.add(topology);
    } else {
      topology = topoCache.get(geomID);
    }

    final double e2B = term2B.partialInteraction(topology, params);

    final double e3B = term3B.partialInteraction(topology, params);

    final double swg = e2B + e3B;

    assert (!Double.isNaN(swg));

    return (Double.isInfinite(swg)) ? FixedValues.NONCONVERGEDENERGY : swg;
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds,
      final double[] grad) {

    Topology topology;
    if (topoCache.size() < geomID + 1) {
      // first the topology
      topology =
          new Topology(
              cartes.getAllAtomTypes(),
              cartes.getAllXYZCoord(),
              bonds,
              cartes.getAllCharges(),
              cartes.getAllSpins(),
              cartes.getAllAtomNumbers(),
              contr13,
              contr14,
              false);
      topoCache.add(topology);
    } else {
      topology = topoCache.get(geomID);
    }

    final double e2B = term2B.partialParamGradient(topology, params, grad);

    final double e3B = term3B.partialParamGradient(topology, params, grad);

    final double swg = e2B + e3B;

    assert (!Double.isNaN(swg));

    return (Double.isInfinite(swg)) ? FixedValues.NONCONVERGEDENERGY : swg;
  }

  @Override
  public double[][] minMaxBordersForParams(final AdaptiveParameters params) {

    final List<double[][]> borders = new ArrayList<>();

    final double[][] borders2B = term2B.bordersForMyParams(params);
    borders.add(borders2B);

    final double[][] borders3B = term3B.bordersForMyParams(params);
    borders.add(borders3B);

    // copy around
    final double[][] allBorders = new double[2][params.getNumberOfParamters()];
    int offset = 0;
    for (final double[][] partBorders : borders) {
      final int noParams = partBorders[0].length;
      System.arraycopy(partBorders[0], 0, allBorders[0], offset, noParams);
      System.arraycopy(partBorders[1], 0, allBorders[1], offset, noParams);
      offset += noParams;
    }

    return allBorders;
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      final ArrayList<CartesianCoordinates> refCartes, final String sMethod) {

    final ArrayList<Topology> tops = new ArrayList<>();
    for (final CartesianCoordinates refCart : refCartes) {
      final BondInfo refBonds = org.ogolem.core.CoordTranslation.checkForBonds(refCart, 1.3);
      final Topology topology =
          new Topology(
              refCart.getAllAtomTypes(),
              refCart.getAllXYZCoord(),
              refBonds,
              refCart.getAllCharges(),
              refCart.getAllSpins(),
              refCart.getAllAtomNumbers(),
              contr13,
              contr14);
      tops.add(topology);
    }

    int paramSum = 0;
    ArrayList<Tuple<String, Integer>> paramsPerKey = new ArrayList<>();

    final Tuple3D<String[], int[], Integer> tupel2B =
        term2B.requiredParams(refCartes, tops, sMethod);
    final String[] keys2B = tupel2B.getObject1();
    final int[] paramsPerKey2B = tupel2B.getObject2();
    for (int i = 0; i < keys2B.length; i++) {
      final Tuple<String, Integer> tup = new Tuple<>(keys2B[i], paramsPerKey2B[i]);
      paramsPerKey.add(tup);
    }
    paramSum += tupel2B.getObject3();

    final Tuple3D<String[], int[], Integer> tupel3B =
        term3B.requiredParams(refCartes, tops, sMethod);
    final String[] keys3B = tupel3B.getObject1();
    final int[] paramsPerKey3B = tupel3B.getObject2();
    for (int i = 0; i < keys3B.length; i++) {
      final Tuple<String, Integer> tup = new Tuple<>(keys3B[i], paramsPerKey3B[i]);
      paramsPerKey.add(tup);
    }
    paramSum += tupel3B.getObject3();

    final String[] allKeys = new String[paramsPerKey.size()];
    final int[] allParamsPerKey = new int[paramsPerKey.size()];
    for (int i = 0; i < allKeys.length; i++) {
      final Tuple<String, Integer> tupel = paramsPerKey.get(i);
      allKeys[i] = tupel.getObject1();
      allParamsPerKey[i] = tupel.getObject2();
    }

    final AdaptiveParameters paramStub =
        new AdaptiveParameters(paramSum, -1, allKeys, allParamsPerKey, sMethod);

    return paramStub;
  }
}
