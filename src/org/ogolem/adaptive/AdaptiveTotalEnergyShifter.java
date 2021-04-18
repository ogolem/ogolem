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

import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.Topology;
import org.ogolem.helpers.Tuple3D;

/**
 * Shifts total energies. Might be helpful sometimes.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class AdaptiveTotalEnergyShifter implements AdaptiveInteractionTerm {

  private static final long serialVersionUID = (long) 20111024;
  private final boolean useCaching;
  private int[] posCache;

  AdaptiveTotalEnergyShifter(final boolean useCache) {
    this.useCaching = useCache;
  }

  private AdaptiveTotalEnergyShifter(final AdaptiveTotalEnergyShifter orig) {
    this.useCaching = orig.useCaching;
    if (orig.posCache != null) this.posCache = orig.posCache.clone();
    else this.posCache = null;
  }

  @Override
  public AdaptiveTotalEnergyShifter copy() {
    return new AdaptiveTotalEnergyShifter(this);
  }

  @Override
  public double partialInteraction(final Topology topology, final AdaptiveParameters params) {

    final String[] atoms = topology.getAtomNames();
    final int noOfAtoms = topology.getNumberOfAtoms();

    if (useCaching && posCache == null) {
      // initialize cache
      posCache = new int[noOfAtoms];
      for (int i = 0; i < noOfAtoms; i++) {
        final String s = "totalenergyshifter:" + atoms[i];
        posCache[i] = params.getStartPointForKey(s);
      }
    }

    double energy = 0.0;
    final double[] p = params.getAllParamters();
    for (int i = 0; i < noOfAtoms; i++) {

      final int pos =
          (useCaching) ? posCache[i] : params.getStartPointForKey("totalenergyshifter:" + atoms[i]);

      energy += p[pos];
    }

    return energy;
  }

  @Override
  public double partialParamGradient(
      final Topology topology, final AdaptiveParameters params, final double[] grad) {

    final String[] atoms = topology.getAtomNames();
    final int noOfAtoms = topology.getNumberOfAtoms();
    final double[] p = params.getAllParamters();

    if (useCaching && posCache == null) {
      // initialize cache
      posCache = new int[noOfAtoms];
      for (int i = 0; i < noOfAtoms; i++) {
        final String s = "totalenergyshifter:" + atoms[i];
        posCache[i] = params.getStartPointForKey(s);
      }
    }

    double energy = 0.0;
    for (int i = 0; i < noOfAtoms; i++) {
      final int pos =
          (useCaching) ? posCache[i] : params.getStartPointForKey("totalenergyshifter:" + atoms[i]);
      grad[pos]++;
      energy += p[pos];
    }

    return energy;
  }

  @Override
  public Gradient partialCartesianGradient(
      final Topology topology, final AdaptiveParameters params) {

    // there is no cartesian gradient, obviously.
    final Gradient g = new Gradient();
    final double[][] grad = new double[3][topology.getNumberOfAtoms()];
    g.setGradientTotal(grad);

    final String[] atoms = topology.getAtomNames();
    final int noOfAtoms = topology.getNumberOfAtoms();

    if (useCaching && posCache == null) {
      // initialize cache
      posCache = new int[noOfAtoms];
      for (int i = 0; i < noOfAtoms; i++) {
        final String s = "totalenergyshifter:" + atoms[i];
        posCache[i] = params.getStartPointForKey(s);
      }
    }

    double energy = 0.0;
    final double[] p = params.getAllParamters();
    for (int i = 0; i < noOfAtoms; i++) {

      final int pos =
          (useCaching) ? posCache[i] : params.getStartPointForKey("totalenergyshifter:" + atoms[i]);

      energy += p[pos];
    }
    g.setTotalEnergy(energy);

    return g;
  }

  @Override
  public Tuple3D<String[], int[], Integer> requiredParams(
      final ArrayList<CartesianCoordinates> cartesians,
      final ArrayList<Topology> topologies,
      final String sMethod) {

    final ArrayList<String> allAtoms = new ArrayList<>();

    for (final CartesianCoordinates cartes : cartesians) {
      final String[] atoms = cartes.getAllAtomTypes();
      for (final String atom : atoms) {
        if (!allAtoms.contains(atom)) allAtoms.add(atom);
      }
    }

    // and back...
    int paramSum = 0;
    final String[] keys = new String[allAtoms.size()];
    final int[] paramsPerKey = new int[allAtoms.size()];
    for (int i = 0; i < allAtoms.size(); i++) {
      keys[i] = "totalenergyshifter:" + allAtoms.get(i);
      paramsPerKey[i] = 1;
      paramSum++;
    }

    final Tuple3D<String[], int[], Integer> result = new Tuple3D<>(keys, paramsPerKey, paramSum);

    return result;
  }

  @Override
  public double[][] bordersForMyParams(final AdaptiveParameters params) {

    // not all of the parameters belong to us
    final List<String> corrKeys = params.getKeysStartingWith("totalenergyshifter:");
    final int noOfKeys = corrKeys.size();

    double[][] borders = new double[2][noOfKeys];

    for (int i = 0; i < noOfKeys; i++) {
      // in hartree seems reasonable.
      borders[0][i] = -100.0;
      borders[1][i] = 0.0;
    }

    return borders;
  }
}
