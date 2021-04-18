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
package org.ogolem.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Chains multiple cartesian full backends after each other.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class ChainedCartesianFullBackend implements CartesianFullBackend {

  private static final long serialVersionUID = (long) 20200720;

  private final List<CartesianFullBackend> backends;
  private double[] gradCache;
  private double[] partsCache;

  ChainedCartesianFullBackend(final List<CartesianFullBackend> backends) {

    assert (backends != null);
    assert (backends.size() > 0);

    this.backends = backends;
  }

  private ChainedCartesianFullBackend(final ChainedCartesianFullBackend orig) {

    this.backends = new ArrayList<>();
    for (final CartesianFullBackend back : orig.backends) {
      backends.add(back.copy());
    }
  }

  @Override
  public ChainedCartesianFullBackend copy() {
    return new ChainedCartesianFullBackend(this);
  }

  @Override
  public String getMethodID() {

    String s = "chained:\n";
    for (final CartesianFullBackend back : backends) {
      s += "\t\t" + back.getMethodID() + "\n";
    }

    return s;
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

    if (gradCache == null || gradCache.length < 3 * iNoOfAtoms) {
      gradCache = new double[3 * iNoOfAtoms];
    } else {
      Arrays.fill(gradCache, 0.0);
    }

    if (partsCache == null || partsCache.length < energyparts.length) {
      partsCache = new double[energyparts.length];
    } else {
      Arrays.fill(partsCache, 0.0);
    }

    int lastOffset = 0;
    for (int i = 0; i < atsPerMol.length - 1; i++) {
      lastOffset += atsPerMol[i];
    }

    final int gradEnd = (hasRigidEnv) ? lastOffset : iNoOfAtoms;

    // just add up
    double e = 0.0;
    for (final CartesianFullBackend back : backends) {
      back.gradientCalculation(
          lID,
          iIteration,
          xyz1D,
          saAtomTypes,
          atomNos,
          atsPerMol,
          energyparts,
          iNoOfAtoms,
          faCharges,
          spins,
          bonds,
          grad,
          hasRigidEnv);

      final double[][] workGradMat = grad.getTotalGradient();
      for (int coord = 0; coord < 3; coord++) {
        for (int at = 0; at < gradEnd; at++) {
          gradCache[coord * iNoOfAtoms + at] += workGradMat[coord][at];
        }
      }

      for (int i = 0; i < energyparts.length; i++) {
        partsCache[i] += energyparts[i];
      }

      e += grad.getFunctionValue();
    }
    final double[][] gradMat = grad.getTotalGradient();
    System.arraycopy(gradCache, 0, gradMat[0], 0, iNoOfAtoms);
    System.arraycopy(gradCache, iNoOfAtoms, gradMat[1], 0, iNoOfAtoms);
    System.arraycopy(gradCache, iNoOfAtoms * 2, gradMat[2], 0, iNoOfAtoms);

    grad.setFunctionValue(e);

    System.arraycopy(partsCache, 0, energyparts, 0, energyparts.length);
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

    if (partsCache == null || partsCache.length < energyparts.length) {
      partsCache = new double[energyparts.length];
    } else {
      Arrays.fill(partsCache, 0.0);
    }

    // just add up
    double e = 0.0;
    for (final CartesianFullBackend back : backends) {
      e +=
          back.energyCalculation(
              lID,
              iIteration,
              xyz1D,
              saAtomTypes,
              atomNos,
              atsPerMol,
              energyparts,
              iNoOfAtoms,
              faCharges,
              spins,
              bonds,
              hasRigidEnv);

      for (int i = 0; i < energyparts.length; i++) {
        partsCache[i] += energyparts[i];
      }
    }

    System.arraycopy(partsCache, 0, energyparts, 0, energyparts.length);

    return e;
  }
}
