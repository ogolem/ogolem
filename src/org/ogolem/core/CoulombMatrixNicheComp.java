/*
 * Copyright (c) 2013, J. M. Dieterich
 *               2015, J. M. Dieterich and B. Hartke
 *               2018, B. Hartke
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted
 * provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list of conditions
 * and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this list of
 * conditions and the following disclaimer in the documentation and/or other materials provided with
 * the distribution.
 *
 * All advertising materials mentioning features or use of this software must display the
 * following acknowledgement:
 *
 * This product includes software of the ogolem.org project developed by J. M. Dieterich and B.
 * Hartke (Christian-Albrechts-University Kiel, Germany) and contributors.
 *
 * Neither the name of the ogolem.org project, the University of Kiel nor the names of its
 * contributors may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package org.ogolem.core;

import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;

/**
 * Nicher based on eigenvalues of geometry Coulomb matrix w/o charges The discrimination power
 * increases with width becoming smaller and with number getting larger.
 *
 * @author Bernd Hartke
 * @version 2021-07-06
 */
class CoulombMatrixNicheComp implements NicheComputer<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20180513;
  private final double width;
  private final int number;

  public static final double DEFAULTWIDTH = 0.1;
  public static final int DEFAULTNUMBER = 1;

  CoulombMatrixNicheComp(final double width, final int number) {
    this.width = width;
    this.number = number;
  }

  private CoulombMatrixNicheComp(final CoulombMatrixNicheComp orig) {
    this.width = orig.width;
    this.number = orig.number;
  }

  @Override
  public NicheComputer<Molecule, Geometry> copy() {
    return new CoulombMatrixNicheComp(this);
  }

  @Override
  public Niche computeNiche(final Geometry g) {
    final double[][] allCOMs = g.getAllCOMs();
    final int noMols = allCOMs[0].length;
    final contrib.jama.Matrix cMatrix = new contrib.jama.Matrix(noMols, noMols);
    setupCoulombMatrix(allCOMs, cMatrix.getArray(), noMols);
    final contrib.jama.EigenvalueDecomposition eigen =
        new contrib.jama.EigenvalueDecomposition(cMatrix, true);
    double[] evals = eigen.getRealEigenvalues();
    double[] percDiffs = new double[noMols - 1];
    // this baseline is arbitrary; evals[0] would be more logical, but this would give div-by-zero
    // below
    final double baseline = evals[0] - 0.1 * Math.abs(evals[0]);
    for (int i = 0; i < noMols - 1; i++) {
      percDiffs[i] = (evals[i + 1] - evals[i]) / (evals[i] - baseline);
    }
    /*
     * now we coarse-grain this eigenvalue vector, to transform it into a niche designator;
     * the idea is to concatenate the top "number" entries from the percDiffs-array (discarding the
     * others), after coarsening them to "width" (both number and width are input options).
     * The coarsening should be done by: (1) divide by width, (2) round-to-integer,
     * (3) multiply by width. However, we skip step (3) here, since it does not really do anything
     * except providing a re-scaling of the (coarsened) numbers, but invites lots of trouble
     * if one wants to avoid including not-needed decimal places into the niche string...
     * Sadly, skipping step (3) destroys the intuitive relationship between input width values
     * and resulting niche designators :-(
     */
    String designator = "cmat-";
    for (int i = noMols - 2;
        i > noMols - 2 - number;
        i--) { // takes "number" values from the end of the percDiffs array, in inverted ordering
      /*
       * the following line could almost be the solution, if there were a foolproof way of getting numDigits, the number of significant
       *digits in a real, despite possible rounding/representation errors:
       *            designator = designator + String.format("%."+numDigits+"f",Math.round(percDiffs[i]/width)*width) + "-";
       * Hence, we drop step (3) and simply do this:
       */
      designator = designator + Math.round(percDiffs[i] / width) + "-";
    }
    designator =
        designator.substring(
            0, designator.length() - 1); // remove trailing dash, just for aesthetics...
    return new Niche(designator);
  }

  private static void setupCoulombMatrix(
      final double[][] allCOMs, final double[][] coulombMatrix, final int noMols) {
    // sets up neutral Coulomb matrix (no charges); full matrix with explizit symmetrization
    for (int i = 0; i < noMols; i++) {
      for (int j = i + 1; j < noMols; j++) {
        final double dx = allCOMs[0][i] - allCOMs[0][j];
        final double dy = allCOMs[1][i] - allCOMs[1][j];
        final double dz = allCOMs[2][i] - allCOMs[2][j];
        final double n = Math.sqrt(dx * dx + dy * dy + dz * dz);
        final double distInv = n == 0 ? 0.0 : 1.0 / n;
        coulombMatrix[i][j] = distInv; // should have been Z_i Z_j / dist_ij in case of charges
        coulombMatrix[j][i] = distInv; // explicit symmetrization
      }
      coulombMatrix[i][i] = 0.0; // should have been some simple function of charge (e.g. 0.5*Z^2.4)
    }
  }
}
