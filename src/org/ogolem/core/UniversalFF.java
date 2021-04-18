/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

import de.gwdg.rmata.atomdroiduff.AtomdroidUFF;

/**
 * Provides energy and gradient calculated using a very simple UFF approach.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class UniversalFF implements CartesianFullBackend {

  // the ID
  private static final long serialVersionUID = (long) 20200830;

  private final int style;
  private final boolean useCaching;
  private int noAtoms;
  private AtomdroidUFF uff;
  private double[][] coords;
  private boolean hasBeenUsed = false;

  public UniversalFF(final int whichStyle, final boolean useCaching) {
    this.style = whichStyle;
    this.useCaching = useCaching;
  }

  UniversalFF(final UniversalFF orig) {
    this.style = orig.style;
    this.useCaching = orig.useCaching;

    if (orig.uff != null) {
      this.uff = new AtomdroidUFF(orig.uff);
    } else {
      this.uff = null;
    }

    this.noAtoms = orig.noAtoms;
    if (orig.coords != null) {
      this.coords = new double[noAtoms][3];
    } else {
      this.coords = null;
    }
  }

  @Override
  public UniversalFF copy() {
    return new UniversalFF(this);
  }

  @Override
  public String getMethodID() {
    return "UniversalFF (Atomdroid implementation)";
  }

  @Override
  public double energyCalculation(
      long lID,
      int iIteration,
      double[] daXYZ1D,
      String[] atoms,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      int noAts,
      float[] charges,
      short[] iaSpins,
      final BondInfo bonds,
      final boolean hasRigidEnv) {

    if (uff == null || !useCaching) {
      initialize(noAts, bonds, atoms, charges, atomNos, atsPerMol, hasRigidEnv);
    }

    copy1d3d(daXYZ1D, coords, noAts, atomNos);
    boolean isFirst = false;
    if (!this.hasBeenUsed) {
      isFirst = true;
      hasBeenUsed = true;
    }
    // final boolean isFirst = (iIteration == -1 || iIteration == 0);
    return uff.energy(style, coords, isFirst) * Constants.KJTOHARTREE;
  }

  @Override
  public void gradientCalculation(
      long lID,
      int iIteration,
      double[] daXYZ1D,
      String[] atoms,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      int noAts,
      float[] charges,
      short[] iaSpins,
      final BondInfo bonds,
      final Gradient gradient,
      final boolean hasRigidEnv) {

    if (uff == null || !useCaching) {
      initialize(noAts, bonds, atoms, charges, atomNos, atsPerMol, hasRigidEnv);
    }

    copy1d3d(daXYZ1D, coords, noAts, atomNos);

    de.gwdg.rmata.atomdroiduff.Gradient g = uff.gradient(style, coords, !hasBeenUsed);
    if (!hasBeenUsed) {
      hasBeenUsed = true; // now it has
    }

    gradient.setTotalEnergy(g.getEnergy() * Constants.KJTOHARTREE);
    final double[][] gradVals = gradient.getTotalGradient();
    copy3d3d(g.getGradient(), noAts, gradVals, atomNos);
  }

  private void initialize(
      final int noAts,
      final BondInfo bonds,
      final String[] names,
      final float[] charges,
      final short[] atomNos,
      final int[] atsPerMol,
      final boolean hasRigidEnv) {

    int noDummy = 0;
    for (int i = 0; i < atomNos.length; i++) {
      if (atomNos[i] == 0) {
        noDummy++;
      }
    }

    noAtoms = noAts - noDummy;
    coords = new double[noAtoms][3];

    // let's work on the bonds first...
    int numBonds = 0;
    final short[][] bonding = new short[noAtoms][noAtoms];
    final double[][] pos = new double[noAtoms][3];
    int c1 = 0;
    for (int i = 0; i < noAts; i++) {
      if (atomNos[i] == 0) {
        continue;
      }
      int c2 = 0;
      for (int j = 0; j < noAts; j++) {
        if (atomNos[j] == 0) {
          continue;
        }
        final short bond = bonds.bondType(i, j);
        if (bond == BondInfo.UNCERTAIN) {
          bonding[c1][c2] = 1;
          numBonds++;
        } else if (bond == BondInfo.VDW) {
          // vdW is no bond for us here
        } else {
          bonding[c1][c2] = bond;
          numBonds++;
        }
        c2++;
      }
      c1++;
    }

    int noFrozen = (hasRigidEnv) ? atsPerMol[atsPerMol.length - 1] : 0;
    int noUnfrozen = noAtoms - noFrozen;

    int noDummies = 0;
    final float[] ch = new float[noAtoms];
    final String[] na = new String[noAtoms];
    final short[] nos = new short[noAtoms];
    c1 = 0;
    for (int i = 0; i < noAts; i++) {
      if (atomNos[i] == 0) {
        if (i >= noUnfrozen) noDummies++;
        continue;
      }
      ch[c1] = charges[i];
      na[c1] = names[i];
      nos[c1] = atomNos[i];
      c1++;
    }

    // subtract dummy atoms that used to be in the frozen part from the number of frozen atoms
    if (hasRigidEnv) noFrozen -= noDummies;

    uff = new AtomdroidUFF(noAtoms, noFrozen, na, bonding, ch, pos, numBonds, nos);
  }

  private static void copy1d3d(
      final double[] x1d, final double[][] x3d, final int noAts, final short[] atomNos) {

    int c = 0;
    for (int i = 0; i < noAts; i++) {
      if (atomNos[i] == 0) {
        continue;
      }
      x3d[c][0] = x1d[i] * Constants.BOHRTOANG;
      x3d[c][1] = x1d[i + noAts] * Constants.BOHRTOANG;
      x3d[c][2] = x1d[i + 2 * noAts] * Constants.BOHRTOANG;
      c++;
    }
  }

  private static void copy3d3d(
      final double[][] o, final int noAts, final double[][] g, final short[] atomNos) {

    assert (g.length == 3);
    assert (g[0].length == noAts);
    assert (g[1].length == noAts);
    assert (g[2].length == noAts);

    int c = 0;
    for (int i = 0; i < noAts; i++) {
      if (atomNos[i] == 0) {
        g[0][i] = 0.0;
        g[1][i] = 0.0;
        g[2][i] = 0.0;
      } else {
        g[0][i] = o[c][0] * Constants.BOHRTOANG * Constants.KJTOHARTREE;
        g[1][i] = o[c][1] * Constants.BOHRTOANG * Constants.KJTOHARTREE;
        g[2][i] = o[c][2] * Constants.BOHRTOANG * Constants.KJTOHARTREE;
        c++;
      }
    }
  }
}
