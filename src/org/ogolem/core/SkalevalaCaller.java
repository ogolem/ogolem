/*
Copyright (c) 2012-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
              2017-2020, J. M. Dieterich and B. Hartke
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

import org.ogolem.skalevala.Configuration;
import org.ogolem.skalevala.EHNDO;
import org.ogolem.skalevala.EHNDOParameters;
import org.ogolem.skalevala.Runot;

/**
 * Calls the sKalevala part of ogolem.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class SkalevalaCaller implements CartesianFullBackend {
  // TODO debug, extend
  // the ID
  private static final long serialVersionUID = (long) 20120103;

  @Override
  public SkalevalaCaller copy() {
    return new SkalevalaCaller();
  }

  @Override
  public String getMethodID() {
    return "sKalevala powered EHNDO";
  }

  @Override
  public double energyCalculation(
      long lID,
      int iIteration,
      double[] daXYZ1D,
      String[] saAtomTypes,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      int iNoOfAtoms,
      float[] faCharges,
      short[] iaSpins,
      final BondInfo bonds,
      final boolean hasRigidEnv) {

    // put the 1D coordinates into 3D ones (layout as sKalevala wants)
    final double[][] daXYZ = new double[iNoOfAtoms][3];
    for (int i = 0; i < iNoOfAtoms; i++) {
      for (int j = 0; j < 3; j++) {
        daXYZ[i][j] = daXYZ1D[j * iNoOfAtoms + i] * Constants.BOHRTOANG;
      }
    }

    float fcharge = 0;
    int spin = 0;
    for (int i = 0; i < faCharges.length; i++) {
      fcharge += faCharges[i];
      spin += iaSpins[i];
    }

    Configuration.printTimings_$eq(false);
    // System.out.println("DEBUG: EHNDO entering energy...");
    final Runot runot = new EHNDO();
    final org.ogolem.skalevala.CartesianCoordinates cartes =
        new org.ogolem.skalevala.CartesianCoordinates(
            daXYZ, saAtomTypes, atomNos, spin, (int) fcharge, bonds.getFullBondMatrix());
    final EHNDOParameters params = new EHNDOParameters();
    params.initializeDefaults();

    return runot.energy(cartes, params) * Constants.EVTOHARTREE;
  }

  @Override
  public void gradientCalculation(
      long lID,
      int iIteration,
      double[] daXYZ1D,
      String[] saAtomTypes,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      int iNoOfAtoms,
      float[] faCharges,
      short[] iaSpins,
      final BondInfo bonds,
      final Gradient gradient,
      final boolean hasRigidEnv) {

    for (int i : atomNos) {
      System.out.println("DEBUG: Have " + i);
    }

    assert (atomNos.length == iNoOfAtoms);
    assert (faCharges.length == iNoOfAtoms);

    // put the 1D coordinates into 3D ones (layout as sKalevala wants)
    final double[][] xyz = new double[iNoOfAtoms][3];
    for (int i = 0; i < iNoOfAtoms; i++) {
      for (int j = 0; j < 3; j++) {
        xyz[i][j] = daXYZ1D[j * iNoOfAtoms + i] * Constants.BOHRTOANG;
      }
    }

    float fcharge = 0;
    int spin = 0;
    for (int i = 0; i < faCharges.length; i++) {
      fcharge += faCharges[i];
      spin += iaSpins[i];
    }

    Configuration.printTimings_$eq(false);
    final Runot runot = new EHNDO();
    final org.ogolem.skalevala.CartesianCoordinates cartes =
        new org.ogolem.skalevala.CartesianCoordinates(
            xyz, saAtomTypes, atomNos, spin, (int) fcharge, bonds.getFullBondMatrix());
    final EHNDOParameters params = new EHNDOParameters();
    params.initializeDefaults();

    // System.out.println("DEBUG: EHNDO entering gradient...");
    final org.ogolem.skalevala.Gradient g = runot.gradient(cartes, params);

    final double[][] gradMat = gradient.getTotalGradient();
    final double[][] origMat = g.gradient().toArray();
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < iNoOfAtoms; j++) {
        gradMat[i][j] = origMat[j][i] * Constants.EVTOHARTREE * Constants.BOHRTOANG;
      }
    }

    // System.out.println("DEBUG: EHNDO exiting gradient...");

    gradient.setTotalEnergy(g.energy() * Constants.EVTOHARTREE);
    System.out.println("Energy was " + g.energy());
    System.out.println("Energy is " + gradient.getTotalEnergy());
  }
}
