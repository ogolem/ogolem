/*
Copyright (c) 2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.microbenchmarks;

import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CoordTranslation;
import org.ogolem.core.UniversalFF;

/**
 * Tests UFF energy cost.
 *
 * @author Johannes Dieterich
 * @version 2020-08-02
 */
class UFFFrozenSurfaceEnergyBenchmark implements SingleMicroBenchmark {

  private final UniversalFF uff;
  private final double[] xyz1D;
  private final String[] atoms;
  private final short[] atomNos;
  private final int[] atsPerMol;
  private final double[] energyparts;
  private final int noAtoms;
  private final float[] charges;
  private final short[] spins;
  private final BondInfo bonds;

  UFFFrozenSurfaceEnergyBenchmark() {

    final CartesianCoordinates cartes = CartesianCoordinatesLibrary.getCO10onPd200UFFMin();

    this.xyz1D = cartes.getAll1DCartes();
    this.atoms = cartes.getAllAtomTypes();
    this.atomNos = cartes.getAllAtomNumbers();
    this.atsPerMol = cartes.getAllAtomsPerMol();
    this.energyparts = new double[cartes.getNoOfMolecules()];
    this.noAtoms = cartes.getNoOfAtoms();
    this.charges = cartes.getAllCharges();
    this.spins = cartes.getAllSpins();
    this.bonds = CoordTranslation.checkForBonds(cartes, 1.0);

    // remove bonds inside the surface
    for (int i = 20; i < 221; i++) {
      for (int j = 20; j < 221; j++) {
        bonds.setBond(i, j, BondInfo.NOBOND);
      }
    }

    // remove bonds between surface and cluster
    for (int i = 0; i < 20; i++) {
      for (int j = 20; j < 221; j++) {
        bonds.setBond(i, j, BondInfo.NOBOND);
        bonds.setBond(j, i, BondInfo.NOBOND);
      }
    }

    this.uff = new UniversalFF(0, true);
  }

  @Override
  public String name() {
    return "CO10 on frozen Pd200 UFF energy";
  }

  @Override
  public double runSingle() throws Exception {

    final double e =
        uff.energyCalculation(
            1,
            1,
            xyz1D,
            atoms,
            atomNos,
            atsPerMol,
            energyparts,
            noAtoms,
            charges,
            spins,
            bonds,
            true);

    return e;
  }
}
