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
import org.ogolem.core.Gradient;
import org.ogolem.core.UniversalFF;

/**
 * Tests UFF for the energy of a large (1428 atoms) molecule.
 *
 * @author Johannes Dieterich
 * @version 2020-08-16
 */
class UFFLargeMoleculeGradientBenchmark implements SingleMicroBenchmark {

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
  private final Gradient gradient;

  UFFLargeMoleculeGradientBenchmark() {

    final CartesianCoordinates cartes = CartesianCoordinatesLibrary.getUFFLargeMolecule();

    this.xyz1D = cartes.getAll1DCartes();
    this.atoms = cartes.getAllAtomTypes();
    this.atomNos = cartes.getAllAtomNumbers();
    this.atsPerMol = cartes.getAllAtomsPerMol();
    this.energyparts = new double[cartes.getNoOfMolecules()];
    this.noAtoms = cartes.getNoOfAtoms();
    this.charges = cartes.getAllCharges();
    this.spins = cartes.getAllSpins();
    this.bonds = CartesianCoordinatesLibrary.getUFFLargeMoleculeBonds();
    this.uff = new UniversalFF(0, true);
    this.gradient = new Gradient(3, cartes.getNoOfAtoms());
  }

  @Override
  public String name() {
    return "large molecule UFF gradient";
  }

  @Override
  public double runSingle() throws Exception {

    uff.gradientCalculation(
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
        gradient,
        false);

    return gradient.getTotalEnergy();
  }
}
