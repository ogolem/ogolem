/*
Copyright (c) 2012-2014, J. M. Dieterich
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

import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.stats.GenericDetailStatistics;

/**
 * A simple adapter for the geometry to gradient/energy provider.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class FullyCartesianCoordinates implements GenericBackend<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20200429;

  public static final double MAXINCREMENTBOUNDS = 10.0; // in bohr

  private final CartesianFullBackend back;

  private Gradient grad;
  private long id;
  private String[] atoms;
  private short[] atomNos;
  private float[] charges;
  private short[] spins;
  private BondInfo bonds;
  private int[] atsPerMol;
  private double[] energyparts;
  private boolean hasRigidEnv;

  public FullyCartesianCoordinates(final CartesianFullBackend back) {
    this.back = back;
  }

  public FullyCartesianCoordinates(final FullyCartesianCoordinates orig) {
    this.back = orig.back.copy();
  }

  @Override
  public FullyCartesianCoordinates copy() {
    return new FullyCartesianCoordinates(this);
  }

  @Override
  public String getMyID() {
    return back.getMethodID();
  }

  @Override
  public int numberOfActiveCoordinates(final Geometry individual) {
    return individual.getNumberOfAtoms() * 3;
  }

  @Override
  public double[] getActiveCoordinates(final Geometry individual) {

    final CartesianCoordinates cartes = individual.getCartesians();

    this.atomNos = cartes.getAllAtomNumbers();
    this.atoms = cartes.getAllAtomTypes();
    this.atsPerMol = cartes.getAllAtomsPerMol();
    this.bonds = individual.getBondInfo();
    this.charges = cartes.getAllCharges();
    this.spins = cartes.getAllSpins();
    this.energyparts = new double[individual.getNumberOfIndieParticles()];
    this.grad = new Gradient(3, atomNos.length);
    this.hasRigidEnv = cartes.containedEnvType() == CartesianCoordinates.ENVTYPE.RIGID;

    return cartes.getAll1DCartes();
  }

  @Override
  public void resetToStable(final double[] coordinates) {
    // nothing
  }

  @Override
  public void updateActiveCoordinates(final Geometry individual, final double[] coordinates) {

    final CartesianCoordinates cartes = individual.getCartesians();
    cartes.setAll1DCartes(coordinates, cartes.getNoOfAtoms());

    CoordTranslation.updateGeometryFromCartesian(cartes, individual);
  }

  @Override
  public double gradient(final double[] currCoords, final double[] gradient, final int iteration) {

    GenericDetailStatistics.incrementFitnessEvals();

    assert (currCoords.length == gradient.length);
    assert (currCoords.length == atoms.length * 3);

    back.gradientCalculation(
        id,
        iteration,
        currCoords,
        atoms,
        atomNos,
        atsPerMol,
        energyparts,
        atomNos.length,
        charges,
        spins,
        bonds,
        grad,
        hasRigidEnv);
    grad.getGradientData(gradient);

    final double e = grad.getTotalEnergy();

    assert (!Double.isInfinite(e));
    assert (!Double.isNaN(e));

    return e;
  }

  @Override
  public double fitness(final double[] currCoords, final int iteration) {
    GenericDetailStatistics.incrementFitnessEvals();
    assert (currCoords.length == atoms.length * 3);

    final double e =
        back.energyCalculation(
            id,
            iteration,
            currCoords,
            atoms,
            atomNos,
            atsPerMol,
            energyparts,
            atoms.length,
            charges,
            spins,
            bonds,
            hasRigidEnv);

    assert (!Double.isInfinite(e));
    assert (!Double.isNaN(e));

    return e;
  }

  @Override
  public Geometry fitness(final Geometry individual, final boolean forceOneEval) {

    GenericDetailStatistics.incrementFitnessEvals();

    final CartesianCoordinates c = individual.getCartesians();

    double[] eparts;
    if (energyparts != null && energyparts.length == c.getNoOfMolecules()) {
      eparts = energyparts;
    } else {
      eparts = new double[c.getNoOfMolecules()];
    }
    final double e =
        back.energyCalculation(
            id,
            42,
            c.getAll1DCartes(),
            c.getAllAtomTypes(),
            c.getAllAtomNumbers(),
            c.getAllAtomsPerMol(),
            eparts,
            c.getNoOfAtoms(),
            c.getAllCharges(),
            c.getAllSpins(),
            individual.getBondInfo(),
            c.containedEnvType() == CartesianCoordinates.ENVTYPE.RIGID);

    assert (!Double.isInfinite(e));
    assert (!Double.isNaN(e));

    individual.setFitness(e);

    return individual;
  }

  CartesianFullBackend getMyBackend() {
    return back;
  }

  @Override
  public BOUNDSTYPE boundariesInRepresentation(final Geometry individual) {
    return BOUNDSTYPE.NONE;
  }

  @Override
  public void bestEstimateBoundaries(
      final double[] currCoords, final double[] low, final double[] high) {

    assert (low.length == high.length);
    assert (low.length == currCoords.length);

    for (int i = 0; i < high.length; i++) {
      high[i] = currCoords[i] + MAXINCREMENTBOUNDS;
      low[i] = currCoords[i] - MAXINCREMENTBOUNDS;
    }
  }
}
