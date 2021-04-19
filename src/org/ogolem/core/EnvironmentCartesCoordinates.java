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

import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.stats.GenericDetailStatistics;
import org.ogolem.helpers.Tuple;

/**
 * A simple adapter for the geometry to gradient/energy provider considering environments. I.e., we
 * keep the environment frozen at all times. Does NOT consider the bounds to the cluster location
 * w.r.t. the environment (yet?).
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class EnvironmentCartesCoordinates implements GenericBackend<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20200525;

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
  private int totNoAtoms;
  private int atomsInEnv;
  private int noAtomsEnvless;
  private double[] fullCoords;
  private boolean hasRigidEnv;

  public EnvironmentCartesCoordinates(final CartesianFullBackend back) {
    this.back = back;
  }

  public EnvironmentCartesCoordinates(final EnvironmentCartesCoordinates orig) {
    this.back = orig.back.copy();
  }

  @Override
  public EnvironmentCartesCoordinates copy() {
    return new EnvironmentCartesCoordinates(this);
  }

  @Override
  public String getMyID() {
    return "environment enabled: " + back.getMethodID();
  }

  @Override
  public int numberOfActiveCoordinates(final Geometry individual) {
    return individual.getNumberOfAtoms() * 3;
  }

  @Override
  public double[] getActiveCoordinates(final Geometry individual) {

    if (!individual.containsEnvironment()) {
      throw new RuntimeException(
          "Geometry does not contain an environment, which is required for EnvironmentCartesCoordinates!");
    }

    final Tuple<CartesianCoordinates, BondInfo> tup =
        individual.getCartesiansAndBondsWithEnvironment();
    final CartesianCoordinates cartes = tup.getObject1();
    final BondInfo bonds = tup.getObject2();

    final Environment env = individual.getEnvironment();
    this.atomsInEnv = env.atomsInEnv();

    this.atomNos = cartes.getAllAtomNumbers();
    this.atoms = cartes.getAllAtomTypes();
    this.atsPerMol = cartes.getAllAtomsPerMol();
    this.bonds = bonds;
    this.charges = cartes.getAllCharges();
    this.spins = cartes.getAllSpins();
    this.hasRigidEnv = cartes.containedEnvType() == CartesianCoordinates.ENVTYPE.RIGID;

    this.energyparts =
        new double[individual.getNumberOfIndieParticles() + this.atomsInEnv]; // for the environment
    this.grad = new Gradient(3, atomNos.length);

    this.totNoAtoms = cartes.getNoOfAtoms();
    this.fullCoords = cartes.getAll1DCartes();

    this.noAtomsEnvless = totNoAtoms - atomsInEnv;

    final double[] envlessCartes = new double[3 * noAtomsEnvless];
    System.arraycopy(fullCoords, 0, envlessCartes, 0, noAtomsEnvless); // x component
    System.arraycopy(fullCoords, totNoAtoms, envlessCartes, noAtomsEnvless, noAtomsEnvless); // y
    System.arraycopy(
        fullCoords, 2 * totNoAtoms, envlessCartes, 2 * noAtomsEnvless, noAtomsEnvless); // z

    return envlessCartes;
  }

  @Override
  public void resetToStable(final double[] coordinates) {
    // nothing
  }

  @Override
  public void updateActiveCoordinates(final Geometry individual, final double[] coordinates) {

    final CartesianCoordinates cartes = individual.getCartesiansWithEnvironment();
    final double[][] xyz = cartes.getAllXYZCoord();

    System.arraycopy(coordinates, 0, xyz[0], 0, noAtomsEnvless); // x component
    System.arraycopy(coordinates, noAtomsEnvless, xyz[1], 0, noAtomsEnvless); // y
    System.arraycopy(coordinates, 2 * noAtomsEnvless, xyz[2], 0, noAtomsEnvless); // z

    CoordTranslation.updateGeometryFromCartesian(cartes, individual);
  }

  @Override
  public double gradient(final double[] currCoords, final double[] gradient, final int iteration) {

    GenericDetailStatistics.incrementFitnessEvals();

    // curr coords must be rewritten (we are reusing the existing env information in fullCoords)
    System.arraycopy(currCoords, 0, fullCoords, 0, noAtomsEnvless); // x component
    System.arraycopy(currCoords, noAtomsEnvless, fullCoords, totNoAtoms, noAtomsEnvless); // y
    System.arraycopy(
        currCoords, 2 * noAtomsEnvless, fullCoords, 2 * totNoAtoms, noAtomsEnvless); // z

    assert (currCoords.length == noAtomsEnvless * 3);

    back.gradientCalculation(
        id,
        iteration,
        fullCoords,
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

    // remove gradient elements belonging to environment
    final double[][] gradDat = grad.getTotalGradient();
    System.arraycopy(gradDat[0], 0, gradient, 0, noAtomsEnvless); // x component
    System.arraycopy(gradDat[1], 0, gradient, noAtomsEnvless, noAtomsEnvless); // y
    System.arraycopy(gradDat[2], 0, gradient, 2 * noAtomsEnvless, noAtomsEnvless); // z

    final double e = grad.getTotalEnergy();

    assert (!Double.isInfinite(e));
    assert (!Double.isNaN(e));

    return e;
  }

  @Override
  public double fitness(final double[] currCoords, final int iteration) {

    GenericDetailStatistics.incrementFitnessEvals();

    // curr coords must be rewritten (we are reusing the existing env information in fullCoords)
    System.arraycopy(currCoords, 0, fullCoords, 0, noAtomsEnvless); // x component
    System.arraycopy(currCoords, noAtomsEnvless, fullCoords, totNoAtoms, noAtomsEnvless); // y
    System.arraycopy(
        currCoords, 2 * noAtomsEnvless, fullCoords, 2 * totNoAtoms, noAtomsEnvless); // z

    assert (currCoords.length == noAtomsEnvless * 3);

    final double e =
        back.energyCalculation(
            id,
            iteration,
            fullCoords,
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

    final CartesianCoordinates c = individual.getCartesiansWithEnvironment();

    double[] eparts;
    if (energyparts != null && energyparts.length == c.getNoOfMolecules()) {
      eparts = energyparts;
    } else {
      eparts =
          new double
              [c.getNoOfMolecules()
                  + individual.getEnvironment().atomsInEnv()]; // for the environment
    }

    // that one will simply work. :-)
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
