/*
Copyright (c) 2010-2013, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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

/**
 * Adds a penalty function to an arbitrary local optimization.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class PenaltyLocOpt implements Newton {

  // the ID
  private static final long serialVersionUID = (long) 20101108;

  private final Newton locopt;

  private final AdditivePenaltyFunction penalty;

  PenaltyLocOpt(final GlobalConfig globConf, final Newton newton) {
    locopt = newton;

    if (globConf.penalty == null) {
      penalty = null;
    } else {
      penalty = globConf.penalty.copy();
    }
  }

  private PenaltyLocOpt(final PenaltyLocOpt orig) {
    locopt = orig.locopt.copy();
    if (orig.penalty == null) {
      penalty = null;
    } else {
      penalty = orig.penalty.copy();
    }
  }

  @Override
  public PenaltyLocOpt copy() {
    return new PenaltyLocOpt(this);
  }

  @Override
  public String myIDandMethod() {
    return "Penalty LocOpt with " + locopt.myIDandMethod();
  }

  @Override
  public Geometry localOptimization(Geometry geom) {

    // optimize
    geom = locopt.localOptimization(geom);

    if (penalty == null) {
      System.err.println("ERROR: No penalty function defined.");
      return geom;
    }

    // penalty
    final double dPenalty = penalty.penalty(geom);

    // put in
    final double dNewFitness = geom.getFitness() + dPenalty;
    geom.setFitness(dNewFitness);

    return geom;
  }

  @Override
  public Molecule localOptimization(Molecule molecule) {

    // optimize
    molecule = locopt.localOptimization(molecule);

    // penalty
    final double dPenalty = penalty.penalty(molecule);

    // put in
    final double dNewFitness = molecule.getEnergy() + dPenalty;
    molecule.setEnergy(dNewFitness);

    return molecule;
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      long id,
      CartesianCoordinates cartes,
      boolean[][] constraints,
      boolean isConstricted,
      final BondInfo bonds)
      throws Exception {

    // optimize
    cartes = locopt.cartesToCartes(id, cartes, constraints, isConstricted, bonds);

    // penalty
    final double penaltyAdd = penalty.penalty(cartes);

    // put in
    final double newFitness = cartes.getEnergy() + penaltyAdd;
    cartes.setEnergy(newFitness);

    return cartes;
  }

  @Override
  public CartesianFullBackend getBackend() {
    return locopt.getBackend();
  }

  @Override
  public long getNumberOfGeomLocalOpts() {
    return locopt.getNumberOfGeomLocalOpts();
  }

  @Override
  public long getNumberOfMolLocalOpts() {
    return locopt.getNumberOfMolLocalOpts();
  }
}
