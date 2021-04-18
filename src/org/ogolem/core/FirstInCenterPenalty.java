/*
Copyright (c) 2010-2012, J. M. Dieterich
              2016-2020, J. M. Dieterich and B. Hartke
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
 * Checks if the first molecule is in the total center of mass.
 *
 * @author Johannes Dieterich
 * @version 2020-12-03
 */
class FirstInCenterPenalty implements AdditivePenaltyFunction {

  // the ID
  private static final long serialVersionUID = (long) 20101108;

  private final double thresh;

  private final double constant;

  private final double power;

  FirstInCenterPenalty(final double threshold, final double factor, final double pow) {
    thresh = threshold;
    constant = factor;
    power = pow;
  }

  @Override
  public AdditivePenaltyFunction copy() {
    return new FirstInCenterPenalty(this.thresh, this.constant, this.power);
  }

  @Override
  public double penalty(final Geometry geomStart) {

    // initialize working copy
    Geometry geom = new Geometry(geomStart);

    CartesianCoordinates cartes;
    if (geom.containsEnvironment()) {
      cartes = geom.getCartesiansWithEnvironment();
    } else {
      cartes = geom.getCartesians();
    }

    Environment refEnv = cartes.getReferenceEnvironmentCopy();
    ZMatrix[] zmatVec = cartes.getZMatrices();

    // first move all coordinates to total COM
    cartes.moveCoordsToCOM();

    // translate back
    cartes.setZMatrices(zmatVec);
    cartes.setRefEnvironment(refEnv);

    // we do not need to adjust the zmatVec since molecules are just moved

    geom =
        new Geometry(
            cartes,
            geom.getID(),
            geom.getNumberOfIndieParticles(),
            cartes.getAllAtomsPerMol(),
            geom.getAllFlexies(),
            geom.getExplicitDoFs(),
            geom.getAllConstraints(false),
            geom.getAllConstraintsXYZ(false),
            geom.getSIDs(),
            geom.getBondInfo());

    geom.setFitness(cartes.getEnergy());
    geom.setFather(geom.getFatherID());
    geom.setMother(geom.getMotherID());

    // check where the first COM is
    final double[] com = geom.getMoleculeAtPosition(0).getExternalCenterOfMass();

    // compute distance to 0/0/0
    final double dist = Math.sqrt(com[0] * com[0] + com[1] * com[1] + com[2] * com[2]);

    if (dist > thresh) {
      // compute penalty
      return (constant * Math.pow((dist - thresh), power));
    } else {
      // all good
      return 0.0;
    }
  }

  @Override
  public double penalty(CartesianCoordinates cartes) {

    // molecular cartesian check
    if (cartes.getNoOfMolecules() == 1) return 0.0;

    // working copy
    final CartesianCoordinates c = new CartesianCoordinates(cartes);

    // coords to COM
    c.moveCoordsToCOM();

    // check where the first COM is
    final double[] com = c.calculateCOMOfMol(0);

    // compute distance to 0/0/0
    final double dist = Math.sqrt(com[0] * com[0] + com[1] * com[1] + com[2] * com[2]);

    if (dist > thresh) {
      // compute penalty
      return (constant * Math.pow((dist - thresh), power));
    } else {
      // all good
      return 0.0;
    }
  }

  @Override
  public double penalty(Molecule mol) {
    return 0.0;
  }
}
