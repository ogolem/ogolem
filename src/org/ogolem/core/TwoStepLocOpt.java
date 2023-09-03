/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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
 * Doing not one but two locopts after each other. This is helpful if the convergence pattern of a
 * given system is bad. A first raw locopt using, e.g., an UFF can make the convergence behaviour of
 * the second locopt to be way better.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class TwoStepLocOpt implements Newton {

  // the ID
  private static final long serialVersionUID = (long) 20101108;

  private final Newton primLocOpt;
  private final Newton secLocOpt;

  TwoStepLocOpt(Newton locOpt1, Newton locOpt2) {
    this.primLocOpt = locOpt1;
    this.secLocOpt = locOpt2;
  }

  private TwoStepLocOpt(final TwoStepLocOpt orig) {
    this.primLocOpt = orig.primLocOpt.copy();
    this.secLocOpt = orig.secLocOpt.copy();
  }

  @Override
  public TwoStepLocOpt copy() {
    return new TwoStepLocOpt(this);
  }

  @Override
  public String myIDandMethod() {
    return "Two coupled locopts "
        + primLocOpt.myIDandMethod()
        + " and "
        + secLocOpt.myIDandMethod();
  }

  @Override
  public Geometry localOptimization(final Geometry start) {

    Geometry back = new Geometry(start);
    final Geometry post1 = primLocOpt.localOptimization(start);

    if (post1.getFitness() >= FixedValues.NONCONVERGEDENERGY) {
      back.setFitness(FixedValues.NONCONVERGEDENERGY);
      return back;
    }

    back = new Geometry(post1);
    final Geometry post2 = secLocOpt.localOptimization(post1);

    if (post2.getFitness() >= FixedValues.NONCONVERGEDENERGY) {
      back.setFitness(FixedValues.NONCONVERGEDENERGY);
      return back;
    }

    return post2;
  }

  @Override
  public Molecule localOptimization(final Molecule start) {

    Molecule back = new Molecule(start);
    final Molecule post1 = primLocOpt.localOptimization(start);

    if (post1.getEnergy() >= FixedValues.NONCONVERGEDENERGY) {
      back.setEnergy(FixedValues.NONCONVERGEDENERGY);
      return back;
    }

    back = new Molecule(post1);
    final Molecule post2 = secLocOpt.localOptimization(post1);

    if (post2.getEnergy() >= FixedValues.NONCONVERGEDENERGY) {
      back.setEnergy(FixedValues.NONCONVERGEDENERGY);
      return back;
    }

    return post2;
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      long id,
      CartesianCoordinates cartes,
      boolean[][] constraints,
      boolean isConstricted,
      final BondInfo bonds)
      throws Exception {

    // might throw expections!
    cartes = primLocOpt.cartesToCartes(id, cartes, constraints, isConstricted, bonds);
    cartes = secLocOpt.cartesToCartes(id, cartes, constraints, isConstricted, bonds);

    return cartes;
  }

  /**
   * @return always returns the secondary backend (if available)
   */
  @Override
  public CartesianFullBackend getBackend() {
    return secLocOpt.getBackend();
  }

  @Override
  public long getNumberOfGeomLocalOpts() {
    System.out.println(
        "INFO: TwoStep reports "
            + primLocOpt.getNumberOfGeomLocalOpts()
            + " for locopt "
            + primLocOpt.myIDandMethod()
            + " and "
            + secLocOpt.getNumberOfGeomLocalOpts()
            + " for locopt "
            + secLocOpt.myIDandMethod());
    return primLocOpt.getNumberOfGeomLocalOpts();
  }

  @Override
  public long getNumberOfMolLocalOpts() {
    System.out.println(
        "INFO: TwoStep reports "
            + primLocOpt.getNumberOfMolLocalOpts()
            + " for locopt "
            + primLocOpt.myIDandMethod()
            + " and "
            + secLocOpt.getNumberOfMolLocalOpts()
            + " for locopt "
            + secLocOpt.myIDandMethod());
    return primLocOpt.getNumberOfMolLocalOpts();
  }
}
