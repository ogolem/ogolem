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
 * This uses an external engine (which is specified through the constructor) to do a local
 * optimization of a molecule or a geometry.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class ExternalLocOpt implements Newton {

  // the ID
  private static final long serialVersionUID = (long) 20101108;

  private final FFEngineWrapper ffeng;

  ExternalLocOpt(GlobalConfig globConf, int ffID, double cutE) {
    ffeng = new FFEngineWrapper(globConf, ffID, cutE);
  }

  private ExternalLocOpt(final ExternalLocOpt orig) {
    ffeng = orig.ffeng.copy();
  }

  @Override
  public ExternalLocOpt copy() {
    return new ExternalLocOpt(this);
  }

  @Override
  public String myIDandMethod() {
    return "external LocOpt";
  }

  @Override
  public Molecule localOptimization(Molecule mStartMolecule) {
    return ffeng.localOptimization(mStartMolecule);
  }

  @Override
  public Geometry localOptimization(Geometry gStartGeometry) {
    return ffeng.localOptimization(gStartGeometry);
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      long id,
      CartesianCoordinates cartes,
      boolean[][] constraints,
      boolean isConstricted,
      final BondInfo bonds)
      throws Exception {
    return ffeng.cartesToCartes(id, cartes, constraints, isConstricted, bonds);
  }

  @Override
  public CartesianFullBackend getBackend() {
    return ffeng.getBackend();
  }

  @Override
  public long getNumberOfGeomLocalOpts() {
    return ffeng.getNumberOfGeomLocalOpts();
  }

  @Override
  public long getNumberOfMolLocalOpts() {
    return ffeng.getNumberOfMolLocalOpts();
  }
}
