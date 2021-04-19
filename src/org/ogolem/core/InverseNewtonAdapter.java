/*
Copyright (c) 2015, J. M. Dieterich
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

import java.util.ArrayList;
import java.util.List;
import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.GenericLocOpt;

/**
 * Translates a generic fitness function based local optimization to the Newton interface.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class InverseNewtonAdapter implements Newton {

  private static final long serialVersionUID = (long) 20200429;

  private final GenericLocOpt<Molecule, Geometry> locopt;
  private long countGeomOpts = 0l;
  private long countMolOpts = 0l;

  public InverseNewtonAdapter(final GenericLocOpt<Molecule, Geometry> locopt) {
    this.locopt = locopt;
  }

  InverseNewtonAdapter(final InverseNewtonAdapter orig) {
    this.locopt = orig.locopt.copy();
  }

  @Override
  public Newton copy() {
    return new InverseNewtonAdapter(this);
  }

  @Override
  public String myIDandMethod() {
    return locopt.getMyID();
  }

  @Override
  public CartesianFullBackend getBackend() {

    final GenericBackend<Molecule, Geometry> coordinates = locopt.getBackend();
    if (coordinates instanceof FullyCartesianCoordinates) {
      final FullyCartesianCoordinates fully = (FullyCartesianCoordinates) coordinates;
      final CartesianFullBackend backend = fully.getMyBackend();
      return backend;
    } else if (coordinates instanceof RigidBodyCoordinates) {
      final RigidBodyCoordinates rigid = (RigidBodyCoordinates) coordinates;
      final RigidBodyBackend backend = rigid.getMyBackend();
      // TODO currently not implemented...
      throw new UnsupportedOperationException("Not supported yet.");
    }

    return null;
  }

  @Override
  public Molecule localOptimization(Molecule mStartMolecule) {
    countMolOpts++;
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public Geometry localOptimization(final Geometry gStartGeometry) {
    countGeomOpts++;

    return locopt.fitness(gStartGeometry, false);
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      final long id,
      final CartesianCoordinates cartes,
      final boolean[][] constraints,
      final boolean isConstricted,
      final BondInfo bonds)
      throws Exception {

    // translate to a geometry
    final boolean hasEnv = cartes.containsEnvironment();
    final int noMols = cartes.getNoOfMolecules();
    final String[] sids = new String[noMols];
    final boolean[] molFlexies = new boolean[noMols];
    final boolean[] molConstraints = new boolean[noMols];
    final int[] atsPerMol = cartes.getAllAtomsPerMol();
    final List<boolean[][]> allFlexies = new ArrayList<>();
    int offset = 0;
    for (int mol = 0; mol < noMols; mol++) {
      final CartesianCoordinates molCartes = cartes.giveMolecularCartes(mol, false);
      sids[mol] = molCartes.createDefaultMolecularType();
      molFlexies[mol] = (molCartes.getMolecularRefZMatrix(mol) != null);
      if (molFlexies[mol]) {
        // actually it does not matter what we fill in here. we will translate it to cartesians
        // again. so anything does
        final int noAtoms = molCartes.getNoOfAtoms();
        final boolean[][] myFlexies = new boolean[noAtoms][3];
        allFlexies.add(myFlexies);
      } else {
        allFlexies.add(null);
      }
      boolean hasConstr = false;
      for (int i = 0; i < 3; i++) {
        for (int at = offset; at < offset + atsPerMol[mol]; at++) {
          if (constraints[i][at]) {
            hasConstr = false;
          }
        }
      }
      molConstraints[mol] = hasConstr;
      offset += atsPerMol[mol];
    }

    final Geometry work =
        CoordTranslation.cartesianToGeometry(
            cartes,
            noMols,
            atsPerMol,
            molFlexies,
            allFlexies,
            molConstraints,
            constraints,
            sids,
            bonds);

    final Geometry opt = locopt.fitness(work, false);

    return (hasEnv) ? opt.getCartesiansWithEnvironment() : opt.getCartesians();
  }

  @Override
  public long getNumberOfGeomLocalOpts() {
    return countGeomOpts;
  }

  @Override
  public long getNumberOfMolLocalOpts() {
    return countMolOpts;
  }
}
