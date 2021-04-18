/*
Copyright (c) 2014, J. M. Dieterich
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

import java.util.ArrayList;
import java.util.List;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.RandomUtils;

/**
 * n-Point genotype crossover for molecules.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class PortugalMoleculeXOver implements GenericCrossover<Double, Molecule> {

  private static final long serialVersionUID = (long) 20140327;
  private final int noCuts;

  PortugalMoleculeXOver(final int noCuts) {
    this.noCuts = noCuts;
  }

  PortugalMoleculeXOver(final PortugalMoleculeXOver orig) {
    this.noCuts = orig.noCuts;
  }

  @Override
  public PortugalMoleculeXOver copy() {
    return new PortugalMoleculeXOver(this);
  }

  @Override
  public String getMyID() {
    return "PORTUGAL\nn-point geometry molecule crossover\n\t #crosses: " + noCuts;
  }

  @Override
  public Tuple<Molecule, Molecule> crossover(
      Molecule mother, Molecule father, final long futureID) {

    final int atoms = mother.getNumberOfAtoms();

    // no of dofs
    final boolean[][] dofs = mother.getDegreesOfFreedom();
    int nodofs = 0;
    final List<int[]> dofpos = new ArrayList<>();
    for (int i = 0; i < dofs.length; i++) {
      for (int j = 0; j < dofs[0].length; j++) {
        if (dofs[i][j]) {
          dofpos.add(new int[] {i, j});
          nodofs++;
        }
      }
    }

    if (noCuts >= nodofs) {
      throw new RuntimeException(
          "WARNING: Too many cuts (" + noCuts + ") for DoFs (" + nodofs + ").");
    }

    final List<Integer> cutPoints = RandomUtils.listOfPoints(noCuts, nodofs);

    // the molecule configurations
    final MoleculeConfig mcMother = mother.returnMyConfig();
    final MoleculeConfig mcFather = father.returnMyConfig();

    final MoleculeConfig mcChildOne = new MoleculeConfig(mcMother);
    final MoleculeConfig mcChildTwo = new MoleculeConfig(mcFather);

    // zmatrices
    final ZMatrix zmatMother = mother.getZMatrix();
    final ZMatrix zmatFather = father.getZMatrix();
    final ZMatrix zmatOne = mcChildOne.zmat;
    final ZMatrix zmatTwo = mcChildTwo.zmat;

    boolean mothFirst = true;
    for (int i = 0; i < dofpos.size(); i++) {

      if (!cutPoints.contains(i)) {
        continue;
      }

      final int[] pos = dofpos.get(i);
      final int atom = pos[0];
      final int dof = pos[1];

      final int[] nextpos = (i < dofpos.size() - 1) ? dofpos.get(i + 1) : null;
      final int nextAtom = (nextpos == null) ? atoms : nextpos[0];
      final int nextDOF = (nextpos == null) ? 3 : nextpos[1];

      // we can just ignore the fact that the zmat contains different amounts of dofs per atom
      AtomLoop:
      for (int at = atom; at < nextAtom; at++) {
        int start = 0;
        int end = 3;
        if (at == atom) {
          start = dof;
        }
        if (at == nextAtom) {
          end = nextDOF;
        }

        DOFLoop:
        for (int df = start; df < end; df++) {

          if (mothFirst) {
            zmatOne.setDOF(at, dof, zmatMother.getDOF(at, dof));
            zmatTwo.setDOF(at, dof, zmatFather.getDOF(at, dof));
          } else {
            zmatTwo.setDOF(at, dof, zmatMother.getDOF(at, dof));
            zmatOne.setDOF(at, dof, zmatFather.getDOF(at, dof));
          }

          if (at == nextAtom - 1 && df == end - 1) {
            break AtomLoop;
          }
        }
      }

      mothFirst = !mothFirst;
    }

    // translate the new zmat to cartesian and create a refCartesian
    final CartesianCoordinates refCartesOne =
        zmatOne.translateToCartesianAndAlign(mother.getCartesians());
    final CartesianCoordinates refCartesTwo =
        zmatTwo.translateToCartesianAndAlign(father.getCartesians());

    // set the zmatrices and the reference cartesians
    mcChildOne.zmat = zmatOne;
    mcChildTwo.zmat = zmatTwo;

    mcChildOne.refXYZ = refCartesOne.getAllXYZCoord();
    mcChildTwo.refXYZ = refCartesTwo.getAllXYZCoord();

    // create the two children
    Molecule childOne = new Molecule(mcChildOne);
    Molecule childTwo = new Molecule(mcChildTwo);

    // TODO test and debug
    return new Tuple<>(childOne, childTwo);
  }

  @Override
  public short hasPriority() {
    return -1;
  }
}
