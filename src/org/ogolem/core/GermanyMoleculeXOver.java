/*
Copyright (c) 2014, J. M. Dieterich
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
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.RandomUtils;

/**
 * The most fundamental genotype operator in the molecule space.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class GermanyMoleculeXOver implements GenericCrossover<Double, Molecule> {

  private static final long serialVersionUID = (long) 20140401;
  private final boolean doOldStyle;
  private final double gaussWidth;

  GermanyMoleculeXOver(final boolean doOldStyle, final double gaussWidth) {
    this.doOldStyle = doOldStyle;
    this.gaussWidth = gaussWidth;
  }

  GermanyMoleculeXOver(final GermanyMoleculeXOver orig) {
    this.doOldStyle = orig.doOldStyle;
    this.gaussWidth = orig.gaussWidth;
  }

  @Override
  public GermanyMoleculeXOver copy() {
    return new GermanyMoleculeXOver(this);
  }

  @Override
  public String getMyID() {
    return "GERMANY\nmolecule genotype crossover";
  }

  @Override
  public Tuple<Molecule, Molecule> crossover(
      Molecule mother, Molecule father, final long futureID) {
    if (doOldStyle) {
      return crossoverOld(mother, father, futureID);
    } else {
      return crossoverNew(mother, father, futureID);
    }
  }

  private Tuple<Molecule, Molecule> crossoverNew(
      Molecule mother, Molecule father, final long futureID) {

    // the children molecules
    Molecule childOne;
    Molecule childTwo;

    final int atoms = mother.getNumberOfAtoms();

    // no of dofs
    final boolean[][] dofs = mother.getDegreesOfFreedom();
    int nodofs = 0;
    final ArrayList<int[]> dofpos = new ArrayList<>();
    for (int i = 0; i < dofs.length; i++) {
      for (int j = 0; j < dofs[0].length; j++) {
        if (dofs[i][j]) {
          dofpos.add(new int[] {i, j});
          nodofs++;
        }
      }
    }

    if (nodofs == 0 || nodofs == 1) {
      return crossoverOld(mother, father, futureID);
    }

    final double r = RandomUtils.gaussDouble(-1.0, 1.0, gaussWidth);

    /*
     * now we have a Gaussian distributes random number in between -1.0 and 1.0.
     * mapping of that to the actual cutting position takes place now
     */

    // the point where the cut through the vectors will take place
    int cut = 0;

    // the "size" of each bin
    double dBinSize = 2.0 / (double) nodofs;

    // the initial start is at -1.0
    double offset = -1.0 + dBinSize;

    // loop till we find the right spot to cut
    while (offset <= 1.0) {
      if (r < offset) {
        // the random number is in the bin
        break;
      } else {
        offset += dBinSize;
        cut++;
      }
    }

    // ok, we know now that we want to cut after a certain degree of freedom
    final int[] realcut = dofpos.get(cut);

    // the molecule configurations
    MoleculeConfig mcMother = mother.returnMyConfig();
    MoleculeConfig mcFather = father.returnMyConfig();

    MoleculeConfig mcChildOne = new MoleculeConfig(mcMother);
    MoleculeConfig mcChildTwo = new MoleculeConfig(mcFather);

    // zmatrices
    ZMatrix zmatMother = mother.getZMatrix();
    ZMatrix zmatFather = father.getZMatrix();
    ZMatrix zmatOne = mcChildOne.zmat;
    ZMatrix zmatTwo = mcChildTwo.zmat;

    // the actual chopping of zmat's
    for (int i = 0; i < realcut[0]; i++) {
      zmatOne.setABondLength(i, zmatMother.getABondLength(i));
      zmatTwo.setABondLength(i, zmatFather.getABondLength(i));

      zmatOne.setABondAngle(i, zmatMother.getABondAngle(i));
      zmatTwo.setABondAngle(i, zmatFather.getABondAngle(i));

      zmatOne.setADihedral(i, zmatMother.getADihedral(i));
      zmatTwo.setADihedral(i, zmatFather.getADihedral(i));
    }

    // the intermediate one
    if (realcut[1] == 0) {
      final int i = realcut[0];
      // after the bond length
      zmatOne.setABondLength(i, zmatMother.getABondLength(i));
      zmatTwo.setABondLength(i, zmatFather.getABondLength(i));

      zmatOne.setABondAngle(i, zmatFather.getABondAngle(i));
      zmatTwo.setABondAngle(i, zmatMother.getABondAngle(i));
      zmatOne.setADihedral(i, zmatFather.getADihedral(i));
      zmatTwo.setADihedral(i, zmatMother.getADihedral(i));
    } else if (realcut[1] == 1) {
      final int i = realcut[0];
      // after the bond angle
      zmatOne.setABondLength(i, zmatMother.getABondLength(i));
      zmatTwo.setABondLength(i, zmatFather.getABondLength(i));
      zmatOne.setABondAngle(i, zmatMother.getABondAngle(i));
      zmatTwo.setABondAngle(i, zmatFather.getABondAngle(i));

      zmatOne.setADihedral(i, zmatFather.getADihedral(i));
      zmatTwo.setADihedral(i, zmatMother.getADihedral(i));
    } else if (realcut[1] == 2) {
      final int i = realcut[0];
      // after the dihedral
      zmatOne.setABondLength(i, zmatMother.getABondLength(i));
      zmatTwo.setABondLength(i, zmatFather.getABondLength(i));
      zmatOne.setABondAngle(i, zmatMother.getABondAngle(i));
      zmatTwo.setABondAngle(i, zmatFather.getABondAngle(i));
      zmatOne.setADihedral(i, zmatMother.getADihedral(i));
      zmatTwo.setADihedral(i, zmatFather.getADihedral(i));
    }

    // behind the cut
    for (int i = realcut[0] + 1; i < atoms; i++) {
      zmatOne.setABondLength(i, zmatFather.getABondLength(i));
      zmatTwo.setABondLength(i, zmatMother.getABondLength(i));

      zmatOne.setABondAngle(i, zmatFather.getABondAngle(i));
      zmatTwo.setABondAngle(i, zmatMother.getABondAngle(i));

      zmatOne.setADihedral(i, zmatFather.getADihedral(i));
      zmatTwo.setADihedral(i, zmatMother.getADihedral(i));
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
    childOne = new Molecule(mcChildOne);
    childTwo = new Molecule(mcChildTwo);

    return new Tuple<>(childOne, childTwo);
  }

  private Tuple<Molecule, Molecule> crossoverOld(
      Molecule mother, Molecule father, final long futureID) {

    // the children molecules
    Molecule mChildOne;
    Molecule mChildTwo;

    int iNoOfAtoms = mother.getNumberOfAtoms();

    final double dRanGauss = RandomUtils.gaussDouble(-1.0, 1.0, gaussWidth);

    /*
     * now we have a Gaussian distributes random number in between -1.0 and 1.0.
     * mapping of that to the actual cutting position takes place now
     */

    // the point where the cut through the vectors will take place
    int iCut = 0;

    // the "size" of each bin
    double dBinSize = 2.0 / (double) iNoOfAtoms;
    // System.out.println("GermanyGlobOpt: the binsize:" + dBinSize);
    // the initial start is at -1.0
    double dOffset = -1.0 + dBinSize;

    // loop till we find the right spot to cut
    while (dOffset <= 1.0) {
      if (dRanGauss < dOffset) {
        // the random number is in the bin
        break;
      } else {
        dOffset += dBinSize;
        iCut++;
      }
    }

    // the molecule configurations
    MoleculeConfig mcMother = mother.returnMyConfig();
    MoleculeConfig mcFather = father.returnMyConfig();

    MoleculeConfig mcChildOne = new MoleculeConfig(mcMother);
    MoleculeConfig mcChildTwo = new MoleculeConfig(mcFather);

    // zmatrices
    ZMatrix zmatMother = mother.getZMatrix();
    ZMatrix zmatFather = father.getZMatrix();
    ZMatrix zmatOne = mcChildOne.zmat;
    ZMatrix zmatTwo = mcChildTwo.zmat;

    // the actual chopping of zmat's
    for (int i = 0; i < iCut; i++) {
      zmatOne.setABondLength(i, zmatMother.getABondLength(i));
      zmatTwo.setABondLength(i, zmatFather.getABondLength(i));

      zmatOne.setABondAngle(i, zmatMother.getABondAngle(i));
      zmatTwo.setABondAngle(i, zmatFather.getABondAngle(i));

      zmatOne.setADihedral(i, zmatMother.getADihedral(i));
      zmatTwo.setADihedral(i, zmatFather.getADihedral(i));
    }

    for (int i = iCut; i < iNoOfAtoms; i++) {
      zmatOne.setABondLength(i, zmatFather.getABondLength(i));
      zmatTwo.setABondLength(i, zmatMother.getABondLength(i));

      zmatOne.setABondAngle(i, zmatFather.getABondAngle(i));
      zmatTwo.setABondAngle(i, zmatMother.getABondAngle(i));

      zmatOne.setADihedral(i, zmatFather.getADihedral(i));
      zmatTwo.setADihedral(i, zmatMother.getADihedral(i));
    }

    // translate the new zmat to cartesian and create a refCartesian
    CartesianCoordinates refCartesOne =
        zmatOne.translateToCartesianAndAlign(mother.getCartesians());
    CartesianCoordinates refCartesTwo =
        zmatTwo.translateToCartesianAndAlign(father.getCartesians());

    // set the zmatrices and the reference cartesians
    mcChildOne.zmat = zmatOne;
    mcChildTwo.zmat = zmatTwo;

    mcChildOne.refXYZ = refCartesOne.getAllXYZCoord();
    mcChildTwo.refXYZ = refCartesTwo.getAllXYZCoord();

    // create the two children
    mChildOne = new Molecule(mcChildOne);
    mChildTwo = new Molecule(mcChildTwo);

    return new Tuple<>(mChildOne, mChildTwo);
  }

  @Override
  public short hasPriority() {
    return -1;
  }
}
