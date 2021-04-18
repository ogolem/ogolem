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

import static org.ogolem.core.Constants.*;

import java.io.Serializable;
import org.ogolem.generic.Copyable;

/**
 * An object representing a z-matrix.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class ZMatrix implements Serializable, Copyable {

  // the serial version UID
  private static final long serialVersionUID = (long) 20130921;

  private final int noOfAtoms;
  private String[] atomTypes;
  private double[] bondLengths;
  private int[] bondConnects;
  private double[] bondAngles;
  private int[] angleConnects;
  private double[] dihedrals;
  private int[] dihedralConnects;

  /**
   * The constructor. We do not want a generic constructor here.
   *
   * @param numAtoms The number of atoms in the molecule that the matrix is representing.
   */
  public ZMatrix(final int numAtoms) {
    noOfAtoms = numAtoms;
    atomTypes = new String[numAtoms];
    bondLengths = new double[numAtoms];
    bondConnects = new int[numAtoms];
    bondAngles = new double[numAtoms];
    angleConnects = new int[numAtoms];
    dihedrals = new double[numAtoms];
    dihedralConnects = new int[numAtoms];
  }

  /**
   * A copy constructor.
   *
   * @param original The original z-Matrix.
   */
  public ZMatrix(final ZMatrix original) {
    this.bondAngles = original.bondAngles.clone();
    this.bondLengths = original.bondLengths.clone();
    this.dihedrals = original.dihedrals.clone();
    this.noOfAtoms = original.noOfAtoms;
    this.angleConnects = original.angleConnects.clone();
    this.bondConnects = original.bondConnects.clone();
    this.dihedralConnects = original.dihedralConnects.clone();
    this.atomTypes = original.atomTypes.clone();
  }

  @Override
  public ZMatrix copy() {
    return new ZMatrix(this);
  }

  public int getNoOfAtoms() {
    return noOfAtoms;
  }

  public String[] getAllAtomNames() {
    return atomTypes;
  }

  public double[] getAllBondLengths() {
    return bondLengths;
  }

  public int[] getAllBondConnects() {
    return bondConnects;
  }

  public double[] getAllBondAngles() {
    return bondAngles;
  }

  public int[] getAllAnglesConnects() {
    return angleConnects;
  }

  public double[] getAllDihedrals() {
    return dihedrals;
  }

  public int[] getAllDihedralConnects() {
    return dihedralConnects;
  }

  double getABondLength(final int i) {
    return bondLengths[i];
  }

  double getABondAngle(final int i) {
    return bondAngles[i];
  }

  double getADihedral(final int i) {
    return dihedrals[i];
  }

  public void setAllAtomNames(final String[] atomNames) {
    atomTypes = atomNames;
  }

  public void setAllBondLengths(final double[] bonds) {
    bondLengths = bonds;
  }

  public void setAllBondConnects(final int[] connects) {
    bondConnects = connects;
  }

  public void setAllBondAngles(final double[] angs) {
    bondAngles = angs;
  }

  public void setAllAnglesConnects(final int[] connects) {
    angleConnects = connects;
  }

  public void setAllDihedrals(final double[] dih) {
    dihedrals = dih;
  }

  public void setAllDihedralConnects(final int[] connects) {
    dihedralConnects = connects;
  }

  void setAZMatrixLine(
      int iWhichAtom,
      String sAtom,
      double dBondLength,
      int iBondConnect,
      double dBondAngle,
      int iAngleConnect,
      double dDihedral,
      int iDihedralConnect) {
    atomTypes[iWhichAtom] = sAtom;
    bondLengths[iWhichAtom] = dBondLength;
    bondConnects[iWhichAtom] = iBondConnect;
    bondAngles[iWhichAtom] = dBondAngle;
    angleConnects[iWhichAtom] = iAngleConnect;
    dihedrals[iWhichAtom] = dDihedral;
    dihedralConnects[iWhichAtom] = iDihedralConnect;
  }

  void setABondLength(final int whichAtom, final double bondLength) {
    assert (!Double.isNaN(bondLength) && !Double.isInfinite(bondLength));
    assert (bondLength >= 0.0);
    bondLengths[whichAtom] = bondLength;
  }

  void setABondAngle(final int whichAtom, final double bondAngle) {
    assert (!Double.isNaN(bondAngle) && !Double.isInfinite(bondAngle));
    if (bondAngle > Math.PI) {
      bondAngles[whichAtom] = Math.PI;
    } else if (bondAngle < 0.0) {
      bondAngles[whichAtom] = 0.0;
    } else {
      bondAngles[whichAtom] = bondAngle;
    }
  }

  void setADihedral(final int whichAtom, final double dihedral) {
    assert (!Double.isNaN(dihedral) && !Double.isInfinite(dihedral));
    if (dihedral > Math.PI) {
      dihedrals[whichAtom] = Math.PI;
    } else if (dihedral < -Math.PI) {
      dihedrals[whichAtom] = -Math.PI;
    } else {
      dihedrals[whichAtom] = dihedral;
    }
  }

  /**
   * The caller is responsible for sane values!
   *
   * @param kind
   * @param atom
   * @param value
   */
  void setAValue(final int kind, final int atom, final double value) {
    if (kind == 0) bondLengths[atom] = value;
    else if (kind == 1) bondAngles[atom] = value;
    else if (kind == 2) dihedrals[atom] = value;
  }

  String[] getPrintableZMatrix() {
    final String[] saOutput = new String[noOfAtoms];
    saOutput[0] = atomTypes[0];
    saOutput[1] = atomTypes[1] + "\t" + bondLengths[1] * BOHRTOANG + "\t" + bondConnects[1];
    saOutput[2] =
        atomTypes[2]
            + "\t"
            + bondLengths[2] * BOHRTOANG
            + "\t"
            + bondConnects[2]
            + "\t"
            + bondAngles[2] * 180 / Math.PI
            + "\t"
            + angleConnects[2];
    for (int i = 3; i < noOfAtoms; i++) {
      saOutput[i] =
          atomTypes[i]
              + "\t"
              + bondLengths[i] * BOHRTOANG
              + "\t"
              + bondConnects[i]
              + "\t"
              + bondAngles[i] * 180 / Math.PI
              + "\t"
              + angleConnects[i]
              + "\t"
              + dihedrals[i] * 180 / Math.PI
              + "\t"
              + dihedralConnects[i];
    }
    return saOutput;
  }

  double getDOF(final int atom, final int dof) {
    if (dof == 0) return bondLengths[atom];
    else if (dof == 1) return bondAngles[atom];
    else if (dof == 2) return dihedrals[atom];
    System.err.println("WARNING: No such kind " + dof + " in getDOF()! Contact author(s)!");
    return 0.0;
  }

  public void setDOF(final int atom, final int dof, final double val) {
    if (dof == 0) {
      // manipulate the bond length
      bondLengths[atom] = val;
    } else if (dof == 1) {
      // manipulate the bond angle
      bondAngles[atom] = val;
      if (bondAngles[atom] > Math.PI) {
        bondAngles[atom] = Math.PI;
      } else if (bondAngles[atom] < 0.0) {
        bondAngles[atom] = 0.0;
      }
    } else if (dof == 2) {
      // manipulate the dihedral
      dihedrals[atom] = val;
      if (dihedrals[atom] > Math.PI) {
        dihedrals[atom] = Math.PI;
      } else if (dihedrals[atom] < -Math.PI) {
        dihedrals[atom] = -Math.PI;
      }
    }
  }

  /**
   * Translate this z-Matrix to Cartesian coordinates and optionally align it.
   *
   * @param refCartes if null: no aligning will take place. Otherwise, this is the reference.
   * @return the translated cartesian coordinates
   */
  public CartesianCoordinates translateToCartesianAndAlign(final CartesianCoordinates refCartes) {
    CartesianCoordinates cartes = CoordTranslation.zMatToCartesians(this);
    if (refCartes != null) {
      // aligning
      cartes = CoordTranslation.alignTwoCartes(refCartes, cartes);
    }
    return cartes;
  }

  /**
   * Translate this z-Matrix to Cartesian coordinates and optionally align it.
   *
   * @param refCartes if null: no aligning will take place. Otherwise, this is the reference.
   * @return the translated cartesian coordinates
   */
  public CartesianCoordinates translateToCartesianAndAlign(final double[][] refCartes) {
    CartesianCoordinates cartes = CoordTranslation.zMatToCartesians(this);
    if (refCartes != null) {
      // aligning
      try {
        cartes = CoordTranslation.alignTwoCartesians(refCartes, cartes).getObject1();
      } catch (Exception e) {
        System.err.println(
            "WARNING: Failed to align molecules. "
                + e.toString()
                + "\nReturning non-aligned cartes!");
      }
    }
    return cartes;
  }
}
