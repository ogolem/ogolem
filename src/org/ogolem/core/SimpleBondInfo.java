/*
Copyright (c) 2013, J. M. Dieterich
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
 * A simple, double matrix-backed implementation of the BondInfo interface.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class SimpleBondInfo extends AbstractBondInfo {

  private static final long serialVersionUID = (long) 20160224;
  private final boolean[][] bonds;
  private final short[][] bondTypes;

  public SimpleBondInfo(final int noAtoms) {
    super(noAtoms);
    this.bonds = new boolean[noAtoms][noAtoms];
    for (int i = 0; i < noAtoms; i++) {
      for (int j = 0; j < noAtoms; j++) {
        bonds[i][j] = false;
      }
    }
    this.bondTypes = new short[noAtoms][noAtoms];
  }

  public SimpleBondInfo(final boolean[][] bondMat) {
    super(bondMat.length);
    this.bonds = new boolean[noAtoms][noAtoms];
    this.bondTypes = new short[noAtoms][noAtoms];
    for (int i = 0; i < noAtoms; i++) {
      for (int j = 0; j < noAtoms; j++) {
        bonds[i][j] = bondMat[i][j];
        if (bonds[i][j]) {
          bondTypes[i][j] = BondInfo.UNCERTAIN;
        }
      }
    }
  }

  private SimpleBondInfo(final SimpleBondInfo orig) {
    super(orig);
    this.bondTypes = orig.bondTypes.clone();
    this.bonds = orig.bonds.clone();
  }

  @Override
  public SimpleBondInfo copy() {
    final SimpleBondInfo cop = new SimpleBondInfo(noAtoms);
    for (int i = 0; i < noAtoms; i++) {
      System.arraycopy(this.bondTypes[i], 0, cop.bondTypes[i], 0, noAtoms);
      System.arraycopy(this.bonds[i], 0, cop.bonds[i], 0, noAtoms);
    }
    return cop;
  }

  @Override
  public SimpleBondInfo shallowCopy() {
    return new SimpleBondInfo(this);
  }

  @Override
  public int getNoOfAtoms() {
    return noAtoms;
  }

  @Override
  public boolean hasBond(final int atom1, final int atom2) {
    return this.bonds[atom1][atom2];
  }

  @Override
  public short bondType(final int atom1, final int atom2) {
    return this.bondTypes[atom1][atom2];
  }

  @Override
  public void setBond(final int atom1, final int atom2, final short bondType) {
    this.bondTypes[atom1][atom2] = bondType;
    this.bondTypes[atom2][atom1] = bondType;
    if (bondType != BondInfo.NOBOND) {
      this.bonds[atom1][atom2] = true;
      this.bonds[atom2][atom1] = true;
    }
  }

  @Override
  public boolean bondMatrixFast() {
    return true;
  }

  /**
   * Due to its double-matrix backing, this operation is fast.
   *
   * @return the full bonds
   */
  @Override
  public boolean[][] getFullBondMatrix() {
    return this.bonds;
  }

  @Override
  public boolean bondInfoFast() {
    return true;
  }

  /**
   * Due to its double-matrix backing, this operation is fast.
   *
   * @return the full bond information
   */
  @Override
  public short[][] getBondInformation() {
    return this.bondTypes;
  }
}
