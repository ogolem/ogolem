/*
Copyright (c) 2013, J. M. Dieterich
              2016-2023, J. M. Dieterich and B. Hartke
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

import java.util.Arrays;
import org.ogolem.math.BoolSymmetricMatrixNoDiag;
import org.ogolem.math.ShortSymmetricMatrixNoDiag;

/**
 * A simple, double matrix-backed implementation of the BondInfo interface.
 *
 * @author Johannes Dieterich
 * @version 2023-09-08
 */
public class SimpleBondInfo extends AbstractBondInfo {

  private static final long serialVersionUID = (long) 20201228;
  private final BoolSymmetricMatrixNoDiag bonds;
  private final ShortSymmetricMatrixNoDiag bondTypes;

  public SimpleBondInfo(final int noAtoms) {
    super(noAtoms);
    this.bonds = new BoolSymmetricMatrixNoDiag(noAtoms);
    final boolean[] buff = bonds.underlyingStorageBuffer();
    Arrays.fill(buff, false);
    this.bondTypes = new ShortSymmetricMatrixNoDiag(noAtoms);
    final var buffT = bondTypes.underlyingStorageBuffer();
    Arrays.fill(buffT, BondInfo.NOBOND);
  }

  public SimpleBondInfo(final BoolSymmetricMatrixNoDiag bondMat) {
    super(bondMat.noRows());
    this.bonds = new BoolSymmetricMatrixNoDiag(noAtoms);
    this.bondTypes = new ShortSymmetricMatrixNoDiag(noAtoms);

    final boolean[] buffBool = bonds.underlyingStorageBuffer();
    final short[] buffShort = bondTypes.underlyingStorageBuffer();

    assert (buffBool.length == buffShort.length);

    for (int i = 0; i < buffBool.length; i++) {
      if (buffBool[i]) buffShort[i] = BondInfo.UNCERTAIN;
    }
  }

  private SimpleBondInfo(final SimpleBondInfo orig) {
    super(orig);
    this.bondTypes = orig.bondTypes.copy();
    this.bonds = orig.bonds.copy();
  }

  @Override
  public SimpleBondInfo copy() {
    return new SimpleBondInfo(this);
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
    assert (atom1 != atom2);
    return this.bonds.getElement(atom1, atom2);
  }

  @Override
  public short bondType(final int atom1, final int atom2) {
    assert (atom1 != atom2);
    return this.bondTypes.getElement(atom1, atom2);
  }

  @Override
  public void setBond(final int atom1, final int atom2, final short bondType) {

    assert (atom1 != atom2);
    assert (atom1 < noAtoms);
    assert (atom2 < noAtoms);

    this.bondTypes.setElement(atom1, atom2, bondType);

    final boolean hasBond = (bondType != BondInfo.NOBOND);
    this.bonds.setElement(atom1, atom2, hasBond);
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
  public BoolSymmetricMatrixNoDiag getFullBondMatrix() {
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
  public ShortSymmetricMatrixNoDiag getBondInformation() {
    return this.bondTypes;
  }
}
