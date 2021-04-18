/*
Copyright (c) 2013, J. M. Dieterich
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

import java.io.Serializable;
import java.util.List;
import org.ogolem.generic.Copyable;

/**
 * Interface describing what bond information implementations need to be capable of.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public interface BondInfo extends Serializable, Copyable {

  public static final short NOBOND = 0;
  public static final short SINGLE = 1;
  public static final short DOUBLE = 2;
  public static final short TRIPLE = 3;
  public static final short AROMATIC = 4;
  public static final short VDW = 5;
  public static final short UNCERTAIN = 99;

  /**
   * Clones this BondInfo object.
   *
   * @return a clone.
   */
  @Override
  public BondInfo copy();

  /**
   * Returns a shallow copy of this BondInfo object, reducing the needed memory requirements but
   * increasing the risk for something to go BOOM!
   *
   * @return a shallow copy
   */
  public BondInfo shallowCopy();

  /**
   * Get the number of atoms this bonding object was made for.
   *
   * @return number of atoms
   */
  public int getNoOfAtoms();

  /**
   * Whether or not the two atoms are bonded.
   *
   * @param atom1 Atom 1 in the structure.
   * @param atom2 Atom 2 in the structure.
   * @return true if there is a bond, false otherwise
   */
  public boolean hasBond(final int atom1, final int atom2);

  /**
   * Returns the bond type of the bond. Following codes are supported: * 0: no bond * 1: single bond
   * * 2: double bond * 3: triple bond * 4: aromatic bond * 5: vdW bond * 99: unspecified bond
   *
   * @param atom1 Atom 1 in the structure
   * @param atom2 Atom 2 in the structure
   * @return the above bond type
   */
  public short bondType(final int atom1, final int atom2);

  /**
   * Set a bond for two atoms.
   *
   * @param atom1 Atom 1 in the structure.
   * @param atom2 Atom 2 in the structure.
   * @param bondType one of the bond types above.
   */
  public void setBond(final int atom1, final int atom2, final short bondType);

  /**
   * Whether or not the return of the full bond matrix is fast or not.
   *
   * @return true if fast, false otherwise.
   */
  public boolean bondMatrixFast();

  /**
   * Get the full bond matrix. May be VERY SLOW depending on the implementation.
   *
   * @return the bond matrix (boolean: bond yes/no)
   */
  public boolean[][] getFullBondMatrix();

  /**
   * Whether or not the return of the full bond information is fast or not.
   *
   * @return true if fast, false otherwise.
   */
  public boolean bondInfoFast();

  /**
   * Get the full bond information. May be VERY SLOW depending on the implementation.
   *
   * @return the bond information as the above described shorts.
   */
  public short[][] getBondInformation();

  /**
   * Translates this bond information to the input formatted one.
   *
   * @return a list of strings ready for .ogo input.
   */
  public List<String> translateToInput();
}
