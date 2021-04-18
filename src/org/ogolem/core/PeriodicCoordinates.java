/*
Copyright (c) 2016, J. M. Dieterich
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

      This product includes software developed at the Universities of
      Kiel, Goettingen (Germany) and Princeton University (USA) by its
      contributors: J. M. Dieterich and B. Hartke.

    * Neither the name of the University of Kiel, the University of Goettingen,
      Princeton University nor the names of its contributors may be used to
      endorse or promote products derived from this software without specific
      prior written permission.

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

/**
 * An object to hold periodic coordinates (fractional). As usual, everything is, where applicable
 * defined in a.u. (except weights).
 *
 * @author Johannes M Dieterich
 * @version 2020-12-30
 */
public class PeriodicCoordinates implements StructuralData {

  private static final long serialVersionUID = (long) 20160721;

  private final int noAtoms;
  private final int noMolecules;
  private final int[] atomsPerMol;
  private final String[] atomNames;
  private final short[] atomNos;
  private final double[][] cellVectors;
  private final double[][] fractionalCoords;
  private final List<ZMatrix> zmats;
  private final short[] spins;
  private final float[] charges;

  public PeriodicCoordinates(
      final int noAtoms,
      final int noMolecules,
      final int[] atomsPerMol,
      final String[] atomNames,
      final double[][] cellVectors,
      final double[][] fractionalCoords,
      final List<ZMatrix> zmats) {
    assert (noAtoms > 0);
    assert (noMolecules > 0);
    assert (atomsPerMol != null);
    assert (atomsPerMol.length == noAtoms);
    assert (zmats != null);
    assert (zmats.size() == noMolecules);
    assert (atomNames != null);
    assert (atomNames.length == noAtoms);
    assert (cellVectors != null);
    assert (cellVectors.length == 3);
    assert (cellVectors[0].length == 3);
    assert (cellVectors[1].length == 3);
    assert (cellVectors[2].length == 3);
    assert (fractionalCoords != null);
    assert (fractionalCoords.length == 3);
    assert (fractionalCoords[0].length == noAtoms);
    assert (fractionalCoords[1].length == noAtoms);
    assert (fractionalCoords[2].length == noAtoms);

    this.noAtoms = noAtoms;
    this.noMolecules = noMolecules;
    this.atomsPerMol = atomsPerMol;
    this.zmats = zmats;
    this.atomNames = atomNames;
    this.atomNos = new short[noAtoms];
    for (int i = 0; i < noAtoms; i++) {
      this.atomNos[i] = AtomicProperties.giveAtomicNumber(atomNames[i]);
    }
    this.cellVectors = cellVectors;
    this.fractionalCoords = fractionalCoords;
    this.spins = new short[noAtoms];
    this.charges = new float[noAtoms];
  }

  private PeriodicCoordinates(final PeriodicCoordinates orig) {
    assert (orig.noMolecules > 0);
    assert (orig.noAtoms > 0);
    assert (orig.atomsPerMol.length == orig.noMolecules);
    assert (orig.cellVectors != null);
    assert (orig.cellVectors.length == 3);
    assert (orig.cellVectors[0].length == 3);

    this.noAtoms = orig.noAtoms;
    this.noMolecules = orig.noMolecules;

    this.cellVectors = new double[3][3];
    cellVectors[0][0] = orig.cellVectors[0][0];
    cellVectors[0][1] = orig.cellVectors[0][1];
    cellVectors[0][2] = orig.cellVectors[0][2];
    cellVectors[1][0] = orig.cellVectors[1][0];
    cellVectors[1][1] = orig.cellVectors[1][1];
    cellVectors[1][2] = orig.cellVectors[1][2];
    cellVectors[2][0] = orig.cellVectors[2][0];
    cellVectors[2][1] = orig.cellVectors[2][1];
    cellVectors[2][2] = orig.cellVectors[2][2];

    assert (orig.fractionalCoords != null);
    assert (orig.fractionalCoords.length == 3);
    assert (orig.fractionalCoords[0].length == noAtoms);
    assert (orig.fractionalCoords[1].length == noAtoms);
    assert (orig.fractionalCoords[2].length == noAtoms);
    this.fractionalCoords = new double[3][];
    fractionalCoords[0] = orig.fractionalCoords[0].clone();
    fractionalCoords[1] = orig.fractionalCoords[1].clone();
    fractionalCoords[2] = orig.fractionalCoords[2].clone();

    assert (orig.atomNos != null);
    assert (orig.atomNos.length == noAtoms);
    this.atomNos = orig.atomNos.clone();

    assert (orig.atomNames != null);
    assert (orig.atomNames.length == noAtoms);
    this.atomNames = orig.atomNames.clone();

    assert (orig.atomsPerMol != null);
    assert (orig.atomsPerMol.length == noMolecules);
    this.atomsPerMol = orig.atomsPerMol.clone();

    this.zmats = new ArrayList<>();
    for (int i = 0; i < orig.zmats.size(); i++) {
      final ZMatrix zmat = orig.zmats.get(i);
      if (zmat != null) {
        this.zmats.add(zmat.copy());
      } else {
        zmats.add(null);
      }
    }
    this.spins = orig.spins.clone();
    this.charges = orig.charges.clone();
  }

  @Override
  public PeriodicCoordinates copy() {
    return new PeriodicCoordinates(this);
  }

  @Override
  public CartesianCoordinates getCartesianCoordinates() {

    final CartesianCoordinates cartes = null; // new CartesianCoordinates();
    System.err.println("PLacehodler for getCartesianCoordinates for periodic, IMPLEMENT ME!!!!");
    // XXX important
    return cartes;
  }

  @Override
  public int getNoOfAtoms() {
    assert (noAtoms > 0);
    return noAtoms;
  }

  @Override
  public int getNoOfMolecules() {
    assert (noMolecules > 0);
    return noMolecules;
  }

  public int[] getAtomsPerMolecule() {
    assert (atomsPerMol != null);
    assert (atomsPerMol.length == noMolecules);
    return atomsPerMol;
  }

  public String[] getAtomNames() {
    assert (atomNames != null);
    assert (atomNames.length == noAtoms);
    return atomNames;
  }

  @Override
  public short[] getAllAtomNumbers() {
    assert (atomNos != null);
    assert (atomNos.length == noAtoms);
    return atomNos;
  }

  public double[][] getCellVectors() {
    assert (cellVectors != null);
    assert (cellVectors.length == 3);
    assert (cellVectors[0].length == 3);
    assert (cellVectors[1].length == 3);
    assert (cellVectors[2].length == 3);
    return cellVectors;
  }

  public double[][] getFractionalCoordinates() {
    assert (fractionalCoords != null);
    assert (fractionalCoords.length == 3);
    assert (fractionalCoords[0].length == noAtoms);
    assert (fractionalCoords[1].length == noAtoms);
    assert (fractionalCoords[2].length == noAtoms);
    return fractionalCoords;
  }

  @Override
  public short[] getAllSpins() {
    return this.spins;
  }

  @Override
  public float[] getAllCharges() {
    return this.charges;
  }

  @Override
  public void setAllCharges(final float[] newCharges) {
    assert (newCharges != null);
    assert (newCharges.length == noAtoms);
  }

  @Override
  public void setAllSpins(final short[] newSpins) {
    assert (newSpins != null);
    assert (newSpins.length == noAtoms);
  }
}
