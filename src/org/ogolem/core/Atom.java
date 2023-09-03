/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.lang.reflect.Field;
import org.ogolem.generic.Copyable;

/**
 * This defines an atom.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
class Atom implements Serializable, Copyable {

  private static final long serialVersionUID = (long) 20170130;
  // which name is it supposed to have?
  private final String sID;
  // which atom number
  private final short atomNo;
  // assign another ID in the molecule, integer should be enough
  private final int iID;
  // energy part of the atom
  private double energyPart;
  // the position in an arbitrary cartesian coordinate system
  private double[] position;

  /** Constructor filling all the fields but the derived ones... */
  Atom(final AtomConfig ac) {
    this.sID = ac.sID;
    this.iID = ac.iID;
    this.atomNo = ac.atomNo;
    assert (ac.atomNo == AtomicProperties.giveAtomicNumber(sID));
    this.energyPart = ac.energypart;
    this.position = ac.position;
  }

  /** A copy constructor. */
  Atom(final Atom original) {
    this.energyPart = original.energyPart;
    this.sID = original.sID;
    this.iID = original.iID;
    this.atomNo = original.atomNo;
    assert (original.atomNo == AtomicProperties.giveAtomicNumber(sID));
    this.position = original.position.clone();
  }

  @Override
  public Atom copy() {
    return new Atom(this);
  }

  /*
   * Getters and Setters
   */
  String getSID() {
    return sID;
  }

  double getWeight() {
    return AtomicProperties.giveWeight(atomNo);
  }

  int getID() {
    return iID;
  }

  double getEnergyPart() {
    return energyPart;
  }

  double[] getPosition() {
    return position;
  }

  void setEnergyPart(double energypart) {
    this.energyPart = energypart;
  }

  double getRadius() {
    return AtomicProperties.giveRadius(atomNo);
  }

  /*
   * Methods
   */
  /**
   * Returns its own atomic configuration as an {@code AtomConfig} object.
   *
   * @return ac
   */
  AtomConfig returnMyConfig() {
    AtomConfig ac = new AtomConfig();
    ac.iID = iID;
    ac.sID = sID;
    return ac;
  }

  private void writeObject(ObjectOutputStream oos) throws IOException {
    oos.writeDouble(energyPart);
    oos.writeInt(iID);
    oos.writeShort(atomNo);
    oos.writeUTF(sID);
    oos.writeDouble(position[0]);
    oos.writeDouble(position[1]);
    oos.writeDouble(position[2]);
  }

  private void readObject(ObjectInputStream ois)
      throws IOException,
          ClassNotFoundException,
          ClassCastException,
          IllegalAccessException,
          NoSuchFieldException {

    energyPart = ois.readDouble();

    final Class<? extends Atom> cl = this.getClass();
    final Field f = cl.getDeclaredField("iID");
    final int tempID = ois.readInt();
    f.setAccessible(true);
    f.set(this, tempID);

    final Field f1 = cl.getDeclaredField("atomNo");
    final short tempNo = ois.readShort();
    f1.setAccessible(true);
    f1.set(this, tempNo);

    final Field f2 = cl.getDeclaredField("sID");
    final String tempSID = ois.readUTF();
    f2.setAccessible(true);
    f2.set(this, tempSID);

    position = new double[3];
    position[0] = ois.readDouble();
    position[1] = ois.readDouble();
    position[2] = ois.readDouble();
  }
}
