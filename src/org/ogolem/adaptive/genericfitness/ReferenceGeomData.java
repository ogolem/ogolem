/*
Copyright (c) 2015-2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive.genericfitness;

import org.ogolem.core.BondInfo;
import org.ogolem.core.StructuralData;
import org.ogolem.properties.Property;

/**
 * Reference geometrical data, never hurts to have, right? :-)
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class ReferenceGeomData<T extends Property, V extends StructuralData>
    implements ReferenceInputData<T> {

  private static final long serialVersionUID = (long) 20160716;

  public final V c;
  public final BondInfo bonds;
  public final int id;

  public ReferenceGeomData(final V c, final BondInfo bonds, final int id) {
    this.c = c;
    this.bonds = bonds;
    this.id = id;
  }

  @SuppressWarnings("unchecked") // no generic clone
  private ReferenceGeomData(final ReferenceGeomData<T, V> orig) {
    this.c = (V) orig.c.copy();
    this.bonds = orig.bonds.copy();
    this.id = orig.id;
  }

  @Override
  public ReferenceGeomData<T, V> copy() {
    return new ReferenceGeomData<>(this);
  }

  @Override
  public int belongsToReferencePoint() {
    return id;
  }
}
