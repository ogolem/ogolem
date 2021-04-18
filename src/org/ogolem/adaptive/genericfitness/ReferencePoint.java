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

import org.ogolem.properties.Property;

/**
 * A straighforward implementation of a generic reference point.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class ReferencePoint<T extends Property, V extends ReferenceInputData<T>>
    implements GenericReferencePoint<T, V> {

  private static final long serialVersionUID = (long) 20151028;

  private final T refProp;
  private final V refInput;
  private final int id;
  private final double weight;
  private final double maxDiff;

  public ReferencePoint(
      final T refProperty,
      final V refInputData,
      final int refID,
      final double refWeight,
      final double refMaxDiff) {

    assert (refProperty != null);
    assert (refInputData != null);
    assert (refWeight > 0.0);
    assert (refMaxDiff > 0.0);

    this.refProp = refProperty;
    this.refInput = refInputData;
    this.id = refID;
    this.weight = refWeight;
    this.maxDiff = refMaxDiff;
  }

  @SuppressWarnings("unchecked")
  ReferencePoint(final ReferencePoint<T, V> orig) {
    this.id = orig.id;
    this.maxDiff = orig.maxDiff;
    this.refInput = (V) orig.refInput.copy();
    this.refProp = (T) orig.refProp.copy();
    this.weight = orig.weight;
  }

  @Override
  public ReferencePoint<T, V> copy() {
    return new ReferencePoint<>(this);
  }

  @Override
  public T getReferenceProperty() {
    return refProp;
  }

  @Override
  public V getReferenceInputData() {
    return refInput;
  }

  @Override
  public int getReferenceID() {
    return id;
  }

  @Override
  public double getRefWeight() {
    return weight;
  }

  @Override
  public double getMaxAllowedDiff() {
    return maxDiff;
  }
}
