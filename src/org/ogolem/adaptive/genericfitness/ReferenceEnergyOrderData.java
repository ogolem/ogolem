/*
Copyright (c) 2016-2020, J. M. Dieterich and B. Hartke
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

import java.util.List;
import org.ogolem.core.StructuralData;
import org.ogolem.helpers.Tuple;
import org.ogolem.properties.EnergyOrder;

/**
 * Data for energy order calculations.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class ReferenceEnergyOrderData<V extends StructuralData>
    implements ReferenceInputData<EnergyOrder> {

  private static final long serialVersionUID = (long) 20160716;

  private final List<Tuple<String, ReferenceGeomData<EnergyOrder, V>>> geoms;
  private final int refPoint;

  public ReferenceEnergyOrderData(
      final int refPoint, final List<Tuple<String, ReferenceGeomData<EnergyOrder, V>>> geoms) {
    assert (geoms != null);
    assert (!geoms.isEmpty());

    this.geoms = geoms;
    this.refPoint = refPoint;
  }

  private ReferenceEnergyOrderData(final ReferenceEnergyOrderData<V> orig) {
    this.geoms = orig.geoms; // just re-use
    this.refPoint = orig.refPoint;
  }

  @Override
  public ReferenceEnergyOrderData<V> copy() {
    return new ReferenceEnergyOrderData<>(this);
  }

  public ReferenceGeomData<EnergyOrder, V> getGeomDataForKey(final String key) {

    for (final Tuple<String, ReferenceGeomData<EnergyOrder, V>> ent : geoms) {
      if (ent.getObject1().equalsIgnoreCase(key)) {
        return ent.getObject2();
      }
    }

    System.err.println(
        "ERROR: No reference data found for key " + key + ". Have fun with the null pointer...");
    return null;
  }

  public List<Tuple<String, ReferenceGeomData<EnergyOrder, V>>> getAllGeomData() {
    return geoms;
  }

  @Override
  public int belongsToReferencePoint() {
    return refPoint;
  }
}
