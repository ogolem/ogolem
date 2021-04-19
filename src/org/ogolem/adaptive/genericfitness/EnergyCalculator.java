/*
Copyright (c) 2012, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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

import org.ogolem.adaptive.Adaptivable;
import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.properties.Energy;

/**
 * Calculates an energy. Adaptor to standard Adaptivable interface.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class EnergyCalculator
    implements PropertyCalculator<Energy, ReferenceGeomData<Energy, CartesianCoordinates>> {
  // we actually want to revamp the adaptivables and remove that interface altogether, I think.
  private static final long serialVersionUID = (long) 20131231;
  private final Adaptivable backend;

  public EnergyCalculator(final Adaptivable adapt) {
    this.backend = adapt;
  }

  private EnergyCalculator(EnergyCalculator orig) {
    this.backend = orig.backend.copy();
  }

  @Override
  public EnergyCalculator copy() {
    return new EnergyCalculator(this);
  }

  @Override
  public Energy calculateProperty(
      final AdaptiveParameters p, final ReferenceGeomData<Energy, CartesianCoordinates> geom) {
    return new Energy(backend.energyOfStructWithParams(geom.c, p, geom.id, geom.bonds));
  }

  @Override
  public Energy calculatePropertyGradient(
      final AdaptiveParameters p,
      final ReferenceGeomData<Energy, CartesianCoordinates> geom,
      final double[] grad) {
    for (int i = 0; i < grad.length; i++) {
      grad[i] = 0.0;
    }
    final double en = backend.gradientOfStructWithParams(geom.c, p, geom.id, geom.bonds, grad);
    return new Energy(en);
  }
}
