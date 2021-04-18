/*
Copyright (c) 2014, J. M. Dieterich
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
package org.ogolem.adaptive;

import org.ogolem.generic.GenericSanityCheck;

/**
 * A parameter sanity check.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class ParameterSanityCheck implements GenericSanityCheck<Double, AdaptiveParameters> {

  private static final long serialVersionUID = (long) 20140502;
  private final boolean checkBorders;
  private final double[] lower;
  private final double[] upper;

  public ParameterSanityCheck(
      final boolean checkBorders, final double[] lower, final double[] upper) {
    assert (upper != null);
    assert (lower != null);
    assert (lower.length == upper.length);
    this.checkBorders = checkBorders;
    this.lower = lower;
    this.upper = upper;
  }

  @Override
  public GenericSanityCheck<Double, AdaptiveParameters> copy() {
    return new ParameterSanityCheck(checkBorders, lower.clone(), upper.clone());
  }

  @Override
  public boolean isSane(final AdaptiveParameters individual) {

    final double[] values = individual.getAllParamters();
    if (values.length != upper.length) {
      return false;
    }
    for (int i = 0; i < values.length; i++) {
      final double d = values[i];
      if (Double.isNaN(d) || Double.isInfinite(d)) {
        return false;
      }
      if (checkBorders) {
        if (d > upper[i] || d < lower[i]) {
          return false;
        }
      }
    }

    return true;
  }
}
