/*
Copyright (c) 2018-2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.properties;

import org.ogolem.core.FixedValues;

/**
 * A Cauchy stress tensor (i.e., a 3x3 matrix).
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class StressTensor extends MatrixProperty {

  private static final long serialVersionUID = (long) 20180123;

  public StressTensor(final double[][] stress) {
    super(stress, false);
    if (stress != null) {
      assert (stress.length == 3);
      assert (stress[0].length == 3);
      assert (stress[1].length == 3);
      assert (stress[2].length == 3);
    }
  }

  private StressTensor(final StressTensor orig) {
    super(orig);
  }

  @Override
  public StressTensor copy() {
    return new StressTensor(this);
  }

  @Override
  protected boolean ensureCorrectProperty(Property p) {
    return (p instanceof StressTensor);
  }

  @Override
  public boolean makeSensible() {

    if (data == null) {
      return false;
    }

    boolean wasTouched = false;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (Double.isInfinite(data[i][j])
            || Double.isNaN(data[i][j])
            || data[i][j] > FixedValues.NONCONVERGEDGRADIENT) {
          data[i][j] = FixedValues.NONCONVERGEDGRADIENT;
          wasTouched = true;
        }
      }
    }

    return wasTouched;
  }

  @Override
  public String printableProperty() {
    if (data == null) {
      return "NULL'D STRESS TENSOR";
    }

    String s = data[0][0] + "\t" + data[0][1] + "\t" + data[0][2];
    s += data[1][0] + "\t" + data[1][1] + "\t" + data[1][2];
    s += data[2][0] + "\t" + data[2][1] + "\t" + data[2][2];

    return s;
  }

  @Override
  public String name() {
    return "STRESS TENSOR";
  }
}
