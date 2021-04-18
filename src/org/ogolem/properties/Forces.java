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
package org.ogolem.properties;

import org.ogolem.core.FixedValues;

/**
 * A forces property.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class Forces extends MatrixProperty {

  private static final long serialVersionUID = (long) 20171215;

  public static final double DEFAULTFORCE = FixedValues.NONCONVERGEDGRADIENT;

  public Forces(final double[][] forces) {
    super(forces, true);
  }

  private Forces(final Forces orig) {
    super(orig);
  }

  @Override
  public Forces copy() {
    return new Forces(this);
  }

  /**
   * This will ALSO return absolute differences!
   *
   * @param p the other property. Must be an instance of Forces.
   * @return the ABSOLUTE difference
   */
  @Override
  public double signedDifference(final Property p) {
    return absoluteDifference(p);
  }

  @Override
  public boolean makeSensible() {

    if (data == null) {
      return false;
    }

    boolean wasTouched = false;
    for (int i = 0; i < data.length; i++) {
      for (int j = 0; j < data[i].length; j++) {
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
      return "NULL'D FORCES";
    }

    String s = "";
    for (int i = 0; i < data.length; i++) {
      s += "Row " + i + ":";
      for (int j = 0; j < data[i].length; j++) {
        s += data[i][j] + "\t";
      }

      s += "\n";
    }

    return s;
  }

  @Override
  public String name() {
    return "FORCES";
  }

  public static Forces getDefaultForces(final int noAtoms) {

    final double[][] forceVals = new double[3][noAtoms];
    for (int coord = 0; coord < 3; coord++) {
      for (int at = 0; at < noAtoms; at++) {
        forceVals[coord][at] = DEFAULTFORCE;
      }
    }

    final Forces forces = new Forces(forceVals);

    return forces;
  }

  @Override
  protected boolean ensureCorrectProperty(final Property p) {
    return (p instanceof Forces);
  }
}
