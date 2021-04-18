/*
Copyright (c) 2012-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
              2017-2020, J. M. Dieterich and B. Hartke
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

/**
 * A degree of freedom.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class DegreeOfFreedom extends ScalarProperty {

  private static final long serialVersionUID = (long) 20171215;
  private final int[] refPoints;

  public DegreeOfFreedom(final double value, final int[] points) {
    super(value);
    this.refPoints = points;
  }

  @Override
  public DegreeOfFreedom copy() {
    return new DegreeOfFreedom(this.getValue(), refPoints.clone());
  }

  @Override
  public boolean makeSensible() {

    // treatment depends on the exact degree of freedom type
    boolean unsensible = false;
    if (Double.isInfinite(this.getValue()) || Double.isNaN(this.getValue())) {
      this.scalar = 0.0;
      unsensible = true;
    }

    if (refPoints.length == 2) {
      // bond
      if (this.scalar <= 0.0) {
        this.scalar = 1.0;
        unsensible = true;
      }
    } else if (refPoints.length == 3) {
      // angle
      // TODO implement
    } else if (refPoints.length == 4) {
      // dihedral
      // TODO implement
    } else {
      // unknown
      System.err.println(
          "ERROR: Degree of freedom is unknown (might be sensible for out-of-plane etc). Contact author(s).");
    }

    return unsensible;
  }

  @Override
  public String printableProperty() {

    String type;
    switch (refPoints.length) {
      case 2:
        type = "bond";
        break;
      case 3:
        type = "angle";
        break;
      case 4:
        type = "dihedral";
        break;
      default:
        type = "unknown!";
        break;
    }

    String points = " atoms: ";
    for (final double p : refPoints) {
      points += p + "  ";
    }

    return "" + this.getValue() + "(" + type + ")" + points;
  }

  @Override
  public String name() {
    return "DEGREE OF FREEDOM";
  }

  @Override
  protected boolean ensureCorrectProperty(Property p) {
    return (p instanceof DegreeOfFreedom);
  }
}
