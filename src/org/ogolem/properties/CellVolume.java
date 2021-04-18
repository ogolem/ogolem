/*
Copyright (c) 2015, J. M. Dieterich and B. Hartke
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
 * An (equilibrium) cell volume property.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class CellVolume extends ScalarProperty {

  private static final long serialVersionUID = (long) 20171215;

  public CellVolume(final double cellVolume) {
    super(cellVolume);
    assert (cellVolume >= 0.0);
  }

  private CellVolume(final CellVolume orig) {
    super(orig.scalar);
  }

  @Override
  public CellVolume copy() {
    return new CellVolume(this);
  }

  @Override
  public boolean makeSensible() {
    if (Double.isInfinite(this.getValue())
        || Double.isNaN(this.getValue())
        || this.getValue() < 0.0) {
      this.scalar = 0.0; // not too sane...
      return true;
    }

    return false;
  }

  @Override
  public String printableProperty() {
    return "" + this.scalar;
  }

  @Override
  public String name() {
    return "CELL VOLUME";
  }

  @Override
  protected boolean ensureCorrectProperty(Property p) {
    return (p instanceof CellVolume);
  }
}
