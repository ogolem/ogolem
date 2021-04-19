/*
Copyright (c) 2017-2020, J. M. Dieterich and B. Hartke
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
 * The base class for a vector property
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public abstract class VectorProperty implements Property {

  private static final long serialVersionUID = (long) 20171215;

  protected final boolean normDifferences;
  protected final double[] data;

  protected VectorProperty(final double[] data, final boolean normDifferences) {
    assert (data != null);
    this.data = data;
    this.normDifferences = normDifferences;
  }

  protected VectorProperty(final VectorProperty orig) {
    this.data = orig.data.clone();
    this.normDifferences = orig.normDifferences;
  }

  @Override
  public abstract VectorProperty copy();

  /**
   * Return the value of a vector property as its vectors norm. If this is not wanted: overriding is
   * necessary.
   *
   * @return the vector norm
   */
  @Override
  public double getValue() {

    if (data == null) {
      return FixedValues.NONCONVERGEDENERGY;
    }

    double sum = 0.0;
    for (int i = 0; i < data.length; i++) {
      sum += data[i] * data[i];
    }

    return Math.sqrt(sum);
  }

  @Override
  public double signedDifference(Property p) {

    if (!ensureCorrectProperty(p)) {
      throw new IllegalArgumentException("Property should be an instance of " + name());
    }

    final VectorProperty vp = (VectorProperty) p;
    if (data == null || vp.data == null) {
      return FixedValues.NONCONVERGEDENERGY;
    }

    if (vp.data.length != this.data.length) {
      throw new RuntimeException(
          "Vector properties differ in lengths: " + vp.data.length + " vs " + this.data.length);
    }

    double diff = 0.0;
    for (int i = 0; i < data.length; i++) {
      diff +=
          (this.data[i]
              - vp.data[
                  i]); // XXX this is suboptimal, but for the time being I have no better definition
      // in my head
    }

    if (this.normDifferences) {
      diff /= data.length;
    }

    return diff;
  }

  @Override
  public double absoluteDifference(Property p) {

    if (!ensureCorrectProperty(p)) {
      throw new IllegalArgumentException("Property should be an instance of " + name());
    }

    final VectorProperty vp = (VectorProperty) p;
    if (data == null || vp.data == null) {
      return FixedValues.NONCONVERGEDENERGY;
    }

    if (vp.data.length != this.data.length) {
      throw new RuntimeException(
          "Vector properties differ in lengths: " + vp.data.length + " vs " + this.data.length);
    }

    double diff = 0.0;
    for (int i = 0; i < data.length; i++) {
      diff += Math.abs(this.data[i] - vp.data[i]);
    }

    if (this.normDifferences) {
      diff /= data.length;
    }

    return diff;
  }

  protected abstract boolean ensureCorrectProperty(final Property p);
}
