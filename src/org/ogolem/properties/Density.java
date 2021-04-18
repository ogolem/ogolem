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
 * A density on a grid (Cartesian, equidistantly spaced) property.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class Density extends TensorProperty {

  private static final long serialVersionUID = (long) 20170925;

  public static final double DEFAULTDENSITY = 0.0;

  public Density(final double[][][] density) {
    super(density, true);
  }

  private Density(final Density orig) {
    super(orig);
  }

  @Override
  public Density copy() {
    return new Density(this);
  }

  /**
   * Will return the integral over the density, except for negative densities!
   *
   * @return integral over the density.
   */
  @Override
  public double getValue() {

    if (data == null) {
      return FixedValues.NONCONVERGEDENERGY;
    }

    double sum = 0.0;
    for (int i = 0; i < data.length; i++) {
      for (int j = 0; j < data[i].length; j++) {
        for (int k = 0; k < data[i][j].length; k++) {
          if (data[i][j][k] < 0.0) continue;
          sum += data[i][j][k];
        }
      }
    }

    return sum;
  }

  /**
   * This will ALSO return absolute differences!
   *
   * @param p the other property. Must be an instance of Density.
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
        for (int k = 0; k < data[i][j].length; k++) {
          if (Double.isInfinite(data[i][j][k])
              || Double.isNaN(data[i][j][k])
              || data[i][j][k] > FixedValues.NONCONVERGEDGRADIENT) {
            data[i][j][k] = 0.0;
            wasTouched = true;
          }
        }
      }
    }

    return wasTouched;
  }

  @Override
  public String printableProperty() {

    if (data == null) {
      return "NULL'D DENSITY";
    }

    String s = "";
    for (int i = 0; i < data.length; i++) {
      for (int j = 0; j < data[i].length; j++) {
        for (int k = 0; k < data[i][j].length; k++) {
          s += data[i][j][k] + "\t";
        }

        s += "\n";
      }

      s += "\n";
    }

    return s;
  }

  @Override
  public String name() {
    return "DENSITY";
  }

  public static Density getDefaultDensity(final int gridX, final int gridY, final int gridZ) {

    final double[][][] densityVals = new double[gridX][gridY][gridZ];
    for (int x = 0; x < gridX; x++) {
      for (int y = 0; y < gridY; y++) {
        for (int z = 0; z < gridZ; z++) {
          densityVals[x][y][z] = 42.0;
        }
      }
    }

    final Density density = new Density(densityVals);

    return density;
  }

  public int getGridDimX() {

    if (data == null) {
      return 0;
    }

    assert (data.length > 0);
    return data.length;
  }

  public int getGridDimY() {

    if (data == null) {
      return 0;
    }

    assert (data[0].length > 0);
    return data[0].length;
  }

  public int getGridDimZ() {

    if (data == null) {
      return 0;
    }

    assert (data[0][0].length > 0);
    return data[0][0].length;
  }

  @Override
  protected boolean ensureCorrectProperty(Property p) {
    return (p instanceof Density);
  }
}
