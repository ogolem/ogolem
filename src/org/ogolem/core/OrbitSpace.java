/*
Copyright (c) 2010, J. M. Dieterich
              2016-2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.core;

import org.ogolem.random.Lottery;

/**
 * Defines an orbit space, i.e., a "shell" of a sphere.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
final class OrbitSpace implements AllowedSpace {

  private static final long serialVersionUID = (long) 20200622;

  private final Lottery random;

  private final double lowOrb;

  private final double highOrb;

  private final double[] center;

  OrbitSpace(final double[] center, final double lowestOrb, final double highestOrb) {
    this.center = center;
    this.lowOrb = lowestOrb;
    this.highOrb = highestOrb;
    this.random = Lottery.getInstance();
  }

  @Override
  public OrbitSpace copy() {

    final double[] tmp = center.clone();

    return new OrbitSpace(tmp, this.lowOrb, this.highOrb);
  }

  @Override
  public double[] getPointInSpace() {

    // we need some randoms
    final double[] spherical = new double[3];
    // radius
    spherical[0] = (highOrb - lowOrb) * random.nextDouble() + lowOrb;

    // phi and omega
    spherical[1] = random.nextDouble() * 2.0 * Math.PI;
    spherical[2] = random.nextDouble() * Math.PI;

    // translate to cartesian
    final double[] cartes = new double[3];
    CoordTranslation.sphericalToCartesianCoord(spherical, cartes);

    // move with respect to center
    for (int i = 0; i < 3; i++) {
      cartes[i] += center[i];
    }

    assert (isPointInSpace(cartes));

    return cartes;
  }

  @Override
  public boolean isPointInSpace(double[] point) {

    // move with respect to middle of sphere
    final double x = point[0] - center[0];
    final double y = point[1] - center[1];
    final double z = point[2] - center[2];

    // calculate only the radius
    final double r = Math.sqrt(x * x + y * y + z * z);

    // check
    return (r >= lowOrb && r <= highOrb);
  }
}
