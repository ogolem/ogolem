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
package org.ogolem.core;

import org.ogolem.random.Lottery;

/**
 * A half-sphere, where the base is always in the x-y plane.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class HalfSphereSpace implements AllowedSpace {

  private static final long serialVersionUID = (long) 20200622;

  private final Lottery random;

  private final double radius;
  private final double[] center;

  HalfSphereSpace(final double[] middle, final double radius) {
    this.center = middle;
    this.radius = radius;
    this.random = Lottery.getInstance();
  }

  private HalfSphereSpace(final HalfSphereSpace orig) {
    this.radius = orig.radius;
    this.center = orig.center;
    this.random = Lottery.getInstance();
  }

  @Override
  public HalfSphereSpace copy() {
    return new HalfSphereSpace(this);
  }

  @Override
  public double[] getPointInSpace() {

    // we need some randoms
    final double[] spherical = new double[3];
    // radius
    spherical[0] = radius * random.nextDouble();

    // phi and omega
    spherical[1] = random.nextDouble() * 2.0 * Math.PI;
    spherical[2] = random.nextDouble() * Math.PI;

    // translate to cartesian
    final double[] point = new double[3];
    CoordTranslation.sphericalToCartesianCoord(spherical, point);

    // move with respect to middle
    for (int i = 0; i < 3; i++) {
      point[i] += center[i];
    }

    // now ensure that it is in the upper half of the sphere
    if (point[2] < center[2]) {
      final double diff = center[2] - point[2];
      point[2] = center[2] + diff;
    }

    assert (isPointInSpace(point));

    return point;
  }

  @Override
  public boolean isPointInSpace(final double[] point) {

    if (point[2] < center[2]) {
      return false; // too low
    }

    // move with respect to middle of sphere
    final double x = point[0] - center[0];
    final double y = point[1] - center[1];
    final double z = point[2] - center[2];

    // calculate only the radius
    final double r = Math.sqrt(x * x + y * y + z * z);

    // check
    return (r <= radius);
  }
}
