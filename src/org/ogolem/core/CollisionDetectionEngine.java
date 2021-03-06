/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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
package org.ogolem.core;

import java.io.Serializable;
import org.ogolem.generic.Copyable;

/**
 * The interface defining what ALL the collision detection engines need to know.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public interface CollisionDetectionEngine extends Serializable, Copyable {

  @Override
  CollisionDetectionEngine copy();

  CollisionInfo checkForCollision(
      final CartesianCoordinates cartesians, final double blowFactor, final BondInfo bonds);

  void checkForCollision(
      final CartesianCoordinates cartesians,
      final double blowFactor,
      final BondInfo bonds,
      final CollisionInfo info);

  boolean checkOnlyForCollision(
      final CartesianCoordinates cartesians, final double blowFactor, final BondInfo bonds);

  /**
   * Check only for collisions in - presumably - a changed part of the Cartesian. Note that this
   * routine is specified to check collisions of the changed part WITH the rest of the Cartesian
   * (i.e., it will check before offset and after endset if applicable)
   *
   * @param cartesians the Cartesian coordinates to be checked
   * @param blowFactor the blow factor to be used to scale atomic radii
   * @param bonds the known bonds in this Cartesian set
   * @param offset the offset from which collisions should be checked (inclusive)
   * @param endset the endset towards which collisions should be checked (exclusive)
   * @return whether a collision was found or not
   */
  boolean checkOnlyForCollision(
      final CartesianCoordinates cartesians,
      final double blowFactor,
      final BondInfo bonds,
      final int offset,
      final int endset);
}
