/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

import java.util.List;
import org.ogolem.core.CollisionInfo.Collision;
import org.ogolem.core.Environment.ENVIRONMENTTYPE;
import org.ogolem.random.RandomUtils;

/**
 * Holds as the father of all environments for the cluster.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public final class SimpleSurface extends SimpleEnvironment {

  private static final long serialVersionUID = (long) 20200719;

  /**
   * The constructor for the initial creation.
   *
   * @param environment
   * @param flexy
   * @param bonds
   */
  SimpleSurface(
      final CartesianCoordinates environment,
      final boolean flexy,
      final BondInfo bonds,
      final List<Integer> secondaryAtoms,
      final CollisionDetection.CDTYPE whichColl,
      final double collBlow,
      final Atom[] referenceAtoms,
      final AllowedSpace allowedSpace,
      final FITMODE fit) {
    super(
        environment,
        flexy,
        bonds,
        secondaryAtoms,
        whichColl,
        collBlow,
        referenceAtoms,
        allowedSpace,
        fit,
        ENVIRONMENTTYPE.SURFACE);
  }

  /**
   * A standard copy constructor.
   *
   * @param original
   */
  private SimpleSurface(final SimpleSurface original) {
    super(original);
  }

  @Override
  public SimpleSurface copy() {
    return new SimpleSurface(this);
  }

  @Override
  public void initializeConnections(
      final CartesianCoordinates clusterCartes, final BondInfo clusterBonds) {

    final CollisionDetection collDetect = new CollisionDetection(whichCollDetect);

    int emergCounter = 0;
    boolean hasColl;
    do {

      if (DEBUG) {
        System.out.println(
            "DEBUG: Trying to initialize connections for " + (emergCounter + 1) + " time");
      }

      if (mode != FITMODE.LAYERONLY) {
        // first the eulers
        RandomUtils.randomEulers(eulersToCluster);
      }

      // now the distances
      distanceCOMs = space.getPointInSpace();

      // check for collisions. DO NOT check for collisions within cluster!
      final CollisionInfo collInfo =
          collDetect.checkForCollision(
              marryThem(clusterCartes, clusterBonds).getObject1(),
              blowColl,
              marryBonds(
                  clusterCartes.getNoOfAtoms(),
                  (clusterCartes.getNoOfAtoms() + cartEnv.getNoOfAtoms())));
      hasColl = collInfo.hasCollision();
      emergCounter++;

      if (hasColl && DEBUG) {
        // get more info
        final List<Collision> colls = collInfo.getCollisions();
        for (final Collision coll : colls) {
          System.out.println(
              "DEBUG: Collision between " + coll.getAtomOne() + " and " + coll.getAtomTwo());
        }
      }

      if (hasColl) {
        // well, easy enough: we need to move up and since we define the surface to be in z, it's
        // easy enough
        // let's start with moving 50 au up
        double crash = distanceCOMs[2];
        distanceCOMs[2] += 50;
        final CollisionInfo ci =
            collDetect.checkForCollision(
                marryThem(clusterCartes, clusterBonds).getObject1(),
                blowColl,
                marryBonds(
                    clusterCartes.getNoOfAtoms(),
                    (clusterCartes.getNoOfAtoms() + cartEnv.getNoOfAtoms())));
        final boolean co = ci.hasCollision();
        if (co) {
          System.err.println(
              "ERROR: Moving points upwards by 50 a.u. not sufficient. That's weird...");
          continue;
        }
        double noCrash = distanceCOMs[2];

        boolean ct = false;
        while (ct && emergCounter < FixedValues.MAXTOEMERGENCY) {

          final double nextGuess = (noCrash - crash) / 2;
          distanceCOMs[2] = nextGuess;
          final CollisionInfo cix =
              collDetect.checkForCollision(
                  marryThem(clusterCartes, clusterBonds).getObject1(),
                  blowColl,
                  marryBonds(
                      clusterCartes.getNoOfAtoms(),
                      (clusterCartes.getNoOfAtoms() + cartEnv.getNoOfAtoms())));
          final boolean cox = cix.hasCollision();

          if (cox) {
            crash = nextGuess;
          } else {
            ct = (noCrash - nextGuess > 0); // rather arbitrary but below 1 a.u. seems good enough?
            noCrash = nextGuess;
          }

          emergCounter++;
        }

        if (!ct) {
          // we found something that works
          hasColl = false;
        }
      }

      // check if the point is still in space
      final boolean isInSpace = space.isPointInSpace(distanceCOMs);
      if (!isInSpace) {
        System.err.println(
            "ERROR: We found a collision free location of the cluster on the surface but it's not in the allowed space. Maybe it's too small in z?");
        hasColl = true;
      }

    } while (hasColl && emergCounter < FixedValues.MAXTOEMERGENCY);

    assert (space.isPointInSpace(distanceCOMs));

    if (DEBUG) {
      System.out.println(
          "Took " + emergCounter + " attempts to initialize connections of cluster with env.");
    }
  }
}
