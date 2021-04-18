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

import java.util.List;
import org.ogolem.core.CollisionInfo.Collision;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.Lottery;

/**
 * Xovers with an environment: Vinland.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
class VinlandGeometryXOver implements GenericCrossover<Molecule, Geometry> {

  static enum MOVEMODE {
    FULLRANDOM,
    SURFACESTYLE
  };

  static final double DEFAULTINCRBOHR = 0.2;
  static final int DEFAULTMAXMOVETRIES = 20; // max: 4 bohr move!
  static final int DEFAULTMAXTRIES = 200;

  private static final long serialVersionUID = (long) 20200719;
  private static final boolean DEBUG = false;

  private final GenericCrossover<Molecule, Geometry> xover;
  private final CollisionDetectionEngine cd;
  private final double blowFac;
  private final Lottery rnd;
  private final MOVEMODE mode;
  private final double movePerStep;
  private final int maxMoveTries;
  private final int maxTries;

  VinlandGeometryXOver(
      final GenericCrossover<Molecule, Geometry> xover,
      final double blowFac,
      final MOVEMODE mode,
      final double movePerStep,
      final int maxMoveTries,
      final int maxTries) {
    assert (maxMoveTries > 0);
    assert (maxTries > 0);
    this.xover = xover;
    this.cd =
        (mode == MOVEMODE.SURFACESTYLE)
            ? new CollisionDetection(CollisionDetection.CDTYPE.SIMPLEPAIRWISE)
            : new CollisionDetection(CollisionDetection.CDTYPE.ADVANCEDPAIRWISE);
    this.blowFac = blowFac;
    this.rnd = Lottery.getInstance();
    this.mode = mode;
    this.movePerStep = movePerStep;
    this.maxMoveTries = maxMoveTries;
    this.maxTries = maxTries;
  }

  VinlandGeometryXOver(final VinlandGeometryXOver orig) {
    this.cd = orig.cd.copy();
    this.blowFac = orig.blowFac;
    this.xover = orig.xover;
    this.rnd = Lottery.getInstance();
    this.mode = orig.mode;
    this.movePerStep = orig.movePerStep;
    this.maxMoveTries = orig.maxMoveTries;
    this.maxTries = orig.maxTries;
  }

  @Override
  public VinlandGeometryXOver copy() {
    return new VinlandGeometryXOver(this);
  }

  @Override
  public String getMyID() {
    return "VINLAND GEOMETRY XOVER:\n" + xover.getMyID();
  }

  @Override
  public Tuple<Geometry, Geometry> crossover(
      final Geometry mother, final Geometry father, final long futureID) {

    boolean geomHasCollision = true;
    int tries = 0;

    Geometry g1 = null;
    Geometry g2 = null;

    do {

      final Tuple<Geometry, Geometry> children = xover.crossover(mother, father, futureID);
      final Geometry gTmp1 = children.getObject1();
      final Geometry gTmp2 = children.getObject2();

      if (gTmp1 != null) {
        final CartesianCoordinates c1 = gTmp1.getCartesians();
        final BondInfo bonds = gTmp1.getBondInfo();
        final boolean hasColl1 = cd.checkOnlyForCollision(c1, blowFac, bonds);
        if (!hasColl1) {
          if (g1 == null) {
            g1 = gTmp1;
          } else {
            g2 = gTmp1;
          }
        }
      }

      if (gTmp2 != null) {
        final CartesianCoordinates c2 = gTmp2.getCartesians();
        final BondInfo bonds = gTmp2.getBondInfo();
        final boolean hasColl2 = cd.checkOnlyForCollision(c2, blowFac, bonds);
        if (!hasColl2) {
          if (g2 == null) {
            g2 = gTmp2;
          } else if (g1 == null) {
            g1 = gTmp2;
          }
        }
      }

      if (g1 != null && g2 != null) {
        // we still could have dissociation, but we do not care...
        geomHasCollision = false; // got two now
      }

      tries++;
    } while (geomHasCollision && tries < maxTries);

    if (tries >= maxTries) {
      return new Tuple<>(g1, g2); // they may be partially null
    }

    if (!g1.containsEnvironment()) return new Tuple<>(g1, g2);

    final CollisionInfo cI = new MultiCollisionInfo();
    if (g1 != null) {
      environmentFitting(g1, cI, tries);
    }

    if (g2 != null) {
      environmentFitting(g2, cI, tries);
    }

    return new Tuple<>(g1, g2);
  }

  @Override
  public short hasPriority() {
    return -1;
  }

  private void environmentFitting(Geometry g, CollisionInfo cI, int tries) {

    boolean hasEnvColl = true;
    int lastColls = 0;
    int dirMove = rnd.nextInt(3);
    boolean positiveMove = true;
    int myTries = 0;
    // only get the bonds once - no need to repeatedly bang on that
    final Tuple<CartesianCoordinates, BondInfo> tup = g.getCartesiansAndBondsWithEnvironment();
    final CartesianCoordinates c = tup.getObject1();
    final double[][] xyz = c.getAllXYZCoord();
    final BondInfo bonds = tup.getObject2();
    final int atomsInEnv =
        c.getAtomsPerMol(
            c.getNoOfMolecules()
                - 1); // we know that we have an environment and we know the last "molecule" is the
    // environment
    final int totAtomsWOEnv = c.getNoOfAtoms() - atomsInEnv;

    final double[] moveVec = new double[3];

    do {

      if (mode == MOVEMODE.SURFACESTYLE) {
        // use offsets to check only for collisions between environment and cluster. not inside
        // cluster/env
        final boolean hasColl =
            cd.checkOnlyForCollision(c, blowFac, bonds, totAtomsWOEnv, totAtomsWOEnv);
        if (!hasColl) break; // we still could have dissociation but we do not care here.

        // move always in z, always positive
        for (int i = 0; i < totAtomsWOEnv; i++) {
          xyz[2][i] += movePerStep;
        }

        moveVec[2] += movePerStep;
      } else {

        cd.checkForCollision(c, blowFac, bonds, cI);

        final boolean hasColl = cI.hasCollision();
        if (!hasColl) {
          // we still could have dissociation, but we do not care...
          break;
        }

        final List<Collision> colls = cI.getCollisions();
        final int thisColls = colls.size();

        if (thisColls > lastColls) {
          // we do need a new direction
          dirMove = rnd.nextInt(3);
          positiveMove = rnd.nextBoolean();
        } // if not: we are doing fine -> continue in that direction

        // let's move the whole cluster a bit
        final double move = (positiveMove) ? movePerStep : -movePerStep;

        for (int i = 0; i < totAtomsWOEnv; i++) {
          xyz[dirMove][i] += move;
        }

        moveVec[dirMove] += move;
        lastColls = thisColls;
      }

      tries++;
      myTries++;
    } while (hasEnvColl && tries < maxTries && myTries < maxMoveTries);

    g.getEnvironment().moveCluster(0, moveVec[0]);
    g.getEnvironment().moveCluster(1, moveVec[1]);
    g.getEnvironment().moveCluster(2, moveVec[2]);
  }
}
