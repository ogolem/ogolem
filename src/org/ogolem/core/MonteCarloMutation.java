/*
Copyright (c) 2012-2014, J. M. Dieterich
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

import java.io.Serializable;
import java.util.List;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * Implements a mutation operator which in essence is only a MC step.
 *
 * @author Johannes Dieterich
 * @version 2020-07-19
 */
public class MonteCarloMutation implements Serializable {

  private static final long serialVersionUID = (long) 20140327;
  private static final boolean VERBOSE = false;

  private static final Lottery random = Lottery.getInstance();

  public static enum MOVEMODE {
    ALL,
    ONE,
    SOME,
    GAUSSIAN,
    ALLMOL,
    ONEMOL,
    SOMEMOL,
    GAUSSIANMOL
  };

  public static Molecule mutate(
      final Molecule mol,
      final MOVEMODE mode,
      final double maxMove,
      final int gaussMax,
      final double gaussWidth) {

    assert (maxMove > 0.0);
    assert (gaussMax >= 0);
    assert (gaussWidth > 0.0);
    assert (mode != null);

    final CartesianCoordinates cartes = mol.getCartesians();
    final double[][] xyz = cartes.getAllXYZCoord();
    final boolean[][] constr = mol.getConstraints();

    switch (mode) {
      case ONE:
        moveOne(xyz, constr, maxMove);
        break;
      case ALL:
        moveAll(xyz, constr, maxMove);
        break;
      case SOME:
        moveSome(xyz, constr, maxMove);
        break;
      case GAUSSIAN:
        moveGauss(xyz, constr, maxMove, gaussMax, gaussWidth);
        break;
      default:
        System.err.println(
            "ERROR: Unknown case in molecular MonteCarloMutation. Contact author(s)!");
        break;
    }

    final Molecule molNew =
        new Molecule(
            cartes,
            mol.getMolPosition(),
            mol.getSID(),
            mol.getFlexy(),
            mol.getDegreesOfFreedom(),
            mol.isConstricted(),
            mol.getConstraints());

    return molNew;
  }

  public static Geometry mutate(
      final Geometry g,
      final MOVEMODE mode,
      final double maxMove,
      final int gaussMax,
      final double gaussWidth) {

    assert (maxMove > 0.0);
    assert (gaussMax > 0);
    assert (gaussWidth > 0.0);

    if (mode == MOVEMODE.ALL
        || mode == MOVEMODE.ONE
        || mode == MOVEMODE.SOME
        || mode == MOVEMODE.GAUSSIAN) {
      final CartesianCoordinates cartes = new CartesianCoordinates(g.getCartesians());
      final double[][] xyz = cartes.getAllXYZCoord();
      final boolean[][] constr = g.getAllConstraintsXYZ(false);

      switch (mode) {
        case ONE:
          moveOne(xyz, constr, maxMove);
          break;
        case ALL:
          moveAll(xyz, constr, maxMove);
          break;
        case SOME:
          moveSome(xyz, constr, maxMove);
          break;
        case GAUSSIAN:
          moveGauss(xyz, constr, maxMove, gaussMax, gaussWidth);
          break;
        default:
          System.err.println(
              "ERROR: Unknown case in geometric MonteCarloMutation. Contact author(s)!");
          break;
      }

      final Geometry gNew = new Geometry(g);
      // we use this routine since it allows us to get above only the cartesians w/o environment and
      // here update the coordinates without needing to consider environment
      gNew.updateAllCoordinates(xyz);

      return gNew;
    } else {

      final Geometry gNew = new Geometry(g);
      final boolean[] constr = g.getAllConstraints(false);

      switch (mode) {
        case ONEMOL:
          moveOne(gNew, constr, maxMove);
          break;
        case ALLMOL:
          moveAll(gNew, constr, maxMove);
          break;
        case SOMEMOL:
          moveSome(gNew, constr, maxMove);
          break;
        case GAUSSIANMOL:
          moveGauss(gNew, constr, maxMove, gaussMax, gaussWidth);
          break;
        default:
          System.err.println(
              "ERROR: Unknown case in geometric MonteCarloMutation2. Contact author(s)!");
          break;
      }

      return gNew;
    }
  }

  private static void moveOne(
      final double[][] xyz, final boolean[][] constr, final double maxMove) {

    final double[] r = new double[3];

    int which;
    boolean cont = true;
    do {
      which = random.nextInt(xyz[0].length);
      if (!constr[0][which] || !constr[1][which] || !constr[2][which]) cont = false;
    } while (cont);

    // move
    randomizer(xyz, constr, maxMove, which, r);
  }

  private static void moveAll(
      final double[][] xyz, final boolean[][] constr, final double maxMove) {

    final double[] r = new double[3];

    for (int i = 0; i < xyz[0].length; i++) {
      if (!constr[0][i] || !constr[1][i] || !constr[2][i]) {
        randomizer(xyz, constr, maxMove, i, r);
      }
    }
  }

  private static void moveSome(
      final double[][] xyz, final boolean[][] constr, final double maxMove) {

    final double target = random.nextDouble();
    final double[] r = new double[3];

    int countMoved = 0;
    for (int i = 0; i < xyz[0].length; i++) {
      if ((!constr[0][i] || !constr[1][i] || !constr[2][i]) && (target > random.nextDouble())) {
        randomizer(xyz, constr, maxMove, i, r);
        countMoved++;
      }
    }

    if (VERBOSE) {
      System.out.println("DEBUG: Moved " + countMoved + " atoms.");
    }
  }

  private static void moveGauss(
      final double[][] xyz,
      final boolean[][] constr,
      final double maxMove,
      final int gaussMax,
      final double gaussWidth) {

    int nonConstr = 0;
    for (int i = 0; i < xyz[0].length; i++) {
      if (!constr[0][i] || !constr[1][i] || !constr[2][i]) nonConstr++;
    }

    final int noMove = (int) RandomUtils.gaussDoubleAroundVal(1, nonConstr, gaussWidth, gaussMax);
    final List<Integer> move = RandomUtils.listOfPoints(noMove, nonConstr);

    int idx = 0;
    final double[] r = new double[3];
    for (int i = 0; i < xyz[0].length; i++) {
      if (!constr[0][i] || !constr[1][i] || !constr[2][i] && move.contains(idx)) {
        randomizer(xyz, constr, maxMove, i, r);
        idx++;
      }
    }

    if (VERBOSE) {
      System.out.println("DEBUG: Moved " + noMove + " atoms.");
    }
  }

  private static void moveOne(final Geometry g, final boolean[] constr, final double maxMove) {

    int which;
    boolean cont = true;
    do {
      which = random.nextInt(g.getNumberOfIndieParticles());
      if (!constr[which]) cont = false;
    } while (cont);

    // move
    final Molecule newMol =
        mutate(
            g.getMoleculeAtPosition(which),
            MOVEMODE.SOME,
            maxMove,
            1,
            1.0); // hardcoding is OK as we are not using the Gaussian mutation
    g.setMoleculeAtPosition(which, newMol);
  }

  private static void moveAll(final Geometry g, final boolean[] constr, final double maxMove) {

    for (int i = 0; i < g.getNumberOfIndieParticles(); i++) {
      if (!constr[i]) {
        final Molecule newMol = mutate(g.getMoleculeAtPosition(i), MOVEMODE.SOME, maxMove, 1, 1.0);
        g.setMoleculeAtPosition(i, newMol);
      }
    }
  }

  private static void moveSome(final Geometry g, final boolean[] constr, final double maxMove) {

    final double target = random.nextDouble();

    int countMoved = 0;
    for (int i = 0; i < g.getNumberOfIndieParticles(); i++) {
      if (!constr[i] && (target > random.nextDouble())) {
        final Molecule newMol = mutate(g.getMoleculeAtPosition(i), MOVEMODE.SOME, maxMove, 1, 1.0);
        g.setMoleculeAtPosition(i, newMol);
        countMoved++;
      }
    }

    if (VERBOSE) {
      System.out.println("DEBUG: Moved " + countMoved + " molecules.");
    }
  }

  private static void moveGauss(
      final Geometry g,
      final boolean[] constr,
      final double maxMove,
      final int gaussMax,
      final double gaussWidth) {

    int nonConstr = 0;
    for (int i = 0; i < constr.length; i++) {
      if (!constr[i]) nonConstr++;
    }

    final int noMove = (int) RandomUtils.gaussDoubleAroundVal(0, nonConstr, gaussWidth, gaussMax);
    final List<Integer> move = RandomUtils.listOfPoints(noMove, nonConstr);

    int idx = 0;
    for (int i = 0; i < g.getNumberOfIndieParticles(); i++) {
      if (!constr[i] && move.contains(idx)) {
        final Molecule newMol = mutate(g.getMoleculeAtPosition(i), MOVEMODE.SOME, maxMove, 1, 1.0);
        g.setMoleculeAtPosition(i, newMol);
        idx++;
      }
    }

    if (VERBOSE) {
      System.out.println("DEBUG: Moved " + noMove + " molecules.");
    }
  }

  private static void randomizer(
      final double[][] xyz,
      final boolean[][] constr,
      final double maxMove,
      final int which,
      final double[] r) {
    RandomUtils.randomVector(r);
    final boolean b = random.nextBoolean();
    final double scaling = (b) ? maxMove * random.nextDouble() : -maxMove * random.nextDouble();
    if (!constr[0][which]) xyz[0][which] += scaling * r[0];
    if (!constr[1][which]) xyz[1][which] += scaling * r[1];
    if (!constr[2][which]) xyz[2][which] += scaling * r[2];
  }
}
