/*
Copyright (c) 2014-2015, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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

import org.apache.commons.math3.util.FastMath;
import org.ogolem.helpers.Tuple;
import org.ogolem.math.LookupedFunction;

/**
 * A collection of helpful snippets for TIPnP type water force fields.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class TIPnPHelpers {

  /**
   * Sanitizes a water molecule to have certain internal degrees of freedom.
   *
   * @param mol must be a water molecule. NOT BEING CHECKED!
   * @param rOH the wished for distance of O to H in a.u.
   * @param angHOH the wished for angle H-O-H in deg
   */
  static void sanitizeWater(final Molecule mol, final double rOH, final double angHOH) {

    // XXX first check if the current one is fine
    // final double[][] currXYZ = mol.getReferenceCartesians();
    // final double distOH1Sq;

    final CartesianCoordinates cMol = mol.getCartesians();

    // this may be not the fastest, but it is certainly the simplest
    final CartesianCoordinates cCorr = cMol.copy();
    final double[][] xyzCorrect = cCorr.getAllXYZCoord();
    // O: in 0/0/0
    xyzCorrect[0][0] = 0.0;
    xyzCorrect[1][0] = 0.0;
    xyzCorrect[2][0] = 0.0;
    // H1: in rOH/0/0
    xyzCorrect[0][1] = rOH;
    xyzCorrect[1][1] = 0.0;
    xyzCorrect[2][1] = 0.0;
    // H2: slightly more involved...
    final double realRot = Math.toRadians(angHOH);
    final double sinR = FastMath.sin(realRot);
    final double cosR = FastMath.cos(realRot);

    // basically roate H1 around the z-axis
    xyzCorrect[0][2] = cosR * rOH;
    xyzCorrect[1][2] = sinR * rOH;
    xyzCorrect[2][2] = 0.0;

    // take care of COM moves
    final double[] comDevMol = cMol.calculateTheCOM();
    cMol.moveCoordsToCOM();
    cCorr.moveCoordsToCOM();

    // align the "correct" with the reference
    final Tuple<CartesianCoordinates, Double> aligned;
    try {
      aligned = CoordTranslation.alignTwoCartesians(cMol, cCorr);
    } catch (Exception e) {
      System.err.println(
          "WARNING: Aligning of rigidified and original failed. Using horrible rotated rigidified one!");
      mol.setReferenceCartesian(cCorr.getAllXYZCoord());
      return;
    }
    aligned.getObject1().moveCoordsToPoint(comDevMol);
    mol.setReferenceCartesian(aligned.getObject1().getAllXYZCoord());
  }

  /**
   * An anti-cold fusion function to be used in TIPnP evaluations as going back to the original
   * phenix work by Bernd Hartke.
   */
  public static class AntiColdFusionFunction implements LookupedFunction {

    private static final long serialVersionUID = (long) 20141221;

    @Override
    public AntiColdFusionFunction copy() {
      return new AntiColdFusionFunction();
    }

    @Override
    public double func(final double radius) {
      return Math.exp(-radius) * Math.pow(radius, -7);
    }

    double grad(final double radius) {
      // TODO
      return 0.0;
    }
  }
}
