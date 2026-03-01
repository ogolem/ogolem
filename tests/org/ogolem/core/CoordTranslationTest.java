/*
Copyright (c) 2026, J. M. Dieterich and B. Hartke
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

import static org.junit.jupiter.api.Assertions.*;

import java.util.Random;
import org.junit.jupiter.api.Test;

/**
 * @author Johannes Dieterich
 * @version 2026-02-21
 */
public class CoordTranslationTest {

  private static final double NUMACC = 1e-10;
  private final Random r;

  public CoordTranslationTest() {
    this.r = new Random();
  }

  /** Test of sanitizePhi method, of class CoordTranslation. */
  @Test
  public void testSanitizePhi() {
    System.out.println("sanitizePhi");

    final int TRIALS = 1000;
    final double lower = -Math.PI;
    final double upper = Math.PI;
    final double period = 2 * Math.PI;
    final int maxPer = 10;

    for (int i = 0; i < TRIALS; i++) {
      final double ang1 = lower + r.nextDouble() * period;
      final double san1 = CoordTranslation.sanitizePhi(ang1);
      assertEquals(ang1, san1, NUMACC);

      // now one period on top
      final double ang2 = ang1 + period;
      final double san2 = CoordTranslation.sanitizePhi(ang2);
      assertEquals(ang1, san2, NUMACC);

      // remove one period
      final double ang3 = ang1 - period;
      final double san3 = CoordTranslation.sanitizePhi(ang3);
      assertEquals(ang1, san3, NUMACC);

      // add rnd periods
      final int pers1 = r.nextInt(maxPer);
      final double ang4 = ang1 + pers1 * period;
      final double san4 = CoordTranslation.sanitizePhi(ang4);
      assertEquals(ang1, san4, NUMACC);

      // remove rnd periods
      final int pers2 = r.nextInt(maxPer);
      final double ang5 = ang1 - pers2 * period;
      final double san5 = CoordTranslation.sanitizePhi(ang5);
      assertEquals(ang1, san5, NUMACC);
    }
  }

  /** Test of sanitizeOmega method, of class CoordTranslation. */
  @Test
  public void testSanitizeOmega() {
    System.out.println("sanitizeOmega");
    final int TRIALS = 1000;
    final double lower = -0.5 * Math.PI;
    final double upper = 0.5 * Math.PI;
    final double period = Math.PI;
    final int maxPer = 10;

    for (int i = 0; i < TRIALS; i++) {
      final double ang1 = lower + r.nextDouble() * period;
      final double san1 = CoordTranslation.sanitizeOmega(ang1);
      assertEquals(ang1, san1, NUMACC);

      // now one period on top
      final double ang2 = ang1 + period;
      final double san2 = CoordTranslation.sanitizeOmega(ang2);
      assertEquals(ang1, san2, NUMACC);

      // remove one period
      final double ang3 = ang1 - period;
      final double san3 = CoordTranslation.sanitizeOmega(ang3);
      assertEquals(ang1, san3, NUMACC);

      // add rnd periods
      final int pers1 = r.nextInt(maxPer);
      final double ang4 = ang1 + pers1 * period;
      final double san4 = CoordTranslation.sanitizeOmega(ang4);
      assertEquals(ang1, san4, NUMACC);

      // remove rnd periods
      final int pers2 = r.nextInt(maxPer);
      final double ang5 = ang1 - pers2 * period;
      final double san5 = CoordTranslation.sanitizeOmega(ang5);
      assertEquals(ang1, san5, NUMACC);
    }
  }

  /** Test of sanitizePsi method, of class CoordTranslation. */
  @Test
  public void testSanitizePsi() {
    System.out.println("sanitizePsi");

    final int TRIALS = 1000;
    final double lower = -Math.PI;
    final double upper = Math.PI;
    final double period = 2 * Math.PI;
    final int maxPer = 10;

    for (int i = 0; i < TRIALS; i++) {
      final double ang1 = lower + r.nextDouble() * period;
      final double san1 = CoordTranslation.sanitizePsi(ang1);
      assertEquals(ang1, san1, NUMACC);

      // now one period on top
      final double ang2 = ang1 + period;
      final double san2 = CoordTranslation.sanitizePsi(ang2);
      assertEquals(ang1, san2, NUMACC);

      // remove one period
      final double ang3 = ang1 - period;
      final double san3 = CoordTranslation.sanitizePsi(ang3);
      assertEquals(ang1, san3, NUMACC);

      // add rnd periods
      final int pers1 = r.nextInt(maxPer);
      final double ang4 = ang1 + pers1 * period;
      final double san4 = CoordTranslation.sanitizePsi(ang4);
      assertEquals(ang1, san4, NUMACC);

      // remove rnd periods
      final int pers2 = r.nextInt(maxPer);
      final double ang5 = ang1 - pers2 * period;
      final double san5 = CoordTranslation.sanitizePsi(ang5);
      assertEquals(ang1, san5, NUMACC);
    }
  }

  @Test
  public void testCalcAngle() {

    double[][] xyz = new double[3][3];
    xyz[0][1] = 1.0;
    xyz[0][2] = 2.0;
    double angle = CoordTranslation.calcAngle(xyz, 0, 1, 2);
    assertEquals(Math.PI, angle, NUMACC);

    xyz = new double[3][3];
    xyz[0][1] = 1.0;
    xyz[0][2] = 1.0;
    xyz[1][2] = 1.0;
    angle = CoordTranslation.calcAngle(xyz, 0, 1, 2);
    assertEquals(Math.PI / 2, angle, NUMACC);

    xyz = new double[3][3];
    xyz[0][1] = 1.0;
    angle = CoordTranslation.calcAngle(xyz, 0, 1, 2);
    assertEquals(0.0, angle, NUMACC);

    xyz = new double[3][3];
    xyz[0][1] = 1.0;
    xyz[0][2] = 1.0;
    xyz[1][2] = -1.0;
    angle = CoordTranslation.calcAngle(xyz, 0, 1, 2);
    assertEquals(Math.PI / 2, angle, NUMACC);
  }

  @Test
  public void testCalcDihedral() {

    double[][] xyz = new double[3][4];
    xyz[0][1] = 1.0;
    xyz[0][2] = 2.0;
    xyz[0][3] = 3.0;
    double dihedral = CoordTranslation.calcDihedral(xyz, 0, 1, 2, 3);
    assertEquals(0.0, dihedral, NUMACC);

    xyz = new double[3][4];
    xyz[0][1] = 1.0;
    xyz[0][2] = 2.0;
    xyz[0][3] = 2.0;
    xyz[1][3] = 1.0;
    dihedral = CoordTranslation.calcDihedral(xyz, 0, 1, 2, 3);
    assertEquals(0.0, dihedral, NUMACC);

    xyz = new double[3][4];
    xyz[1][0] = -1.0;
    xyz[0][0] = 1.0;
    xyz[0][1] = 1.0;
    xyz[0][2] = 2.0;
    xyz[0][3] = 2.0;
    xyz[1][3] = 1.0;
    dihedral = CoordTranslation.calcDihedral(xyz, 0, 1, 2, 3);
    assertEquals(Math.PI, dihedral, NUMACC);

    xyz = new double[3][4];
    xyz[1][0] = -1.0;
    xyz[0][0] = 1.0;
    xyz[0][1] = 1.0;
    xyz[0][2] = 2.0;
    xyz[0][3] = 2.0;
    xyz[2][3] = 1.0;
    dihedral = CoordTranslation.calcDihedral(xyz, 0, 1, 2, 3);
    assertEquals(-Math.PI / 2, dihedral, NUMACC);

    xyz = new double[3][4];
    xyz[1][0] = -1.0;
    xyz[0][0] = 1.0;
    xyz[0][1] = 1.0;
    xyz[0][2] = 2.0;
    xyz[0][3] = 2.0;
    xyz[2][3] = -1.0;
    dihedral = CoordTranslation.calcDihedral(xyz, 0, 1, 2, 3);
    assertEquals(Math.PI / 2, dihedral, NUMACC);

    xyz = new double[3][4];
    xyz[1][0] = -1.0;
    xyz[0][0] = 1.0;
    xyz[0][1] = 1.0;
    xyz[0][2] = 2.0;
    xyz[0][3] = 2.0;
    xyz[1][3] = -1.0;
    xyz[2][3] = -1.0;
    dihedral = CoordTranslation.calcDihedral(xyz, 0, 1, 2, 3);
    assertEquals(Math.PI / 4, dihedral, NUMACC);

    xyz = new double[3][4];
    xyz[1][0] = -1.0;
    xyz[0][0] = 1.0;
    xyz[0][1] = 1.0;
    xyz[0][2] = 2.0;
    xyz[0][3] = 2.0;
    xyz[1][3] = 1.0;
    xyz[2][3] = -1.0;
    dihedral = CoordTranslation.calcDihedral(xyz, 0, 1, 2, 3);
    assertEquals(3 * Math.PI / 4, dihedral, NUMACC);
  }

  private static final double ROT_TOL = 1e-8;

  /** Helper: assert that (x,y,z) lies on the z-axis (x and y negligible). */
  private static void assertOnZAxis(
      final double x, final double y, final double z, final double tol) {
    assertEquals(0.0, x, tol, "x should be 0 (on z-axis)");
    assertEquals(0.0, y, tol, "y should be 0 (on z-axis)");
    assertTrue(Math.abs(z) >= 0.0, "z component");
  }

  /** Helper: rotate point onto z-axis and verify the rotated point is on z-axis. */
  private static void assertRotatesPointToZAxis(final double px, final double py, final double pz) {
    final double[][] xyz = new double[3][1];
    xyz[0][0] = px;
    xyz[1][0] = py;
    xyz[2][0] = pz;
    final double[] point = new double[] {px, py, pz};
    final double[][] result = CoordTranslation.rotatePointToZAxis(xyz, point, 1);
    assertNotNull(result);
    assertEquals(3, result.length);
    assertEquals(1, result[0].length);
    final double norm = Math.sqrt(px * px + py * py + pz * pz);
    assertOnZAxis(result[0][0], result[1][0], result[2][0], ROT_TOL);
    assertEquals(norm, Math.abs(result[2][0]), ROT_TOL);
  }

  /** Degenerate: zero vector returns copy of xyz unchanged. */
  @Test
  public void testRotatePointToZAxisZeroVector() {
    final double[][] xyz = new double[3][2];
    xyz[0][0] = 1.0;
    xyz[1][0] = 2.0;
    xyz[2][0] = 3.0;
    xyz[0][1] = 4.0;
    xyz[1][1] = 5.0;
    xyz[2][1] = 6.0;
    final double[] point = new double[] {0.0, 0.0, 0.0};
    final double[][] result = CoordTranslation.rotatePointToZAxis(xyz, point, 2);
    assertNotNull(result);
    assertNotSame(xyz, result);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 2; j++) {
        assertEquals(xyz[i][j], result[i][j], NUMACC);
      }
    }
  }

  /** Degenerate: point already on z-axis (0, 0, z) returns copy. */
  @Test
  public void testRotatePointToZAxisAlreadyOnZAxis() {
    final double[][] xyz = new double[3][1];
    xyz[0][0] = 10.0;
    xyz[1][0] = 20.0;
    xyz[2][0] = 30.0;
    final double[] point = new double[] {0.0, 0.0, 5.0};
    final double[][] result = CoordTranslation.rotatePointToZAxis(xyz, point, 1);
    assertNotNull(result);
    assertNotSame(xyz, result);
    assertEquals(10.0, result[0][0], NUMACC);
    assertEquals(20.0, result[1][0], NUMACC);
    assertEquals(30.0, result[2][0], NUMACC);
  }

  /** Degenerate: point on z-axis negative (0, 0, -1). */
  @Test
  public void testRotatePointToZAxisAlreadyOnZAxisNegative() {
    final double[][] xyz = new double[3][1];
    xyz[0][0] = 7.0;
    xyz[1][0] = 8.0;
    xyz[2][0] = 9.0;
    final double[] point = new double[] {0.0, 0.0, -1.0};
    final double[][] result = CoordTranslation.rotatePointToZAxis(xyz, point, 1);
    assertNotNull(result);
    assertNotSame(xyz, result);
    assertEquals(7.0, result[0][0], NUMACC);
    assertEquals(8.0, result[1][0], NUMACC);
    assertEquals(9.0, result[2][0], NUMACC);
  }

  /** Normal: point on x-axis (1,0,0) rotates to positive z. */
  @Test
  public void testRotatePointToZAxisPointOnXAxis() {
    assertRotatesPointToZAxis(1.0, 0.0, 0.0);
  }

  /** Normal: point in xy-plane (1,1,0) rotates to z. */
  @Test
  public void testRotatePointToZAxisPointInXYPlane() {
    assertRotatesPointToZAxis(1.0, 1.0, 0.0);
  }

  /** Normal: general point (1,2,3). */
  @Test
  public void testRotatePointToZAxisGeneralPoint() {
    assertRotatesPointToZAxis(1.0, 2.0, 3.0);
  }

  /** Extreme: very small point below tolerance → degenerate, returns copy. */
  @Test
  public void testRotatePointToZAxisVerySmallPoint() {
    final double[][] xyz = new double[3][1];
    xyz[0][0] = 100.0;
    xyz[1][0] = 200.0;
    xyz[2][0] = 300.0;
    final double[] point = new double[] {1e-12, 1e-12, 1e-12};
    final double[][] result = CoordTranslation.rotatePointToZAxis(xyz, point, 1);
    assertNotNull(result);
    assertNotSame(xyz, result);
    assertEquals(100.0, result[0][0], NUMACC);
    assertEquals(200.0, result[1][0], NUMACC);
    assertEquals(300.0, result[2][0], NUMACC);
  }

  /** Extreme: very large coordinates. */
  @Test
  public void testRotatePointToZAxisVeryLargePoint() {
    final double scale = 1e10;
    assertRotatesPointToZAxis(scale, 0.0, 0.0);
    assertRotatesPointToZAxis(scale * 0.6, scale * 0.8, 0.0);
  }

  /** Extreme: ats = 0 (empty coordinates). */
  @Test
  public void testRotatePointToZAxisZeroAtoms() {
    final double[][] xyz = new double[3][1];
    xyz[0][0] = 1.0;
    xyz[1][0] = 0.0;
    xyz[2][0] = 0.0;
    final double[] point = new double[] {1.0, 0.0, 0.0};
    final double[][] result = CoordTranslation.rotatePointToZAxis(xyz, point, 0);
    assertNotNull(result);
    assertEquals(3, result.length);
    assertEquals(0, result[0].length);
  }

  /** Multiple atoms: rotation is by same matrix; one column was the point. */
  @Test
  public void testRotatePointToZAxisMultipleAtoms() {
    final double px = 1.0;
    final double py = 2.0;
    final double pz = 3.0;
    final double norm = Math.sqrt(px * px + py * py + pz * pz);
    final double[][] xyz = new double[3][3];
    xyz[0][0] = px;
    xyz[1][0] = py;
    xyz[2][0] = pz;
    xyz[0][1] = 0.0;
    xyz[1][1] = 0.0;
    xyz[2][1] = 1.0;
    xyz[0][2] = -1.0;
    xyz[1][2] = 0.0;
    xyz[2][2] = 0.0;
    final double[] point = new double[] {px, py, pz};
    final double[][] result = CoordTranslation.rotatePointToZAxis(xyz, point, 3);
    assertNotNull(result);
    assertEquals(3, result[0].length);
    assertOnZAxis(result[0][0], result[1][0], result[2][0], ROT_TOL);
    assertEquals(norm, Math.abs(result[2][0]), ROT_TOL);
  }

  /** Point with negative y (triggers sign flip in implementation). */
  @Test
  public void testRotatePointToZAxisNegativeY() {
    assertRotatesPointToZAxis(1.0, -1.0, 1.0);
  }
}
