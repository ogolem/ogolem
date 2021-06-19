/*
Copyright (c) 2021, D. Behrens, J. M. Dieterich, B. Hartke
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

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

import java.util.ArrayList;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.ogolem.generic.genericpool.Niche;

/**
 * Test Class for Coulomb-Matrix
 *
 * @author Dominik Behrens
 * @version 2021-07-06
 */
public class CoulombMatrixTest {
  private CoulombMatrixNicheComp comp;

  @Before
  public void setUp() {
    comp =
        new CoulombMatrixNicheComp(
            CoulombMatrixNicheComp.DEFAULTWIDTH, CoulombMatrixNicheComp.DEFAULTNUMBER);
  }

  @Test
  public void testSetupAllSameSpot() {
    Geometry sample = constructGoldCube(0, 0);
    Niche niche = comp.computeNiche(sample);
    Assert.assertEquals("cmat-0", niche.getID());
  }

  @Test
  public void testSetupWidthLower() {
    Geometry sample = constructGoldCube(5, -0.219);
    Niche niche = comp.computeNiche(sample);
    Assert.assertEquals("cmat-45", niche.getID());
  }

  @Test
  public void testSetupWidthUpper() {
    Geometry sample = constructGoldCube(5, 0.228);
    Niche niche = comp.computeNiche(sample);
    Assert.assertEquals("cmat-45", niche.getID());
  }

  @Test
  public void testSetupJustOutWidthLower() {
    Geometry sample = constructGoldCube(5, -0.220);
    Niche niche = comp.computeNiche(sample);
    Assert.assertEquals("cmat-44", niche.getID());
  }

  @Test
  public void testSetupJustOutWidthUpper() {
    Geometry sample = constructGoldCube(5, 0.229);
    CartesianCoordinates dfv = sample.getCartesians();
    Niche niche = comp.computeNiche(sample);
    Assert.assertEquals("cmat-44", niche.getID());
  }

  private Geometry constructGoldCube(double a, double b) {
    MoleculeConfig org = new MoleculeConfig(false);
    org.atomNumbers = new short[] {79};
    org.atomTypes = new String[] {"Au"};
    org.charges = new float[] {0};
    org.noOfAtoms = 1;
    org.spins = new short[] {0};
    org.refXYZ = new double[3][1];
    double[][] xyzs =
        new double[][] {
          {0, 0, 0}, {0, a, 0}, {0, 0, a}, {0, a, a},
          {a, 0, 0}, {a, a, 0}, {a, 0, a}, {a, a, a + b}
        };
    ArrayList<MoleculeConfig> mcs = new ArrayList<>(4);
    for (int i = 0; i < 8; i++) {
      MoleculeConfig mc = org.copy();
      mc.externalCOM = xyzs[i];
      mcs.add(mc);
    }
    GeometryConfig gc = new GeometryConfig();
    gc.geomMCs = mcs;
    gc.bonds = new SimpleBondInfo(8);
    gc.noOfParticles = 8;
    return new Geometry(gc);
  }
}
