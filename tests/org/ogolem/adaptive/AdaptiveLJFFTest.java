/*
Copyright (c) 2020-2021, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive;

import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import org.junit.jupiter.api.Test;
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.LennardJonesFFTest;
import org.ogolem.core.SimpleBondInfo;

/**
 * @author Johannes Dieterich
 * @version 2021-07-17
 */
public class AdaptiveLJFFTest {

  private static final double ENERGYLJ38 = -0.0659856619926505;
  private static final double ENERGYLJ55 = -0.10594240151944964;
  private static final double NUMACC = 1.0e-8;

  private final AdaptiveLJFF adapLJFF;

  private final int[] atsPerMol38;
  private final int[] atsPerMol55;

  private final CartesianCoordinates lj38;
  private final CartesianCoordinates lj55;

  public AdaptiveLJFFTest() {

    this.atsPerMol38 = new int[38];
    for (int i = 0; i < 38; i++) {
      atsPerMol38[i] = 1;
    }

    this.atsPerMol55 = new int[55];
    for (int i = 0; i < 55; i++) {
      atsPerMol55[i] = 1;
    }

    CartesianCoordinates lj38 = null;
    CartesianCoordinates lj55 = null;

    try {
      lj38 =
          org.ogolem.core.Input.parseCartesFromFileData(
              LennardJonesFFTest.getLJ38MIN().split("\n"),
              38,
              atsPerMol38,
              new short[38],
              new float[38]);

      lj55 =
          org.ogolem.core.Input.parseCartesFromFileData(
              LennardJonesFFTest.getLJ55MIN().split("\n"),
              55,
              atsPerMol55,
              new short[55],
              new float[55]);
    } catch (Exception e) {
      fail(e.toString());
    }

    assert (lj38 != null);
    assert (lj55 != null);

    this.lj38 = lj38;
    this.lj55 = lj55;

    /*
    parameters to constructor:
        final boolean isInAdaptive, final boolean easyLBMix,
        final int startPot, final int endPot, final int potIncr,
        final AdaptiveParameters parameters, final boolean use3Body,
        final double closeCutBlow, final double distCutBlow,
        final boolean useCache, final boolean discard2Body
    */

    final AdaptiveLJFF dummy =
        new AdaptiveLJFF(
            true /* pretend to be in adaptive*/,
            false,
            6,
            12,
            6,
            null /* only then can we set params to null */,
            false,
            1.2,
            20.0,
            false /* no cache, dummy */,
            false);
    final ArrayList<CartesianCoordinates> cartes = new ArrayList<>();
    cartes.add(lj38);

    final AdaptiveParameters params = dummy.createInitialParameterStub(cartes, dummy.getMethodID());

    // now set the parameters to our regular LJ parameters for Ar
    final double sigma = AtomicProperties.giveLennardJonesSigma("Ar");
    final double epsilon = AtomicProperties.giveLennardJonesEpsilon("Ar");

    final String[] keys = params.getAllKeysCopy();
    assert (keys.length == 1); // there should only be one

    final int noParams = params.getAmountOfParametersForKey(keys[0]);
    assert (noParams == 3); // must be one epsilon and two sigma (same sigma)

    final double[] paramVals = params.getAllParamters();
    paramVals[0] = epsilon;
    paramVals[1] = -sigma; // to make it 12-6
    paramVals[2] = sigma;

    adapLJFF =
        new AdaptiveLJFF(
            false /* pretend to be in core*/,
            false,
            6,
            12,
            6,
            params /* real params */,
            false,
            1.2,
            20.0,
            false /* no cache, different LJ cluster sizes */,
            false);
  }

  @Test
  public void testEnergyCalculationAdapLJMins() {

    final BondInfo bonds38 = new SimpleBondInfo(38);
    final BondInfo bonds55 = new SimpleBondInfo(55);

    final double e38 =
        adapLJFF.energyCalculation(
            -1,
            0,
            lj38.getAll1DCartes(),
            lj38.getAllAtomTypes(),
            lj38.getAllAtomNumbers(),
            atsPerMol38,
            new double[38],
            lj38.getNoOfAtoms(),
            lj38.getAllCharges(),
            lj38.getAllSpins(),
            bonds38,
            false);

    assertEquals(ENERGYLJ38, e38, NUMACC);

    final double e55 =
        adapLJFF.energyCalculation(
            -1,
            0,
            lj55.getAll1DCartes(),
            lj55.getAllAtomTypes(),
            lj55.getAllAtomNumbers(),
            atsPerMol55,
            new double[55],
            lj55.getNoOfAtoms(),
            lj55.getAllCharges(),
            lj55.getAllSpins(),
            bonds55,
            false);

    assertEquals(ENERGYLJ55, e55, NUMACC);
  }

  @Test
  public void testGradientCalculationAdapLJMins() {

    final BondInfo bonds38 = new SimpleBondInfo(38);
    final BondInfo bonds55 = new SimpleBondInfo(55);

    final Gradient grad38 = new Gradient(3, 38);
    adapLJFF.gradientCalculation(
        -1,
        0,
        lj38.getAll1DCartes(),
        lj38.getAllAtomTypes(),
        lj38.getAllAtomNumbers(),
        atsPerMol38,
        new double[38],
        lj38.getNoOfAtoms(),
        lj38.getAllCharges(),
        lj38.getAllSpins(),
        bonds38,
        grad38,
        false);

    assertEquals(ENERGYLJ38, grad38.getTotalEnergy(), NUMACC);

    final double[][] gradient38 = grad38.getTotalGradient();

    double gradTot38 = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 38; j++) gradTot38 += gradient38[i][j];
    }

    assertEquals(0.0, gradTot38, NUMACC);

    final Gradient grad55 = new Gradient(3, 55);
    adapLJFF.gradientCalculation(
        -1,
        0,
        lj55.getAll1DCartes(),
        lj55.getAllAtomTypes(),
        lj55.getAllAtomNumbers(),
        atsPerMol55,
        new double[55],
        lj55.getNoOfAtoms(),
        lj55.getAllCharges(),
        lj55.getAllSpins(),
        bonds55,
        grad55,
        false);

    assertEquals(ENERGYLJ55, grad55.getTotalEnergy(), NUMACC);

    final double[][] gradient55 = grad55.getTotalGradient();

    double gradTot55 = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 55; j++) gradTot55 += gradient55[i][j];
    }

    assertEquals(0.0, gradTot55, NUMACC);
  }
}
