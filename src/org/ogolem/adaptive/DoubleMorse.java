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
package org.ogolem.adaptive;

import java.util.ArrayList;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.helpers.Machine;

/**
 * A dirty hack for Mr Oxford supplying him with a double Morse potential. ATTENTION: This uses fake
 * cartesian coordinates. It expects two atoms (dummies), where the "x-value" of the cartesians is
 * one of the two degrees of freedom.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
final class DoubleMorse extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20150727;

  // this is currently not used!
  private final boolean bSymmetric;

  private final double dMax;

  private final double dMachinePrecision;

  private static final boolean bDebugGradient = false;

  DoubleMorse(final boolean bSymPotential, final double dDMax) {
    this.bSymmetric = bSymPotential;
    this.dMax = dDMax;
    this.dMachinePrecision = Machine.calcMachinePrecision();
  }

  @Override
  public DoubleMorse copy() {
    return new DoubleMorse(bSymmetric, dMax);
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {

    final double[][] daXYZ = cartes.getAllXYZCoordsCopy();

    final double[] daParams = params.getAllParamters();

    // ATTENTION: this is hard-coded and ugly!
    final double dRho = daXYZ[0][0];
    final double dDelta = daXYZ[0][1];

    // now calculate the functional form
    final double dF1 =
        (1.0 + (daParams[0] * Math.pow(dRho, daParams[1]) * Math.exp(daParams[2] * dRho)))
            * daParams[3];
    final double dF2 =
        (1.0 + (daParams[17] * Math.pow(dRho, daParams[18]) * Math.exp(daParams[19] * dRho)))
            * daParams[20];

    final double dG1 = daParams[4] + daParams[5] * dRho + daParams[6] * Math.pow(dRho, 2);
    final double dG2 = daParams[21] + daParams[22] * dRho + daParams[23] * Math.pow(dRho, 2);

    final double dH1 =
        ((daParams[7] + daParams[8] * dRho + daParams[9] * Math.pow(dRho, 2))
                - Math.log(daParams[10] + daParams[11] * dRho + daParams[12] * Math.pow(dRho, 2)))
            / (daParams[13] + daParams[14] * dRho + daParams[15] * Math.pow(dRho, 2));

    final double dH2 =
        ((daParams[24] + daParams[25] * dRho + daParams[26] * Math.pow(dRho, 2))
                - Math.log(daParams[27] + daParams[28] * dRho + daParams[29] * Math.pow(dRho, 2)))
            / (daParams[30] + daParams[31] * dRho + daParams[32] * Math.pow(dRho, 2));

    final double dC1 = daParams[16];
    final double dC2 = daParams[33];

    final double dShift = daParams[34];

    // put it together
    final double dEnergy =
        dF1 * (Math.pow((1.0 - Math.exp(dG1 * (dH1 - dDelta))), 2) + dC1)
            + dF2 * (Math.pow((1.0 - Math.exp(dG2 * (dDelta - (dMax - dH2)))), 2) + dC2)
            + dShift;

    return dEnergy;
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds,
      final double[] daGrad) {

    final double[][] daXYZ = cartes.getAllXYZCoordsCopy();

    // ATTENTION: this is hard-coded and ugly!
    final double dRho = daXYZ[0][0];
    final double dDelta = daXYZ[0][1];

    final int iNoOfParams = params.getNumberOfParamters();

    final double[] daParams = params.getAllParamters();

    final double dF1 =
        (1.0 + (daParams[0] * Math.pow(dRho, daParams[1]) * Math.exp(daParams[2] * dRho)))
            * daParams[3];
    final double dF2 =
        (1.0 + (daParams[17] * Math.pow(dRho, daParams[18]) * Math.exp(daParams[19] * dRho)))
            * daParams[20];

    final double dG1 = daParams[4] + daParams[5] * dRho + daParams[6] * Math.pow(dRho, 2);
    final double dG2 = daParams[21] + daParams[22] * dRho + daParams[23] * Math.pow(dRho, 2);

    final double dH1 =
        ((daParams[7] + daParams[8] * dRho + daParams[9] * Math.pow(dRho, 2))
                - Math.log(daParams[10] + daParams[11] * dRho + daParams[12] * Math.pow(dRho, 2)))
            / (daParams[13] + daParams[14] * dRho + daParams[15] * Math.pow(dRho, 2));

    final double dH2 =
        ((daParams[24] + daParams[25] * dRho + daParams[26] * Math.pow(dRho, 2))
                - Math.log(daParams[27] + daParams[28] * dRho + daParams[29] * Math.pow(dRho, 2)))
            / (daParams[30] + daParams[31] * dRho + daParams[32] * Math.pow(dRho, 2));

    final double dC1 = daParams[16];
    final double dC2 = daParams[33];

    final double dShift = daParams[34];

    // put it together
    final double dEnergy =
        dF1 * (Math.pow((1.0 - Math.exp(dG1 * (dH1 - dDelta))), 2) + dC1)
            + dF2 * (Math.pow((1.0 - Math.exp(dG2 * (dDelta - (dMax - dH2)))), 2) + dC2)
            + dShift;

    // some more helpers
    final double dW1 = 1.0 - Math.exp(dG1 * dH1 - dDelta * dG1);
    final double dW2 = 1.0 - Math.exp(dG2 * dH2 + (dDelta - dMax) * dG2);

    final double dV1 = Math.pow(dW1, 2) + daParams[16];
    final double dV2 = Math.pow(dW2, 2) + daParams[33];

    // do not include the shifting derivative in the loop
    for (int i = 0; i < iNoOfParams - 1; i++) {

      double dDF1;
      switch (i) {
        case 0:
          dDF1 = daParams[3] * Math.pow(dRho, daParams[1]) * Math.exp(daParams[2] * dRho);
          break;
        case 1:
          dDF1 =
              daParams[0]
                  * daParams[3]
                  * Math.exp(daParams[2] * dRho)
                  * Math.log(dRho)
                  * Math.exp(daParams[1] * Math.log(dRho));
          break;
        case 2:
          dDF1 =
              daParams[0]
                  * daParams[3]
                  * Math.pow(dRho, (daParams[1] + 1.0))
                  * Math.exp(daParams[2] * dRho);
          break;
        case 3:
          dDF1 = 1.0 + daParams[0] * Math.pow(dRho, daParams[1]) * Math.exp(daParams[2] * dRho);
          break;
        default:
          dDF1 = 0.0;
          break;
      }

      double dDF2;
      switch (i) {
        case 17:
          dDF2 = daParams[20] * Math.pow(dRho, daParams[18]) * Math.exp(daParams[19] * dRho);
          break;
        case 18:
          dDF2 =
              daParams[17]
                  * daParams[20]
                  * Math.exp(daParams[19] * dRho)
                  * Math.log(dRho)
                  * Math.exp(daParams[18] * Math.log(dRho));
          break;
        case 19:
          dDF2 =
              daParams[17]
                  * daParams[20]
                  * Math.pow(dRho, (daParams[18] + 1.0))
                  * Math.exp(daParams[19] * dRho);
          break;
        case 20:
          dDF2 = 1.0 + daParams[17] * Math.pow(dRho, daParams[18]) * Math.exp(daParams[19] * dRho);
          break;
        default:
          dDF2 = 0.0;
          break;
      }

      double dDG1;
      switch (i) {
        case 4:
          dDG1 = 1.0;
          break;
        case 5:
          dDG1 = dRho;
          break;
        case 6:
          dDG1 = Math.pow(dRho, 2);
          break;
        default:
          dDG1 = 0.0;
          break;
      }

      double dDG2;
      switch (i) {
        case 21:
          dDG2 = 1.0;
          break;
        case 22:
          dDG2 = dRho;
          break;
        case 23:
          dDG2 = Math.pow(dRho, 2);
          break;
        default:
          dDG2 = 0.0;
          break;
      }

      // some intermediate terms
      final double dA1 =
          1.0
              / ((daParams[10] + daParams[11] * dRho + daParams[12] * Math.pow(dRho, 2))
                  * (daParams[13] + daParams[14] * dRho + daParams[15] * Math.pow(dRho, 2)));
      final double dA2 =
          1.0
              / ((daParams[27] + daParams[28] * dRho + daParams[29] * Math.pow(dRho, 2))
                  * (daParams[30] + daParams[31] * dRho + daParams[32] * Math.pow(dRho, 2)));

      final double dB1 =
          (Math.log(daParams[10] + daParams[11] * dRho + daParams[12] * Math.pow(dRho, 2))
                  - (daParams[7] + daParams[8] * dRho + daParams[9] * Math.pow(dRho, 2)))
              / (daParams[13] + daParams[14] * dRho + daParams[15] * Math.pow(dRho, 2));
      final double dB2 =
          (Math.log(daParams[27] + daParams[28] * dRho + daParams[29] * Math.pow(dRho, 2))
                  - (daParams[24] + daParams[25] * dRho + daParams[26] * Math.pow(dRho, 2)))
              / (daParams[30] + daParams[31] * dRho + daParams[32] * Math.pow(dRho, 2));

      double dDH1;
      switch (i) {
        case 7:
          dDH1 = 1.0 / (daParams[13] + daParams[14] * dRho + daParams[15] * Math.pow(dRho, 2));
          break;
        case 8:
          dDH1 = dRho / (daParams[13] + daParams[14] * dRho + daParams[15] * Math.pow(dRho, 2));
          break;
        case 9:
          dDH1 =
              Math.pow(dRho, 2)
                  / (daParams[13] + daParams[14] * dRho + daParams[15] * Math.pow(dRho, 2));
          break;
        case 10:
          dDH1 = -dA1;
          break;
        case 11:
          dDH1 = -dA1 * dRho;
          break;
        case 12:
          dDH1 = -dA1 * Math.pow(dRho, 2);
          break;
        case 13:
          dDH1 = dB1;
          break;
        case 14:
          dDH1 = dB1 * dRho;
          break;
        case 15:
          dDH1 = dB1 * Math.pow(dRho, 2);
          break;
        default:
          dDH1 = 0.0;
          break;
      }

      double dDH2;
      switch (i) {
        case 24:
          dDH2 = 1.0 / (daParams[30] + daParams[31] * dRho + daParams[32] * Math.pow(dRho, 2));
          break;
        case 25:
          dDH2 = dRho / (daParams[30] + daParams[31] * dRho + daParams[32] * Math.pow(dRho, 2));
          break;
        case 26:
          dDH2 =
              Math.pow(dRho, 2)
                  / (daParams[30] + daParams[31] * dRho + daParams[32] * Math.pow(dRho, 2));
          break;
        case 27:
          dDH2 = -dA2;
          break;
        case 28:
          dDH2 = -dA2 * dRho;
          break;
        case 29:
          dDH2 = -dA2 * Math.pow(dRho, 2);
          break;
        case 30:
          dDH2 = dB2;
          break;
        case 31:
          dDH2 = dB2 * dRho;
          break;
        case 32:
          dDH2 = dB2 * Math.pow(dRho, 2);
          break;
        default:
          dDH2 = 0.0;
          break;
      }

      final double dDW1 =
          -(dDG1 * dH1 + dG1 * dDH1 - dDelta * dDG1) * Math.exp(dG1 * dH1 - dDelta * dG1);
      final double dDW2 =
          -(dDG2 * dH2 + dG2 * dDH2 + (dDelta - dMax) * dDG2)
              * Math.exp(dG2 * dH2 + (dDelta - dMax) * dG2);

      double dDV1 = 2.0 * dW1 * dDW1;
      double dDV2 = 2.0 * dW2 * dDW2;

      if (i == 16) {
        dDV1 += 1.0;
      }

      if (i == 33) {
        dDV2 += 1.0;
      }

      final double dGradVal = (dDF1 * dV1 + dDV1 * dF1) + (dDF2 * dV2 + dDV2 * dF2);
      daGrad[i] = dGradVal;
    }

    // now also the shifting derivative
    daGrad[34] = 1.0;

    if (bDebugGradient) {
      final double[] numGrad = daGrad.clone();
      final double numE =
          NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, numGrad);

      if (numE != dEnergy)
        System.err.println(
            "DEBUG: Energy mismatch between routine and numerical gradient of " + (dEnergy - numE));

      // compare it with the analytical one
      for (int i = 0; i < numGrad.length; i++) {
        if (Math.abs(numGrad[i] - daGrad[i]) > 1e-6) {
          System.out.println(
              "DEBUG: difference in analytical vs numerical parameter gradient: "
                  + daGrad[i]
                  + " and "
                  + numGrad[i]
                  + ". This is a difference of "
                  + (Math.abs(numGrad[i] - daGrad[i]))
                  + ".");
        }
      }
    }

    return dEnergy;
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      final ArrayList<CartesianCoordinates> refCartes, final String sMethod) {

    final String[] saAtoms = {"XX"};
    final int[] iaParamsPerAt = {35};
    final AdaptiveParameters paramStub =
        new AdaptiveParameters(35, -1, saAtoms, iaParamsPerAt, sMethod);

    return paramStub;
  }

  @Override
  public double[][] minMaxBordersForParams(final AdaptiveParameters params) {

    final int iNoOfParams = params.getNumberOfParamters();

    if (iNoOfParams != 35) {
      System.err.println(
          "ERROR: Something went wrong! Your amount of parameters for the "
              + "double morse potential is unequal 35: "
              + iNoOfParams
              + ". FIX THIS!");
    }

    final double[][] daBorders = new double[2][35];

    daBorders[0][0] = 0.0;
    daBorders[1][0] = 5.0e6;

    daBorders[0][1] = -10.0;
    daBorders[1][1] = 0.0;

    daBorders[0][2] = -1.0;
    daBorders[1][2] = 0.0;

    daBorders[0][3] = -1.0;
    daBorders[1][3] = 1.0;

    daBorders[0][4] = -100.0;
    daBorders[1][4] = 100.0;

    daBorders[0][5] = -100.0;
    daBorders[1][5] = 100.0;

    daBorders[0][6] = -10.0;
    daBorders[1][6] = 10.0;

    daBorders[0][7] = -100.0;
    daBorders[1][7] = 100.0;

    daBorders[0][8] = -10.0;
    daBorders[1][8] = 10.0;

    daBorders[0][9] = -100.0;
    daBorders[1][9] = 100.0;

    daBorders[0][10] = -10.0;
    daBorders[1][10] = 10.0;

    daBorders[0][11] = -1.0e4;
    daBorders[1][11] = 1.0e4;

    daBorders[0][12] = -1.0e3;
    daBorders[1][12] = 1.0e3;

    daBorders[0][13] = -1.0e3;
    daBorders[1][13] = 1.0e3;

    daBorders[0][14] = -1.0e3;
    daBorders[1][14] = 1.0e3;

    daBorders[0][15] = -1.0e3;
    daBorders[1][15] = 1.0e3;

    daBorders[0][16] = 1.0e-3;
    daBorders[1][16] = 1.0e3;

    daBorders[0][17] = 1.0e-2;
    daBorders[1][17] = 1.0e2;

    daBorders[0][18] = -100.0;
    daBorders[1][18] = 100.0;

    daBorders[0][19] = 0.0;
    daBorders[1][19] = 5.0e6;

    daBorders[0][20] = -10.0;
    daBorders[1][20] = 0.0;

    daBorders[0][21] = -1.0;
    daBorders[1][21] = 0.0;

    daBorders[0][22] = -1.0;
    daBorders[1][22] = 1.0;

    daBorders[0][23] = -100.0;
    daBorders[1][23] = 100.0;

    daBorders[0][24] = -100.0;
    daBorders[1][24] = 100.0;

    daBorders[0][25] = -10.0;
    daBorders[1][25] = 10.0;

    daBorders[0][26] = -100.0;
    daBorders[1][26] = 100.0;

    daBorders[0][27] = -1.0e4;
    daBorders[1][27] = 1.0e4;

    daBorders[0][28] = -1.0e4;
    daBorders[1][28] = 1.0e4;

    daBorders[0][29] = -1.0e3;
    daBorders[1][29] = 1.0e3;

    daBorders[0][30] = -1.0e3;
    daBorders[1][30] = 1.0e3;

    daBorders[0][31] = 1.0e-3;
    daBorders[1][31] = 1.0e3;

    daBorders[0][32] = 1.0e-2;
    daBorders[1][32] = 1.0e2;

    daBorders[0][33] = -100.0;
    daBorders[1][33] = 100.0;

    daBorders[0][34] = -141.0;
    daBorders[1][34] = -139.0;

    /*        daBorders[0][0] = 98.5;
    daBorders[1][0] = 98.6;

    daBorders[0][1] = -0.8;
    daBorders[1][1] = -0.7;

    daBorders[0][2] = -0.4;
    daBorders[1][2] = -0.3;

    daBorders[0][3] = 0.2;
    daBorders[1][3] = 0.3;

    daBorders[0][4] = 4.0;
    daBorders[1][4] = 5.0;

    daBorders[0][5] = 1.3;
    daBorders[1][5] = 1.4;

    daBorders[0][6] = 0.0;
    daBorders[1][6] = 0.1;

    daBorders[0][7] = 13.0;
    daBorders[1][7] = 13.1;

    daBorders[0][8] = -0.3;
    daBorders[1][8] = -0.2;

    daBorders[0][9] = 0.0;
    daBorders[1][9] = 0.1;

    daBorders[0][10] = 3249.2;
    daBorders[1][10] = 3249.3;

    daBorders[0][11] = -575.7;
    daBorders[1][11] = -575.6;

    daBorders[0][12] = 57.4;
    daBorders[1][12] = 57.5;

    daBorders[0][13] = 28.8;
    daBorders[1][13] = 28.9;

    daBorders[0][14] = 1.5;
    daBorders[1][14] = 1.6;

    daBorders[0][15] = -0.2;
    daBorders[1][15] = -0.1;

    daBorders[0][16] = -4.4;
    daBorders[1][16] = -4.3;

    daBorders[0][17] = 241.8;
    daBorders[1][17] = 241.9;

    daBorders[0][18] = -1.2;
    daBorders[1][18] = -1.1;

    daBorders[0][19] = -0.4;
    daBorders[1][19] = -0.3;

    daBorders[0][20] = 0.1;
    daBorders[1][20] = 0.2;

    daBorders[0][21] = 3.2;
    daBorders[1][21] = 3.3;

    daBorders[0][22] = 1.3;
    daBorders[1][22] = 1.4;

    daBorders[0][23] = 0.0;
    daBorders[1][23] = 0.1;

    daBorders[0][24] = 8.0;
    daBorders[1][24] = 8.1;

    daBorders[0][25] = 0.3;
    daBorders[1][25] = 0.4;

    daBorders[0][26] = -0.1;
    daBorders[1][26] = 0.0;

    daBorders[0][27] = 5849.4;
    daBorders[1][27] = 5849.5;

    daBorders[0][28] = -2026.3;
    daBorders[1][28] = -2026.2;

    daBorders[0][29] = 610.1;
    daBorders[1][29] = 610.2;

    daBorders[0][30] = 6.2;
    daBorders[1][30] = 6.3;

    daBorders[0][31] = -1.5;
    daBorders[1][31] = -1.4;

    daBorders[0][32] = 0.0;
    daBorders[1][32] = 0.1;

    daBorders[0][33] = 4.3;
    daBorders[1][33] = 4.4;

    daBorders[0][34] = -140.2;
    daBorders[1][34] = -140.1;*/

    return daBorders;
  }
}
