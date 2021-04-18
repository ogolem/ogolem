/*
Copyright (c) 2010     , J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive;

import static java.lang.Math.sqrt;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import org.apache.commons.math3.util.FastMath;
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.helpers.Machine;

/**
 * Implements a GUPTA force field as used by Johnston et al. for cluster optimization.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public final class AdaptiveGUPTA extends AbstractAdaptiveBackend {

  // the ID
  private static final long serialVersionUID = (long) 20150727;

  private final int whichExp;
  private AdaptiveParameters params;
  private final double dMachinePrecision;
  private final double blowDistCut;
  private final boolean useCaching;
  private int[] paramOffsetCache;
  private double[] partsCache; // can be overwritten...

  private static final boolean DEBUG = false;

  public AdaptiveGUPTA(
      final boolean bIsInAdaptive,
      final AdaptiveParameters parameters,
      final double blowDistCutoff,
      final boolean useCache,
      final int expF) {
    this.dMachinePrecision = Machine.calcMachinePrecision();
    this.blowDistCut = blowDistCutoff;
    this.useCaching = useCache;
    this.whichExp = expF;

    if (!bIsInAdaptive) {
      if (parameters != null) {
        params = parameters;
      } else {
        System.out.println(
            "INFO: Performance penalty encountered, by "
                + "needing to read the parameters in again and again. Consider changing. :-)");
        final String sFile = "adaptive-ljff.param";
        String[] saData;
        try {
          saData = Input.ReadFile(sFile);
        } catch (Exception e) {
          System.err.println("ERROR: Couldn't read parameter file in: " + e.toString());
          saData = new String[0];
        }
        int iID = 0;
        if (DEBUG) {
          AdaptiveParameters paramsTemp = null;
          try {
            paramsTemp = new AdaptiveParameters(saData, iID);
          } catch (Exception e) {
            System.err.println("ERROR: Exception occured in creating parameter object.");
            e.printStackTrace(System.err);
          }
          params = paramsTemp;
        } else {
          params = new AdaptiveParameters(saData, iID);
        }
      }
    } else {
      // we simply do not need them
      params = null;
    }

    // explicitly initialize the caches
    this.paramOffsetCache = null;
  }

  private AdaptiveGUPTA(final AdaptiveGUPTA orig) {
    this.dMachinePrecision = orig.dMachinePrecision;
    this.blowDistCut = orig.blowDistCut;
    this.useCaching = orig.useCaching;
    this.whichExp = orig.whichExp;

    if (orig.params != null) {
      this.params = new AdaptiveParameters(orig.params);
    } else {
      this.params = null;
    }

    if (orig.paramOffsetCache != null) {
      this.paramOffsetCache = orig.paramOffsetCache.clone();
    } else {
      this.paramOffsetCache = null;
    }
  }

  @Override
  public AdaptiveGUPTA copy() {
    return new AdaptiveGUPTA(this);
  }

  @Override
  public String getMethodID() {
    return "Adaptive GUPTA";
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {

    this.params = params;

    if (partsCache == null) {
      partsCache = new double[cartes.getNoOfMolecules()];
    }

    assert (partsCache.length == cartes.getNoOfAtoms());

    final double d =
        energyCalculation(
            (long) -1,
            1,
            cartes.getAll1DCartes(),
            cartes.getAllAtomTypes(),
            cartes.getAllAtomNumbers(),
            cartes.getAllAtomsPerMol(),
            partsCache,
            cartes.getNoOfAtoms(),
            cartes.getAllCharges(),
            cartes.getAllSpins(),
            bonds,
            cartes.containedEnvType() == CartesianCoordinates.ENVTYPE.RIGID);

    this.params = null;

    return d;
  }

  @Override
  public double energyCalculation(
      long lID,
      int iIteration,
      double[] daXYZ1D,
      String[] saAtomTypes,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      int iNoOfAtoms,
      float[] faCharges,
      short[] iaSpins,
      final BondInfo bonds,
      final boolean hasRigidEnv) {

    // zero the partial contributions
    for (int i = 0; i < energyparts.length; i++) {
      energyparts[i] = 0.0;
    }

    final String[] saAtoms = saAtomTypes;

    if (useCaching && paramOffsetCache == null) {
      // initialize caches
      this.paramOffsetCache = new int[(iNoOfAtoms * iNoOfAtoms - iNoOfAtoms) / 2];
      int counter = 0;
      for (int i = 0; i < iNoOfAtoms - 1; i++) {
        for (int j = i + 1; j < iNoOfAtoms; j++) {
          final String s = buildKey(saAtoms, i, j, params);
          paramOffsetCache[counter] = params.getStartPointForKey(s);
          counter++;
        }
      }
    }

    final double[] xyz1 = new double[3];
    double[] daParams = params.getAllParamters();
    int offset = 0;
    double dEnergy = 0.0;
    int counter = -1;

    int lastOffset = 0;
    for (int i = 0; i < atsPerMol.length - 1; i++) {
      lastOffset += atsPerMol[i];
    }

    final int firstLoopAtomNo = (hasRigidEnv) ? lastOffset - 1 : iNoOfAtoms - 1;
    OuterLoop:
    for (int i = 0; i < firstLoopAtomNo; i++) {

      final double rad1 = AtomicProperties.giveRadius(atomNos[i]);
      xyz1[0] = daXYZ1D[i];
      xyz1[1] = daXYZ1D[i + iNoOfAtoms];
      xyz1[2] = daXYZ1D[i + 2 * iNoOfAtoms];

      for (int j = i + 1; j < iNoOfAtoms; j++) {
        counter++;

        final double rad2 = AtomicProperties.giveRadius(atomNos[j]);

        // calculate pair distance
        final double dDistX = xyz1[0] - daXYZ1D[j];
        final double dDistY = xyz1[1] - daXYZ1D[j + iNoOfAtoms];
        final double dDistZ = xyz1[2] - daXYZ1D[j + 2 * iNoOfAtoms];
        final double dDist = sqrt(dDistX * dDistX + dDistY * dDistY + dDistZ * dDistZ);

        // check for cutoff
        if (dDist > blowDistCut * (rad1 + rad2)) {
          continue;
        }

        // get the params
        if (!useCaching) {
          daParams = getParamsForPair(saAtoms, i, j, params);
        } else {
          offset = paramOffsetCache[counter];
        }

        if (daParams == null) {
          continue;
        }
        if (offset == -1) {
          continue;
        }

        final double dExpInside = (dDist / daParams[offset + 2]) - 1.0;

        /*
         * calculate the repulsive pair term
         * A(a,b)   = params[0]
         * p(a,b)   = params[1]
         * r0(a,b)  = params[2]
         */
        final double dRepulsivePair = daParams[offset] * exp(-daParams[offset + 1] * dExpInside);

        /*
         * calculate the attractive "effective many-body" term
         * chi(a,b) = params[3]
         * q(a,b)   = params[4]
         */
        final double dAttrManyBody =
            sqrt(
                daParams[offset + 3]
                    * daParams[offset + 3]
                    * exp(-2.0 * daParams[offset + 4] * dExpInside));

        // add them up
        final double tmp = dRepulsivePair - dAttrManyBody;
        dEnergy += tmp;

        energyparts[i] += tmp;
        energyparts[j] += tmp;

        if (DEBUG) {
          System.out.println(
              " DEBUG: GUPTA energy "
                  + i
                  + " \t"
                  + j
                  + " \t"
                  + (dRepulsivePair - dAttrManyBody)
                  + " \t"
                  + dDist);
        }
      }
    }

    return dEnergy;
  }

  @Override
  public void gradientCalculation(
      long lID,
      int iIteration,
      double[] daXYZ1D,
      String[] saAtomTypes,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      int iNoOfAtoms,
      float[] faCharges,
      short[] iaSpins,
      final BondInfo bonds,
      final Gradient gradient,
      final boolean hasRigidEnv) {

    // zero the partial contributions
    for (int i = 0; i < energyparts.length; i++) {
      energyparts[i] = 0.0;
    }

    if (useCaching && paramOffsetCache == null) {
      // initialize caches
      this.paramOffsetCache = new int[(iNoOfAtoms * iNoOfAtoms - iNoOfAtoms) / 2];
      int counter = 0;
      for (int i = 0; i < iNoOfAtoms - 1; i++) {
        for (int j = i + 1; j < iNoOfAtoms; j++) {
          final String s = buildKey(saAtomTypes, i, j, params);
          paramOffsetCache[counter] = params.getStartPointForKey(s);
          counter++;
        }
      }
    }

    gradient.zeroGradient();
    final double[][] daGradientMat = gradient.getTotalGradient();

    final double[] xyz1 = new double[3];
    double[] daParams = params.getAllParamters();
    int counter = -1;
    double energy = 0.0;
    int offset = 0;

    int lastOffset = 0;
    for (int i = 0; i < atsPerMol.length - 1; i++) {
      lastOffset += atsPerMol[i];
    }

    final int firstLoopAtomNo = (hasRigidEnv) ? lastOffset - 1 : iNoOfAtoms - 1;
    for (int i = 0; i < firstLoopAtomNo; i++) {

      final double rad1 = AtomicProperties.giveRadius(atomNos[i]);
      xyz1[0] = daXYZ1D[i];
      xyz1[1] = daXYZ1D[i + iNoOfAtoms];
      xyz1[2] = daXYZ1D[i + 2 * iNoOfAtoms];

      for (int j = i + 1; j < iNoOfAtoms; j++) {
        counter++;

        final double rad2 = AtomicProperties.giveRadius(atomNos[j]);

        // calculate pair distance
        final double dDistX = xyz1[0] - daXYZ1D[j];
        final double dDistY = xyz1[1] - daXYZ1D[j + iNoOfAtoms];
        final double dDistZ = xyz1[2] - daXYZ1D[j + 2 * iNoOfAtoms];
        final double dDist = sqrt(dDistX * dDistX + dDistY * dDistY + dDistZ * dDistZ);

        // check for cutoff
        if (dDist > blowDistCut * (rad1 + rad2)) {
          continue;
        }

        // get the params
        if (!useCaching) {
          daParams = getParamsForPair(saAtomTypes, i, j, params);
        } else {
          offset = paramOffsetCache[counter];
        }

        if (daParams == null) {
          continue;
        }
        if (offset == -1) {
          continue;
        }

        final double dExpInside = (dDist / daParams[offset + 2]) - 1.0;
        final double dExpT = exp(-daParams[offset + 1] * dExpInside);
        final double dNum = -daParams[offset] * daParams[offset + 1] * dExpT;
        final double dDenom = daParams[offset + 2] * dDist;
        final double dFac = dNum / dDenom;

        // derivatives of the repulsive contribution
        final double dGradRepulX = dFac * dDistX;
        final double dGradRepulY = dFac * dDistY;
        final double dGradRepulZ = dFac * dDistZ;

        final double dPSq = daParams[offset + 3] * daParams[offset + 3];
        final double dExpT2 = exp(-2 * daParams[offset + 4] * dExpInside);
        final double dNum2 = -dPSq * daParams[offset + 4] * dExpT2;
        final double dDenom2 = daParams[offset + 2] * dDist * sqrt(dPSq * dExpT2);
        final double dFac2 = dNum2 / dDenom2;

        // derivatives of the attractive contribution
        final double dGradAttrX = dFac2 * dDistX;
        final double dGradAttrY = dFac2 * dDistY;
        final double dGradAttrZ = dFac2 * dDistZ;

        final double gradX = (dGradRepulX - dGradAttrX);
        final double gradY = (dGradRepulY - dGradAttrY);
        final double gradZ = (dGradRepulZ - dGradAttrZ);

        daGradientMat[0][i] += gradX;
        daGradientMat[0][j] -= gradX;
        daGradientMat[1][i] += gradY;
        daGradientMat[1][j] -= gradY;
        daGradientMat[2][i] += gradZ;
        daGradientMat[2][j] -= gradZ;

        // energy
        final double dRepulsivePair = daParams[offset + 0] * dExpT;
        final double dAttrManyBody = sqrt(dPSq * dExpT2);

        // add them up
        final double tmp = dRepulsivePair - dAttrManyBody;
        energy += tmp;

        energyparts[i] += tmp;
        energyparts[j] += tmp;
      }
    }

    gradient.setTotalEnergy(energy);

    if (DEBUG) {
      // compare analytical and numerical gradient
      final Gradient numericalGrad =
          NumericalGradients.numericalGradient(
              lID,
              iIteration,
              daXYZ1D,
              saAtomTypes,
              atomNos,
              atsPerMol,
              energyparts,
              iNoOfAtoms,
              faCharges,
              iaSpins,
              bonds,
              this,
              hasRigidEnv);

      final double[][] daNumGrad = numericalGrad.getTotalGradient();

      for (int i = 0; i < daNumGrad.length; i++) {
        for (int j = 0; j < daNumGrad[0].length; j++) {
          if (Math.abs(daNumGrad[i][j] - daGradientMat[i][j]) >= 1E-7) {
            // we have a noticable difference
            System.err.println(
                "DEBUG: Analytical vs numerical gradient is: "
                    + daGradientMat[i][j]
                    + " vs "
                    + daNumGrad[i][j]);
          } else {
            System.out.println("DEBUG: Analytical vs numerical gradient was fine. :-)");
          }
        }
      }

      numericalGrad.setGradientTotal(daNumGrad);
      numericalGrad.setTotalEnergy(gradient.getTotalEnergy());

      // we return the numerical gradient
      gradient.copyDataIn(numericalGrad);
    }
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds,
      final double[] grad) {

    final int noOfAtoms = cartes.getNoOfAtoms();
    final int noOfParams = params.getNumberOfParamters();
    final double[][] xyz = cartes.getAllXYZCoord();
    final String[] atoms = cartes.getAllAtomTypes();
    final double[] daParams = params.getAllParamters();
    final short[] atomNos = cartes.getAllAtomNumbers();

    if (useCaching && paramOffsetCache == null) {
      // initialize caches
      this.paramOffsetCache = new int[(noOfAtoms * noOfAtoms - noOfAtoms) / 2];
      int counter = 0;
      for (int i = 0; i < noOfAtoms - 1; i++) {
        for (int j = i + 1; j < noOfAtoms; j++) {
          final String s = buildKey(atoms, i, j, params);
          paramOffsetCache[counter] = params.getStartPointForKey(s);
          counter++;
        }
      }
    }

    double energy = 0.0;
    int counter = -1;
    for (int i = 0; i < noOfAtoms - 1; i++) {

      final double rad1 = AtomicProperties.giveRadius(atomNos[i]);

      for (int j = i + 1; j < noOfAtoms; j++) {
        counter++;

        final double rad2 = AtomicProperties.giveRadius(atomNos[j]);

        // calculate pair distance
        final double dDistX = xyz[0][i] - xyz[0][j];
        final double dDistY = xyz[1][i] - xyz[1][j];
        final double dDistZ = xyz[2][i] - xyz[2][j];
        final double dDist = sqrt(dDistX * dDistX + dDistY * dDistY + dDistZ * dDistZ);

        // check for cutoff
        if (dDist > blowDistCut * (rad1 + rad2)) {
          continue;
        }

        // get the params
        int offset;
        if (!useCaching) {
          final String s = buildKey(atoms, i, j, params);
          offset = params.getStartPointForKey(s);
        } else {
          offset = paramOffsetCache[counter];
        }

        if (offset == -1) {
          continue;
        }

        /*
         * calculate the repulsive pair term
         * A(a,b)   = params[0]
         * p(a,b)   = params[1]
         * r0(a,b)  = params[2]
         */
        final double r0Sq = daParams[offset + 2] * daParams[offset + 2];
        final double dRepulsivePair =
            daParams[offset] * exp(-daParams[offset + 1] * ((dDist / daParams[offset + 2]) - 1.0));

        /*
         * calculate the attractive "effective many-body" term
         * chi(a,b) = params[3]
         * q(a,b)   = params[4]
         */
        final double chiSq = daParams[offset + 3] * daParams[offset + 3];

        final double dAttrManyBody =
            Math.sqrt(
                chiSq * exp(-2.0 * daParams[offset + 4] * ((dDist / daParams[offset + 2]) - 1.0)));

        // derivative of the parameters
        final double expT1 = exp(daParams[offset + 1] * (dDist / daParams[offset + 2] - 1));
        final double expT2 = exp(2 * daParams[offset + 4] * (dDist / daParams[offset + 2] - 1));
        final double sig1 = sqrt(chiSq / expT2);

        grad[offset] += 1 / expT1;
        grad[offset + 1] += -daParams[offset] * (dDist / daParams[offset + 2] - 1) / expT1;
        grad[offset + 2] +=
            dDist * daParams[offset] * daParams[offset + 1] / (r0Sq * expT1)
                - dDist * chiSq * daParams[offset + 4] / (r0Sq * expT2 * sig1);
        grad[offset + 3] += -daParams[offset + 3] / (expT2 * sig1);
        grad[offset + 4] += chiSq * (2 * dDist / daParams[offset + 2] - 2) / (2 * expT2 * sig1);

        // add them up
        energy += dRepulsivePair;
        energy -= dAttrManyBody;

        if (DEBUG) {
          System.out.println(" " + i + " \t" + j + " \t" + (dRepulsivePair - dAttrManyBody));
        }
      }
    }

    return energy;
  }

  @Override
  public double[][] minMaxBordersForParams(final AdaptiveParameters params) {

    final int iNoOfParams = params.getNumberOfParamters();
    final double[][] daBorders = new double[2][iNoOfParams];

    int i = 0;
    while (i < iNoOfParams) {

      // A(a,b) = params[0]
      daBorders[0][i] = -2;
      daBorders[1][i] = 2;
      i++;

      // p(a,b) = params[1]
      daBorders[0][i] = -20;
      daBorders[1][i] = 20;
      i++;

      // r0(a,b)= params[2]
      daBorders[0][i] = 0;
      daBorders[1][i] = 20;
      i++;

      // chi(a,b)= params[3]
      daBorders[0][i] = -2;
      daBorders[1][i] = 2;
      i++;

      // q(a,b)  = params[4]
      daBorders[0][i] = -20;
      daBorders[1][i] = 20;
      i++;
    }

    return daBorders;
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      final ArrayList<CartesianCoordinates> refCartes, final String sMethod) {

    // loop through all reference geometries and figure a list of non-redundant atom tyes out
    final LinkedList<String> llAtoms = new LinkedList<>();

    final Iterator<CartesianCoordinates> itRefGeoms = refCartes.iterator();

    while (itRefGeoms.hasNext()) {
      final CartesianCoordinates cartesTemp = itRefGeoms.next();
      final String[] saAtoms = cartesTemp.getAllAtomTypes();

      for (final String saAtom : saAtoms) {
        if (!llAtoms.contains(saAtom)) {
          // add it to the list
          llAtoms.add(saAtom);
        }
      }
    }

    int iNoOfPairs = 0;
    for (int i = 0; i < llAtoms.size(); i++) {
      for (int j = 0; j <= i; j++) {
        iNoOfPairs++;
      }
    }

    final String[] saPairs = new String[iNoOfPairs];
    final int[] iaParamsPerAt = new int[iNoOfPairs];
    final int iNoOfAtoms = llAtoms.size();
    int iParamSum = 0;

    int iCounter = 0;
    for (int i = 0; i < iNoOfAtoms; i++) {
      for (int j = i; j < iNoOfAtoms; j++) {
        final String sPair = "adaptivegupta:" + llAtoms.get(i) + llAtoms.get(j);
        saPairs[iCounter] = sPair;
        // always five
        iaParamsPerAt[iCounter] = 5;
        iParamSum += iaParamsPerAt[iCounter];
        iCounter++;
      }
    }

    final AdaptiveParameters paramStub =
        new AdaptiveParameters(iParamSum, -1, saPairs, iaParamsPerAt, sMethod);

    return paramStub;
  }

  private static double[] getParamsForPair(
      final String[] atoms, final int i, final int j, final AdaptiveParameters parameters) {

    double[] params;

    final StringBuilder sb = new StringBuilder(20);
    sb.append("adaptivegupta:");
    sb.append(atoms[i]);
    sb.append(atoms[j]);
    params = parameters.getParametersForKey(sb.toString());

    if (params != null) {
      return params;
    }

    final StringBuilder sb2 = new StringBuilder(20);
    sb2.append("adaptivegupta:");
    sb2.append(atoms[j]);
    sb2.append(atoms[i]);
    params = parameters.getParametersForKey(sb2.toString());

    if (params != null) {
      return params;
    } else {
      System.err.println("WARNING: No parameters for key " + sb.toString() + ". Returning null.");
      return null;
    }
  }

  private static String buildKey(
      final String[] atoms, final int i, final int j, final AdaptiveParameters parameters) {

    final StringBuilder sb = new StringBuilder(20);
    sb.append("adaptivegupta:");
    sb.append(atoms[i]);
    sb.append(atoms[j]);
    final double[] params = parameters.getParametersForKey(sb.toString());

    if (params != null) {
      return sb.toString();
    }

    final StringBuilder sb2 = new StringBuilder(25);
    sb2.append("adaptivegupta:");
    sb2.append(atoms[j]);
    sb2.append(atoms[i]);
    // needs to be this one now. in good and in bad...
    return sb2.toString();
  }

  private double exp(final double x) {
    switch (whichExp) {
      case 0:
        return Math.exp(x);
      case 1:
        return org.ogolem.math.FastFunctions.fastExp(x);
      case 2:
        return org.ogolem.math.FastFunctions.fastExp2(x);
      case 3:
        return org.ogolem.math.FastFunctions.fastCorrExp(x);
      case 4:
        return org.ogolem.math.FastFunctions.lim128Exp(x);
      case 5:
        return org.ogolem.math.FastFunctions.lim256Exp(x);
      case 6:
        return org.ogolem.math.FastFunctions.lim512Exp(x);
      case 7:
        return org.ogolem.math.FastFunctions.lim1024Exp(x);
      case 8:
        return org.ogolem.math.FastFunctions.lim2048Exp(x);
      case 9:
        return org.ogolem.math.FastFunctions.lim4096Exp(x);
      case 10:
        return FastMath.exp(x);
      default:
        return Math.exp(x);
    }
  }
}
