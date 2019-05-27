/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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
import java.util.Iterator;
import java.util.LinkedList;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.helpers.Machine;

/**
 * Provides a simple, yet adaptivable universal force field using electrostatic, Lennard-Jones and
 * Axilrod-Teller terms. It contains an additional energy shift to better reproduce ab-initio
 * reference calculations.
 * @author Johannes Dieterich
 * @version 2015-07-27
 */
public final class AdaptiveUFF extends AbstractAdaptiveBackend{

    // the ID
    private static final long serialVersionUID = (long) 20150727;

    private final double dMachinePrecision;

    private final AdaptiveParameters params;

    /**
     * Constructor for the use as an adaptive backend in the adaptive subpart.
     */
    public AdaptiveUFF(final boolean bIsInAdaptive){
        dMachinePrecision = Machine.calcMachinePrecision();

        if(!bIsInAdaptive){
            final String sFile = "adaptive-uff.param";
            String[] saData = null;
            try {
                saData = Input.ReadFile(sFile);
            } catch (Exception e) {
                System.err.println("ERROR: Couldn't read parameter file in: " + e.toString());
                saData = new String[0];
            }
            int iID = 0;
            params = new AdaptiveParameters(saData, iID);
        } else {
            // we simply do not need them
            params = null;
        }
    }

    @Override
    public AdaptiveUFF clone(){
        return new AdaptiveUFF( (params == null) );
    }

    @Override
    public String getMethodID(){
        return "Adaptive UFF";
    }

    @Override
    public void gradientCalculation(final long lID, final int iIteration, final double[] daXYZ1D,
            final String[] saAtomTypes, final short[] atomNos, final int[] atsPerMol, final double[] energyparts, final int iNoOfAtoms,
            final float[] faCharges, final short[] iaSpins, final BondInfo bonds, final Gradient gradient){

        final Gradient analyticalGrad = new Gradient();

        // the matrix for the gradient
        final double[][] daGradientMat = analyticalGrad.getTotalGradient();

        // calculate all gradient parts

        // the part of LJ can be more or less taken as such, but one should perhaps optimize the variable access a little
        // the difficult part will be the derivatives of the axilrod-teller term
        // the electrostatic one should then be again fair
        
        // TODO analytical gradient, for the time being numerical
        final Gradient numericalGrad = NumericalGradients.numericalGradient(lID, iIteration, daXYZ1D, saAtomTypes, atomNos, atsPerMol, energyparts, iNoOfAtoms, faCharges, iaSpins, bonds, this);
        gradient.copyDataIn(numericalGrad);
    }

    @Override
    public double energyCalculation(final long lID, final int iIteration, final double[] daXYZ1D,
            final String[] saAtomTypes, final short[] atomNos, final int[] atsPerMol, final double[] energyparts,
            final int iNoOfAtoms, final float[] faCharges, final short[] iaSpins, final BondInfo bonds){

        // this is a fake and just allowed since it is intermediate
        final int[] iaAtsPerMol = {iNoOfAtoms};
        final CartesianCoordinates cartes = new CartesianCoordinates(iNoOfAtoms, 1, iaAtsPerMol);
        cartes.setAllCharges(faCharges);
        cartes.setAllSpins(iaSpins);

        final double dEnergy = energyOfStructWithParams(cartes, params, (int) lID, bonds);

        return dEnergy;
    }

    @Override
    public double energyOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds){

        final int iNoOfAtoms = cartes.getNoOfAtoms();

        final double[][] daXYZ = cartes.getAllXYZCoordsCopy();

        final double[][] daDistXYZ = new double[iNoOfAtoms][iNoOfAtoms];

        final String[] saAtoms = cartes.getAllAtomTypes();

        // first: calculate all pair distances, since they are used more than once
        for(int i = 0; i < iNoOfAtoms; i++){
            for(int j = i; j < iNoOfAtoms; j++){
                daDistXYZ[i][j] =  Math.sqrt(Math.pow((daXYZ[0][i] - daXYZ[0][j]), 2.0)
                        + Math.pow((daXYZ[1][i] - daXYZ[1][j]), 2.0)
                        + Math.pow((daXYZ[2][i] - daXYZ[2][j]), 2.0));
                daDistXYZ[j][i] = daDistXYZ[i][j];
            }
        }

        /*
         * the (most likely) biggest contribution: electrostatic potential
         */
        final float[] faCharges = cartes.getAllCharges();

        double dElectroPot = 0.0;

        for(int i = 0; i < iNoOfAtoms - 1; i++){
            for (int j = i + 1; j < iNoOfAtoms; j++) {
                dElectroPot += (faCharges[i]*faCharges[j])/daDistXYZ[i][j];
            }
        }

        /*
         * three-body terms, axilrod-teller potential
         */
        double dAxilrodTellerPot = 0;

        for(int i = 0; i < iNoOfAtoms; i++){
            final double dCParam1 = params.getParametersForKey(saAtoms[i])[0];
            for(int j = 0; j < i-1; j++){
                final double dCParam2 = params.getParametersForKey(saAtoms[j])[0];
                for(int k = 0; k < j-1; k++){
                    final double dCParam3 = params.getParametersForKey(saAtoms[k])[0];
                    // arithmetric mean is used
                    final double dCParam = (dCParam1 + dCParam2 + dCParam3) / 3.0;

                    final double dDistIJ = daDistXYZ[i][j];
                    final double dDistJK = daDistXYZ[j][k];
                    final double dDistIK = daDistXYZ[i][k];

                    final double dCosAngle1 = (Math.pow(dDistIJ, 2) + Math.pow(dDistJK, 2) - Math.pow(dDistIK, 2))
                            / (2.0 * dDistIJ * dDistJK);
                    final double dCosAngle2 = (Math.pow(dDistIJ, 2) + Math.pow(dDistIK, 2) - Math.pow(dDistJK, 2))
                            / (2.0 * dDistIJ * dDistIK);
                    final double dCosAngle3 = (Math.pow(dDistJK, 2) + Math.pow(dDistIK, 2) - Math.pow(dDistIJ, 2))
                            / (2.0 * dDistJK * dDistIK);

                    dAxilrodTellerPot += dCParam * (1 + 3 * dCosAngle1 * dCosAngle2 * dCosAngle3) /
                            Math.pow((dDistIJ * dDistJK * dDistIK), 3);
                }
            }
        }

        /*
         * classical Lennard-Jones potential
         */
        double dLJPot = 0.0;
        double dPotEnergy = 0.0;
        for (int i = 0; i < iNoOfAtoms - 1; i++) {

            // epsilon and sigma
            final double dEpsilon1 = params.getParametersForKey(saAtoms[i])[1];
            final double dSigma1 = params.getParametersForKey(saAtoms[i])[2];

            for (int j = i + 1; j < iNoOfAtoms; j++) {

                // the Lennard-Jones parameters
                final double dEpsilon2 = params.getParametersForKey(saAtoms[j])[1];
                final double dSigma2 = params.getParametersForKey(saAtoms[j])[2];

                // do lorentz-berthelot mixing
                final double dEpsilon = Math.sqrt(dEpsilon1 * dEpsilon2);
                final double dSigma = 0.5 * (dSigma1 + dSigma2);

                // the cutoff distance
                final double dSeam = 0.64 * dSigma;

                // more constants... needed for cutting of the potential
                final double dConst1 = (4.0 * dEpsilon * (Math.pow((dSigma / dSeam), 12.0) - Math.pow((dSigma / dSeam), 6.0)) - 10000.0) / dSeam;
                final double dConst2 = 10000.0;

                final double dDist = daDistXYZ[i][j];

                if (true){//dDist > dSeam) {
                    final double dInvR = dSigma / dDist;
                    final double dInvRPow6 = Math.pow(dInvR, 6);
                    final double dInvRPow12 = Math.pow(dInvRPow6, 2);

                    /*
                     * The used Lennard-Jones potential is of the form
                     * V(r_{ij})= 4\cdot\epsilon\cdot \left[\left(\dfrac{\sigma}{r_{ij}}\right)^{12}-\left(\dfrac{\sigma}{r_{ij}}^{6}\right)\right].
                     */

                    dPotEnergy = 4.0 * dEpsilon * (dInvRPow12 - dInvRPow6);

                } else {
                    System.out.println("WARNING: Atoms are rather close together. Using cutoff potential.");
                    dPotEnergy = dConst1 * dDist + dConst2;
                }
                dLJPot += dPotEnergy;
            }
        }

        // a last contribution: a constant represented by the last parameter
        double dEnergyShift = 0.0;
        for (int i = 0; i < iNoOfAtoms; i++) {
            dEnergyShift += params.getParametersForKey(saAtoms[i])[3];
        }

        // add the three contributions up
        double dTotEnergy = dElectroPot + dAxilrodTellerPot + dLJPot + dEnergyShift;

        return dTotEnergy;
    }

    @Override
    public double gradientOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds,
            final double[] daGrad){

        final int iNoOfParams = params.getNumberOfParamters();

        final int iNoOfAtoms = cartes.getNoOfAtoms();

        final String[] saAtoms = cartes.getAllAtomTypes();

        final String[] saForAtoms = params.getForWhichAtoms();

        /*
         * anyway we first need all pairwise distances
         */
        final double[][] daXYZ = cartes.getAllXYZCoordsCopy();

        final double[][] daDistXYZ = new double[iNoOfAtoms][iNoOfAtoms];

        for(int i = 0; i < iNoOfAtoms; i++){
            for(int j = i; j < iNoOfAtoms; j++){
                daDistXYZ[i][j] =  Math.sqrt(Math.pow((daXYZ[0][i] - daXYZ[0][j]), 2.0)
                        + Math.pow((daXYZ[1][i] - daXYZ[1][j]), 2.0)
                        + Math.pow((daXYZ[2][i] - daXYZ[2][j]), 2.0));
                daDistXYZ[j][i] = daDistXYZ[i][j];
            }
        }

        int iParamCounter = 0;
        
        for(int iForAtom = 0; iForAtom < saForAtoms.length; iForAtom++){

            final String sCurrentAtom = saForAtoms[iForAtom];

            /*
             * first parameter: axilrod teller, this should be fairly easy
             * we need to keep in mind that we JUST want to have those three-bodies
             * where there is AT LEAST one atom of the kind we are looking at ATM
             */
            double dAxilrodTellerGrad = 0.0;

            for (int i = 0; i < iNoOfAtoms; i++) {
                for (int j = 0; j < i - 1; j++) {
                    for (int k = 0; k < j - 1; k++) {

                        final boolean bMatchIC = saAtoms[i].equalsIgnoreCase(sCurrentAtom);
                        final boolean bMatchJC = saAtoms[j].equalsIgnoreCase(sCurrentAtom);
                        final boolean bMatchKC = saAtoms[k].equalsIgnoreCase(sCurrentAtom);

                        double dFront = 0.0;

                        if (bMatchIC && bMatchJC && bMatchKC) {
                            // all atoms are the same AND of the correct kind
                            dFront = 1.0;
                        } else if (bMatchIC && bMatchJC && !bMatchKC) {
                            // i and j are of the right kind
                            dFront = 2. / 3.;
                        } else if (bMatchIC && bMatchKC && !bMatchJC) {
                            // i and k are of the right kind
                            dFront = 2. / 3.;
                        } else if (bMatchJC && bMatchKC && !bMatchIC) {
                            // j and k are of the right kind
                            dFront = 2. / 3.;
                        } else if (bMatchIC && !bMatchJC && !bMatchKC) {
                            // i is of the right kind
                            dFront = 1. / 3.;
                        } else if (!bMatchIC && bMatchJC && !bMatchKC) {
                            // j is of the right kind
                            dFront = 1. / 3.;
                        } else if (!bMatchIC && !bMatchJC && bMatchKC) {
                            // k is of the right kind
                            dFront = 1. / 3.;
                        } else {
                            // nothing correct
                            dFront = 0.0;
                        }

                        // multiply it with the second part
                        final double dDistIJ = daDistXYZ[i][j];
                        final double dDistJK = daDistXYZ[j][k];
                        final double dDistIK = daDistXYZ[i][k];

                        final double dCosAngle1 = (Math.pow(dDistIJ, 2) + Math.pow(dDistJK, 2) - Math.pow(dDistIK, 2)) / (2.0 * dDistIJ * dDistJK);
                        final double dCosAngle2 = (Math.pow(dDistIJ, 2) + Math.pow(dDistIK, 2) - Math.pow(dDistJK, 2)) / (2.0 * dDistIJ * dDistIK);
                        final double dCosAngle3 = (Math.pow(dDistJK, 2) + Math.pow(dDistIK, 2) - Math.pow(dDistIJ, 2)) / (2.0 * dDistJK * dDistIK);

                        dAxilrodTellerGrad += dFront * (1 + 3 * dCosAngle1 * dCosAngle2 * dCosAngle3) /
                                Math.pow((dDistIJ * dDistJK * dDistIK), 3);
                    }
                }
            }

            // put it to the correct spot
            daGrad[iParamCounter] = dAxilrodTellerGrad;
            iParamCounter++;

            /*
             * then the lennard-jones epsilon and sigma: again, consider the atom types
             * 2 parameters
             */
            double dLJEpsilonGrad = 0.0;
            double dLJSigmaGrad = 0.0;
            for(int i = 0; i < iNoOfAtoms -1; i++){
                final double dEpsParam1 = params.getParametersForKey(saAtoms[i])[1];
                final double dSigmaParam1 = params.getParametersForKey(saAtoms[i])[2];

                for(int j = i + 1; j < iNoOfAtoms; j++){
                    final double dEpsParam2 = params.getParametersForKey(saAtoms[j])[1];
                    final double dSigmaParam2 = params.getParametersForKey(saAtoms[i])[2];

                    final double dSigma = 0.5 * (dSigmaParam1 + dSigmaParam2);
                    final double dEpsilon = Math.sqrt(dEpsParam1 * dEpsParam2);

                    final boolean bMatchIC = saAtoms[i].equalsIgnoreCase(sCurrentAtom);
                    final boolean bMatchJC = saAtoms[j].equalsIgnoreCase(sCurrentAtom);

                    double dFrontEps = 0.0;
                    double dSigmaForGrad = 0.0;
                    double dSigmaFront = 0.0;

                    if(bMatchIC && bMatchJC){
                        dFrontEps = 1.0;
                        dSigmaForGrad = dSigmaParam1;
                        dSigmaFront = 1.0;
                    } else if(bMatchIC && !bMatchJC){
                        dFrontEps = Math.sqrt(dEpsParam2) / 2 * Math.sqrt(dEpsParam1);
                        dSigmaForGrad = (dSigmaParam2 + dSigmaParam1) / 2.0;
                        dSigmaFront = 0.5;
                    } else if(!bMatchIC && bMatchJC){
                        dFrontEps = Math.sqrt(dEpsParam1) / 2 * Math.sqrt(dEpsParam2);
                        dSigmaForGrad = (dSigmaParam1 + dSigmaParam2) / 2.0;
                        dSigmaFront = 0.5;
                    } else{
                        dFrontEps = 0.0;
                        dSigmaForGrad = 0.0;
                        dSigmaFront = 0.0;
                    }

                    // multiply it with the second half
                    final double dDist = daDistXYZ[i][j];
                    final double dInvR = dSigma / dDist;
                    final double dInvRPow6 = Math.pow(dInvR, 6);
                    final double dInvRPow12 = Math.pow(dInvRPow6, 2);

                    dLJEpsilonGrad += 4.0 * dFrontEps * (dInvRPow12 - dInvRPow6);

                    final double dFirst = dSigmaFront * 12 * Math.pow(dSigmaForGrad, 11) / Math.pow(dDist, 12);
                    final double dSecond = dSigmaFront * 6 * Math.pow(dSigmaForGrad, 5) / Math.pow(dDist, 6);

                    dLJSigmaGrad += 4.0 * dEpsilon * (dFirst - dSecond);
                }
            }

            // put it to the correct spot
            daGrad[iParamCounter] = dLJEpsilonGrad;
            iParamCounter++;
            
            daGrad[iParamCounter] = dLJSigmaGrad;
            iParamCounter++;

            /*
             * now the deviation of our shifting constant: 1.0
             * times the number of atoms
             */
            int iNoOfCorrectAtoms = 0;
            for(int i = 0; i < iNoOfAtoms; i++){
                if(saAtoms[i].equalsIgnoreCase(sCurrentAtom)){
                    iNoOfCorrectAtoms++;
                }
            }
            daGrad[iParamCounter] = 1.0 * iNoOfCorrectAtoms;
            iParamCounter++;
        }

        // TODO disable numerical gradient calculation
        final double[] daNumericalGrad = new double[iNoOfParams];
        final double numE = NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, daNumericalGrad);

        // compare analytical and numerical, for debugging purposes
        for(int i = 0; i < daGrad.length; i++){
            System.out.println("DEBUG: analytical vs numerical gradient: " + daNumericalGrad[i] + " resp " + daGrad[i]
                    + " at position " + i);
            if(Math.abs(daGrad[i] - daNumericalGrad[i]) > 1E-6){
                System.err.println("DEBUG: Problem in analytical vs numerical: " +
                        daGrad[i] + " resp " + daNumericalGrad[i] +
                        " at position " + i);
            }
        }

        //XXX
        final double dEnergy = energyOfStructWithParams(cartes, params, geomID, bonds);

        if(Math.abs(dEnergy - numE) > 1E-6){
            System.err.println("DEBUG: Problem in analytical vs numerical fitness: " +
                    dEnergy + " versus " + numE);
        }

        //TODO we are NOT sure about the analytical gradient!!!
        return dEnergy;
    }

    @Override
    public AdaptiveParameters createInitialParameterStub(final ArrayList<CartesianCoordinates> refCartes,
            final String sMethod){

        // loop through all reference geometries and figure a list of non-redundant atom tyes out
        final LinkedList<String> llAtoms = new LinkedList<>();

        final Iterator<CartesianCoordinates> itRefGeoms = refCartes.iterator();

        while(itRefGeoms.hasNext()){
            final CartesianCoordinates cartesTemp = itRefGeoms.next();
            final String[] saAtoms = cartesTemp.getAllAtomTypes();

            for (int i = 0; i < saAtoms.length; i++) {
                if (!llAtoms.contains(saAtoms[i])) {
                    // add it to the list
                    llAtoms.add(saAtoms[i]);
                }
            }
        }

        int iParamSum = 0;

        String[] saAtoms = new String[llAtoms.size()];
        int[] iaParamsPerAt = new int[llAtoms.size()];

        // for each atom we need a set of parameters
        for (int i = 0; i < saAtoms.length; i++) {
            String sTempAtom = llAtoms.get(i);
            iaParamsPerAt[i] = numberParamsRequiredForAtom(sTempAtom, sMethod);
            iParamSum += iaParamsPerAt[i];
            saAtoms[i] = sTempAtom;
        }

        final AdaptiveParameters paramStub = new AdaptiveParameters(iParamSum, -1, saAtoms, iaParamsPerAt, sMethod);

        return paramStub;

    }

    private int numberParamsRequiredForAtom(final String sAtomID, final String sMethod){
        /*
         * 1 for axilrod teller
         * 1 lj epsilon
         * 1 lj sigma
         * 1 energy shift
         * independent of the atom type
         */
        return 4;
    }

    @Override
    public double[][] minMaxBordersForParams(final AdaptiveParameters params){

        final double[][] daBorders = new double[2][params.getNumberOfParamters()];

        final String[] saAtoms = params.getForWhichAtoms();

        int iCounter = 0;

        for (int i = 0; i < saAtoms.length; i++) {
            // Axilrod-Teller C: uber-speculative (TM)
            daBorders[0][iCounter] = 0.00;
            daBorders[1][iCounter] = 100.00;
            iCounter++;

            // LJ epsilon: speculative (TM)
            daBorders[0][iCounter] = 0.00001;
            daBorders[1][iCounter] = 0.00090;
            iCounter++;

            // LJ sigma: speculative (TM)
            daBorders[0][iCounter] = 1.00;
            daBorders[1][iCounter] = 15.0;
            iCounter++;

            // energy shift: uber-speculative(TM)
            // basis are some calculations of neutral and charged single atoms with DFT, I do realize
            // that this cannot be considered as heavy testing but for the time being this is the
            // best that we can get
            daBorders[0][iCounter] = -10000.00;
            daBorders[1][iCounter] = 0.0;
            iCounter++;
        }
        //TODO C borders are difficult: adjust
        //TODO energy shift is difficult: adjust?
        return daBorders;
    }
}
