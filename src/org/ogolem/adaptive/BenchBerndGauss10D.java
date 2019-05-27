/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;

/**
 * Bernds gaussian benchmark function in 10D.
 * @author Johannes Dieterich
 * @version 2015-07-27
 */
final class BenchBerndGauss10D extends AbstractAdaptivable{

    private static final long serialVersionUID = (long) 20150727;
    private final double[][] daGaussCoords;

    private final double[] daGaussWidth;

    private final double[] daGaussWeight;

    private static final boolean bDebug = false;

    BenchBerndGauss10D(){

        String[] saGaussians;

        try{
            saGaussians = Input.ReadFile("gaussians.in");
        } catch(Exception e){
            System.err.println("ERROR: Couldn't read in gaussian aux file. " + e.toString());
            saGaussians = null;
        }

        this.daGaussCoords = new double[10][2000];
        this.daGaussWidth  = new double[2000];
        this.daGaussWeight = new double[2000];

        if (saGaussians != null) {
            int iLine = 0;
            for (int i = 0; i < 2000; i++) {

                String sTemp = saGaussians[iLine].trim();

                try {
                    // first weight, then width
                    final int iIndex = sTemp.indexOf(" ");
                    final String sWeight = sTemp.substring(0, iIndex);
                    daGaussWeight[i] = Double.parseDouble(sWeight);

                    final String sWidth = sTemp.substring(iIndex).trim();
                    daGaussWidth[i] = Double.parseDouble(sWidth);

                } catch (Exception e) {
                    System.err.println("ERROR: Some error occured during parsing gaussian width and/or weight. " + e.toString());
                }

                iLine++;

                // now ten coordinates
                for (int j = 0; j < 10; j++) {
                    String sTemp2 = saGaussians[iLine].trim();
                    try {
                        double d = Double.parseDouble(sTemp2);
                        daGaussCoords[j][i] = d;
                    } catch (Exception e) {
                        System.err.println("ERROR: Couldn't parse gaussian coordinate" + e.toString());
                    }
                    iLine++;
                }
            }
        } else {
            System.err.println("ERROR: We DO NOT have ANY values for the benchmark " +
                    "due to a lacking gaussian.in file. This run makes NO sense.");
        }
    }

    @Override
    public BenchBerndGauss10D clone(){
        return new BenchBerndGauss10D();
    }

    @Override
    public double energyOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds){

        final double[] daParams = params.getAllParamters();

        double dValue = 0.0;

        for(int i = 0; i < 2000; i++){

            double r = 0.0;
            for(int j = 0; j < 10; j++){
                final double d = daParams[j] - daGaussCoords[j][i];
                r += d*d;
            }

            dValue += daGaussWeight[i] * Math.exp(-daGaussWidth[i] * r);

            for(int j = 0; j < 10; j++){
                if(daParams[j] > 10.0 || daParams[j] < 0.0){
                    // take the cutoff potential
                    dValue += daParams[j]*daParams[j];
                }
            }
        }

        return dValue;
    }

    @Override
    public double gradientOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds,
            final double[] daGrad){
        // we of course do not actually use the cartesian coordinates, just the
        // parameters

        final int iNoOfParams = params.getNumberOfParamters();

        final double[] daParams = params.getAllParamters();

        double totEnergy = 0.0;
        for (int i = 0; i < 2000; i++) {

            double r = 0.0;
            for (int j = 0; j < 10; j++) {
                final double d = daParams[j] - daGaussCoords[j][i];
                r += d*d;
            }

            final double t1 = daGaussWeight[i] * Math.exp(-daGaussWidth[i] * r);
            for (int j = 0; j < 10; j++) {
                daGrad[j] += -2.0 * daGaussWidth[i] * (daParams[j] - daGaussCoords[j][i]) * t1;
            }

            for (int j = 0; j < 10; j++) {
                if (daParams[j] > 10.0 || daParams[j] < 0.0) {
                    // take the cutoff potential
                    daGrad[j] += daParams[j]*daParams[j];
                }
            }
            
            totEnergy += t1;
        }

        if(bDebug){
            // calculate a numerical gradient and compare it to the analytical one
            final double[] numGrad = new double[iNoOfParams];
            final double numE = NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, numGrad);

            for(int i = 0; i < numGrad.length; i++){
                if(Math.abs(numGrad[i]-daGrad[i]) >= 1E-7){
                    // we have a noticable difference
                    System.err.println("DEBUG: Analytical vs numerical gradient is: "
                            + daGrad[i] + " vs " + numGrad[i]);
                }
            }
        }

        return totEnergy;
    }

    @Override
    public AdaptiveParameters createInitialParameterStub(final ArrayList<CartesianCoordinates> refCartes,
            final String sMethod){

        final String[] saAtoms = {"XX"};
        final int[] iaParamsPerAt = {10};
        final AdaptiveParameters paramStub = new AdaptiveParameters(10, -1, saAtoms, iaParamsPerAt, sMethod);

        return paramStub;
    }

    @Override
    public double[][] minMaxBordersForParams(final AdaptiveParameters params){

        final double[][] daMinMax = new double[2][10];

        for(int i = 0; i < 10; i++){
            daMinMax[0][i] = 0.0;
            daMinMax[1][i] = 10.0;
        }

        return daMinMax;
    }
}
