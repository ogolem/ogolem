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
 * Rastrigins benchmark function in nD. Global minimum: f(x)=0; xi=0.
 * @author Johannes Dieterich
 * @version 2015-07-27
 */
final class BenchRastrigin extends AbstractAdaptivable {

    private static final long serialVersionUID = (long) 20150727;
    public static final boolean DEBUG = false;

    private final int iDims;

    BenchRastrigin(final int dims){
        this.iDims = dims;
    }

    @Override
    public BenchRastrigin clone(){
        return new BenchRastrigin(iDims);
    }

    @Override
    public double energyOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds){

        final double[] daParams = params.getAllParamters();

        double dEnergy = 10.0*iDims;

        for(int i = 0; i < daParams.length; i++){
            if(daParams[i] > 5.12 || daParams[i] < -5.12){
                // take the cutoff potential
                dEnergy += 10.0*daParams[i]*daParams[i];
            } else {
                dEnergy += daParams[i]*daParams[i]-10.0*Math.cos(2.0*Math.PI*daParams[i]);
            }
        }

        return dEnergy;
    }

    @Override
    public double gradientOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds,
            final double[] daGrad){

        final double[] daParams = params.getAllParamters();

        final int iNoOfParams = params.getNumberOfParamters();

        double dEnergy = 10.0*iDims;
        for(int i = 0; i < iNoOfParams; i++){
            if(daParams[i] > 5.12 || daParams[i] < -5.12){
                // take the deviation of the cutoff potential
                daGrad[i] = 10.0 * daParams[i];
                dEnergy += 10.0*daParams[i]*daParams[i];
            } else {
                /*
                 * f'(x(i)) = 2x_i + 20pi*sin(2pi x_i)
                 */
                final double t1 = 2.0*Math.PI*daParams[i];
                daGrad[i] = 2.0*daParams[i] + 20.0 * Math.PI*Math.sin(t1);
                dEnergy += daParams[i]*daParams[i]-10.0*Math.cos(t1);
            }
        }

        if(DEBUG){
            // calculate a numerical gradient and compare it to the analytical one
            final double[] daNumGrad = new double[iNoOfParams];
            final double numE = NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, daNumGrad);

            for(int i = 0; i < daNumGrad.length; i++){
                if(Math.abs(daNumGrad[i]-daGrad[i]) >= 1E-7){
                    // we have a noticable difference
                    System.err.println("DEBUG: Analytical vs numerical gradient is: "
                            + daGrad[i] + " vs " + daNumGrad[i]);
                } else {
                    System.out.println("DEBUG: Analytical vs numerical gradient was fine. :-)");
                }
            }
        }

        return dEnergy;
    }

    @Override
    public AdaptiveParameters createInitialParameterStub(final ArrayList<CartesianCoordinates> refCartes,
            final String sMethod){

        final String[] saAtoms = {"XX"};
        final int[] iaParamsPerAt = {iDims};
        final AdaptiveParameters paramStub = new AdaptiveParameters(iDims, -1, saAtoms, iaParamsPerAt, sMethod);

        return paramStub;
    }

    @Override
    public double[][] minMaxBordersForParams(final AdaptiveParameters params){

        final double[][] daMinMax = new double[2][iDims];

        for(int i = 0; i < iDims; i++){
            daMinMax[0][i] = -5.12;
            daMinMax[1][i] = 5.12;
        }

        return daMinMax;
    }
}
