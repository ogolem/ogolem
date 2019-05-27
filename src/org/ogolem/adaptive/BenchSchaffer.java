/**
Copyright (c) 2010     , J. M. Dieterich and B. Hartke
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
 * Schaffers F7 benchmark function in N-D. Global minimum: f(x,y)=0: x=0,y=0.
 * @author Johannes Dieterich
 * @version 2015-07-27
 */
final class BenchSchaffer extends AbstractAdaptivable {

    private static final long serialVersionUID = (long) 20150727;
    private final static boolean DEBUG = false;
    private final int iDims;

    BenchSchaffer(final int dims){
        this.iDims = dims;
    }

    @Override
    public BenchSchaffer clone(){
        return new BenchSchaffer(iDims);
    }

    @Override
    public double energyOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds){

        final double[] daParams = params.getAllParamters();

        double dEnergy = 0.0;
        for(int i = 0; i < iDims - 1; i++){
            final double pi = daParams[i];
            final double piSq = pi*pi;
            final double piN = daParams[i+1];
            final double piNSq = piN*piN;
            final double si = Math.sqrt(piSq + piNSq);
            final double siR = Math.sqrt(si);
            final double t1 = Math.sin(50*Math.pow(si,1.0/5.0));
            dEnergy += siR *(t1*t1+1);
        }

        dEnergy /= (iDims-1.0);
        dEnergy *= dEnergy;

        return dEnergy;
    }

    @Override
    public double gradientOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds,
            final double[] grad){
        
        final double[] daParams = params.getAllParamters();

        final double norm = 1/(iDims-1.0);
        double dEnergy = 0.0;
        for(int i = 0; i < iDims - 1; i++){
            
            final double pi = daParams[i];
            final double piSq = pi*pi;
            final double piN = daParams[i+1];
            final double piNSq = piN*piN;
            final double si = Math.sqrt(piSq + piNSq);
            final double siR = Math.sqrt(si);
            final double t1 = 50*Math.pow(si,1.0/5.0);
            final double t2s = Math.sin(t1);
            final double t2sSq = t2s*t2s;
            
            final double t3 = siR *(t2sSq+1);
            
            final double si34 = siR*siR*siR;
            final double secDiv = 1/(2*si34);
            final double secX = pi*secDiv;
            final double secY = piN*secDiv;
            
            final double third = t2sSq*secDiv;
            final double thirdX = pi*third;
            final double thirdY = piN*third;
            
            final double si1320 = Math.pow((piSq+piNSq),13.0/20.0);
            final double t2c = Math.cos(t1);
            final double fourth = 20*t2s*t2c/si1320;
            final double fourthX = pi*fourth;
            final double fourthY = piN*fourth;
            
            grad[i  ] += (secX+thirdX+fourthX);
            grad[i+1] += (secY+thirdY+fourthY);
            dEnergy += t3;
        }
        
        dEnergy *= norm;
        for(int i = 0; i < iDims; i++) grad[i] *= dEnergy*2*norm;

        dEnergy *= dEnergy;
        
        if(DEBUG){
            // calculate a numerical gradient and compare it to the analytical one
            final double[] daNumGrad = new double[iDims];
            final double numE = NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, daNumGrad);

            for(int i = 0; i < daNumGrad.length; i++){
                if(Math.abs(daNumGrad[i]-grad[i]) >= 1E-7){
                    // we have a noticable difference
                    System.err.println("DEBUG: Analytical vs numerical gradient is: "
                            + grad[i] + " vs " + daNumGrad[i]);
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
            daMinMax[0][i] = -100.0;
            daMinMax[1][i] = 100.0;
        }

        return daMinMax;
    }
}
