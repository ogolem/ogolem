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
 * Schaffers F6 benchmark function in 2D. Global minimum: f(x,y)=0: x=0,y=0.
 * @author Johannes Dieterich
 * @version 2015-07-27
 */
final class BenchSchafferF6 extends AbstractAdaptivable{

    private static final long serialVersionUID = (long) 20150727;
    private static final boolean DEBUG = false;
    
    BenchSchafferF6(){
    }

    @Override
    public BenchSchafferF6 clone(){
        return new BenchSchafferF6();
    }

    @Override
    public double energyOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds){
        
        final double[] daParams = params.getAllParamters();
        final double p0 = daParams[0];
        final double p0Sq= p0*p0;
        final double p1 = daParams[1];
        final double p1Sq= p1*p1;

        final double t1 = (Math.sin(Math.sqrt(p0Sq + p1Sq)));
        final double dNom = t1*t1 - 0.5;

        final double t2 = 1.0 + 0.001*(p0Sq + p1Sq);
        final double dDenom = t2*t2;

        // we make no check whether the params are within the boundaries, seems
        // to work anyways. :-)

        return (0.5 + dNom/dDenom);
    }

    @Override
    public double gradientOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds,
            final double[] grad){
        
        final double[] daParams = params.getAllParamters();
        final double p0 = daParams[0];
        final double p0Sq= p0*p0;
        final double p1 = daParams[1];
        final double p1Sq= p1*p1;

        final double t0 = p0Sq + p1Sq;
        final double t00 = Math.sqrt(t0);
        final double t1 = Math.sin(t00);
        final double t11 = Math.sin(2*t00);
        final double dNom = t1*t1 - 0.5;

        final double t2 = 1.0 + 0.001*t0;
        final double dDenom = 1.0/(t2*t2);
        
        final double func = (0.5 + dNom*dDenom);
        
        final double gradT = (-0.004*t1*t1+(t2*t11)/t00+0.002)*dDenom/t2;
        grad[0] = p0*gradT;
        grad[1] = p1*gradT;

        if(DEBUG){
            final double[] numGrad = new double[grad.length];

            final double numFunc = NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, numGrad);
            System.out.println("DEBUG OUTPUT FOR NUMGRAD");
            System.out.println("Func " + numFunc + "\t" + func + "\t" + (numFunc-func));
            for(int i = 0; i < grad.length; i++){
                System.out.println(numGrad[i]+ "\t" + grad[i] + "\t" + (numGrad[i]-grad[i]));
                grad[i] = numGrad[i];
            }
                
            return numFunc;
        }
        
        return func;
    }

    @Override
    public AdaptiveParameters createInitialParameterStub(final ArrayList<CartesianCoordinates> refCartes,
            final String sMethod){

        final String[] saAtoms = {"XX"};
        final int[] iaParamsPerAt = {2};
        final AdaptiveParameters paramStub = new AdaptiveParameters(2, -1, saAtoms, iaParamsPerAt, sMethod);

        return paramStub;
    }

    @Override
    public double[][] minMaxBordersForParams(final AdaptiveParameters params){

        final double[][] daMinMax = new double[2][2];

        for(int i = 0; i < 2; i++){
            daMinMax[0][i] = -100.0;
            daMinMax[1][i] = 100.0;
        }

        return daMinMax;
    }
}
