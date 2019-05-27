/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2012     , J. M. Dieterich
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

import org.ogolem.generic.genericpool.GenericPool;

/**
 * Providing some static analyzation functions.
 * @author Johannes Dieterich
 * @version 2012-02-20
 */
class Analytics {

    /**
     * Calculates the mean of all paramters of a given population.
     * @param pop The population.
     * @return An array with the mean values.
     */
    static double[] meanOfParameters(final GenericPool<Double,AdaptiveParameters> pop){

        int iNoOfParams = pop.getIndividualAtPosition(0).getNumberOfParamters();
        double[] daResults = new double[pop.getPoolSize()];

        for(int i = 0; i < pop.getPoolSize(); i++){

            double[] daCurrentParams = pop.getIndividualAtPosition(i).getAllParamters();

            for(int j = 0; j < iNoOfParams; j++){
                daResults[j] = daResults[j] + daCurrentParams[j];
            }
        }

        // normalize it
        for(int i = 0; i < iNoOfParams; i++){
            daResults[i] = daResults[i] / pop.getPoolSize();
        }
        
        return daResults;
    }

    /**
     * Calculates the variety of a given population.
     * @param pop The population.
     * @return An array with the variety.
     */
    static double[] varOfParameters(final GenericPool<Double,AdaptiveParameters> pop){

        int iNoOfParams = pop.getIndividualAtPosition(0).getNumberOfParamters();
        
        double[] daMeanValues = meanOfParameters(pop);

        double[] daVarValues = new double[iNoOfParams];

        for(int i = 0; i < pop.getPoolSize(); i++){

            double[] daCurrentParams = pop.getIndividualAtPosition(i).getAllParamters();

            for(int j = 0; j < iNoOfParams; j++){
                daVarValues[j] = daVarValues[j] + Math.pow((daCurrentParams[j] - daMeanValues[j]),2);
            }
        }

        // normalize it
        for(int i = 0; i < iNoOfParams; i++){
            daVarValues[i] = daVarValues[i] / pop.getPoolSize();
        }

        return daVarValues;
    }

    /**
     * Calculates a covariance matrix out of a given population.
     * @param pop The population.
     * @return The covariance matrix.
     */
    static double[][] covarianceOfParameters(final GenericPool<Double,AdaptiveParameters> pop){

        int iNoOfParams = pop.getIndividualAtPosition(0).getNumberOfParamters();

        int iSizeOfPop = pop.getPoolSize();

        double[] daMeanValues = meanOfParameters(pop);

        double[][] daCovarianceMatrix = new double[iNoOfParams][iNoOfParams];

        double[][] daAllParameters = new double[iSizeOfPop][iNoOfParams];

        for (int i = 0; i < iSizeOfPop; i++) {
            daAllParameters[i] = pop.getIndividualAtPosition(i).getAllParamters();
        }

        for(int i = 0; i < iNoOfParams; i++){
            for(int j = i; j < iNoOfParams; j++){
                for(int k = 0; k < iSizeOfPop; k++){
                    daCovarianceMatrix[i][j] = daCovarianceMatrix[i][j] +
                            (daAllParameters[k][i] - daMeanValues[i]) * (daAllParameters[k][j] - daMeanValues[j]);
                }
                daCovarianceMatrix[j][i] = daCovarianceMatrix[i][j];
                
            }
        }

        // normalize it
        for(int i = 0; i < daCovarianceMatrix.length; i++){
            for(int j = 0; j < daCovarianceMatrix[0].length; j++){
                daCovarianceMatrix[i][j] = daCovarianceMatrix[i][j] / iSizeOfPop;
            }
        }

        return daCovarianceMatrix;
    }
}
