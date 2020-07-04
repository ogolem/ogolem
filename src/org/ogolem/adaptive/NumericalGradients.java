/**
Copyright (c) 2010-2013, J. M. Dieterich
              2015-2016, J. M. Dieterich and B. Hartke
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

import org.ogolem.adaptive.genericfitness.PropertyCalculator;
import org.ogolem.adaptive.genericfitness.ReferenceInputData;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.Topology;
import org.ogolem.helpers.Machine;
import org.ogolem.properties.Property;

/**
 * Computes numerical parameter gradients.
 * @author Johannes Dieterich
 * @version 2016-02-24
 */
public class NumericalGradients extends org.ogolem.core.NumericalGradients {
    
    static double calculateParamGrad(final CartesianCoordinates cartes,
            final AdaptiveParameters params, Adaptivable adaptivable, final int geomID,
            final BondInfo bonds, final double[] daGrad){
        return calculateParamGrad(cartes, params, adaptivable, geomID, bonds, daGrad, false);
    }

    static double calculateParamGrad(final CartesianCoordinates cartes,
            final AdaptiveParameters params, Adaptivable adaptivable, final int geomID,
            final BondInfo bonds, final double[] daGrad, final boolean twoPoint){

        final double dMachinePrecision = Machine.calcMachinePrecision();
        final int iNoOfParams = params.getNumberOfParamters();
        final double[] daOrigParamsCopy = params.getAllParamtersCopy();

        final double dTotEnergy = adaptivable.energyOfStructWithParams(cartes, params, geomID, bonds);
        for(int i = 0; i < iNoOfParams; i++){

            double dH = 1000 * Math.sqrt(dMachinePrecision) * daOrigParamsCopy[i];
            if(dH == 0.0) {dH = Math.sqrt(dMachinePrecision);}

            final double dVariableVal1 = daOrigParamsCopy[i] + dH;
            final double dVariableVal2 = daOrigParamsCopy[i] - dH;

            // put them in
            params.setParameterAtPos(dVariableVal1, i);
            final double dEnergy1 = adaptivable.energyOfStructWithParams(cartes, params, geomID, bonds);
            
            if(twoPoint){
                // put the original value back in
                params.setParameterAtPos(daOrigParamsCopy[i], i);
                // calculate the two point stencil
                daGrad[i] = (dEnergy1 - dTotEnergy) / dH;
                continue;
            }

            params.setParameterAtPos(dVariableVal2,i);
            final double dEnergy2 = adaptivable.energyOfStructWithParams(cartes, params, geomID, bonds);

            // put the original value back in
            params.setParameterAtPos(daOrigParamsCopy[i], i);

            // calculate the three point stencil
            daGrad[i] = (dEnergy1 - dEnergy2) / (2.0 * dH);
        }

        return dTotEnergy;
    }

    static ParameterGradient calculateParamGrad(final Topology topology,
            final AdaptiveParameters params, AdaptiveInteractionTerm term,
            int start, int end){

        final double dMachinePrecision = Machine.calcMachinePrecision();

        final int iNoOfParams = params.getNumberOfParamters();

        final ParameterGradient grad = new ParameterGradient(iNoOfParams);

        final double[] daGrad = new double[iNoOfParams];

        final double[] daOrigParamsCopy = params.getAllParamtersCopy();

        for(int i = start; i <= end; i++){

            double dH = 1000 * Math.sqrt(dMachinePrecision) * daOrigParamsCopy[i];
            if(dH == 0.0) {dH = Math.sqrt(dMachinePrecision);}

            final double dVariableVal1 = daOrigParamsCopy[i] + dH;
            final double dVariableVal2 = daOrigParamsCopy[i] - dH;

            // put them in
            params.setParameterAtPos(dVariableVal1, i);
            final double dEnergy1 = term.partialInteraction(topology, params);

            params.setParameterAtPos(dVariableVal2,i);
            final double dEnergy2 = term.partialInteraction(topology, params);

            // put the original value back in
            params.setParameterAtPos(daOrigParamsCopy[i], i);

            // calculate the three point stencil
            daGrad[i] = (dEnergy1 - dEnergy2) / (2.0 * dH);
        }

        // now one last energy calculation
        final double dTotEnergy = term.partialInteraction(topology, params);

        grad.setTotalEnergy(dTotEnergy);
        grad.setGradient(daGrad);

        return grad;
    }

    static Gradient calculateGrad(AdaptiveParameters params,
            AdaptiveInteractionTerm term, Topology refTopo){
        
        System.err.println("WARNING: SOMETHING OFF IN THIS METHOD");

        final double dMachinePrecision = Machine.calcMachinePrecision();
        final double machPrecSqrt = Math.sqrt(dMachinePrecision);

        final int iNoOfAtoms = refTopo.getNumberOfAtoms();
        final double[][] daGradMatrix = new double[3][iNoOfAtoms];
        final double[][] daXYZ = refTopo.getPositions();
        final String[] saAtoms = refTopo.getAtomNames();
        
        // an initial energy calculation
        Topology topo = new Topology(saAtoms,daXYZ,refTopo.getBonds(), refTopo.getCharges(), refTopo.getSpins(), refTopo.getAtomicNumbers(),
                        refTopo.get13Contributions(), refTopo.get14Contributions());
        final double dEnergy = term.partialInteraction(topo, params);

        for(int i = 0; i < 3; i++){
            for(int j = 0; j < iNoOfAtoms; j++){

                double dH = 1000 * machPrecSqrt * daXYZ[i][j];
                if(dH < machPrecSqrt) {dH = machPrecSqrt;}
             
                System.out.println(" comp " + i + " atom " + j + " dH " + dH + " " + machPrecSqrt);
             
                final double orig = daXYZ[i][j];
                
                // add it
                daXYZ[i][j] = orig + dH;

                topo = new Topology(saAtoms,daXYZ,refTopo.getBonds(), refTopo.getCharges(), refTopo.getSpins(), refTopo.getAtomicNumbers(),
                        refTopo.get13Contributions(), refTopo.get14Contributions());
                final double dEnergy1 = term.partialInteraction(topo, params);

                // substract it
                daXYZ[i][j] =  orig - dH;
                topo = new Topology(saAtoms,daXYZ,refTopo.getBonds(), refTopo.getCharges(), refTopo.getSpins(), refTopo.getAtomicNumbers(),
                        refTopo.get13Contributions(), refTopo.get14Contributions());
                final double dEnergy2 = term.partialInteraction(topo, params);

                // correct it
                daXYZ[i][j] = orig;

                // compute stencil
                daGradMatrix[i][j] = (dEnergy1 - dEnergy2) / (2.0 * dH);
            }
        }

        // create Gradient
        final Gradient gradient = new Gradient();
        gradient.setGradientTotal(daGradMatrix);
        gradient.setTotalEnergy(dEnergy);

        return gradient;
    }
    
    /**
     * Calculates the numerical gradient for a certain property. Please note that this relies on getValue() for this property, which may not make sense!
     * @param <T> the type of property
     * @param calculator the property calculator
     * @param p the paramters
     * @param geom the geometry
     * @param grad the gradient array, changed on exit
     * @return the property for the adaptive parameters
     */
    static <T extends Property, V extends ReferenceInputData<T>> T numericalPropertyGradient(final PropertyCalculator<T,V> calculator, final AdaptiveParameters p, final V data, final double[] grad){
        
        final double[] params = p.getAllParamters();
        final int noOfParams = params.length;
        assert(params.length == grad.length);
        
        final double machinePrecision = Machine.calcMachinePrecision();

        for(int i = 0; i < noOfParams; i++){

            double dH = 1000 * Math.sqrt(machinePrecision) * params[i];
            if(dH == 0.0) {dH = Math.sqrt(machinePrecision);}

            final double origParam = params[i];
            final double plusParam = origParam + dH;
            final double minusParam = origParam - dH;
            
            // plus evaluation
            params[i] = plusParam;
            final T propPlus = calculator.calculateProperty(p, data);
            final double resPlus = propPlus.getValue();
            
            // minus evaluation
            params[i] = minusParam;
            final T propMinus = calculator.calculateProperty(p, data);
            final double resMinus = propMinus.getValue();
            
            // reset and compute gradient
            params[i] = origParam;
            grad[i] = (resPlus - resMinus) / (2.0 * dH);
        }
        
        // one last evaluation
        final T prop = calculator.calculateProperty(p, data);
        
        return prop;
    }
}
