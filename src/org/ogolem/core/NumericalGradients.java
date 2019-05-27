/**
Copyright (c) 2011-2015, J. M. Dieterich
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
package org.ogolem.core;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.helpers.Machine;

/**
 * Gradients using numerical algorithms.
 * @author Johannes Dieterich
 * @version 2015-01-01
 */
public class NumericalGradients {

    private static final double numericalPrecision = Math.sqrt(Machine.calcMachinePrecision());

    public static Gradient numericalGradient(final long id, final int iteration, final double[] xyz,
            final String[] atoms, final short[] atomNos, final int[] atsPerMol, final double[] energyparts,
            final int noOfAtoms, final float[] charges, final short[] spins,
            final BondInfo bonds, CartesianFullBackend backend){

        // the matrix for the gradient
        final double[][] gradientMat = new double[3][noOfAtoms];
        // a 1D version
        final double[] gradient1D = new double[3 * noOfAtoms];
        // working copy of the coordinates
        final double[] xyzCopy = xyz.clone();

        for(int i = 0; i < gradient1D.length; i++){

            // increment
            final double h = (xyzCopy[i] == 0.0) ? 100*numericalPrecision : 100 * numericalPrecision * xyzCopy[i];

            // first plus it
            xyzCopy[i] += h;
            final double energy1 = backend.energyCalculation(id, iteration, xyzCopy,
                    atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds);

            // then minus it (there was the dH left from the step before)
            xyzCopy[i] -=  2*h;
            final double energy2 = backend.energyCalculation(id, iteration, xyzCopy,
                    atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds);
            
            xyzCopy[i] = xyz[i];

            // also checks for NaN and other cases cases
            if(energy1 > FixedValues.NONCONVERGEDENERGY || energy2 > FixedValues.NONCONVERGEDENERGY
                    || Double.isNaN(energy1) || Double.isNaN(energy2)){
                // we have at least one non converged result
                System.err.println("WARNING: Numerical gradient fails at "
                        + i +" for id " + id);
                gradient1D[i] = 0.0;
            } else {
                // calculate the three point stencil
                gradient1D[i] = (energy1 - energy2) / (2 * h);
            }
        }

        // now one last energy calculation
        final double totEnergy = backend.energyCalculation(id, iteration, xyzCopy,
                    atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds);
        
        // copy the arrays
        System.arraycopy(gradient1D, 0            , gradientMat[0], 0, noOfAtoms);
        System.arraycopy(gradient1D, noOfAtoms    , gradientMat[1], 0, noOfAtoms);
        System.arraycopy(gradient1D, noOfAtoms * 2, gradientMat[2], 0, noOfAtoms);

        // put the information in
        final Gradient gradient = new Gradient();
        gradient.setGradientTotal(gradientMat);
        gradient.setTotalEnergy(totEnergy);

        // return it
        return gradient;
    }
    
    public static double gradientForRigidBody(final Geometry ref, final List<CartesianCoordinates> origCartes,
            final double[][] coms, final double[][] eulers, final double[] gradient, final RigidBodyBackend backend, final int counter){
        
        // get a copy of the cartesians for work
        final List<CartesianCoordinates> copyCartes = new ArrayList<>();
        for(int mol = 0; mol < origCartes.size(); mol++){
            final CartesianCoordinates clone = origCartes.get(mol).clone();
            copyCartes.add(clone);
        }
        
        // zero: get one energy (need to rotate stuff for that)
        int iter = counter;
        // rotate once
        rotateCartesians(origCartes,copyCartes,eulers);
        final double energy = backend.energy(ref, copyCartes, coms, iter);
        
        final double hDef = 100*numericalPrecision;
        iter++;
        
        // first take care of all the COM gradient elements
        for(int mol = 0; mol < origCartes.size(); mol++){
            
            for(int coord = 0; coord < 3; coord++){
                // increment
                final double h = (coms[mol][coord] == 0.0) ? hDef : hDef * coms[mol][coord];
                // first plus it
                final double comOld = coms[mol][coord];
                coms[mol][coord] += h;
                final double eCoordP = backend.energy(ref, copyCartes, coms, iter);
                iter++;
                // then minus it
                coms[mol][coord] = comOld - h;
                final double eCoordM = backend.energy(ref, copyCartes, coms, iter);
                iter++;
                // then reset and calculate element
                coms[mol][coord] = comOld;
                gradient[mol*6+coord] = (eCoordP - eCoordM) / (2 * h);
                
                if(false){
                    System.out.println("DEBUG: Individual " + ref.getID() + " Molecule " + mol + " COM coordinate " + coord + " gradient " + gradient[mol*6+coord]);
                }
            
            }
            
        }
        
        for(int mol = 0; mol < origCartes.size(); mol++){
            
            // now the eulers, slightly more involved
            
            for(int euler = 0; euler < 3; euler++){
            
                // increment Euler and rotate
                final double eulerBak = eulers[mol][euler];
                eulers[mol][euler] += hDef;
                rotateCartesians(origCartes,copyCartes,eulers);
                final double eP = backend.energy(ref, copyCartes, coms, iter);
                iter++;
            
                // decrement Euler and rotate
                eulers[mol][euler] -= 2*hDef;
                rotateCartesians(origCartes,copyCartes,eulers);
                final double eM = backend.energy(ref, copyCartes, coms, iter);
                iter++;
                
                // reset and calculate gradient
                eulers[mol][euler] = eulerBak;
                gradient[mol*6+3+euler] = (eP - eM) / (2 * hDef);
                
                if(false){
                    System.out.println("DEBUG: Individual " + ref.getID() + " Molecule " + mol + " Euler orientation " + euler + " gradient " + gradient[mol*6+3+euler]);
                }
            }
        }
        
        return energy;
    }
    
    private static void rotateCartesians(final List<CartesianCoordinates> origCartes, final List<CartesianCoordinates> rotated,
            final double[][] eulers){
        
        assert(origCartes.size() == rotated.size());
        assert(eulers.length == origCartes.size());
        
        final double[][] rotMat = new double[3][3];
        
        for(int mol = 0; mol < origCartes.size(); mol++){
        
            final double[][] rotP = CoordTranslation.rotateXYZ(origCartes.get(mol).getAllXYZCoord(), eulers[mol], rotMat);
            rotated.get(mol).setAllXYZ(rotP);
        }
    }
}
