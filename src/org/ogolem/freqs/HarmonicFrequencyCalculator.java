/**
Copyright (c) 2010-2013, J. M. Dieterich
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
package org.ogolem.freqs;

import contrib.jama.*;
import java.io.PrintWriter;
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Constants;

/**
 * Calculates harmonic frequencies. Rather primitively, to be honest. Partly,
 * this is due to Jama (which seemed like a nice and easy choice) not supporting
 * any fancinesses like symmetric matrices. Also, the NumericalHessianCalculator
 * does currently not exploit the symmetry of the Hessian matrix.
 * @author Johannes Dieterich
 * @version 2013-09-23
 */
public final class HarmonicFrequencyCalculator {
    
    private final static boolean DEBUG = false;

    public static Frequencies calculateFrequencies(final CartesianCoordinates cartes,
            final FrequencyMethod method, final BondInfo bonds){

        // first move to COM, just in case.
        cartes.moveCoordsToCOM();

        final int noOfAtoms = cartes.getNoOfAtoms();
        final String[] atoms = cartes.getAllAtomTypes();

        // calculate hessian
        final double[][] hessian = method.calculateHessian(-1, cartes, bonds);

        // calculate mass-weighting matrix
        final double[] masses = new double[noOfAtoms*3];
        for(int i = 0; i < noOfAtoms; i++){
            final double invMass = 1.0/Math.sqrt(AtomicProperties.giveWeight(atoms[i])*Constants.U_AU);
            if(DEBUG){
                System.out.println("DEBUG: Mass of atom " + atoms[i] + " is " + AtomicProperties.giveWeight(atoms[i])*Constants.U_AU
                    + " inv: " + invMass);
            }
            masses[i            ] = invMass;
            masses[i+  noOfAtoms] = invMass;
            masses[i+2*noOfAtoms] = invMass;
        }

        final Matrix hessianMat = new Matrix(hessian);
        if(DEBUG){
            System.out.println("DEBUG: Hessian matrix coming.");
            hessianMat.print(new PrintWriter(System.out,true),10,8);
            System.out.println("");
            
            // symmetry check
            System.out.println("DEBUG: Starting Hessian symmetry check");
            for(int i = 0; i < hessian.length-1; i++){
                for(int j = i+1; j < hessian.length; j++){
                    if(Math.abs(hessian[i][j]-hessian[j][i]) > 1e-6){
                        System.out.println("DEBUG: Hessian seems non-symmetric in position "
                                + i + "/" + j + "  " + hessian[i][j] + " vs " + hessian[j][i]);
                    }
                }
            }
            System.out.println("DEBUG: Finished Hessian symmetry check");
        }
        
        /*for(int i = 0; i < hessian.length; i++){
            for(int j = 0; j < hessian[i].length; j++){
                hessian[i][j] *= Constants.HARTREETOKCAL*Constants.BOHRTOANG*Constants.BOHRTOANG;
            }
        }*/

        // setup mass-weighted hessian
        final Matrix massHessianMat = new Matrix(hessian.length, hessian.length);
        final double[][] massHessian = massHessianMat.getArray();
        for(int i = 0; i < massHessian.length; i++){
            for(int j = 0; j < massHessian[i].length; j++){
                massHessian[i][j] = hessian[i][j]*masses[i]*masses[j];
                if(DEBUG){
                    System.out.println("DEBUG: mass weighted element " + i + " " + j + " " + massHessian[i][j] + " was " + hessian[i][j] + " " + masses[i] + " " + masses[j]);
                }
            }
        }

        if(DEBUG){
            System.out.println("DEBUG: Mass weighted Hessian matrix coming.");
            massHessianMat.print(new PrintWriter(System.out,true),10,8);
            System.out.println("");
            // symmetry check
            System.out.println("DEBUG: Starting mass weighted Hessian symmetry check");
            final double[][] massHess = massHessianMat.getArray();
            for(int i = 0; i < massHess.length-1; i++){
                for(int j = i+1; j < massHess.length; j++){
                    if(Math.abs(massHess[i][j]-massHess[j][i]) > 1e-6){
                        System.out.println("DEBUG: Mass weighted Hessian seems non-symmetric in position "
                                + i + "/" + j + "  " + massHess[i][j] + " vs " + massHess[j][i]);
                    }
                }
            }
            System.out.println("DEBUG: Finished mass-weighted Hessian symmetry check");
        }
        
        // diagonalize hessian
        final EigenvalueDecomposition eigenH = new EigenvalueDecomposition(hessianMat);
        // get eigenvalues of hessian
        final double[] imagEigenH = eigenH.getImagEigenvalues();
        final double[] realEigenH = eigenH.getRealEigenvalues();
        
        if(DEBUG){
            System.out.println("DEBUG: Imaginary eigenvalues of hessian");
            for(int i = 0; i < imagEigenH.length; i++){
                System.out.println(i + "\t" + imagEigenH[i]*Constants.HARTREETOKCAL);
            }
            System.out.println("DEBUG: Real eigenvalues of hessian");
            for(int i = 0; i < realEigenH.length; i++){
                System.out.println(i + "\t" + realEigenH[i]*Constants.HARTREETOKCAL);
            }
        }


        // diagonalize mass-weighted hessian
        final EigenvalueDecomposition eigen = new EigenvalueDecomposition(massHessianMat);

        // get eigenvalues of mass-weighted hessian
        final double[] imagEigen = eigen.getImagEigenvalues();
        final double[] realEigen = eigen.getRealEigenvalues();
        
        final double[] freqsReal = new double[realEigen.length];
        for(int i = 0; i < realEigen.length; i++){
            freqsReal[i] = Math.sqrt(Math.abs(realEigen[i]))/2/Math.PI;
        }
        
        if(DEBUG){
            System.out.println("DEBUG: Imaginary eigenvalues");
            for(int i = 0; i < imagEigen.length; i++){
                System.out.println(i + "\t" + imagEigen[i]);
            }
            System.out.println("DEBUG: Real eigenvalues");
            for(int i = 0; i < realEigen.length; i++){
                System.out.println(i + "\t" + realEigen[i]);
            }
        }
        
        // get eigenvectors
        final Matrix eigenV = eigen.getV();
        
        if(DEBUG){
            System.out.println("DEBUG: eigen vectors");
            eigenV.print(16, 8);
        }

        // XXX intesities unknown since no dipol moments (and no way to obtain them) available for anything ATM
        final Frequencies freqs = new Frequencies(imagEigen, freqsReal, null, eigenV.getArray());

        return freqs;
    }    
}
