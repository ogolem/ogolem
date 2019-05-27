/**
Copyright (c) 2010-2014, J. M. Dieterich
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
package org.ogolem.freqs;

import org.ogolem.core.CartesianFullBackend;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;

/**
 * Numeric calculation of the Hessian.
 * @author Johannes Dieterich
 * @version 2015-03-27
 */
public final class NumericalHessianCalculator {
    
    private static final boolean DEBUG = false;
    private static final double NUMERICALINCREMENT = 1e-4;
    
    public static double[][] numericalHessian (final long id, final CartesianCoordinates cartes,
            final CartesianFullBackend backend, final BondInfo bonds){
    
        final double[][] hessGrads = numericalHessianFromGrads(id, cartes,backend,bonds);
        
        if(DEBUG){
            final double[][] hessPoints = numericalHessianFromPoints(id, cartes,backend,bonds);
            for(int i = 0; i < hessPoints.length; i++){
                for(int j = 0; j < hessPoints.length; j++){
                    final double diff = Math.abs(hessPoints[i][j] - hessGrads[i][j]);
                    if(diff >= 1E-7){System.err.println("DEBUG: Mismatch at " + i + "\t" + j + " diff: " + diff
                            + "\t" + hessPoints[i][j] + "\t" + hessGrads[i][j]);}
                    else {System.err.println("DEBUG: No mismatch at " + i + "\t" + j + " diff: " + diff
                            + "\t" + hessPoints[i][j]);}
                }
            }
            
            return hessGrads;
        }
        
        return hessGrads;
    }
    
    public static double[][] numericalHessianFromPoints(final long id, final CartesianCoordinates cartes,
            final CartesianFullBackend backend, final BondInfo bonds){

        final int noOfAtoms = cartes.getNoOfAtoms();
        final double[] xyz = cartes.getAll1DCartes();
        final String[] atoms = cartes.getAllAtomTypes();
        final short[] spins = cartes.getAllSpins();
        final float[] charges = cartes.getAllCharges();
        final short[] atomNos = cartes.getAllAtomNumbers();
        final int[] atsPerMol = cartes.getAllAtomsPerMol();
        final double[] energyparts = new double[cartes.getAllAtomsPerMol().length];

        final double[][] hessian = new double[noOfAtoms*3][noOfAtoms*3];
        
        final double e = backend.energyCalculation(id, 0, 
                xyz, atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds);
        
        int iter = 1;
        for(int i = 0; i < noOfAtoms*3; i++){
            
            final double incI = (xyz[i] == 0.0) ? NUMERICALINCREMENT : xyz[i]*NUMERICALINCREMENT;
            final double coordI = xyz[i];    
            
            // simple for i == j
            xyz[i] = coordI + incI;
            final double eP = backend.energyCalculation(id, iter,
                    xyz, atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds);
            
            xyz[i] = coordI - incI;
            final double eM = backend.energyCalculation(id, iter,
                    xyz, atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds);
            
            if(DEBUG){
                System.out.println("DEBUG: " + i + " e  " + eP);
                System.out.println("DEBUG: " + i + " eP  " + eP);
                System.out.println("DEBUG: " + i + " eM  " + eM);
                System.out.println("DEBUG: " + i + " incI " + incI);
                System.out.println("DEBUG: " + i + " div  " + (2*incI));
                System.out.println("DEBUG: " + i + " num  " + (eP - eM));
            }
            
            hessian[i][i] = (eP + eM -2*e) / (incI*incI);
            
            for(int j = i+1; j < noOfAtoms*3; j++){
                
                final double incJ = (xyz[j] == 0.0) ? NUMERICALINCREMENT : xyz[j]*NUMERICALINCREMENT;
                
                final double coordJ = xyz[j];
                xyz[i] = coordI + incI;
                xyz[j] = coordJ + incJ;
                final double ePP = backend.energyCalculation(id, iter,
                    xyz, atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds);
                iter++;
                
                xyz[i] = coordI - incI;
                final double eMP = backend.energyCalculation(id, iter,
                    xyz, atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds);
                iter++;
                
                xyz[j] = coordJ - incJ;
                final double eMM = backend.energyCalculation(id, iter,
                    xyz, atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds);
                iter++;
                
                xyz[i] = coordI + incI;
                final double ePM = backend.energyCalculation(id, iter,
                    xyz, atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds);
                iter++;
                
                if(DEBUG){
                    System.out.println("DEBUG: " + i + " " + j + " ePP  " + ePP);
                    System.out.println("DEBUG: " + i + " " + j + " ePM  " + ePM);
                    System.out.println("DEBUG: " + i + " " + j + " eMP  " + eMP);
                    System.out.println("DEBUG: " + i + " " + j + " eMM  " + eMM);
                    System.out.println("DEBUG: " + i + " " + j + " incI " + incI);
                    System.out.println("DEBUG: " + i + " " + j + " incJ " + incJ);
                    System.out.println("DEBUG: " + i + " " + j + " div  " + (4*incI*incJ));
                    System.out.println("DEBUG: " + i + " " + j + " num  " + (ePP + eMM - ePM - eMP));
                }
                
                hessian[i][j] = (ePP + eMM - ePM - eMP) / (4*incI*incJ);
                hessian[j][i] = hessian[i][j];
                
                xyz[i] = coordI;
                xyz[j] = coordJ;
            }
        }
        
        return hessian;
   }

    public static double[][] numericalHessianFromGrads(final long id, final CartesianCoordinates cartes,
            final CartesianFullBackend backend, final BondInfo bonds){
        
        final int noOfAtoms = cartes.getNoOfAtoms();
        final double[] xyz = cartes.getAll1DCartes();
        final String[] atoms = cartes.getAllAtomTypes();
        final short[] spins = cartes.getAllSpins();
        final float[] charges = cartes.getAllCharges();
        final short[] atomNos = cartes.getAllAtomNumbers();
        final int[] atsPerMol = cartes.getAllAtomsPerMol();
        final double[] energyparts = new double[cartes.getAllAtomsPerMol().length];

        final double[][] hessian = new double[noOfAtoms*3][noOfAtoms*3];
        
        final double e = backend.energyCalculation(id, 0, 
                xyz, atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds);
        
        final Gradient gPlus = new Gradient(3,noOfAtoms);
        final Gradient gMinus = new Gradient(3,noOfAtoms);

        int iter = 1;
        for(int i = 0; i < noOfAtoms*3; i++){
            
            final double incr = (xyz[i] == 0.0) ? NUMERICALINCREMENT : xyz[i]*NUMERICALINCREMENT;
            final double invIncr = 1.0/(2*incr);
           
            final double coord = xyz[i];
            xyz[i] += incr;
            backend.gradientCalculation(id, iter,
                    xyz, atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds,
                    gPlus);
            final double[] gradPlus = gPlus.getGradient();
            final double eP = gPlus.getTotalEnergy();
            iter++;
            
            xyz[i] = coord - incr;
            backend.gradientCalculation(id, iter,
                    xyz, atoms, atomNos, atsPerMol, energyparts, noOfAtoms, charges, spins, bonds,
                    gMinus);
            final double[] gradMinus = gMinus.getGradient();
            final double eM = gMinus.getTotalEnergy();
            iter++;

            // calculate the second derivative
            for (int j = 0; j < noOfAtoms * 3; j++) {
                if (DEBUG) {
                    System.out.println("DEBUG: j " + j + " gradPlus " + gradPlus[j] + " gradMinus " + gradMinus[j]);
                }
                if(i == j){
                    hessian[i][i] = (eP + eM - 2*e) / (incr*incr);
                } else {
                    hessian[i][j] = (gradPlus[j] - gradMinus[j]) * invIncr;
                }
            }

            // put this coord back in place
            xyz[i] = coord;
        }

        return hessian;
    }
}
