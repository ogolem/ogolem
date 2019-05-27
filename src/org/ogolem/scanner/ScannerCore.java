/**
Copyright (c) 2013, J. M. Dieterich
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
package org.ogolem.scanner;

import org.ogolem.core.CartesianFullBackend;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Newton;
import org.ogolem.core.ZMatrix;
import org.ogolem.io.OutputPrimitives;

/**
 * The core of the scanner. Essentially a recursively working algorithm.
 * @author Johannes Dieterich
 * @version 2013-09-23
 */
class ScannerCore {
    
    private ScannerCore(){}; // no instantiation!// no instantiation!
    
    static void scan(final ScannerConfig scanConf, final CartesianCoordinates cartes, final Newton locopt, final String outPrefix,
            final boolean[][] constraints, final boolean isConstricted, final BondInfo bonds){
        
        // while we are there: write out our first geometry
        cartes.updateCartesians();
        if (scanConf.singlePoints) {
            final CartesianFullBackend back = locopt.getBackend();
            if (back == null) {
                System.err.println("WARNING: Backend is null although you want single point energies. Well, you will not get any.");
            } else {
                final double e = back.energyCalculation(0, 0, cartes.getAll1DCartes(),
                        cartes.getAllAtomTypes(), cartes.getAllAtomNumbers(),
                        cartes.getAllAtomsPerMol(), new double[cartes.getNoOfAtoms()],
                        cartes.getNoOfAtoms(), cartes.getAllCharges(), cartes.getAllSpins(), bonds);
                cartes.setEnergy(e);
            }
        }
        final String[] father = cartes.createPrintableCartesians();
        try {
            OutputPrimitives.writeOut(outPrefix + "-orig-nonopt.xyz", father, false);
        } catch (Exception e) {
            System.err.println("WARNING: Failure to write original snon-optimized geometry.");
            e.printStackTrace(System.err);
        }
        
        if (scanConf.optimize) {
            final CartesianCoordinates copy = new CartesianCoordinates(cartes);
            try {
                final CartesianCoordinates opt = locopt.cartesToCartes(-1, copy, constraints, isConstricted, bonds);
                final String[] optCoords = opt.createPrintableCartesians();
                OutputPrimitives.writeOut(outPrefix + "-orig-opt.xyz", optCoords, false);
            } catch (Exception e) {
                System.err.println("WARNING: Failure to optimize or write original optimized geometry.");
                e.printStackTrace(System.err);
            }
        }

        
        final State state = new State();
        scanner(state,0,scanConf,cartes,locopt,outPrefix,constraints,isConstricted,bonds);
        
    }
    
    private static void scanner(final State state, final int myLevel, final ScannerConfig scanConf, final CartesianCoordinates cartes, final Newton locopt, final String outPrefix,
            final boolean[][] constraints, final boolean isConstricted, final BondInfo bonds){
        
        if(myLevel >= scanConf.scanDoFs.size()) {
            // well, it seems we are at the bottom. so write out the geometry and (if wanted) optimize it
            cartes.updateCartesians();
            state.geomID++;
            
            if(scanConf.singlePoints){
                final CartesianFullBackend back  = locopt.getBackend();
                if(back == null){
                    System.err.println("WARNING: Backend is null although you want single point energies. Well, you will not get any.");
                } else {
                    final double e = back.energyCalculation(state.geomID, 0, cartes.getAll1DCartes(),
                            cartes.getAllAtomTypes(), cartes.getAllAtomNumbers(),
                            cartes.getAllAtomsPerMol(), new double[cartes.getNoOfAtoms()],
                            cartes.getNoOfAtoms(), cartes.getAllCharges(), cartes.getAllSpins(), bonds);
                    cartes.setEnergy(e);
                }
            }
            
            final String[] unopt = cartes.createPrintableCartesians();
            try{
                OutputPrimitives.writeOut(outPrefix + state.geomID + "-nonopt.xyz", unopt, false);
            } catch (Exception e){
                System.err.println("WARNING: Failure to write non-optimized geometry " + state.geomID);
                e.printStackTrace(System.err);
            }
            
            if(scanConf.optimize){
                final CartesianCoordinates copy = new CartesianCoordinates(cartes);
                try{
                    final CartesianCoordinates opt = locopt.cartesToCartes(state.geomID, copy, constraints, isConstricted, bonds);
                    final String[] optCoords = opt.createPrintableCartesians();
                    OutputPrimitives.writeOut(outPrefix + state.geomID + "-opt.xyz", optCoords, false);
                } catch(Exception e){
                    System.err.println("WARNING: Failure to optimize or write optimized geometry " + state.geomID);
                    e.printStackTrace(System.err);
                }
            }
            
            return;
        }
                
        final ScannerDoF myDoF = scanConf.scanDoFs.get(myLevel);
        final ZMatrix zmat = cartes.getMolecularRefZMatrix(myDoF.getMolecule());
        for(final double currVal : myDoF){
            zmat.setDOF(myDoF.getAtom(), myDoF.getDoF(), currVal);
            // now go one down
            scanner(state, myLevel+1,scanConf,cartes,locopt,outPrefix,constraints,isConstricted,bonds);
        }
    }
    
    private static class State {
        
        long geomID = (long) 0;
    }
}
