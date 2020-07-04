/**
Copyright (c) 2019, J. M. Dieterich and B. Hartke
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
package org.ogolem.microbenchmarks;

import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.ScaTTM3FBackend;

/**
 * Tests TTM3F gradient cost with a random water 25 geometry.
 * @author Johannes Dieterich
 * @version 2020-05-25
 */
class WaterTTM3FLargeBenchmark implements SingleMicroBenchmark {
    
    private final ScaTTM3FBackend back;
    private final double[] xyz1D;
    private final String[] atoms;
    private final short[] atomNos;
    private final int[] atsPerMol;
    private final double[] energyparts;
    private final int noAtoms;
    private final float[] charges;
    private final short[] spins;
    private final BondInfo bonds;
    private final Gradient gradient;
    
    WaterTTM3FLargeBenchmark(){
        
        final CartesianCoordinates cartes = CartesianCoordinatesLibrary.getWater25Random();
        
        this.xyz1D = cartes.getAll1DCartes();
        this.atoms = cartes.getAllAtomTypes();
        this.atomNos = cartes.getAllAtomNumbers();
        this.atsPerMol = cartes.getAllAtomsPerMol();
        this.energyparts = new double[cartes.getNoOfMolecules()];
        this.noAtoms = cartes.getNoOfAtoms();
        this.charges = cartes.getAllCharges();
        this.spins = cartes.getAllSpins();
        this.bonds = null; // can be null since we really don't use it
        this.gradient = new Gradient(3, cartes.getNoOfAtoms());
        
        this.back = new ScaTTM3FBackend(cartes.getNoOfMolecules(), 0.5, false, "scattm3fout", false);
    }

    @Override
    public double runSingle() throws Exception {
        
        back.gradientCalculation(1, 1, xyz1D, atoms, atomNos, atsPerMol, energyparts, noAtoms, charges, spins, bonds, gradient, false);
        
        return gradient.getTotalEnergy();
    }

    @Override
    public String name() {
        return "water25 TTM3F gradient";
    }
}
