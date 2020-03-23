/**
Copyright (c) 2020, J. M. Dieterich and B. Hartke
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
import org.ogolem.core.MixedLJForceField;
import org.ogolem.core.SimpleBondInfo;

/**
 * Benchmarks energy calculation in the mixed Lennard-Jones FF for Ar55.
 * @author Johannes Dieterich
 * @version 2020-02-01
 */
class MixedLJFFEnergyBench implements SingleMicroBenchmark {

    private final CartesianCoordinates lj55;
    private final BondInfo bonds55;
    private final MixedLJForceField ljFF;
    private final double[] energyparts;
    private final double[] xyz1D;
    
    MixedLJFFEnergyBench(){
        
        this.lj55 = CartesianCoordinatesLibrary.getAr55GlobMin();
        this.bonds55 = new SimpleBondInfo(55);
        this.ljFF = new MixedLJForceField();
        this.energyparts = new double[55];
        this.xyz1D = lj55.getAll1DCartes();
    }
    
    @Override
    public String name() {
        return "Ar55 mixed LJ energy bench";
    }
    
    @Override
    public double runSingle() throws Exception {
        
        final double e55 = ljFF.energyCalculation(-1, 0, xyz1D,
                lj55.getAllAtomTypes(), lj55.getAllAtomNumbers(), lj55.getAllAtomsPerMol(),
                energyparts, lj55.getNoOfAtoms(), lj55.getAllCharges(), lj55.getAllSpins(),
                bonds55);
        
        return e55;
    }
}
