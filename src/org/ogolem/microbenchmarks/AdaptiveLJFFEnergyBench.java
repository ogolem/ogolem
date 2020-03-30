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

import java.util.ArrayList;
import org.ogolem.adaptive.AdaptiveLJFF;
import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.SimpleBondInfo;

/**
 * Benchmarks energy calculation of the adaptive LJ FF for a homogeneous LJ55
 * configured as if it were a homogeneous LJ-6-12-6 FF.
 * @author Johannes Dieterich
 * @version 2020-03-21
 */
class AdaptiveLJFFEnergyBench implements SingleMicroBenchmark {
    
    private final CartesianCoordinates lj55;
    private final double[] xyz1D;
    private final BondInfo bonds55;
    private final double[] energyparts;
    private final AdaptiveLJFF ljFF;
    
    AdaptiveLJFFEnergyBench(){
        
        this.lj55 = CartesianCoordinatesLibrary.getAr55GlobMin();
        this.bonds55 = new SimpleBondInfo(55);
        this.energyparts = new double[55];
        this.xyz1D = lj55.getAll1DCartes();
        
        final AdaptiveLJFF dummy = new AdaptiveLJFF(true /* pretend to be in adaptive*/, false,
                6, 12, 6, null /* only then can we set params to null */, false, 1.2, 20.0,
                false /* no cache, dummy */, false);
        final ArrayList<CartesianCoordinates> cartes = new ArrayList<>();
        cartes.add(lj55);
        
        final AdaptiveParameters params = dummy.createInitialParameterStub(cartes, dummy.getMethodID());
        
        // now set the parameters to our regular LJ parameters for Ar
        final double sigma = AtomicProperties.giveLennardJonesSigma("Ar");
        final double epsilon = AtomicProperties.giveLennardJonesEpsilon("Ar");
        
        final String[] keys = params.getAllKeysCopy();
        assert(keys.length == 1); // there should only be one
        
        final int noParams = params.getAmountOfParametersForKey(keys[0]);
        assert(noParams == 3); // must be one epsilon and two sigma (same sigma)
        
        final double[] paramVals = params.getAllParamters();
        paramVals[0] = epsilon;
        paramVals[1] = -sigma; // to make it 12-6
        paramVals[2] = sigma;
        
        this.ljFF = new AdaptiveLJFF(false /* pretend to be in core*/, false,
                6, 12, 6, params /* real params */, false, 1.2, 20.0,
                true /* cache */, false);
        
    }
    
    @Override
    public String name() {
        return "Ar55 adaptive LJ energy bench";
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
