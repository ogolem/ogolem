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

import org.ogolem.core.AdvancedPairWise;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CoordTranslation;
import org.ogolem.core.DummyCollisionStrengthComputer;

/**
 * Benchmark the advanced pairwise CD which is our work horse for CD. Check only
 * for collisions w/o calculating explicit distances or storing them.
 * @author Johannes Dieterich
 * @version 2020-02-01
 */
class AdvPairwiseCDCheckOnlyBenchmark implements SingleMicroBenchmark {
    
    private static final double BLOWFAC = 1.4;
    private final CartesianCoordinates cartesians;
    private final AdvancedPairWise cd;
    private final BondInfo bonds;
    
    AdvPairwiseCDCheckOnlyBenchmark(){
        this.cartesians = CartesianCoordinatesLibrary.getKanamycinAPM3Opt();
        this.cd = new AdvancedPairWise(false, new DummyCollisionStrengthComputer());
        this.bonds = CoordTranslation.checkForBonds(cartesians, BLOWFAC);
    }
    
    @Override
    public String name() {
        return "advanced pairwise collision detection (check only) bench";
    }
    
    @Override
    public double runSingle() throws Exception {
        
        final boolean coll = cd.checkOnlyForCollision(cartesians, BLOWFAC, bonds);
        
        return (coll) ? 0 : 1;
    }
}
