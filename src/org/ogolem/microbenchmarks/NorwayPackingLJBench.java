/**
Copyright (c) 2019-2020, J. M. Dieterich and B. Hartke
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

import java.util.List;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CollisionDetection;
import org.ogolem.core.DissociationDetection;
import org.ogolem.core.Geometry;
import org.ogolem.core.NorwayGeometryMutation;
import org.ogolem.core.SimpleBondInfo;
import org.ogolem.random.Lottery;
import org.ogolem.random.StandardRNG;

/**
 * Benchmark the Norway packing mutation/init for a LJ cluster.
 * @author Johannes Dieterich
 * @version 2020-08-09
 */
class NorwayPackingLJBench implements SingleMicroBenchmark {

    private static final int RNGSEED = 42;
    private static final double BLOWCOLL = 1.1;
    private static final double BLOWDISS = 2.2;
    private static final String LJATOM = "Xe";
    
    private final int noLJAtoms;
    private final NorwayGeometryMutation norway;
    private final Geometry geom;
    
    NorwayPackingLJBench(final int noLJAtoms){
        
        Lottery.setGenerator(new StandardRNG(RNGSEED));

        this.norway = new NorwayGeometryMutation(NorwayGeometryMutation.PACKDIM.THREED, CollisionDetection.CDTYPE.SIMPLEPAIRWISE, BLOWCOLL, BLOWDISS,
            DissociationDetection.DEFAULTDD, NorwayGeometryMutation.MUTMODE.ASCENDING);

        assert(noLJAtoms > 0);
        this.noLJAtoms = noLJAtoms;

        final BondInfo bonds = new SimpleBondInfo(noLJAtoms);

        final int[] atsPerMol = new int[noLJAtoms];
        final boolean[] molFlexies = new boolean[noLJAtoms];
        final List<boolean[][]> allFlexies = null;
        final boolean[] molConstraints = new boolean[noLJAtoms];
        final boolean[][] constrXYZ = new boolean[3][noLJAtoms];
        final String[] sids = new String[noLJAtoms];
        
        for(int i = 0; i < noLJAtoms; i++){
            atsPerMol[i] = 1;
            molFlexies[i] = false;
            molConstraints[i] = false;
            sids[i] = LJATOM;
            constrXYZ[0][i] = false;
            constrXYZ[1][i] = false;
            constrXYZ[2][i] = false;
        }
        
        final CartesianCoordinates cartes = new CartesianCoordinates(noLJAtoms, noLJAtoms, atsPerMol);
        final String[] atoms = cartes.getAllAtomTypes();
        
        for(int i = 0; i < noLJAtoms; i++){
            atoms[i] = LJATOM;
        }
        
        cartes.recalcAtomNumbersForced();
        
        this.geom = new Geometry(cartes, 0, noLJAtoms, atsPerMol, molFlexies, allFlexies, molConstraints,
            constrXYZ, sids, bonds);
    }
    
    @Override
    public String name() {
        return "Norway packing LJ bench for " + this.noLJAtoms + " LJ atoms";
    }
    
    @Override
    public double runSingle() throws Exception {
        
        final Geometry packed = norway.mutate(geom);
        return packed.getFitness();
    }
}
