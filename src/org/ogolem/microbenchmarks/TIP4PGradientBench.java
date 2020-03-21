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
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CoordTranslation;
import org.ogolem.core.Geometry;
import org.ogolem.core.RigidBodyBackend;
import org.ogolem.core.RigidBodyCoordinates;
import org.ogolem.core.SimpleBondInfo;
import org.ogolem.core.TIP4PForceField;

/**
 * Benchmark for TIP4P gradient on a medium-sized water cluster.
 * @author Johannes Dieterich
 * @version 2020-02-08
 */
class TIP4PGradientBench implements SingleMicroBenchmark {
     
    private final Geometry w50;
    private final double[] rigidCoords;
    private final double[] grad;
    private final RigidBodyCoordinates rCoords;
    
    TIP4PGradientBench(){
        
        final RigidBodyBackend rBack = new TIP4PForceField();
        this.rCoords = new RigidBodyCoordinates(rBack);
        
        final int[] atsPerMol = new int[50];
        final boolean[] molFlexies = new boolean[50];
        final boolean[] molConstraints = new boolean[50];
        final String[] sids = new String[50];
        final ArrayList<boolean[][]> degreesOfFreedom = new ArrayList<>(50);
        for(int i = 0; i < 50; i++){
            atsPerMol[i] = 3;
            molFlexies[i] = false;
            molConstraints[i] = false;
            degreesOfFreedom.add(null);
            sids[i] = "water";
        }
        
        final boolean[][] constraintsXYZ = new boolean[3][150];
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 150; j++) constraintsXYZ[i][j] = false;
        }
        
        final BondInfo bonds = new SimpleBondInfo(150);
        for(int i = 0; i < 50; i++){
            final int o = i*3+0;
            final int h1 = i*3+1;
            final int h2 = i*3+1;
            bonds.setBond(o, h1, (short)1);
            bonds.setBond(o, h2, (short)1);            
        }
        
        final CartesianCoordinates cart = CartesianCoordinatesLibrary.getWater50LocMin();
        this.w50 = CoordTranslation.cartesianToGeometry(cart, 50, atsPerMol, molFlexies,
                degreesOfFreedom, molConstraints, constraintsXYZ, sids, bonds);
        
        this.rigidCoords = rCoords.getActiveCoordinates(w50);
        this.grad = new double[rigidCoords.length];
    }
    
    @Override
    public double runSingle() throws Exception {
        
        final double energy = rCoords.gradient(rigidCoords, grad, -1);
        
        return energy;
    }

    @Override
    public String name() {
        return "TIP4P water 50 gradient";
    }
}
