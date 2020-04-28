/**
Copyright (c) 2015, J. M. Dieterich
              2016, J. M. Dieterich and B. Hartke
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

import java.util.List;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.RandomUtils;

/**
 * A genotype crossover only taking the molecular orientations into account.
 * @author Johannes Dieterich
 * @version 2016-12-18
 */
public class FoehrGeometryXOver implements GenericCrossover<Molecule,Geometry>{
    
    private static final long serialVersionUID = (long) 20150217;
    
    private final int noCuts;
    
    FoehrGeometryXOver(final int noCuts){
        this.noCuts = noCuts;
    }
    
    FoehrGeometryXOver(final FoehrGeometryXOver orig){
        this.noCuts = orig.noCuts;
    }

    @Override
    public FoehrGeometryXOver clone() {
        return new FoehrGeometryXOver(this);
    }

    @Override
    public String getMyID() {
        
        String s = "FOEHR GEOMETRY X-OVER\n";
        s += "\t no. of cuts " + noCuts;
        
        return s;
    }

    @Override
    public Tuple<Geometry, Geometry> crossover(final Geometry mother, final Geometry father, final long futureID) {
        
        final int noMols = mother.getNumberOfIndieParticles();
        if(noCuts >= noMols){
            System.err.println("ERROR: Too many cutting points for the number of mols (" + noCuts + " vs " + noMols + ").");
            return new Tuple<>(null,null);
        }

        // get cutting points
        final List<Integer> cutPoints = RandomUtils.listOfPoints(noCuts, mother.getNumberOfIndieParticles());
        
        final Geometry gChild1 = mother.clone();
        final Geometry gChild2 = father.clone();
        gChild1.setID(futureID);
        gChild2.setID(futureID);
                
        int cutter = 0;
        boolean swap = false;
        for(int mol = 0; mol < noMols; mol++){
            if(cutter < noCuts && cutPoints.get(cutter) == mol){
                swap = !swap;
                cutter++;
            }
            
            if(swap){
                final double[] eulers1 = gChild1.getMoleculeAtPosition(mol).getOrientation().clone();
                final double[] eulers2 = gChild2.getMoleculeAtPosition(mol).getOrientation().clone();
                
                gChild1.getMoleculeAtPosition(mol).setOrientation(eulers2);
                gChild2.getMoleculeAtPosition(mol).setOrientation(eulers1);
            }
        }
        
        return new Tuple<>(gChild1, gChild2);
    }

    @Override
    public short hasPriority() {
        return 0;
    }
}
