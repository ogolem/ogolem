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
package org.ogolem.core;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 * 
 * @author Johannes Dieterich
 * @version 2020-02-01
 */
public class AdvancedPairWiseTest {
   
    @Test
    public void testCheckOnlyForCollision() {
        
        final int noLJAtoms = 4;
        final int[] atsPerMol = new int[noLJAtoms];
        
        for(int i = 0; i < noLJAtoms; i++){
            atsPerMol[i] = 1;
        }
        
        final CartesianCoordinates cartes = new CartesianCoordinates(noLJAtoms, noLJAtoms, atsPerMol);
        final String[] atoms = cartes.getAllAtomTypes();
        final double[][] xyz = cartes.getAllXYZCoord();
        
        for(int i = 0; i < noLJAtoms; i++){
            atoms[i] = "Ar";
            xyz[0][i] = 4*i;
        }
        cartes.recalcAtomNumbersForced();
        
        final BondInfo bonds = new SimpleBondInfo(4);
        
        final AdvancedPairWise cd = new AdvancedPairWise(true, new DummyCollisionStrengthComputer());
        
        System.out.println("check only for collision: no collision");
        final boolean res1 = cd.checkOnlyForCollision(cartes, 1.0, bonds, 2, 3);
        assertTrue(!res1); // no collision

        System.out.println("check only for collision: no collision in window");
        xyz[0][1] = 0.0;
        final boolean res2 = cd.checkOnlyForCollision(cartes, 1.0, bonds, 2, 3);
        xyz[0][1] = 4.0;
        assertTrue(!res2); // no collision
        
        System.out.println("check only for collision: collision with outside window");
        xyz[0][3] = 4*2.0;
        final boolean res3 = cd.checkOnlyForCollision(cartes, 1.0, bonds, 2, 3);
        xyz[0][3] = 4*3.0;
        assertTrue(res3); // collision
        
        System.out.println("check only for collision: collision with outside window 2");
        xyz[0][0] = 8.0;
        final boolean res4 = cd.checkOnlyForCollision(cartes, 1.0, bonds, 2, 3);
        xyz[0][0] = 0.0;
        assertTrue(res4); // collision
        
        System.out.println("check only for collision: collision inside window");
        xyz[0][2] = 4*3.0;
        final boolean res5 = cd.checkOnlyForCollision(cartes, 1.0, bonds, 2, 3);
        xyz[0][2] = 4*2.0;
        assertTrue(res5); // collision
        
        System.out.println("check only for collision: collision inside window but bond");
        xyz[0][2] = 4*3.0;
        bonds.setBond(2,3,(short)1);
        final boolean res6 = cd.checkOnlyForCollision(cartes, 1.0, bonds, 2, 3);
        xyz[0][2] = 4*2.0;
        assertTrue(!res6); // no collision
    }
}
