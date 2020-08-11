/**
Copyright (c) 2016-2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive;

import org.junit.Test;
import static org.junit.Assert.*;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CollisionDetection;
import org.ogolem.core.DissociationDetection;
import org.ogolem.core.Geometry;

import static org.ogolem.core.Constants.ANGTOBOHR;
import static org.ogolem.core.Constants.EVTOHARTREE;

import java.util.ArrayList;
import java.util.List;

import org.ogolem.core.Gradient;
import org.ogolem.core.NorwayGeometryMutation;
import org.ogolem.core.SimpleBondInfo;
import org.ogolem.core.Topology;
import org.ogolem.random.Lottery;
import org.ogolem.random.StandardRNG;

/**
*
* @author Johannes Dieterich
* @version 2020-08-09
*/
public class AdaptiveSWGFFTest {

    @Test
    public void testCartesianGradient() {
        
        System.out.println("CartesianGradient Si6 with caching");
        final AdaptiveParameters params = setupParam();
        final AdaptiveSWGFF instance = new AdaptiveSWGFF(false, params, true, 1.4);
        
        // assemble a random set of 6 silicon atoms
        final int noAtoms = 6;
        final BondInfo bonds = new SimpleBondInfo(noAtoms);
        final int[] atsPerMol = new int[noAtoms];
        final boolean[] molFlexies = new boolean[noAtoms];
        final List<boolean[][]> allFlexies = new ArrayList<>();
        final boolean[] molConstraints = new boolean[noAtoms];
        final boolean[][] constrXYZ = new boolean[3][noAtoms];
        final String[] sids = new String[noAtoms];
        
        for(int i = 0; i < noAtoms; i++){
            atsPerMol[i] = 1;
            molFlexies[i] = false;
            molConstraints[i] = false;
            sids[i] = "Si";
            constrXYZ[0][i] = false;
            constrXYZ[1][i] = false;
            constrXYZ[2][i] = false;
            allFlexies.add(null);
        }
        
        CartesianCoordinates cartes = new CartesianCoordinates(noAtoms, noAtoms, atsPerMol);
        final String[] atoms = cartes.getAllAtomTypes();
        
        for(int i = 0; i < noAtoms; i++){
            atoms[i] = "Si";
        }
        
        cartes.recalcAtomNumbersForced();

        Lottery.setGenerator(new StandardRNG(142));
        final NorwayGeometryMutation norway = new NorwayGeometryMutation(NorwayGeometryMutation.PACKDIM.THREED, CollisionDetection.CDTYPE.SIMPLEPAIRWISE, 1.4, 4.0,
                DissociationDetection.DEFAULTDD, NorwayGeometryMutation.MUTMODE.ASCENDING);
        
        final Geometry geom = new Geometry(cartes, 0, noAtoms, atsPerMol, molFlexies, allFlexies, molConstraints,
                constrXYZ, sids, bonds);
        final Geometry packed = norway.mutate(geom);
        
        cartes = packed.getCartesians();
        
        final double[] energyparts = new double[noAtoms];
        final Gradient g = new Gradient(3, noAtoms);
        instance.gradientCalculation(0, 0, cartes.getAll1DCartes(), atoms, cartes.getAllAtomNumbers(),
        		atsPerMol, energyparts, noAtoms, cartes.getAllCharges(), cartes.getAllSpins(), bonds, g, false);
        
        final Gradient numG = NumericalGradients.numericalGradient(0, 0, cartes.getAll1DCartes(), atoms, cartes.getAllAtomNumbers(),
        		atsPerMol, energyparts, noAtoms, cartes.getAllCharges(), cartes.getAllSpins(), bonds, instance, false);
        
        System.out.println("E reference " + g.getTotalEnergy());
        System.out.println("E numerical " + numG.getTotalEnergy());
        
        assertEquals(g.getTotalEnergy(), numG.getTotalEnergy(), 1e-8);
        
        final double[] grad1D = g.getGradient();
        for(int i = 0; i < grad1D.length; i++) {
        	System.out.println("Analytical flattened gradient element " + i + " " + grad1D[i]);
        }
        
        final boolean same = g.compare(numG, 1e-6, true);
        assertTrue(same);
    }
    
    @Test
    public void testCartesianGradient42() {
        
        System.out.println("CartesianGradient Si42 with caching");
        final AdaptiveParameters params = setupParam();
        final AdaptiveSWGFF instance = new AdaptiveSWGFF(false, params, true, 1.4);
        
        // assemble a random set of 42 silicon atoms
        final int noAtoms = 42;
        final BondInfo bonds = new SimpleBondInfo(noAtoms);
        final int[] atsPerMol = new int[noAtoms];
        final boolean[] molFlexies = new boolean[noAtoms];
        final List<boolean[][]> allFlexies = new ArrayList<>();
        final boolean[] molConstraints = new boolean[noAtoms];
        final boolean[][] constrXYZ = new boolean[3][noAtoms];
        final String[] sids = new String[noAtoms];
        
        for(int i = 0; i < noAtoms; i++){
            atsPerMol[i] = 1;
            molFlexies[i] = false;
            molConstraints[i] = false;
            sids[i] = "Si";
            constrXYZ[0][i] = false;
            constrXYZ[1][i] = false;
            constrXYZ[2][i] = false;
            allFlexies.add(null);
        }
        
        CartesianCoordinates cartes = new CartesianCoordinates(noAtoms, noAtoms, atsPerMol);
        final String[] atoms = cartes.getAllAtomTypes();
        
        for(int i = 0; i < noAtoms; i++){
            atoms[i] = "Si";
        }
        
        cartes.recalcAtomNumbersForced();

        Lottery.setGenerator(new StandardRNG(141));
        final NorwayGeometryMutation norway = new NorwayGeometryMutation(NorwayGeometryMutation.PACKDIM.THREED, CollisionDetection.CDTYPE.SIMPLEPAIRWISE, 1.4, 4.0,
                DissociationDetection.DEFAULTDD, NorwayGeometryMutation.MUTMODE.ASCENDING);
        
        final Geometry geom = new Geometry(cartes, 0, noAtoms, atsPerMol, molFlexies, allFlexies, molConstraints,
                constrXYZ, sids, bonds);
        final Geometry packed = norway.mutate(geom);
        
        cartes = packed.getCartesians();
        
        final double[] energyparts = new double[noAtoms];
        final Gradient g = new Gradient(3, noAtoms);
        instance.gradientCalculation(0, 0, cartes.getAll1DCartes(), atoms, cartes.getAllAtomNumbers(),
        		atsPerMol, energyparts, noAtoms, cartes.getAllCharges(), cartes.getAllSpins(), bonds, g, false);
        
        final Gradient numG = NumericalGradients.numericalGradient(0, 0, cartes.getAll1DCartes(), atoms, cartes.getAllAtomNumbers(),
        		atsPerMol, energyparts, noAtoms, cartes.getAllCharges(), cartes.getAllSpins(), bonds, instance, false);
        
        System.out.println("E reference " + g.getTotalEnergy());
        System.out.println("E numerical " + numG.getTotalEnergy());
        
        assertEquals(g.getTotalEnergy(), numG.getTotalEnergy(), 1e-8);
        
        final double[] grad1D = g.getGradient();
        for(int i = 0; i < grad1D.length; i++) {
        	System.out.println("Analytical flattened gradient element " + i + " " + grad1D[i]);
        }
        
        final boolean same = g.compare(numG, 1e-6, true);
        assertTrue(same);
    }
	
    private AdaptiveParameters setupParam(){
        
        final int numberOfParameters = 9;
        final int id = -1;
        
        final String[] keys = new String[3];
        keys[0] = "adaptiveswg3b2b:SiSi";
        keys[1] = "adaptiveswg3b:SiSiSi";
        keys[2] = "adaptiveswg2b:SiSi";
        
        final int[] noOfParamsPerKey = new int[]{2,3,4};
        
        final String paramsForMethod = "ADAPTIVESWGTESTONLY";
        
        final AdaptiveParameters params = new AdaptiveParameters(numberOfParameters, id, keys,
            noOfParamsPerKey, paramsForMethod);
        
        final int offset1 = params.getStartPointForKey(keys[0]);
        final double[] p = params.getAllParamters();
        
        p[offset1] = 3.77*ANGTOBOHR;
        p[offset1+1] = 2.4*ANGTOBOHR;
        
        final int offset2 = params.getStartPointForKey(keys[1]);
        p[offset2] = -0.5;
        p[offset2+1] = 0.15;
        p[offset2+2] = 4.0*EVTOHARTREE;
        
        final int offset3 = params.getStartPointForKey(keys[2]);
        p[offset3] = 3.77*ANGTOBOHR;
        p[offset3+1] = 16.3*EVTOHARTREE;
        p[offset3+2] = 11.581*ANGTOBOHR*ANGTOBOHR*ANGTOBOHR*ANGTOBOHR;
        p[offset3+3] = 2.095*ANGTOBOHR;
        
        return params;
    }
}
