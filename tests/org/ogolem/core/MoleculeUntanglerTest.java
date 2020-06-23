/**
 Copyright (c) 2020, D. Behrens, J. M. Dieterich, B. Hartke
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 * Neither the name of the copyright holder nor the names of its
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

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
import org.ogolem.locopt.LBFGSLocOpt;

import java.util.ArrayList;

/**
 * A simple test for the Untangler. For now, it simply runs through the methods and uses the Untangler in a very
 * simple L-BFGS Locopt, just to make sure everything runs smoothly. No benchmarking here.
 */
public class MoleculeUntanglerTest {
    @Test
    public void testTerms() {
        MoleculeUntangler instance = new MoleculeUntangler(1.0, true);

        double[][] distances = new double[3][];
        // Three corner cases using 1-atom molecules
        //  * The atoms are exceptionally close (same coordinates in double precision)
        distances[0] = new double[]{0, 0, 0, 0, 0, 0};
        //  * The atoms are in between lower and upper limits
        distances[1] = new double[]{0, 5, 0, 0, 0, 0};
        //  * The atoms are within summedRadius distance
        distances[2] = new double[]{0, 5.14 * Constants.ANGTOBOHR, 0, 0, 0, 0};

        // Expected Values
        double[] energies = new double[3];
        double[] gradx = new double[3];
        energies[0] = 1.9999999999221;
        energies[1] = 1.1922724326284E-4;
        energies[2] = 0;
        gradx[0] = 7.78210116701E-6;
        gradx[1] = 2.3195961724288E-4;
        gradx[2] = 0;

        // Basic stuff we'll need to call the methods
        int noOfAtoms = 2;
        int[] atsPerMol = new int[]{1, 1};
        // Two gold atoms
        short[] atomNos = new short[]{(short) 79, (short) 79};

        for (int i = 0; i < 3; i++) {
            Gradient grad = new Gradient(3, noOfAtoms);
            instance.gradientCalculation(0, 42, distances[i], null, atomNos, atsPerMol, null, noOfAtoms,
                    null, null, null, grad);
            assertEquals(energies[i], grad.getTotalEnergy(), 1e-12);
            assertEquals(gradx[i], grad.getGradient()[0], 1e-12);
        }
    }

    @Test
    public void testLocOpt() {
        // Build an actual example using an OGOLEM-style locopt cycle.
        MoleculeUntangler untangler = new MoleculeUntangler(1.0, true);

        CartesianCoordinates coords = new CartesianCoordinates(1, 1, new int[]{1});
        coords.setAllAtomTypes(new String[]{"Au"});
        coords.setAtomNumbers(new short[]{(short) 79});

        Molecule mol1 = new Molecule(coords, 0, "au0");
        Molecule mol2 = new Molecule(coords, 1, "au1");
        Molecule mol3 = new Molecule(coords, 2, "au2");
        mol3.setExternalCenterOfMass(0, 5.14 * Constants.ANGTOBOHR, 0);
        GeometryConfig gc = new GeometryConfig();
        gc.geomMCs = new ArrayList<>();
        gc.geomMCs.add(mol1.returnMyConfig());
        gc.geomMCs.add(mol2.returnMyConfig());
        gc.geomMCs.add(mol3.returnMyConfig());
        gc.bonds = new SimpleBondInfo(3);
        gc.noOfParticles = 3;
        Geometry geom = new Geometry(gc);

        FullyCartesianCoordinates xyz = new FullyCartesianCoordinates(untangler);
        LBFGSLocOpt<Molecule, Geometry> LocOpt = new LBFGSLocOpt<>(xyz, 7, 400, 1e-15, 0.9, 50, 0);
        Geometry opt = LocOpt.fitness(geom, false);

        // Just make sure the L-BFGS converged fine.
        assert (opt.getFitness() < 1e-14);
    }
}
