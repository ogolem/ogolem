/**
Copyright (c) 2014, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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

import org.ogolem.generic.GenericAbstractDarwin;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericMutation;
import org.ogolem.generic.GenericSanityCheck;
import org.ogolem.generic.IndividualWriter;
import org.ogolem.helpers.Tuple;

/**
 * A generified (or de-generified?) Darwin implementation based off GenericAbstractDarwin.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class GenericGeometryDarwin extends GenericAbstractDarwin<Molecule,Geometry> {
    
    private static final long serialVersionUID = (long) 20200429;
    private final GenericCrossover<Double,Molecule> molXOver;
    private final GenericMutation<Double,Molecule> molMutation;
    private final double molXOverProb;
    private final double molMutProb;

    GenericGeometryDarwin(final GenericCrossover<Molecule,Geometry> cross, final GenericMutation<Molecule,Geometry> mut,
            final GenericSanityCheck<Molecule,Geometry> sanity, final GenericFitnessFunction<Molecule,Geometry> fitness,
            final IndividualWriter<Geometry> writer, final double crossPoss, final double mutPoss,
            final boolean printBeforeFitness, final int noOfTries, final double molXOverProb, final double molMutProb,
            final GenericCrossover<Double,Molecule> molXOver, final GenericMutation<Double,Molecule> molMutation){
        super(cross, mut, sanity, fitness, writer, crossPoss, mutPoss,
            printBeforeFitness, noOfTries);
        assert(molXOverProb >= 0.0 && molXOverProb <= 1.0);
        assert(molMutProb >= 0.0 && molMutProb <= 1.0);
        this.molMutProb = molMutProb;
        this.molXOverProb = molXOverProb;
        this.molMutation = molMutation;
        this.molXOver = molXOver;
    }
    
    GenericGeometryDarwin(final GenericGeometryDarwin orig){
        super(orig);
        this.molMutProb = orig.molMutProb;
        this.molXOverProb = orig.molXOverProb;
        this.molMutation = orig.molMutation.clone();
        this.molXOver = orig.molXOver.clone();
    }
    
    @Override
    public GenericGeometryDarwin copy() {
        return new GenericGeometryDarwin(this);
    }

    @Override
    public String getMyID() {
        return "GENERIC GEOMETRY DARWIN IMPLEMENTATION"
                + "\nusing (for geometry): "
                + "\n\txover    " + xover.getMyID() + " which probability " + crossPoss*100 + "%"
                + "\n\tmutation " + mutation.getMyID() + " which probability " + mutPoss*100 + "%"
                + "\nusing (for molecule): "
                + "\n\txover    " + molXOver.getMyID() + " which probability " + molXOverProb*100 + "%"
                + "\n\tmutation " + molMutation.getMyID() + " which probability " + molMutProb*100 + "%"
                + "\nfitness  " + fitness.getMyID();
    }

    @Override
    protected void postXOver(final Geometry individual1, final Geometry individual2, final long futureID) {
        assert(individual1.isThereAFlexy() == individual2.isThereAFlexy());
        if(individual1.isThereAFlexy()){
            // we need to try molecular crossover
            final int noIndies = individual1.getNumberOfIndieParticles();
            assert(noIndies == individual2.getNumberOfIndieParticles());
            for(int i = 0; i < individual1.getNumberOfIndieParticles(); i++){
                
                final Molecule mol1 = individual1.getMoleculeAtPosition(i);
                final Molecule mol2 = individual2.getMoleculeAtPosition(i);
                assert(mol1.getFlexy() == mol2.getFlexy());
                if(mol1.getFlexy() && r.nextDouble() <= molXOverProb){
                    
                    final Tuple<Molecule,Molecule> res = molXOver.crossover(mol1, mol2, futureID);
                    individual1.setMoleculeAtPosition(i, res.getObject2());
                    individual2.setMoleculeAtPosition(i, res.getObject1());
                }
            }
        }
    }

    @Override
    protected void postMutation(final Geometry individual) {
        if(individual.isThereAFlexy()){
            // we need to try molecular mutation
            for(int i = 0; i < individual.getNumberOfIndieParticles(); i++){
                
                final Molecule mol = individual.getMoleculeAtPosition(i);
                if(mol.getFlexy() && r.nextDouble() <= molMutProb){
                    final Molecule mutated = molMutation.mutate(mol);
                    individual.setMoleculeAtPosition(i, mutated);
                }
            }
        }
    }

    @Override
    protected void runAfterEachTry() {
        // does not need to be implemented
    }
}
