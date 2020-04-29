/**
Copyright (c) 2014, J. M. Dieterich
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

import java.util.ArrayList;
import java.util.List;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.RandomUtils;

/**
 * n-Point genotype crossover for geometries.
 * @author Johannes Dieterich
 * @version 2014-03-27
 */
public class PortugalGeometryXOver implements GenericCrossover<Molecule,Geometry> {
    
    private static final long serialVersionUID = (long) 20140327;
    private final int noCuts;
    
    PortugalGeometryXOver(final int noCuts){
        this.noCuts = noCuts;
    }
    
    PortugalGeometryXOver(final PortugalGeometryXOver orig){
        this.noCuts = orig.noCuts;
    }

    @Override
    public PortugalGeometryXOver clone() {
        return new PortugalGeometryXOver(this);
    }

    @Override
    public String getMyID() {
        return "PORTUGAL\nn-point geometry genotype crossover\n\t #crosses: " + noCuts;
    }

    @Override
    public Tuple<Geometry, Geometry> crossover(Geometry mother, Geometry father, final long futureID) {
        
        final int noOfIndies = mother.getNumberOfIndieParticles();

        if(noOfIndies <= noCuts){
            throw new RuntimeException("Too many cuts for number of molecules in geometry. " + noCuts + " vs " + noOfIndies);
        }
        
        final List<Integer> cutPoints = RandomUtils.listOfPoints(noCuts, noOfIndies);
                
        // the geometry configurations
        final GeometryConfig gcChildOne = new GeometryConfig();
        final GeometryConfig gcChildTwo = new GeometryConfig();

        final GeometryConfig gcMother = mother.returnMyConfig();
        final GeometryConfig gcFather = father.returnMyConfig();

        //set some small values
        gcChildOne.motherID = gcMother.lID;
        gcChildTwo.motherID = gcMother.lID;
        gcChildOne.fatherID = gcFather.lID;
        gcChildTwo.fatherID = gcFather.lID;
        gcChildOne.noOfParticles = noOfIndies;
        gcChildTwo.noOfParticles = noOfIndies;

        // the actual chopping of MC's...
        final List<MoleculeConfig> alMCOne = new ArrayList<>(noOfIndies);
        final List<MoleculeConfig> alMCTwo = new ArrayList<>(noOfIndies);
        
        boolean mothFirst = true;
        for(int mol = 0; mol < noOfIndies; mol++){
            
            if(mothFirst){
                alMCOne.add(gcMother.geomMCs.get(mol));
                alMCTwo.add(gcFather.geomMCs.get(mol));
            } else{
                alMCTwo.add(gcMother.geomMCs.get(mol));
                alMCOne.add(gcFather.geomMCs.get(mol));
            }
                        
            if(cutPoints.contains(mol)) {
                mothFirst = !mothFirst;
            }
        }
        
        // put the MCs into the GC
        gcChildOne.geomMCs = alMCOne;
        gcChildTwo.geomMCs = alMCTwo;

        gcChildOne.bonds = mother.getBondInfo().clone();
        gcChildTwo.bonds = father.getBondInfo().clone();

        if(mother.containsEnvironment()){
            // the environment
            final List<Environment> envChilds = gcChildOne.env.createOffspring(gcChildTwo.env);
            gcChildOne.env = envChilds.get(0);
            gcChildTwo.env = envChilds.get(1);
        }
        
        // create the two children
        final Geometry gChildOne = new Geometry(gcChildOne);
        final Geometry gChildTwo = new Geometry(gcChildTwo);
        
        //TODO test and debug
        return new Tuple<>(gChildOne,gChildTwo);
    }

    @Override
    public short hasPriority() {
        return -1;
    }
}
