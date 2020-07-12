/**
Copyright (c) 2014, J. M. Dieterich
              2014-2016, J. M. Dieterich and B. Hartke
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
import static org.ogolem.core.GlobOptAtomics.randomRotation;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;

/**
 * This is a straightforward phenotype algorithm parting the geometry (and molecule)
 * into two distinct parts and putting it back together. It returns the fitter of
 * the two resulting children after mutating them.
 * Any rumors concerning this algorithm running better on alcohol are of course
 * exaggerated! But they might not be completely false... ;-)
 * Named in honor of Dr Mats Eriksson, the hardest rocker of them all.
 * @author Johannes Dieterich
 * @version 2016-08-29
 */
public class SwedenGeometryXOver implements GenericCrossover<Molecule,Geometry> {
    
    private static final long serialVersionUID = (long) 20160403;
    private static final boolean DEBUG = false;
    private final PhenotypeGeometryXOver pheno;
    private final GlobOptAtomics.CUTTINGMODE whichGlobOpt;
    
    SwedenGeometryXOver(final GlobOptAtomics.CUTTINGMODE whichGlobOpt){
        this.pheno = new PhenotypeGeometryXOver(whichGlobOpt);
        this.whichGlobOpt = whichGlobOpt;
    }
    
    SwedenGeometryXOver(final SwedenGeometryXOver orig){
        this.pheno = orig.pheno.clone();
        this.whichGlobOpt = orig.whichGlobOpt;
    }
    
    @Override
    public SwedenGeometryXOver clone() {
        return new SwedenGeometryXOver(this);
    }

    @Override
    public String getMyID() {
        return "SWEDEN GEOMETRY XOVER:\n" + pheno.getMyID();
    }

    @Override
    public Tuple<Geometry, Geometry> crossover(final Geometry mother, final Geometry father, final long futureID) {
        
        if(mother.getNumberOfIndieParticles() == 1){

            final Geometry gChild1 = new Geometry(mother);
            final Geometry gChild2 = new Geometry(father);

            return new Tuple<>(gChild1,gChild2);
        }

        // create the children
        Geometry gChild1 = new Geometry(father);
        Geometry gChild2 = new Geometry(mother);

        final CartesianCoordinates cartesMother = gChild2.getCartesians();
        final CartesianCoordinates cartesFather = gChild1.getCartesians();

        cartesMother.moveCoordsToCOM();
        cartesFather.moveCoordsToCOM();
        
        if(DEBUG){
            final String[] printCart1 = cartesMother.createPrintableCartesians();
            final String[] printCart2 = cartesFather.createPrintableCartesians();
            System.out.println("DEBUG: Cartesians prior to rotation: mother");
            for(final String s : printCart1){
                System.out.println(s);
            }
            System.out.println("DEBUG: Cartesians prior to rotation: father");
            for(final String s : printCart2){
                System.out.println(s);
            }
        }

        // rotate the two
        randomRotation(cartesMother);
        randomRotation(cartesFather);
        
        if(DEBUG){
            final String[] printCart1 = cartesMother.createPrintableCartesians();
            final String[] printCart2 = cartesFather.createPrintableCartesians();
            System.out.println("DEBUG: Cartesians after rotation: mother");
            for(final String s : printCart1){
                System.out.println(s);
            }
            System.out.println("DEBUG: Cartesians after rotation: father");
            for(final String s : printCart2){
                System.out.println(s);
            }
        }

         // translate it to geometries again
        gChild1.updateAllCoordinates(cartesFather.getAllXYZCoord());
        gChild2.updateAllCoordinates(cartesMother.getAllXYZCoord());

         // set the IDs of mother and father
         gChild1.setMother(mother.getID());
         gChild1.setFather(father.getID());
         gChild2.setMother(mother.getID());
         gChild2.setFather(father.getID());


        // in case of environments, add it in at this point
        if(mother.containsEnvironment()){
            // the environment
            final List<Environment> envChilds = mother.getEnvironment().createOffspring(father.getEnvironment());
            gChild1.setEnvironment(envChilds.get(0));
            gChild2.setEnvironment(envChilds.get(1));
        }

        return pheno.crossover(gChild1, gChild2, futureID);
    }
    
    @Override
    public short hasPriority() {
        return -1;
    }
}
