/**
Copyright (c) 2012-2014, J. M. Dieterich
              2015-2016, J. M. Dieterich and B. Hartke
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

/**
 * A partial energy based offspring of Sweden. No political statement intended.
 * Tries to optimize the cutting plane wrt the partial energy contributions.
 * Uses a center of quality approach (see BXH papers for that).
 * @author Johannes Dieterich
 * @version 2016-08-29
 */
public class AlandGeometryXOver implements GenericCrossover<Molecule,Geometry>{
    
    private static final long serialVersionUID = (long) 20140328;
    private static final boolean DEBUG = false;
    private final boolean doMergingGlobOpt;
    private final PhenotypeGeometryXOver pheno;
    private final MergingPhenoXOver merger;
    private final CartesianFullBackend backend;
    
    AlandGeometryXOver(final boolean doMerging, final CartesianFullBackend back,
            final MergingPhenoXOver.MergingPhenoConfig config, final double blowColl,
            final double blowDiss, final CollisionDetectionEngine collDetect, 
            final int whichGlobOpt, final boolean doCDDD, final boolean doRandomTrans){
        if(!doMerging){
            pheno = new PhenotypeGeometryXOver(whichGlobOpt);
            merger = null;
            backend = back;
        } else{
            pheno = null;
            backend = back;
            merger = new MergingPhenoXOver(back, blowColl, blowDiss,collDetect,config,doCDDD,doRandomTrans);
        }
        this.doMergingGlobOpt = doMerging;
    }
    
    AlandGeometryXOver(final AlandGeometryXOver orig){
        this.doMergingGlobOpt = orig.doMergingGlobOpt;
        this.merger = (orig.merger != null) ? orig.merger.clone() : null;
        this.pheno = (orig.pheno != null) ? orig.pheno.clone() : null;
        this.backend = (orig.backend != null) ? orig.backend.clone() : null;
    }
    
    @Override
    public GenericCrossover<Molecule, Geometry> clone() {
        return new AlandGeometryXOver(this);
    }

    @Override
    public String getMyID() {
        return "FITNESS BASED PHENOTYPE\n\t";
    }

    @Override
    public Tuple<Geometry, Geometry> crossover(Geometry mother, Geometry father, final long futureID) {        
        
        /*if(false){
            // shortcut...
            final CartesianCoordinates c = father.getCartesians();
            final double[][] coords = c.getAllXYZCoord();
            final double[] com = c.calculateTheCOM();
            System.out.println("DEBUG: COM " + com[0] + "\t" + com[1] + "\t" + com[2]);
            
            final String[] sa = c.createPrintableCartesians();
            for(final String s : sa) System.out.println(s);
            
            final double[][] rot = CoordTranslation.rotatePointToZAxis(coords, com, c.getNoOfAtoms());
            
            c.setAllXYZ(rot);
            
            final String[] sa2 = c.createPrintableCartesians();
            for(final String s : sa2) System.out.println(s);
            
            final double[] com2 = c.calculateTheCOM();
            System.out.println("DEBUG: COM2 " + com2[0] + "\t" + com2[1] + "\t" + com2[2]);
            
            return new Tuple<>(null,null);
        }*/
        
        
        if(mother.getNumberOfIndieParticles() == 1){

            final Geometry gChild1 = new Geometry(mother);
            final Geometry gChild2 = new Geometry(father);
            
            return new Tuple<>(gChild1,gChild2);
        }
        
        // create the children
        Geometry gChild1 = new Geometry(father);
        gChild1.setID(futureID);
        Geometry gChild2 = new Geometry(mother);
        gChild2.setID(futureID);

        final CartesianCoordinates cFather = gChild1.getCartesians();
        final CartesianCoordinates cMother = gChild2.getCartesians();

        cFather.moveCoordsToCOM();
        cMother.moveCoordsToCOM();
        
        if(DEBUG){
            final double[] comF = cFather.calculateTheCOM();
            final double[] comM = cMother.calculateTheCOM();
            System.out.println("DEBUG: CoM father after move: " + comF[0] + "\t" + comF[1] + "\t" + comF[2]);
            System.out.println("DEBUG: CoM mother after move: " + comM[0] + "\t" + comM[1] + "\t" + comM[2]);
        }
        
        CoordTranslation.updateGeometryFromCartesian(cFather, gChild1);
        CoordTranslation.updateGeometryFromCartesian(cMother, gChild2);
                
        final int noMols = cMother.getNoOfMolecules();
        final int noAts = cMother.getNoOfAtoms();
        assert(noMols == cFather.getNoOfMolecules());
        assert(noAts == cFather.getNoOfAtoms());
        final double[] ePFather = new double[noMols];
        final double[] ePMother = new double[noMols];
        assert(ePFather.length == ePMother.length);
        
        if(backend == null){throw new RuntimeException("In merging globopt w/ no backend. This is just wrong.");}
        final double eF = backend.energyCalculation(father.getID(), 0, cFather.getAll1DCartes(),
                cFather.getAllAtomTypes(), cFather.getAllAtomNumbers(), cFather.getAllAtomsPerMol(),
                ePFather, noAts, cFather.getAllCharges(), cFather.getAllSpins(),
                father.getBondInfo());
        final double eM = backend.energyCalculation(mother.getID(), 0, cMother.getAll1DCartes(),
                cMother.getAllAtomTypes(), cMother.getAllAtomNumbers(), cMother.getAllAtomsPerMol(),
                ePMother, noAts, cMother.getAllCharges(), cMother.getAllSpins(),
                mother.getBondInfo());
        
        // translate the partial contributions to positive values (only then the CoQ is defined)
        //TODO this might not hold true for all cases. CoQ should be located, where the energy contributions are highest. In the context of negative
        // interaction energies, this translates to "where they are most negative". Due to the shifting, I think that the CoQ is inverted.
        final double[] ePFTmp = ePFather.clone();
        final double[] ePMTmp = ePMother.clone();
        double mF = Double.MAX_VALUE;
        double mM = Double.MAX_VALUE;
        for(int i = 0; i < noMols; i++){
            mF = Math.min(mF,ePFather[i]);
            mM = Math.min(mM,ePMother[i]);
            if(DEBUG){System.out.println("DEBUG: Partial energy for mol " + i + " of father " + ePFather[i] + " and mother " + ePMother[i]);}
        }
        
        mF = 1.01*Math.abs(mF);
        mM = 1.01*Math.abs(mM);
        
        for(int i = 0; i < noMols; i++){
            ePFTmp[i] += mF;
            ePMTmp[i] += mM;
        }
        
        if(DEBUG) printCartes("Mother BEFORE CoQ rotation", cMother, ePMTmp);
        if(DEBUG) printCartes("Father BEFORE CoQ rotation", cFather, ePFTmp);
        
        // calculate the centres of quality
        final double[] coqF = new double[3];
        final double[] coqM = new double[3];
        for(int i = 0; i < noMols; i++){
            
            final double[] comF = gChild1.getCOM(i);
            final double[] comM = gChild2.getCOM(i);
            
            coqF[0] += comF[0]* ePFTmp[i];
            coqF[1] += comF[1]* ePFTmp[i];
            coqF[2] += comF[2]* ePFTmp[i];
            
            coqM[0] += comM[0]* ePMTmp[i];
            coqM[1] += comM[1]* ePMTmp[i];
            coqM[2] += comM[2]* ePMTmp[i];
            
        }
        
        if(DEBUG){
            System.out.println("DEBUG: CoQ father before rotation: " + coqF[0] + "\t" + coqF[1] + "\t" + coqF[2]);
            System.out.println("DEBUG: CoQ mother before rotation: " + coqM[0] + "\t" + coqM[1] + "\t" + coqM[2]);
        }
        
        // rotate the coqs onto the z-axis
        final double[][] rotF = CoordTranslation.rotatePointToZAxis(cFather.getAllXYZCoord(), coqF, noAts);
        final double[][] rotM = CoordTranslation.rotatePointToZAxis(cMother.getAllXYZCoord(), coqM, noAts);     
        
        cFather.setAllXYZ(rotF);
        cMother.setAllXYZ(rotM);
        
        if(DEBUG) printCartes("Mother AFTER CoQ rotation", cMother, ePMTmp);
        if(DEBUG) printCartes("Father AFTER CoQ rotation", cFather, ePFTmp);
                

        // translate it to geometries again
        CoordTranslation.updateGeometryFromCartesian(cFather, gChild1);
        CoordTranslation.updateGeometryFromCartesian(cMother, gChild2);
                
        // check the CoQs
        final double[] coqF2 = new double[3];
        final double[] coqM2 = new double[3];
        for(int i = 0; i < noMols; i++){
            
            final double[] comF = gChild1.getCOM(i);
            final double[] comM = gChild2.getCOM(i);
            
            coqF2[0] += comF[0]* ePFTmp[i];
            coqF2[1] += comF[1]* ePFTmp[i];
            coqF2[2] += comF[2]* ePFTmp[i];
            
            coqM2[0] += comM[0]* ePMTmp[i];
            coqM2[1] += comM[1]* ePMTmp[i];
            coqM2[2] += comM[2]* ePMTmp[i];
            
        }
                
        if(DEBUG){
            System.out.println("DEBUG: CoQ father after rotation: " + coqF2[0] + "\t" + coqF2[1] + "\t" + coqF2[2]);
            System.out.println("DEBUG: CoQ mother after rotation: " + coqM2[0] + "\t" + coqM2[1] + "\t" + coqM2[2]);
        }
        
        if(coqF2[2] > 0.0){
            cFather.mirrorCoordinates(2);
            CoordTranslation.updateGeometryFromCartesian(cFather, gChild1);
            if(DEBUG) System.out.println("DEBUG: Father had positive z component of CoQ, mirrored.");
        }
        
        if(coqM2[2] < 0.0){
            cMother.mirrorCoordinates(2);
            CoordTranslation.updateGeometryFromCartesian(cMother, gChild2);
            if(DEBUG) System.out.println("DEBUG: Mother had negative z component of CoQ, mirrored.");
        }

        // set the IDs of mother and father
        gChild1.setMother(mother.getID());
        gChild1.setFather(father.getID());
        gChild2.setMother(mother.getID());
        gChild2.setFather(father.getID());


        // in case of environments, add it in at this point
        if (mother.containsEnvironment()) {
            // the environment
            final List<Environment> envChilds = mother.getEnvironment().createOffspring(father.getEnvironment());
            gChild1.setEnvironment(envChilds.get(0));
            gChild2.setEnvironment(envChilds.get(1));
        }

        if(doMergingGlobOpt){
            return merger.crossover(gChild1, gChild2, futureID);
        } else{
            return pheno.crossover(gChild1, gChild2, futureID);
        }
    }
    
    @Override
    public short hasPriority() {
        return 0;
    }
    
    private static void printCartes(final String id, final CartesianCoordinates c, final double[] partial){
        
        if(partial.length != c.getNoOfAtoms()) return;
        
        final String[] sa = c.createPrintableCartesians();
        System.out.println("####################################################");
        System.out.println("####################################################");
        System.out.println("" + id  + " coming now.... ");
        System.out.println("####################################################");
        System.out.println(sa[0]);
        System.out.println(sa[1]);
        for(int i = 0; i < sa.length-2; i++){
            System.out.println(sa[i+2] + " \t" + partial[i]);
        }
        System.out.println("####################################################");
        System.out.println("####################################################");
    }
}
