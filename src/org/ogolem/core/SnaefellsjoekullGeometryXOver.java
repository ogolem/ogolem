/**
Copyright (c) 2014, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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

import contrib.bobyqa.AbstractBOBYQAMethod;
import contrib.bobyqa.BOBYQAOptimizer;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.apache.commons.math3.util.FastMath;
import static org.ogolem.core.GlobOptAtomics.assignMolecularTypes;
import static org.ogolem.core.GlobOptAtomics.findOptimalSphereRadius;
import static org.ogolem.core.GlobOptAtomics.randomRotation;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.RandomUtils;
import org.ogolem.helpers.Tuple;

/**
 * A sphere-based phenotype cut.
 * @author Johannes Dieterich
 * @version 2020-02-09
 */
public class SnaefellsjoekullGeometryXOver implements GenericCrossover<Molecule,Geometry> {

    public static enum CUTTYPE {FULLRANDOM, GAUSSDISTRIBUTED, INVERTEDGAUSSDISTRIBUTED};
    
    private static final long serialVersionUID = (long) 20200209;
    private static final boolean DEBUG = false;
    private final boolean intermediateSanityChecks;
    private final Random r;
    private final boolean adjustRadius;
    private final CUTTYPE cut;
    private final double gaussWidth;
    private final boolean doRigidOpt;
    private final boolean checkForCollAndInflate;
    private final boolean doCD;
    private final boolean doDD;
    private final double blowFacCD;
    private final double blowFacDD;
    private final CollisionDetectionEngine cd;
    private final CartesianFullBackend back;
    private final OptConfig config;
    
    SnaefellsjoekullGeometryXOver(final boolean adjustRadius, final CUTTYPE cut, final double gaussWidth,
            final boolean doRigidOpt, final boolean checkForCollAndInflate, final boolean doCD,
            final boolean doDD, final double blowFacCD, final double blowFacDD, 
            final CollisionDetectionEngine cd, final CartesianFullBackend back, final OptConfig config,
            final boolean intermediateSanityChecks){
        this.r = new Random();
        this.adjustRadius = adjustRadius;
        this.cut = cut;
        this.gaussWidth = gaussWidth;
        this.doRigidOpt = doRigidOpt;
        this.checkForCollAndInflate = checkForCollAndInflate;
        this.config = config;
        this.doCD = doCD;
        this.doDD = doDD;
        this.blowFacCD = blowFacCD;
        this.blowFacDD = blowFacDD;
        this.cd = cd;
        this.back = back;
        this.intermediateSanityChecks = intermediateSanityChecks;
    }
    
    SnaefellsjoekullGeometryXOver(final SnaefellsjoekullGeometryXOver orig){
        this.r = new Random();
        this.adjustRadius = orig.adjustRadius;
        this.cut = orig.cut;
        this.gaussWidth = orig.gaussWidth;
        this.doRigidOpt = orig.doRigidOpt;
        this.checkForCollAndInflate = orig.checkForCollAndInflate;
        this.config = new OptConfig(orig.config);
        this.doCD = orig.doCD;
        this.doDD = orig.doDD;
        this.blowFacCD = orig.blowFacCD;
        this.blowFacDD = orig.blowFacDD;
        this.cd = orig.cd.clone();
        this.back = orig.back.clone();
        this.intermediateSanityChecks = orig.intermediateSanityChecks;
    }
    
    @Override
    public GenericCrossover<Molecule, Geometry> clone() {
        return new SnaefellsjoekullGeometryXOver(this);
    }

    @Override
    public String getMyID() {
        
        String s = "SNAEFELLSJOEKULL GEOMETRY X OVER";
        s += "\n\t intermediate sanity checks " + intermediateSanityChecks;
        s += "\n\t adjusting radius? " + adjustRadius;
        s += "\n\t which cut? " + cut.name();
        s += "\n\t gaussian width " + gaussWidth;
        s += "\n\t adjust radius " + adjustRadius;
        s += "\n\t check for collision and inflate " + checkForCollAndInflate;
        if(doRigidOpt || checkForCollAndInflate){
            s += "\n\t do rigid opt " + doRigidOpt;
            s += "\n\t backend is " + back.getMethodID();
            s += "\n\t blow factor CD " + blowFacCD;
            s += "\n\t blow factor DD " + blowFacDD;
            s += "\n\t do CD " + doCD;
            s += "\n\t do DD " + doDD;
            s += config.toString();
        }
        
        return s;
    }

    @Override
    public Tuple<Geometry, Geometry> crossover(final Geometry mother, final Geometry father, final long futureID) {
        
        if(mother.getNumberOfIndieParticles() == 1){

            final Geometry gChild1 = new Geometry(mother);
            final Geometry gChild2 = new Geometry(father);

            return new Tuple<>(gChild1,gChild2);
        }
        
        /*
         * create work geometries and rotate them
         */
        final Geometry gChild2 = mother.clone();
        final Geometry gChild1 = father.clone();
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
        CoordTranslation.updateGeometryFromCartesian(cartesFather, gChild1);
        CoordTranslation.updateGeometryFromCartesian(cartesMother, gChild2);

        // set the IDs of mother and father
        gChild1.setMother(mother.getMotherID());
        gChild1.setFather(father.getFatherID());
        gChild2.setMother(mother.getMotherID());
        gChild2.setFather(father.getFatherID());
        gChild1.setID(futureID);
        gChild2.setID(futureID);
        
        final int noIndies = gChild1.getNumberOfIndieParticles();
        assert(noIndies == gChild2.getNumberOfIndieParticles());
        
        // in case of environments, add it in at this point
        if(mother.containsEnvironment()){
            // the environment
            final List<Environment> envChilds = mother.getEnvironment().createOffspring(father.getEnvironment());
            gChild1.setEnvironment(envChilds.get(0));
            gChild2.setEnvironment(envChilds.get(1));
        }
                
        // first get all COMs and put them into matrices/double arrays
        final double[][] motherCOMs = new double[3][noIndies];
        final double[][] fatherCOMs = new double[3][noIndies];
        final double[] motherCOMDists = new double[noIndies];
        final double[] fatherCOMDists = new double[noIndies];

        double maxCOMDistMom = 0.0;
        
        for (int i = 0; i < noIndies; i++) {
            final double[] tmpM = gChild2.getMoleculeAtPosition(i).getExternalCenterOfMass();
            final double[] tmpF = gChild1.getMoleculeAtPosition(i).getExternalCenterOfMass();
            for (int j = 0; j < 3; j++) {
                motherCOMs[j][i] = tmpM[j];
                fatherCOMs[j][i] = tmpF[j];
            }
            
            motherCOMDists[i] = FastMath.sqrt(tmpM[0]*tmpM[0]+tmpM[1]*tmpM[1]+tmpM[2]*tmpM[2]);
            fatherCOMDists[i] = FastMath.sqrt(tmpF[0]*tmpF[0]+tmpF[1]*tmpF[1]+tmpF[2]*tmpF[2]);
            maxCOMDistMom = Math.max(maxCOMDistMom, motherCOMDists[i]);
        }
        
        final List<Integer> motherInside = new ArrayList<>(noIndies);
        final List<Integer> motherOutside = new ArrayList<>(noIndies);
        double cutDistMom = 0.0;
        boolean cont = true;
        int ec = 0;
        do {

            motherInside.clear();
            motherOutside.clear();
            
            switch(cut){
                case FULLRANDOM: cutDistMom = maxCOMDistMom*r.nextDouble(); break;
                case GAUSSDISTRIBUTED: cutDistMom = maxCOMDistMom*RandomUtils.halfgaussDouble(0, 1, gaussWidth); break;
                case INVERTEDGAUSSDISTRIBUTED: cutDistMom = maxCOMDistMom*(1-RandomUtils.halfgaussDouble(0, 1, gaussWidth)); break;
            }
        
            // figure out how many are inside and outside this for mom
            for(int i = 0; i < noIndies; i++){
                if(motherCOMDists[i] < cutDistMom){
                    motherInside.add(i);
                } else {
                    motherOutside.add(i);
                }
            }
            
            // make sure that there are indeed molecules inside and outside the spherical cut
            if(!motherInside.isEmpty() && !motherOutside.isEmpty()){
                cont = false;
            } else {
                ec++;
            }
            
        } while(cont && ec < FixedValues.MAXTOEMERGENCY);
        
        final int noMotherInside = motherInside.size();
        final int noMotherOutside = motherOutside.size();
        
        if(noMotherInside == 0 || noMotherOutside == 0){
            System.err.println("WARNING: Snaefellsjoekull empty sphere!");
            return new Tuple<>(null,null);
        }
        
        // now figure out how we must set cutDistDad so that the same number of molecules is inside and outside as for mom
        final double cutDistDad = findOptimalSphereRadius(fatherCOMDists,noMotherInside);
        if(Double.isNaN(cutDistDad)){
            // two COMs were too close together, return null'd geometries
            System.err.println("WARNING: Couldn't find cutting sphere for father, two COMs too close.");
            return new Tuple<>(null,null);
        }
        
        final List<Integer> fatherInside = new ArrayList<>(noIndies);
        final List<Integer> fatherOutside = new ArrayList<>(noIndies);

        for (int ii = 0; ii < noIndies; ii++) {
            if (fatherCOMDists[ii] <= cutDistDad) {
                fatherInside.add(ii);
            } else {
                fatherOutside.add(ii);
            }
        }

        if(DEBUG){
            System.out.println("DEBUG: Mother with outside ones marked as dummy coming.");
            final CartesianCoordinates cM = gChild2.getCartesians();
            final String[] atomsM = cM.getAllAtomTypes();
            motherOutside.forEach((i) -> {
                atomsM[i] = "XX";
            });
            final String[] cMP = cM.createPrintableCartesians();
            for(final String s : cMP){System.out.println(s);}
            
            System.out.println("DEBUG: Father with outside ones marked as dummy coming.");
            final CartesianCoordinates cF = gChild1.getCartesians();
            final String[] atomsF = cF.getAllAtomTypes();
            fatherOutside.forEach((i) -> {
                atomsF[i] = "XX";
            });
            final String[] cFP = cF.createPrintableCartesians();
            for(final String s : cFP){System.out.println(s);}
        }
        
        /*
         * now we want to see which molecular types exist above and underneath for both geometries
         * first: the father
         */
        final List<Molecule> fatherMols = new ArrayList<>(noIndies);
        for(int i = 0; i < noIndies; i++){
            fatherMols.add(gChild2.getMoleculeAtPosition(i));
        }

        final int[] molTypesFather = assignMolecularTypes(fatherMols);

        // now we want to figure out how many of each are above and underneath

        //find the "highest" molecular type
        int highestMolType = 0;
        for(int i = 0; i < noIndies; i++){
            highestMolType = Math.max(highestMolType, molTypesFather[i]);
        }


        final int[] totTypeCounter = new int[highestMolType+1];
        final int[] fatherOutsideCounter = new int[highestMolType+1];
        final int[] motherOutsideCounter = new int[highestMolType+1];


        for(int i = 0; i < noIndies; i++){
            // total
            totTypeCounter[molTypesFather[i]]++;
        }
        
        for(int i = 0; i < fatherOutside.size(); i++){
            //father
            fatherOutsideCounter[molTypesFather[fatherOutside.get(i)]]++;

            //mother, has p.d. the same overall molecular types as father
            motherOutsideCounter[molTypesFather[motherOutside.get(i)]]++;
        }


        final int[] discrepancy = new int[fatherOutsideCounter.length];
        for(int i = 0; i<discrepancy.length; i++){
            /*
              * this means, that if both are same the result is 0, if mother has more it is negative
              * and if father has more it is positive
              */
            discrepancy[i] = fatherOutsideCounter[i] - motherOutsideCounter[i];
            if(DEBUG){System.out.println("DEBUG: Discrepancy " + i + "\t" + discrepancy[i]);}
        }

        for(int i = 0; i < discrepancy.length; i++){
            if(discrepancy[i] == 0){
                continue;
            } else if(discrepancy[i] > 0){
                // from the mother needs to come outside
                while(discrepancy[i] > 0){
                    // which one of the inside goes outside?
                    final int whichOutside = r.nextInt(totTypeCounter[i] - motherOutsideCounter[i]);

                    // which type to change with?
                    int noOfPoss = 0;
                    for(final int dis : discrepancy){
                        if(dis < 0){
                            noOfPoss++;
                        }
                    }

                    final int whichKindInside = r.nextInt(noOfPoss);

                    // map the random to the type
                    int counter = -1;
                    int type = -1;
                    for(int j = 0; j < discrepancy.length; j++){
                        if(discrepancy[j] < 0) {
                            counter++;
                        }

                        if(counter == whichKindInside) {
                            type = j;
                            break;
                        }
                    }

                    // which molecule of type "type"?
                    final int whichMolOfType = r.nextInt(motherOutsideCounter[type]);

                    // now we need to have the ID for the one coming outside
                    int molIDOut  = -1;
                    int molCounter =-1;
                    for(int iMol = 0; iMol < molTypesFather.length; iMol++){
                        if(molTypesFather[iMol] == i && !motherOutside.contains(iMol)){
                            molCounter++;
                        }

                        if(molCounter == whichOutside){
                            molIDOut = iMol;
                            break;
                        }
                    }

                    // and for the one going down
                    int molIDInner = -1;
                    molCounter = -1;
                    for(int iMol = 0; iMol < molTypesFather.length; iMol++){
                        if(molTypesFather[iMol] == type && motherOutside.contains(iMol)){
                            molCounter++;
                        }

                        if(molCounter == whichMolOfType){
                            molIDInner = iMol;
                            break;
                        }
                    }

                    // now exchange the coordinates and angles of the two
                    final double[] tmpCOM1 = gChild2.getMoleculeAtPosition(molIDOut).getExternalCenterOfMass().clone();
                    final double[] tmpOrient1 = gChild2.getMoleculeAtPosition(molIDOut).getOrientation().clone();

                    final double[] tmpCOM2 = gChild2.getMoleculeAtPosition(molIDInner).getExternalCenterOfMass().clone();
                    final double[] tmpOrient2 = gChild2.getMoleculeAtPosition(molIDInner).getOrientation().clone();

                    gChild2.getMoleculeAtPosition(molIDOut).setExternalCenterOfMass(tmpCOM2);
                    gChild2.getMoleculeAtPosition(molIDOut).setOrientation(tmpOrient2);

                    gChild2.getMoleculeAtPosition(molIDInner).setExternalCenterOfMass(tmpCOM1);
                    gChild2.getMoleculeAtPosition(molIDInner).setOrientation(tmpOrient1);

                    // adjust things a little
                    discrepancy[i]--;
                    discrepancy[type]++;
                    motherOutsideCounter[i]++;
                    motherOutsideCounter[type]--;

                    // and also the arraylists
                    motherOutside.remove(Integer.valueOf(molIDInner));
                    motherOutside.add(molIDOut);
                    motherInside.remove(Integer.valueOf(molIDOut));
                    motherInside.add(molIDInner);
                }
            }

            /*
             * we do not need to explicitly handle the discrepancy < 0 case, the > 0 does that for us :-)
             */
        }

        if(intermediateSanityChecks){
            /*
             * a quick check whether everything is ok
             */
            for(final int i : discrepancy){
                if(i != 0) {
                    System.out.println("ERROR: Discrepancy is still existing: " + i + ".");
                }
            }

            for (int i = 0; i < motherOutsideCounter.length; i++) {
                if (motherOutsideCounter[i] != fatherOutsideCounter[i]) {
                    System.out.println("ERROR: Mismatch father mother: "
                        + motherOutsideCounter[i] + "   " + fatherOutsideCounter[i]);
                }
            }
        }


        /*
         * now do the real exchanging
         */

        /* still we need to exchange parts though
         * remember that we need to adjust with the cutting height (so
         * substracting it from the z component of the mother COMs)
         */
         // let's adjust the mother COMs again, since we swapped them round
         for(int i = 0; i < noIndies; i++){
            final double[] tmp1 = gChild2.getMoleculeAtPosition(i).getExternalCenterOfMass();
            for(int j = 0; j < 3; j++){
                motherCOMs[j][i] = tmp1[j];
            }
            motherCOMDists[i] = FastMath.sqrt(tmp1[0]*tmp1[0]+tmp1[1]*tmp1[1]+tmp1[2]*tmp1[2]);
        }

        // first always null the ArrayLists
        fatherInside.clear();
        fatherOutside.clear();
        motherInside.clear();
        motherOutside.clear();
        for (int kk = 0; kk < noIndies; kk++) {
            if (motherCOMDists[kk] < cutDistMom) {
                motherInside.add(kk);
            } else {
                motherOutside.add(kk);
            }
        }
        for (int kk = 0; kk < noIndies; kk++) {
            if (fatherCOMDists[kk] < cutDistDad) {
                fatherInside.add(kk);
            } else {
                fatherOutside.add(kk);
            }
        }
        
        if(intermediateSanityChecks){
            
            final int iMotherCOMs = motherInside.size();
            if (iMotherCOMs != fatherInside.size()){
                // we f***** it up
                System.err.println("ERROR: We f***** it up in Snaefellsjoekull. Mother: " + iMotherCOMs + " father: " + fatherInside.size());
                return new Tuple<>(null,null);
            }

            // and again we need to figure out the molecular types above
            final int[] iaFatherAboveCounter2 = new int[highestMolType + 1];
            final int[] iaMotherAboveCounter2 = new int[highestMolType + 1];
            for(int i = 0; i < fatherOutside.size(); i++){
                iaFatherAboveCounter2[molTypesFather[fatherOutside.get(i)]]++;
                // mother has p.d. the same overall molecular types as the father
                iaMotherAboveCounter2[molTypesFather[motherOutside.get(i)]]++;
            }

            // check whether the exchange worked
            for(int i = 0; i < iaMotherAboveCounter2.length; i++){
                if (iaMotherAboveCounter2[i] != iaFatherAboveCounter2[i]){
                    System.err.println("ERROR: Exchange didn't work out. At position " + i + " mismatch " + iaMotherAboveCounter2[i]
                        + " to " + iaFatherAboveCounter2[i] + ".");
                    
                    if(DEBUG){
                        System.out.println("Fucked up geometry ");
                    }
                    
                    return new Tuple<>(null,null);
                }
            }

            // check whether we messed something up badly
            final int[] typesChild1 = assignMolecularTypes(gChild1.getMolecules());
            final int[] typesChild2 = assignMolecularTypes(gChild2.getMolecules());
            for(int i = 0; i < typesChild1.length; i++){
                if(typesChild1[i] != molTypesFather[i]
                        || typesChild2[i] != molTypesFather[i]){
                    // darn, didn't work
                    System.err.println("ERROR: Problem in Snaefellsjoekull after internal exchange " +
                        "at position " + i + " father " + molTypesFather[i]
                        + " child 1 " + typesChild1[i] + " child 2 "
                        + typesChild2[i]);
                    return new Tuple<>(null,null);
                }
            }
        }
        
        if(DEBUG){
            final GeometryConfig gFaInner = new GeometryConfig();
            gFaInner.noOfParticles = fatherInside.size();
            gFaInner.geomMCs = new ArrayList<>();
            for(final int under : fatherInside){
                final Molecule m = new Molecule(gChild1.getMoleculeAtPosition(under));
                gFaInner.geomMCs.add(m.returnMyConfig());
            }
            
            final GeometryConfig gFaOuter = new GeometryConfig();
            gFaOuter.noOfParticles = fatherOutside.size();
            gFaOuter.geomMCs = new ArrayList<>();
            for(final int above : fatherOutside){
                final Molecule m = new Molecule(gChild1.getMoleculeAtPosition(above));
                gFaOuter.geomMCs.add(m.returnMyConfig());
            }
            
            final GeometryConfig gMaInner = new GeometryConfig();
            gMaInner.noOfParticles = motherInside.size();
            gMaInner.geomMCs = new ArrayList<>();
            for(final int under : motherInside){
                final Molecule m = new Molecule(gChild2.getMoleculeAtPosition(under));
                gMaInner.geomMCs.add(m.returnMyConfig());
            }
            
            final GeometryConfig gMaOuter = new GeometryConfig();
            gMaOuter.noOfParticles = motherOutside.size();
            gMaOuter.geomMCs = new ArrayList<>();
            for(final int above : motherOutside){
                final Molecule m = new Molecule(gChild2.getMoleculeAtPosition(above));
                gMaOuter.geomMCs.add(m.returnMyConfig());
            }
            
            final Geometry g1 = new Geometry(gFaOuter);
            final Geometry g2 = new Geometry(gFaInner);
            final Geometry g3 = new Geometry(gMaOuter);
            final Geometry g4 = new Geometry(gMaInner);
            
            System.out.println("DEBUG: Father above");
            for(final String s : g1.makePrintableAbsoluteCoord(true)){System.out.println(s);}
            System.out.println("DEBUG: Father under");
            for(final String s : g2.makePrintableAbsoluteCoord(true)){System.out.println(s);}
            System.out.println("DEBUG: Mother above");
            for(final String s : g3.makePrintableAbsoluteCoord(true)){System.out.println(s);}
            System.out.println("DEBUG: Mother under");
            for(final String s : g4.makePrintableAbsoluteCoord(true)){System.out.println(s);}
        }

        // exchange the molecules but keep in mind that the order needs to
        // be conserved!
        final double radiusRatio = cutDistMom/cutDistDad;
        final int[] molTypesMother = molTypesFather.clone();
        final List<Molecule> allMols1 = gChild1.getMolecules();
        final List<Molecule> allMols2 = gChild2.getMolecules();
        for (final int whichMol : fatherOutside) {
            final int whichType = molTypesFather[whichMol];
            for (int j = 0; j < motherOutside.size(); j++) {
                final int tmpMol = motherOutside.get(j);
                final int tmpType = molTypesMother[tmpMol];
                if (tmpType == whichType) {
                    // exchange
                    final Molecule m1 = allMols1.set(whichMol, null); // replace w/ null temporarily
                    final Molecule m2 = allMols2.set(tmpMol, null);

                    // adjust the COM heights
                    final double[] tmpCOM1 = m1.getExternalCenterOfMass();
                    final double[] tmpCOM2 = m2.getExternalCenterOfMass();

                    if(adjustRadius){
                        for(int i = 0; i < 3; i++){
                            // COM1 is mother, COM2 is father
                            tmpCOM1[i] /= radiusRatio;
                            tmpCOM2[i] *= radiusRatio;
                        }
                    }
                    
                    gChild1.setMoleculeAtPosition(whichMol, m2);
                    gChild2.setMoleculeAtPosition(tmpMol, m1);
                    molTypesMother[tmpMol] = -10000;

                    break;
                }
                if (j == motherOutside.size() - 1) {
                    // if we end here, we messed it up
                    System.err.println("ERROR: We messed it up in Snaefellsjoekull.");
                    return new Tuple<>(null,null);
                }
            }
        }

        if(intermediateSanityChecks){
            // again check whether we messed something up badly
            final int[] typesChild12 = assignMolecularTypes(gChild1.getMolecules());
            final int[] typesChild22 = assignMolecularTypes(gChild2.getMolecules());
            for(int i = 0; i < typesChild12.length; i++){
                if(typesChild12[i] != molTypesFather[i]
                        || typesChild22[i] != molTypesFather[i]){
                    // darn, didn't work
                   System.err.println("ERROR: Problem in Snaefellsjoekull after internal exchange " +
                            "at position " + i + " father " + molTypesFather[i]
                            + " child 1 " + typesChild12[i] + " child 2 "
                            + typesChild22[i]);
                    return new Tuple<>(null,null);
                }
            }
        }

        // ensure that the molecular IDs are again correct
        for(int i = 0; i < noIndies; i++){
            gChild1.getMoleculeAtPosition(i).setID(i);
            gChild2.getMoleculeAtPosition(i).setID(i);
        }

        // adding huge fitness values, just in case something goes really wrong
        gChild1.setFitness(Double.MAX_VALUE);
        gChild2.setFitness(Double.MAX_VALUE);

        if(DEBUG){
            System.out.println("DEBUG: First child coming:");
            for(final String s: gChild1.makePrintableAbsoluteCoord(true)){
                System.out.println(s);
            }

            System.out.println("DEBUG: Second child coming:");
            for(final String s: gChild2.makePrintableAbsoluteCoord(true)){
                System.out.println(s);
            }
        }
        
        if(DEBUG){
            final boolean hasCollision1 = cd.checkOnlyForCollision(gChild1.getCartesians(), blowFacCD, gChild1.getBondInfo());
            if(DEBUG){System.out.println("DEBUG: Child 1 has collision? " + hasCollision1);}
            final boolean hasCollision2 = cd.checkOnlyForCollision(gChild2.getCartesians(), blowFacCD, gChild2.getBondInfo());
            if(DEBUG){System.out.println("DEBUG: Child 2 has collision? " + hasCollision2);}
        }
        
        if(intermediateSanityChecks){
            final boolean child1Fine = GlobOptAtomics.areGeometriesStillCorrect(mother,gChild1);
            final boolean child2Fine = GlobOptAtomics.areGeometriesStillCorrect(mother,gChild2);
            if(!child1Fine || !child2Fine){
                System.err.println("ERROR: After crossover, child1 is " + child1Fine + " and child2 is " + child2Fine);
                return new Tuple<>(null,null);
            }
        }
        
        if(doRigidOpt || checkForCollAndInflate){
            // do it for both children
            final Geometry gChildOpt1 = rigidOptimization(doRigidOpt, gChild1, fatherOutside, cd, checkForCollAndInflate,
                    doCD, doDD, blowFacCD, blowFacDD, config.maxInflate, config.incrInflate,
                    back, config.fracMinCut, config.fracMaxCut, config);
            
            final Geometry gChildOpt2 = rigidOptimization(doRigidOpt, gChild2, motherOutside, cd, checkForCollAndInflate,
                    doCD, doDD, blowFacCD, blowFacDD, config.maxInflate, config.incrInflate,
                    back, config.fracMinCut, config.fracMaxCut, config);
            
            return new Tuple<>(gChildOpt1,gChildOpt2);
        }
        
        return new Tuple<>(gChild1,gChild2);
    }

    @Override
    public short hasPriority() {
        return -1;
    }
    
    private static Geometry rigidOptimization(final boolean doOpt, final Geometry g, final List<Integer> outer,
            final CollisionDetectionEngine cd, final boolean checkForCollAndInflate,
            final boolean doCD, final boolean doDD, final double blowFacCD, final double blowFacDD,
            final double maxInflate, final double incrInflate, final CartesianFullBackend back,
            final double fracMinCut, final double fracMaxCut, final OptConfig config){
        
        final double[][] outerCOMs = new double[3][outer.size()];
        
        int c = 0;
        for(final int out : outer){
            final double[] com = g.getCOM(out);
            outerCOMs[0][c] = com[0];
            outerCOMs[1][c] = com[1];
            outerCOMs[2][c] = com[2];
            c++;
        }
        
        // check for collisions
        if(checkForCollAndInflate){
            int ec = 0;
            double inflate = 1.0;
            do {
                final boolean hasCollision = cd.checkOnlyForCollision(g.getCartesians(), blowFacCD, g.getBondInfo());
                if(DEBUG){System.out.println("DEBUG: At inflation " + inflate + " collision? " + hasCollision);}
                if(!hasCollision){
                    // cool, no collision
                    break;
                }
            
                // inflate the outer sphere part a bit
                inflate += incrInflate;
                int cx = 0;
                for(final int out : outer){
                    final double[] com = g.getCOM(out);
                    com[0] = inflate*outerCOMs[0][cx];
                    com[1] = inflate*outerCOMs[1][cx];
                    com[2] = inflate*outerCOMs[2][cx];
                    cx++;
                }
            
            } while(inflate <= maxInflate && ec < FixedValues.MAXTOEMERGENCY);
        
            if(ec == FixedValues.MAXTOEMERGENCY){
                System.err.println("WARNING: No collision-free structure found after " + ec + "tries.");
                return g;
            }
            
            // get the new COMs after inflation
            c = 0;
            for(final int out : outer){
                final double[] com = g.getCOM(out);
                outerCOMs[0][c] = com[0];
                outerCOMs[1][c] = com[1];
                outerCOMs[2][c] = com[2];
                c++;
            }
        }
        
        if(!doOpt){
            return g;
        }
        
        // now do the actual local optimization
        // we do a 6D rigid body optimization: 3 eulers for the rotation of the outer shell
        // and 3 expansion cartesians
        final double[][] bounds = new double[2][6];
        bounds[0][0] = -Math.PI;
        bounds[1][0] = Math.PI;
        bounds[0][1] = -0.5*Math.PI;
        bounds[1][1] = 0.5*Math.PI;
        bounds[0][2] = -Math.PI;
        bounds[1][2] = Math.PI;
        bounds[0][3] = fracMinCut;
        bounds[1][3] = fracMaxCut;
        bounds[0][4] = fracMinCut;
        bounds[1][4] = fracMaxCut;
        bounds[0][5] = fracMinCut;
        bounds[1][5] = fracMaxCut;
        
        // "normalize" on the fly for the initial guess
        final double[] guess = new double[]{0.5,0.5,0.5,
            (1.0-fracMinCut)/fracMaxCut,(1.0-fracMinCut)/fracMaxCut,(1.0-fracMinCut)/fracMaxCut
        };
        
        final BOBYQAOptimizer opt = new BOBYQAOptimizer(bounds,
                config.numberOfInterpolationPoints,config.initialTrustRegionRadius,
                config.stoppingTrustRegionRadius,config.noIterations,true);
        
        final OptimizeFunction func = new OptimizeFunction(cd,
                back, g, blowFacCD, blowFacDD, g.getBondInfo(), outer, bounds,
                doCD, doDD, outerCOMs);
                
        // normalize
        final Tuple<Double,double[]> res = opt.doOptimize(guess, func); // guess is already normalized
        if(DEBUG){
            System.out.println("DEBUG: Snaefellsjoekull shell opter: best fitness: " + res.getObject1());
        }
        
        return func.g;
    }
    
    private static class OptimizeFunction extends AbstractBOBYQAMethod {
        
        private final CollisionDetectionEngine cd;
        private final CartesianFullBackend back;
        private Geometry g;
        private final CartesianCoordinates orig;
        private final CartesianCoordinates work;
        private final double[] xyz1D;
        private final double blowBonds;
        private final double blowDissoc;
        private final boolean doCD;
        private final boolean doDD;
        private final BondInfo bonds;
        private final List<Integer> outerIDs;
        private final int noOuters;
        private final double[] eparts;
        private final double[] point;
        private final double[][] origCOMs;
        private final double[][] bounds;
        private final double[][] comOrig;
        private final double[][] transResult;
        private final double[] euler = new double[3];
        private int iter = 0;
        
        OptimizeFunction(final CollisionDetectionEngine cd,
                final CartesianFullBackend back, final Geometry g, final double blowBonds, final double blowDissoc,
                final BondInfo bonds, final List<Integer> outerIDs, final double[][] bounds,
                final boolean doCD, final boolean doDD, final double[][] comOrig){
            this.cd = cd;
            this.back = back;
            this.g = g;
            this.blowBonds = blowBonds;
            this.blowDissoc = blowDissoc;
            this.bonds = bonds;
            this.outerIDs = outerIDs;
            this.noOuters = outerIDs.size();
            this.eparts = new double[g.getNumberOfIndieParticles()];
            this.point = new double[6];
            this.bounds = bounds;
            this.transResult = new double[3][noOuters];
            this.doCD = doCD;
            this.doDD = doDD;
            this.comOrig = comOrig;
            this.orig = g.getCartesians();
            this.work = orig.clone();
            this.xyz1D = new double[work.getNoOfAtoms()*3];
            origCOMs = g.getAllCOMs();
        }
        
        @Override
        public double computeObjectiveValue(final double[] normalized){
            
            // denormalize
            denormalizeToBounds(bounds, normalized, point);
            
            // first expand/contract, then rotate, then set
            for(int i = 0; i < noOuters; i++){
                transResult[0][i] = point[3]*comOrig[0][i];
                transResult[1][i] = point[4]*comOrig[1][i];
                transResult[2][i] = point[5]*comOrig[2][i];
            }
            euler[0] = point[0];
            euler[1] = point[1];
            euler[2] = point[2];
            
            // update the geometry 
            final double[][] rotResult = CoordTranslation.rotateXYZ(transResult, euler);
            for(int i = 0; i < noOuters; i++){
                final int whichMol = outerIDs.get(i);
                final double[] com = g.getCOM(whichMol);
                com[0] = rotResult[0][i];
                com[1] = rotResult[1][i];
                com[2] = rotResult[2][i];
            }
            
            // NOTE: the rotation is NOT entirely "correct". a "correct" rotation should
            // have an influence on the euler angles of the molecules (or the orientation of
            // the individual molecules, same thing). Hence, this is a somewhat of a skewed rotation.

            iter++;
            
            // put the coordinates into the cartesian coordinates
            final double[][] xyzOrig = orig.getAllXYZCoord();
            final double[][] xyz = work.getAllXYZCoord();
            final int[] atsMol = orig.getAllAtomsPerMol();
            final int noMols = g.getNumberOfIndieParticles();
            int off = 0;
            for(int mol = 0; mol < noMols; mol++){
                final int molAts = atsMol[mol];
                final double[] newCOM = g.getMoleculeAtPosition(mol).getExternalCenterOfMass();
                for(int coord = 0; coord < 3; coord++){
                    for(int at = 0; at < molAts; at++){
                        xyz[coord][off+at] = xyzOrig[coord][off+at] + newCOM[coord] - origCOMs[coord][mol];
                    }
                }
                off += molAts;
            }
            
            if(doCD){
                // cd
                final CollisionInfo inf = cd.checkForCollision(work, blowBonds, bonds);
                if(inf.hasCollision()) {return FixedValues.NONCONVERGEDENERGY;}
                
                if(doDD){
                    // dd
                    final boolean diss = DissociationDetection.checkForDissociation(inf.getPairWiseDistances(),
                        work.getAllAtomTypes(), work.getAllAtomNumbers(), blowDissoc, DissociationDetection.DEFAULTDD);
                    // we assume that dissociation is "better" than collision as it proves at least that a collision free state exists
                    if(diss) {return FixedValues.NONCONVERGEDENERGY/2;}
                }
            }
            
            work.getAll1DCartes(xyz1D);            
            final double obj = (back == null) ? 0.0 : back.energyCalculation(g.getID(), iter, xyz1D,
                    work.getAllAtomTypes(), work.getAllAtomNumbers(), work.getAllAtomsPerMol(),
                    eparts, work.getNoOfAtoms(), work.getAllCharges(), work.getAllSpins(), g.getBondInfo());
            
            if(DEBUG){
                System.out.println("DEBUG: At iter " + iter + " COM translation is " + point[3] + " " + point[4] + " " + point[5]);
                System.out.println("DEBUG: At iter " + iter + " Euler rotation is " + point[0] + " " + point[1] + " " + point[2]);
                System.out.println("DEBUG: At iter " + iter + " rigid opt reports " + obj);
            }
            
            return obj;
        }
        
        @Override
        public boolean doesNormalize(){
            return true;
        }
    }
    
    public static class OptConfig implements Serializable {
        private static final long serialVersionUID = (long) 20130327;
        public int numberOfInterpolationPoints = 13; // recommended setting is 2N+1
        public double initialTrustRegionRadius = 1E-1;
        public double stoppingTrustRegionRadius = 1E-5;
        public int noIterations = 250;
        public double maxInflate = 1.3;
        public double incrInflate = 0.05;
        public double fracMinCut = 0.8;
        public double fracMaxCut = 1.5;
        
        public OptConfig(){}
        
        public OptConfig(final OptConfig orig){
            this.initialTrustRegionRadius = orig.initialTrustRegionRadius;
            this.noIterations = orig.noIterations;
            this.numberOfInterpolationPoints = orig.numberOfInterpolationPoints;
            this.stoppingTrustRegionRadius = orig.stoppingTrustRegionRadius;
            this.maxInflate = orig.maxInflate;
            this.incrInflate = orig.incrInflate;
            this.fracMaxCut = orig.fracMaxCut;
            this.fracMinCut = orig.fracMinCut;
        }
        
        @Override
        public String toString(){
            return "\n\t interpoints: " + this.numberOfInterpolationPoints
                    + "\n\t initial trust: " + this.initialTrustRegionRadius
                    + "\n\t stopping trust " + this.stoppingTrustRegionRadius
                    + "\n\t max iterations: " + this.noIterations
                    + "\n\t maximum inflate " + this.maxInflate
                    + "\n\t inflation " + this.incrInflate
                    + "\n\t min border radius opt " + this.fracMinCut
                    + "\n\t max border radius opt " + this.fracMaxCut;
        }
    }
}
