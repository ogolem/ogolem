/**
Copyright (c) 2013-2014, J. M. Dieterich
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
import static org.ogolem.core.Constants.ANGTOBOHR;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * A more intelligent phenotype merging the cluster parts together using another
 * local optimization engine in the translation and rotation dimension. Requires
 * a backend. Original implementation in some version of phenix, this is an independent
 * rewrite though, never seen the original code (so all bugs are mine, all ideas are Bernds).
 * @author Johannes Dieterich
 * @version 2020-07-03
 */
public class MergingPhenoXOver implements GenericCrossover<Molecule,Geometry> {
   
    private static final long serialVersionUID = (long) 20200209;
    private static final boolean DEBUG = false;
    
    private final boolean doCDDD;
    private final boolean randomTrans;    
    private final CartesianFullBackend back;
    private final CollisionDetectionEngine colldetect;
    private final double blowBonds;
    private final double blowDissoc;
    private final MergingPhenoXOver.MergingPhenoConfig config;                
    
    MergingPhenoXOver(final CartesianFullBackend back, final double blowColl, final double blowDiss,
            final CollisionDetectionEngine colldetect, final MergingPhenoXOver.MergingPhenoConfig config,
            final boolean doCDDD, final boolean randomTrans){
        this.back = back;
        this.blowBonds = blowColl;
        this.blowDissoc = blowDiss;
        this.colldetect = colldetect;
        this.config = config;
        this.doCDDD = doCDDD;
        this.randomTrans = randomTrans;
    }
    
    MergingPhenoXOver(final MergingPhenoXOver orig){
        this.back = orig.back.clone();
        this.blowBonds = orig.blowBonds;
        this.blowDissoc = orig.blowDissoc;
        this.config = new MergingPhenoXOver.MergingPhenoConfig(orig.config);
        this.colldetect = orig.colldetect.clone();
        this.doCDDD = orig.doCDDD;
        this.randomTrans = orig.randomTrans;
    }
    
    @Override
    public MergingPhenoXOver clone() {
        return new MergingPhenoXOver(this);
    }

    @Override
    public String getMyID() {
        return "MERGING GLOBOPT\n\tbackend: " + back.getMethodID() +
                "\n\tblow bonds: " + blowBonds + "\n\tblow dissoc: " + blowDissoc
                + "\n\t config: " + config.toString();
    }

    @Override
    public Tuple<Geometry, Geometry> crossover(final Geometry mother, final Geometry father, final long futureID) {
        
        final Geometry gChild1 = new Geometry(father);
        gChild1.setID(futureID);
        final Geometry gChild2 = new Geometry(mother);
        gChild2.setID(futureID);
        
        final Lottery random = Lottery.getInstance();

         // first get all COMs and put them into matrices/double arrays
        final double[][] motherCOMs = new double[3][gChild1.getNumberOfIndieParticles()];
        final double[][] fatherCOMs = new double[3][gChild1.getNumberOfIndieParticles()];

        for(int i = 0; i < gChild1.getNumberOfIndieParticles(); i++){
            final double[] tmp1 = gChild1.getMoleculeAtPosition(i).getExternalCenterOfMass();
            final double[] tmp2 = gChild2.getMoleculeAtPosition(i).getExternalCenterOfMass();
            for(int j = 0; j < 3; j++){
                fatherCOMs[j][i] = tmp1[j];
                motherCOMs[j][i] = tmp2[j];
            }
        }
        
        /*
          * now we want to take a plane cutting through both
          * 1) which height on the z axis?
          */
        final int noMols = gChild1.getNumberOfIndieParticles();
        final List<Integer> fatherAbove = new ArrayList<>(noMols);
        final List<Integer> fatherUnder = new ArrayList<>(noMols);
        double planeHeightFather;
        int ec = 0;
        boolean cont;
        int noFatherUnder;
        do{
            planeHeightFather = GlobOptAtomics.randomCuttingZPlane(config.planeMode, fatherCOMs, random);
            // check which COMs are above and underneath the plane (for the father)
            for (int i = 0; i < noMols; i++) {
                if (fatherCOMs[2][i] <= planeHeightFather) {
                    // underneath
                    fatherUnder.add(i);
                } else {
                    // above
                    fatherAbove.add(i);
                }
            }

            noFatherUnder = fatherUnder.size();
            if (noFatherUnder == 0) {
                // then there was none underneath the plane, this is NOT ok.
                cont = true;
                fatherUnder.clear();
                fatherAbove.clear();
            } else if (noFatherUnder >= noMols) {
                // all are above the plane, also not OK.
                cont = true;
                fatherUnder.clear();
                fatherAbove.clear();
            } else{
                cont = false;
            }
            
            ec++;
        } while(cont && ec < FixedValues.MAXTOEMERGENCY);

        if(cont){
            System.err.println("Too many attemps in finding a first plane height!");
            return new Tuple<>(null,null);
        }
        
        // now move the plane in the mother geometry till enough COMs are above/underneath
        final double planeHeightMum = GlobOptAtomics.findOptimalZPlaneHeight(motherCOMs, noFatherUnder);

        if(Double.isNaN(planeHeightMum)){
            // two COMs were too close together, return null'd geometries
            System.err.println("WARNING: Couldn't find cutting plane for mother, two COMs too close.");
            return new Tuple<>(null,null);
        }
        
        final List<Integer> motherUnder = new ArrayList<>(noMols);
        final List<Integer> motherAbove = new ArrayList<>(noMols);

        for (int ii = 0; ii < noMols; ii++) {
            if (motherCOMs[2][ii] <=planeHeightMum) {
                // underneath
                motherUnder.add(ii);
            } else {
                // above
                motherAbove.add(ii);
            }
        }

        if(noFatherUnder != motherUnder.size() || motherAbove.size() != fatherAbove.size()
                || fatherUnder.size() != motherUnder.size()){
            System.err.println("DEBUG: Plane height is bad. Returning. Mismatch: " + noFatherUnder + " vs " + motherUnder.size()
                + " or " + motherAbove.size() + " vs " + fatherAbove.size());
            return new Tuple<>(null,null);
        }

        /*
         * now we want to see which molecular types exist above and underneath for both geometries
         * first: the father
         */
        final List<Molecule> fatherMols = new ArrayList<>(noMols);

        for(int i = 0; i < noMols; i++){
            fatherMols.add(gChild2.getMoleculeAtPosition(i));
        }

        final int[] molTypesDad = GlobOptAtomics.assignMolecularTypes(fatherMols);

        // now we want to figure out how many of each are above and underneath

        //find the "highest" molecular type
        int highestType = 0;
        for(int i = 0; i < noMols; i++){
            highestType = Math.max(highestType, molTypesDad[i]);
        }


        final int[] totTypeCounter = new int[highestType+1];
        final int[] dadAboveCounter = new int[highestType+1];
        final int[] mumAboveCounter = new int[highestType+1];


        for(int i = 0; i < noMols; i++){
            // total
            totTypeCounter[molTypesDad[i]]++;

        }

        for(int i = 0; i < fatherAbove.size(); i++){
            //father
            dadAboveCounter[molTypesDad[fatherAbove.get(i)]]++;

            //mother, has p.d. the same overall molecular types as father
            mumAboveCounter[molTypesDad[motherAbove.get(i)]]++;
        }


        final int[] discrepancy = new int[dadAboveCounter.length];
        for(int i = 0; i < discrepancy.length; i++){
            /*
              * this means, that if both are same the result is 0, if mother has more it is negative
              * and if father has more it is positive
              */
            discrepancy[i] = dadAboveCounter[i] - mumAboveCounter[i];
        }

        for(int i = 0; i < discrepancy.length; i++){
            if(discrepancy[i] == 0){
                continue;
            } else if(discrepancy[i] > 0){
                // from the mother needs to come up
                while(discrepancy[i] > 0){
                    // which one of the under goes up?
                    final int iWhichUp = random.nextInt( totTypeCounter[i] - mumAboveCounter[i]);

                    // which type to change with?
                    int noPoss = 0;
                    for(int j = 0; j < discrepancy.length; j++){
                        if(discrepancy[j] < 0){
                            noPoss++;
                        }
                    }

                    final int whichKindDown = random.nextInt(noPoss);

                    // map the random to the type
                    int counter = -1;
                    int type = -1;
                    for(int j = 0; j < discrepancy.length; j++){
                        if(discrepancy[j] < 0) {
                            counter++;
                        }

                        if(counter == whichKindDown) {
                            type = j;
                            break;
                        }
                    }

                    // which molecule of type iType?
                    final int whichMolType = random.nextInt(mumAboveCounter[type]);

                    // now we need to have the ID for the one coming up
                    int molIDUp  = -1;
                    int molCounter =-1;
                    for(int mol = 0; mol < molTypesDad.length; mol++){
                        if(molTypesDad[mol] == i && !motherAbove.contains(mol)){
                            molCounter++;
                        }

                        if(molCounter == iWhichUp){
                            molIDUp = mol;
                            break;
                        }
                    }

                    // and for the one going down
                    int molIDDown = -1;
                    molCounter = -1;
                    for(int mol = 0; mol < molTypesDad.length; mol++){
                        if(molTypesDad[mol] == type && motherAbove.contains(mol)){
                            molCounter++;
                        }

                        if(molCounter == whichMolType){
                            molIDDown = mol;
                            break;
                        }
                    }

                    // now exchange the coordinates and angles of the two
                    final double[] tmpCOM1 = gChild2.getMoleculeAtPosition(molIDUp).getExternalCenterOfMass().clone();
                    final double[] tmpOrient1 = gChild2.getMoleculeAtPosition(molIDUp).getOrientation().clone();

                    final double[] tmpCOM2 = gChild2.getMoleculeAtPosition(molIDDown).getExternalCenterOfMass().clone();
                    final double[] daTempOrient2 = gChild2.getMoleculeAtPosition(molIDDown).getOrientation().clone();

                    gChild2.getMoleculeAtPosition(molIDUp).setExternalCenterOfMass(tmpCOM2);
                    gChild2.getMoleculeAtPosition(molIDUp).setOrientation(daTempOrient2);

                    gChild2.getMoleculeAtPosition(molIDDown).setExternalCenterOfMass(tmpCOM1);
                    gChild2.getMoleculeAtPosition(molIDDown).setOrientation(tmpOrient1);

                    // adjust things a little
                    discrepancy[i]--;
                    discrepancy[type]++;
                    mumAboveCounter[i]++;
                    mumAboveCounter[type]--;

                    // and also the arraylists
                    motherAbove.remove(Integer.valueOf(molIDDown));
                    motherAbove.add(molIDUp);
                    motherUnder.remove(Integer.valueOf(molIDUp));
                    motherUnder.add(molIDDown);
                }
            }

            /*
             * we do not need to explicitly handle the discrepancy < 0 case, the >0 does that for us :-)
             */
        }

        /*
         * a quick check whether everything is ok
         */
         for(int j = 0; j < discrepancy.length; j++){
             if(discrepancy[j] !=0) {
                 System.err.println("ERROR: Discrepancy is still existing: " + discrepancy[j] + ".");
                 return new Tuple<>(null,null);
             }
         }

        for (int i = 0; i < mumAboveCounter.length; i++) {
            if (mumAboveCounter[i] != dadAboveCounter[i]) {
                System.err.println("ERROR: Mismatch father mother: "
                        + mumAboveCounter[i] + "   " + dadAboveCounter[i]);
                return new Tuple<>(null,null);
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
         for(int i = 0; i < noMols; i++){
            final double[] daTemp1 = gChild2.getMoleculeAtPosition(i).getExternalCenterOfMass().clone();
            for(int j = 0; j < 3; j++){
                motherCOMs[j][i] = daTemp1[j];
            }
        }
         
        // first always null the ArrayLists
        motherUnder.clear();
        motherAbove.clear();
        for (int kk = 0; kk < noMols; kk++) {
            if (motherCOMs[2][kk] <= planeHeightMum) {
                // underneath
                motherUnder.add(kk);
            } else {
                // above
                motherAbove.add(kk);
            }
        }
        final int noMotherUnder = motherUnder.size();

        if (noMotherUnder != noFatherUnder){
            // we f***** it up
            System.err.println("DEBUG: We f***** it up in merge globopt. Mother: " + noMotherUnder + " father: " + noFatherUnder);
            return new Tuple<>(null,null);
        }

        // and again we need to figure out the molecular types above
        final int[] mumAboveCounter2 = new int[highestType + 1];
        motherAbove.forEach((motherAbove1) -> {
            //mother has p.d. the same overall molecular types as the father
            mumAboveCounter2[molTypesDad[motherAbove1]]++;
        });

        // check whether the exchange worked
        for(int i = 0; i < mumAboveCounter2.length; i++){
            if (mumAboveCounter2[i] != dadAboveCounter[i]){
                System.err.println("DEBUG: Exchange didn't work out. At position " + i + " mismatch " + mumAboveCounter2[i]
                        + " to " + dadAboveCounter[i] + ".");
                return new Tuple<>(null,null);
            }
        }

        // check whether we messed something up badly
        final int[] typesChild1 = GlobOptAtomics.assignMolecularTypes(gChild1.getMolecules());
        final int[] typesChild2 = GlobOptAtomics.assignMolecularTypes(gChild2.getMolecules());
        for(int i = 0; i < typesChild1.length; i++){
            if(typesChild1[i] != molTypesDad[i]
                    || typesChild2[i] != molTypesDad[i]){
                // darn, didn't work
                System.err.println("DEBUG: Problem in merge globopt after internal exchange " +
                        "at position " + i + " father " + molTypesDad[i]
                        + " child 1 " + typesChild1[i] + " child 2 "
                        + typesChild2[i]);
                return new Tuple<>(null,null);
            }
        }
        
        if(DEBUG){
            final GeometryConfig gFaUnder = new GeometryConfig();
            gFaUnder.noOfParticles = fatherUnder.size();
            gFaUnder.geomMCs = new ArrayList<>();
            for(final int under : fatherUnder){
                final Molecule m = new Molecule(gChild1.getMoleculeAtPosition(under));
                gFaUnder.geomMCs.add(m.returnMyConfig());
            }
            
            final GeometryConfig gFaAbove = new GeometryConfig();
            gFaAbove.noOfParticles = fatherAbove.size();
            gFaAbove.geomMCs = new ArrayList<>();
            for(final int above : fatherAbove){
                final Molecule m = new Molecule(gChild1.getMoleculeAtPosition(above));
                gFaAbove.geomMCs.add(m.returnMyConfig());
            }
            
            final GeometryConfig gMaUnder = new GeometryConfig();
            gMaUnder.noOfParticles = motherUnder.size();
            gMaUnder.geomMCs = new ArrayList<>();
            for(final int under : motherUnder){
                final Molecule m = new Molecule(gChild2.getMoleculeAtPosition(under));
                gMaUnder.geomMCs.add(m.returnMyConfig());
            }
            
            final GeometryConfig gMaAbove = new GeometryConfig();
            gMaAbove.noOfParticles = motherAbove.size();
            gMaAbove.geomMCs = new ArrayList<>();
            for(final int above : motherAbove){
                final Molecule m = new Molecule(gChild2.getMoleculeAtPosition(above));
                gMaAbove.geomMCs.add(m.returnMyConfig());
            }
            
            final Geometry g1 = new Geometry(gFaAbove);
            final Geometry g2 = new Geometry(gFaUnder);
            final Geometry g3 = new Geometry(gMaAbove);
            final Geometry g4 = new Geometry(gMaUnder);
            
            System.out.println("DEBUG: Father above");
            for(final String s : g1.makePrintableAbsoluteCoord(true)){System.out.println(s);}
            System.out.println("DEBUG: Father under");
            for(final String s : g2.makePrintableAbsoluteCoord(true)){System.out.println(s);}
            System.out.println("DEBUG: Mother above");
            for(final String s : g3.makePrintableAbsoluteCoord(true)){System.out.println(s);}
            System.out.println("DEBUG: Mother under");
            for(final String s : g4.makePrintableAbsoluteCoord(true)){System.out.println(s);}
        }
        
        /*
         * NOW DO THE REAL EXCHANGE OF MOLECULES
         */
        // exchange the molecules but keep in mind that the order needs to
        // be conserved!
        
        // stores the ids of molecules above the original cutting plane after the Xchange
        final int[] molTypesMum = molTypesDad.clone();
        final List<Molecule> allMols1 = gChild1.getMolecules();
        final List<Molecule> allMols2 = gChild2.getMolecules();
        for (final int fatherAbove1 : fatherAbove) {
            final int whichMol = fatherAbove1;
            final int whichType = molTypesDad[whichMol];
            for (int j = 0; j < motherAbove.size(); j++) {
                final int tempMol = motherAbove.get(j);
                final int tempType = molTypesMum[tempMol];
                if (tempType == whichType) {
                    // exchange
                    final Molecule m1 = allMols1.set(whichMol, null); // replace w/ null temporarily
                    final Molecule m2 = allMols2.set(tempMol, null);

                    //adjust the COM heights
                    final double[] tmpCOM1 = m1.getExternalCenterOfMass();
                    final double[] tmpCOM2 = m2.getExternalCenterOfMass();

                    tmpCOM1[2] = tmpCOM1[2] - planeHeightMum + planeHeightFather;
                    tmpCOM2[2] = tmpCOM2[2] + planeHeightMum - planeHeightFather;

                    gChild1.setMoleculeAtPosition(whichMol, m2);
                    gChild2.setMoleculeAtPosition(tempMol, m1);
                    molTypesMum[tempMol] = -10000;

                    break;
                }
                if (j == motherAbove.size() - 1) {
                    // if we end here, we messed it up
                    System.err.println("DEBUG: We messed it up in MergingPheno.");
                    return new Tuple<>(null,null);
                }
            }
        }

        // again check whether we messed something up badly
        final int[] typesChild12 = GlobOptAtomics.assignMolecularTypes(gChild1.getMolecules());
        final int[] typesChild22 = GlobOptAtomics.assignMolecularTypes(gChild2.getMolecules());
        for(int i = 0; i < typesChild12.length; i++){
            if(typesChild12[i] != molTypesDad[i]
                    || typesChild22[i] != molTypesDad[i]){
                // darn, didn't work
                System.err.println("DEBUG: Problem in MergingPheno after internal exchange " +
                        "at position " + i + " father " + molTypesDad[i]
                        + " child 1 " + typesChild12[i] + " child 2 "
                        + typesChild22[i]);
                return new Tuple<>(null,null);
            }
        }

        // ensure that the molecular IDs are again correct
        for(int i = 0; i < gChild1.getNumberOfIndieParticles(); i++){
            gChild1.getMoleculeAtPosition(i).setID(i);
            gChild2.getMoleculeAtPosition(i).setID(i);
        }

        // adding huge fitness values, just in case something goes really wrong
        gChild1.setFitness(Double.MAX_VALUE);
        gChild2.setFitness(Double.MAX_VALUE);

        if(DEBUG){
            final String[] saCoords1 = gChild1.makePrintableAbsoluteCoord(true);
            final String[] saCoords2 = gChild2.makePrintableAbsoluteCoord(true);

            System.out.println("DEBUG: First child coming:");
            for(final String s: saCoords1){
                System.out.println(s);
            }

            System.out.println("DEBUG: Second child coming:");
            for(final String s: saCoords2){
                System.out.println(s);
            }
        }
        
        // XXX: in an effort to unify, we should be stopping here for sweden and
        // remove the "old" pheno
        
        /*
         * AT THIS POINT OF THE ALGORITHM, WE HAVE TWO PLANE HEIGHTS (FOR MOM AND DAD)
         * AND WE HAVE EQUILLIBRATED THE MOLECULAR TYPES EXISTING ABOVE AND UNDERNEATH
         * SAID PLANE FOR BOTH CLUSTERS. ALSO, THINGS HAVE BEEN EXCHANGED.
         * 
         * NOW, WE WILL "SEPARATE" THE CLUSTER HALFES AND DO OUR MERGING LOCOPT
         * 
         * ACTUALLY, WE WILL BE JUST USING TWO COPIES AND THE SET OF FLAGS TO DECIDE
         * WHICH MOLECULAR TYPE TO USE AS THIS WILL SIGNIFICANTLY SIMPLIFY ALL THE
         * HANDLING CODE.
         */
        
        // prepare the optimization for both candidate structures
        final BOBYQAOptimizer opt = new BOBYQAOptimizer(config.optBounds,
                config.numberOfInterpolationPoints,config.initialTrustRegionRadius,
                config.stoppingTrustRegionRadius,config.noIterations,true);
        
        if(DEBUG){
            System.out.println("DEBUG: bounds[0][0] " + config.optBounds[0][0]);
            System.out.println("DEBUG: bounds[0][1] " + config.optBounds[0][1]);
            System.out.println("DEBUG: bounds[1][0] " + config.optBounds[1][0]);
            System.out.println("DEBUG: bounds[1][1] " + config.optBounds[1][1]);
        }
        
        double[] guess1;
        double[] guess2;
        double[] normG1;
        double[] normG2;
        if(config.use6D){
            guess1 = new double[6];
            normG1 = new double[6];
            if(randomTrans){
                guess1[0] = RandomUtils.halfgaussDouble(config.optBounds[0][0],config.optBounds[1][0],0.1);
                guess1[1] = RandomUtils.halfgaussDouble(config.optBounds[0][1],config.optBounds[1][1],0.1);
                guess1[2] = RandomUtils.halfgaussDouble(config.optBounds[0][2],config.optBounds[1][2],0.1);
            } // else: 0.0 in already
            guess1[3] = config.optBounds[0][3] + random.nextDouble()*(config.optBounds[1][3]-config.optBounds[0][3]);
            guess1[4] = config.optBounds[0][4] + random.nextDouble()*(config.optBounds[1][4]-config.optBounds[0][4]);
            guess1[5] = config.optBounds[0][5] + random.nextDouble()*(config.optBounds[1][5]-config.optBounds[0][5]);
            if(DEBUG){
                System.out.println("DEBUG: Structure 1: initial guess " + guess1[0] + "\t" + guess1[1] + "\t" + guess1[2] + "\t" + guess1[3] + "\t" + guess1[4] + "\t" + guess1[5]);
            }
        
            guess2 = new double[6];
            normG2 = new double[6];
            if(randomTrans){
                guess2[0] = RandomUtils.halfgaussDouble(config.optBounds[0][0],config.optBounds[1][0],0.1);
                guess2[1] = RandomUtils.halfgaussDouble(config.optBounds[0][1],config.optBounds[1][1],0.1);
                guess2[2] = RandomUtils.halfgaussDouble(config.optBounds[0][2],config.optBounds[1][2],0.1);
            } // else: 0.0 in already
            guess2[3] = config.optBounds[0][3] + random.nextDouble()*(config.optBounds[1][3]-config.optBounds[0][3]);
            guess2[4] = config.optBounds[0][4] + random.nextDouble()*(config.optBounds[1][4]-config.optBounds[0][4]);
            guess2[5] = config.optBounds[0][5] + random.nextDouble()*(config.optBounds[1][5]-config.optBounds[0][5]);
            if(DEBUG){
                System.out.println("DEBUG: Structure 2: initial guess " + guess2[0] + "\t" + guess2[1] + "\t" + guess2[2] + "\t" + guess2[3] + "\t" + guess2[4] + "\t" + guess2[5]);
            }
        } else {
            guess1 = new double[2];
            normG1 = new double[2];
            guess1[0] = RandomUtils.halfgaussDouble(config.optBounds[0][0],config.optBounds[1][0]);
            guess1[1] = config.optBounds[0][1] + random.nextDouble()*(config.optBounds[1][1]-config.optBounds[0][1]);
            if(DEBUG){
                System.out.println("DEBUG: Structure 1: initial guess " + guess1[0] + "\t" + guess1[1]);
            }

            guess2 = new double[2];
            normG2 = new double[2];
            guess2[0] = RandomUtils.halfgaussDouble(config.optBounds[0][0],config.optBounds[1][0]);
            guess2[1] = config.optBounds[0][1] + random.nextDouble()*(config.optBounds[1][1]-config.optBounds[0][1]);
            if(DEBUG){
                System.out.println("DEBUG: Structure 2: initial guess " + guess2[0] + "\t" + guess2[1]);
            }
        }
        
        final MergingPhenoXOver.OptimizeFunction func1 = new MergingPhenoXOver.OptimizeFunction(config.use6D, colldetect,
                back, gChild1, gChild1.getCartesians(), blowBonds, blowDissoc, fatherAbove, config.optBounds, doCDDD);
        
        // normalize
        func1.normalizeFromBounds(config.optBounds,guess1,normG1);
        
        final Tuple<Double,double[]> res1 = opt.doOptimize(normG1, func1);
        if(DEBUG){
            System.out.println("DEBUG: Structure 1: best fitness: " + res1.getObject1());
            if(config.use6D){
                System.out.println("DEBUG: Structure 1: coordinates (normalized) " + res1.getObject2()[0] + "\t" + res1.getObject2()[1] + "\t" + res1.getObject2()[2] + "\t" + res1.getObject2()[3] + "\t" + res1.getObject2()[4] + "\t" + res1.getObject2()[5]);
            }
        }
        
        final MergingPhenoXOver.OptimizeFunction func2 = new MergingPhenoXOver.OptimizeFunction(config.use6D, colldetect,
                back, gChild2, gChild2.getCartesians(), blowBonds, blowDissoc, motherAbove, config.optBounds, doCDDD);
        
        // normalize
        func2.normalizeFromBounds(config.optBounds,guess2,normG2);
        
        final Tuple<Double,double[]> res2 = opt.doOptimize(normG2, func2);
        if(DEBUG){
            System.out.println("DEBUG: Structure 2: best fitness: " + res2.getObject1());
            if(config.use6D){
                System.out.println("DEBUG: Structure 2: coordinates (normalized) " + res2.getObject2()[0] + "\t" + res2.getObject2()[1] + "\t" + res2.getObject2()[2] + "\t" + res2.getObject2()[3] + "\t" + res2.getObject2()[4] + "\t" + res2.getObject2()[5]);
            }
        }
        
        // we now have two optimized solutions, put them together again
        // no denormalization necessary :-)
        final CartesianCoordinates c1 = func1.getLastCartes();
        final CartesianCoordinates c2 = func2.getLastCartes();
        
        CoordTranslation.updateGeometryFromCartesian(c1, gChild1);
        gChild1.setFather(father.getID());
        gChild1.setMother(mother.getID());        
        gChild1.setFitness(res1.getObject1());
        
        CoordTranslation.updateGeometryFromCartesian(c2, gChild2);
        gChild2.setFather(father.getID());
        gChild2.setMother(mother.getID());        
        gChild2.setFitness(res2.getObject1());
                
        final boolean child1Fine = GlobOptAtomics.areGeometriesStillCorrect(mother,gChild1);
        final boolean child2Fine = GlobOptAtomics.areGeometriesStillCorrect(mother,gChild2);
        if(!child1Fine || !child2Fine){
            System.err.println("ERROR: After crossover, child1 is " + child1Fine + " and child2 is " + child2Fine);
            return new Tuple<>(null,null);
        }
        
        return new Tuple<>(gChild1,gChild2);
    }

    @Override
    public short hasPriority() {
        return -1;
    }

    private static class OptimizeFunction extends AbstractBOBYQAMethod {
        
        private final CollisionDetectionEngine cd;
        private final CartesianFullBackend back;
        private final Geometry stub;
        private final CartesianCoordinates half1;
        private final double blowBonds;
        private final double blowDissoc;
        private final BondInfo bonds;
        private final List<Integer> upIDs;
        private final double[] eparts;
        private final double[] point;
        private final double[][] bounds;
        private final double[][] rotResult;
        private final double[][] transResult;
        private final CartesianCoordinates lastCartes;
        private final int[] starts;
        private final double[] trans;
        private final double[] rot;
        private final double[] xyz1D;
        private final boolean use6D;
        private final boolean doCDDD;
        
        private int iter = 0;
        
        OptimizeFunction(final boolean use6D, final CollisionDetectionEngine cd,
                final CartesianFullBackend back, final Geometry stubHalf2, final CartesianCoordinates half1,
                final double blowBonds, final double blowDissoc,
                final List<Integer> upIDs, final double[][] bounds, final boolean doCDDD){
            this.cd = cd;
            this.back = back;
            this.stub = stubHalf2;
            this.half1 = half1;
            this.blowBonds = blowBonds;
            this.blowDissoc = blowDissoc;
            this.bonds = stubHalf2.getBondInfo();
            this.upIDs = upIDs;
            this.eparts = new double[stub.getNumberOfIndieParticles()];
            this.point = (use6D) ? new double[6] : new double[2];
            this.bounds = bounds;
            this.lastCartes = stub.getCartesians();
            this.rotResult = new double[3][lastCartes.getNoOfAtoms()];
            this.transResult = new double[3][lastCartes.getNoOfAtoms()];  
            this.starts = new int[lastCartes.getNoOfMolecules()];
            final int[] atsPerMol = lastCartes.getAllAtomsPerMol();
            int off = 0;
            for(int mol = 0; mol < starts.length; mol++){
                starts[mol] = off;
                off += atsPerMol[mol];
            }
            this.rot = new double[3];
            this.trans = new double[3];
            this.xyz1D = new double[3*lastCartes.getNoOfAtoms()];
            this.use6D = use6D;
            this.doCDDD = doCDDD;
        }
        
        @Override
        public double computeObjectiveValue(final double[] normalized){
            
            // denormalize
            denormalizeToBounds(bounds, normalized, point);
            
            if(DEBUG){
                System.out.println("DEBUG: Geom merger lower work cartesian:");
                final String[] sa = new Geometry(stub.returnMyConfig()).getCartesians().createPrintableCartesians();
                for(final String s : sa){
                    System.out.println(s);
                }
            }
            
            if(use6D){
                trans[0] = point[0];
                trans[1] = point[1];
                trans[2] = point[2];
                rot[0] = point[3];
                rot[1] = point[4];
                rot[2] = point[5];
            
                // merge
                geomMerger(half1, trans, rot, upIDs, transResult, rotResult, starts, lastCartes);
            } else {
                final double translation = point[0];
                final double rotation = point[1];
                
                // merge
                geomMerger(half1, translation, rotation, upIDs, transResult, rotResult, starts, lastCartes);
            }
                        
            iter++;
            
            if(doCDDD){
                // cd
                final CollisionInfo inf = cd.checkForCollision(lastCartes, blowBonds, bonds);
                if(inf.hasCollision()) {
                    if(DEBUG){
                        System.out.println("DEBUG: Collision detected.");
                        System.out.println("DEBUG: At iter " + iter + " merging reports " + FixedValues.NONCONVERGEDENERGY);
                    }
                    return FixedValues.NONCONVERGEDENERGY;
                }
            
                // dd
                final boolean diss = DissociationDetection.checkForDissociation(inf.getPairWiseDistances(),
                        lastCartes.getAllAtomTypes(), lastCartes.getAllAtomNumbers(), blowDissoc, DissociationDetection.DEFAULTDD);
                // we assume that dissociation is "better" than collision as it proves at least that a collision free state exists
                if(diss) {
                    if(DEBUG){
                        System.out.println("DEBUG: Dissociation detected.");
                        System.out.println("DEBUG: At iter " + iter + " merging reports " + FixedValues.NONCONVERGEDENERGY/2);
                    }
                    return FixedValues.NONCONVERGEDENERGY/2;
                }
            }
            
            lastCartes.getAll1DCartes(xyz1D);
            final double obj = (back == null) ? 0.0 : back.energyCalculation(stub.getID(), iter, xyz1D,
                    lastCartes.getAllAtomTypes(), lastCartes.getAllAtomNumbers(), lastCartes.getAllAtomsPerMol(),
                    eparts, lastCartes.getNoOfAtoms(), lastCartes.getAllCharges(), lastCartes.getAllSpins(), stub.getBondInfo(),
                    lastCartes.containedEnvType() == CartesianCoordinates.ENVTYPE.RIGID);
            
            if(DEBUG){
                System.out.println("DEBUG: At iter " + iter + " translation is " + point[0]);
                System.out.println("DEBUG: At iter " + iter + " rotation is " + point[1]);
                System.out.println("DEBUG: At iter " + iter + " merging reports " + obj);
            }
            
            return obj;
        }
        
        @Override
        public boolean doesNormalize(){
            return true;
        }
        
        private static void geomMerger(final CartesianCoordinates half1,
                final double[] translations, final double[] rotations, List<Integer> upIDs,
                final double[][] transResult, final double[][] rotResult, final int[] starts,
                final CartesianCoordinates fullCluster){

            
            if(DEBUG){
                final String[] sa = half1.createPrintableCartesians();
                System.out.println("DEBUG: Geom merger upper work cartesian:");
                for(final String s : sa){
                    System.out.println(s);
                }
            }
            
            final double[][] xyz1 = half1.getAllXYZCoord();
            
            // first do the translation
            for(int coord = 0; coord < 3; coord++){
                for(int at = 0; at < xyz1[2].length; at++){
                    transResult[coord][at] = xyz1[coord][at] + translations[coord];
                }
            }
            if(DEBUG){
                final CartesianCoordinates cUp = new CartesianCoordinates(half1);
                cUp.setAllXYZ(transResult);
                final String[] sa = cUp.createPrintableCartesians();
                System.out.println("DEBUG: Geom merger translation was " + translations[0] + " / "
                        + translations[1] + " / " + translations[2]);
                System.out.println("DEBUG: Geom merger cartesian after translation:");
                for(final String s : sa){
                    System.out.println(s);
                }
            }
            
            // then do the rotation
            CoordTranslation.rotateXYZ(xyz1, rotations, rotResult);
            if(DEBUG){
                final CartesianCoordinates cUp = new CartesianCoordinates(half1);
                cUp.setAllXYZ(rotResult);
                final String[] sa = cUp.createPrintableCartesians();
                System.out.println("DEBUG: Geom merger rotation angle was " + rotations[0] + " / " + rotations[1] + " / " + rotations[2]);
                System.out.println("DEBUG: Geom merger cartesian after rotation:");
                for(final String s : sa){
                    System.out.println(s);
                }
            }
            
            // we want to return a set of cartesian coordinates, hence figure out where to copy our parts in
            final double[][] fullXYZ = fullCluster.getAllXYZCoord();
            final int[] lengths = fullCluster.getAllAtomsPerMol();
            for(int i = 0; i < upIDs.size(); i++){
                final int molID = upIDs.get(i);
                if(DEBUG){
                    System.out.println("DEBUG: Changing " + molID + " above the cutting plane");
                }
                // get start and lengths
                final int start = starts[molID];
                final int length = lengths[molID];
                
                // copy
                System.arraycopy(rotResult[0], start, fullXYZ[0], start, length);
                System.arraycopy(rotResult[1], start, fullXYZ[1], start, length);
                System.arraycopy(rotResult[2], start, fullXYZ[2], start, length);
            }
            
            if(DEBUG){
                final String[] cUpStr = fullCluster.createPrintableCartesians();
                System.out.println("DEBUG: Geom merger merged cartesian is:");
                for(final String s : cUpStr){
                    System.out.println(s);
                }
            }
        }
        
        private static void geomMerger(final CartesianCoordinates half1,
                final double translation, final double rotation, List<Integer> upIDs,
                final double[][] transResult, final double[][] rotResult,
                final int[] starts, final CartesianCoordinates fullCluster){


            if(DEBUG){
                final String[] sa = half1.createPrintableCartesians();
                System.out.println("DEBUG: Geom merger upper work cartesian:");
                for(final String s : sa){
                    System.out.println(s);
                }
            }

            final double[][] xyz1 = half1.getAllXYZCoord();

            // first do the translation
            for(int at = 0; at < xyz1[2].length; at++){
                transResult[2][at] = xyz1[2][at] += translation;
            }
            
            if(DEBUG){
                final CartesianCoordinates cUp = new CartesianCoordinates(half1);
                cUp.setAllXYZ(transResult);
                final String[] sa = cUp.createPrintableCartesians();
                System.out.println("DEBUG: Geom merger translation was " + translation);
                System.out.println("DEBUG: Geom merger cartesian after translation:");
                for(final String s : sa){
                    System.out.println(s);
                }
            }

            // then do the rotation
            CoordTranslation.rotateXYZAroundZ(xyz1, rotation, rotResult);
            if(DEBUG){
                final CartesianCoordinates cUp = new CartesianCoordinates(half1);
                cUp.setAllXYZ(rotResult);
                final String[] sa = cUp.createPrintableCartesians();
                System.out.println("DEBUG: Geom merger rotation angle was " + rotation);
                System.out.println("DEBUG: Geom merger cartesian after rotation:");
                for(final String s : sa){
                    System.out.println(s);
                }
            }
            
            // we want to return a set of cartesian coordinates, hence figure out where to copy our parts in
            final double[][] fullXYZ = fullCluster.getAllXYZCoord();
            final int[] lengths = fullCluster.getAllAtomsPerMol();
            for(final int molID : upIDs){
                if(DEBUG){
                    System.out.println("DEBUG: Changing " + molID + " above the cutting plane");
                }
                // get start and lengths
                final int start = starts[molID];
                final int length = lengths[molID];

                // copy
                System.arraycopy(rotResult[0], start, fullXYZ[0], start, length);
                System.arraycopy(rotResult[1], start, fullXYZ[1], start, length);
                System.arraycopy(rotResult[2], start, fullXYZ[2], start, length);
            }

            if(DEBUG){
                final String[] cUpStr = fullCluster.createPrintableCartesians();
                System.out.println("DEBUG: Geom merger merged cartesian is:");
                for(final String s : cUpStr){
                    System.out.println(s);
                }
            }
        }
        
        CartesianCoordinates getLastCartes(){
            return lastCartes;
        }
    }
    
    public static class MergingPhenoConfig implements Serializable {
        
        private static final long serialVersionUID = (long) 20150801;
        
        private final boolean use6D;
        public int numberOfInterpolationPoints;
        public double initialTrustRegionRadius;
        public double stoppingTrustRegionRadius;
        public final double[][] optBounds;
        public int noIterations;
        public GlobOptAtomics.CUTTINGMODE planeMode;
        
        public MergingPhenoConfig(final boolean use6D){
            this.use6D = use6D;
            if(use6D){
                numberOfInterpolationPoints = 13; // recommended setting is 2N+1
                initialTrustRegionRadius = 5E-2;
                stoppingTrustRegionRadius = 1E-4;
                optBounds = new double[][]{
                    {-5.0*ANGTOBOHR,-5.0*ANGTOBOHR,-5.0*ANGTOBOHR,-Math.PI,-0.5*Math.PI,-Math.PI},
                    {10*ANGTOBOHR,10*ANGTOBOHR,10*ANGTOBOHR,Math.PI,0.5*Math.PI,Math.PI}};
                noIterations = 200;
                planeMode = GlobOptAtomics.CUTTINGMODE.ZEROZ;
            } else{
                numberOfInterpolationPoints = 5; // recommended setting is 2N+1
                initialTrustRegionRadius = 1E-1;
                stoppingTrustRegionRadius = 1E-5;
                optBounds = new double[][]{
                    {-1.0*ANGTOBOHR,0.0}, {15*ANGTOBOHR,2*Math.PI}};
                noIterations = 100;
                planeMode = GlobOptAtomics.CUTTINGMODE.ZEROZ;
            }
        }
        
        public MergingPhenoConfig(final MergingPhenoConfig orig){
            this.initialTrustRegionRadius = orig.initialTrustRegionRadius;
            this.noIterations = orig.noIterations;
            this.numberOfInterpolationPoints = orig.numberOfInterpolationPoints;
            this.planeMode = orig.planeMode;
            this.stoppingTrustRegionRadius = orig.stoppingTrustRegionRadius;
            this.use6D = orig.use6D;
            this.optBounds = orig.optBounds.clone(); // shallow copy should be sufficient
        }

        @Override
        public String toString(){
            return "\t use 6D?" + this.use6D + "\n\t interpoints: " + this.numberOfInterpolationPoints + "\n\t initial trust: " + this.initialTrustRegionRadius
                    + "\n\t stopping trust " + this.stoppingTrustRegionRadius + "\n\t max iterations: " + this.noIterations
                    + "\n\t plane modus: " + this.planeMode.toString();
        }
    }
}
