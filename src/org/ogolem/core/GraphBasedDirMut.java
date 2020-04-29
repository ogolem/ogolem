/**
Copyright (c) 2013-2014, J. M. Dieterich and B. Hartke
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

import contrib.bobyqa.AbstractBOBYQAMethod;
import contrib.bobyqa.BOBYQAOptimizer;
import java.awt.Color;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.io.*;
import org.ogolem.generic.GenericMutation;
import org.ogolem.math.TrivialLinearAlgebra;
import static org.ogolem.core.CoordTranslation.distance;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * A directed mutation using a graph-based analysis of the cluster in question.
 * @author Johannes Dieterich
 * @author Bernd Hartke
 * @version 2016-12-18
 */
public class GraphBasedDirMut implements GenericMutation<Molecule,Geometry>{
    
    private static final long serialVersionUID = (long) 20140525;
    
    public static final double DEFAULTEULERINCR = Math.toRadians(30);
    public static final double DEFAULTGRIDINCR = 0.5*Constants.ANGTOBOHR;
    public static final double DEFAULTGRIDSIZE = 5.0*Constants.ANGTOBOHR;
    public static final double JUMPDETECTFACTOR = 3.0;
    public static final double DISTSCALEFAC = 1.042;
    
    private static final boolean DEBUGCDDD = true;
    
    private final boolean DEBUG;
    private static final boolean VERBOSEDEBUG = false;
    private final double blowBonds;
    private final CollisionDetectionEngine collDetect;
    private final double blowColl;
    private final Newton newton;
    private final CartesianFullBackend back;
    private final double gridHalfLength;
    private final double gridIncr;
    private final double eulerIncr;
    private final boolean useGridSearch;
    private final boolean useOptSearch;
    private final boolean doLocOpt;
    private final boolean doCDFirst;
    private final double realCD;
    private final double realDD;
    private final GDMBOBYQAConfig bobConf;
    private final boolean fullyRelaxed;
    
    /**
     * 
     * @param blowBonds the blow factor for the connectivity detection
     * @param collDetect the collision detection engine to be used
     * @param blowColl the blow factor for the collision detection
     * @param newton the local optimization engine (including a backend!) to be used
     * @param gridHalfLength the half length of the grid to be used
     * @param gridIncr the grid increment to be used
     * @param eulerIncr the increment for the Euler grid to be used
     * @param useGridApproach should we use a full grid approach or the water specific extension
     * @param useOptApproach should we use a optimizing approach in the grid space?
     * @param realCDBlow the CD blow factor used in the rest of the program, for debugging purposes
     * @param realDDBlow the DD blow factor used in the rest of the program, for debugging purposes
     * @param doLocOptFirst do a local optimization first, prior to the GDM
     * @param doCDFirst do a CD first, AFTER the local optimization (if defined) and the GDM
     * @param fullyRelaxed do a fully relaxed locopt (only used if useOptApproach == true).
     */
    public GraphBasedDirMut(final double blowBonds, final CollisionDetectionEngine collDetect,
            final double blowColl, final Newton newton, final double gridHalfLength,
            final double gridIncr, final double eulerIncr, final boolean useGridApproach,
            final boolean useOptApproach, final double realCDBlow, final double realDDBlow,
            final boolean doLocOptFirst, final boolean fullyRelaxed,
            final boolean doCDFirst){
        this.DEBUG = (GlobalConfig.DEBUGLEVEL > 1);
        this.blowBonds = blowBonds;
        this.collDetect = collDetect;
        this.blowColl = blowColl;
        this.doLocOpt = doLocOptFirst;
        this.doCDFirst = doCDFirst;
        if(newton.getBackend() == null){
            throw new RuntimeException("GDM must have a non-null Backend from Newton!");
        }
        this.newton = newton;
        this.back = newton.getBackend();
        this.gridHalfLength = gridHalfLength;
        this.gridIncr = gridIncr;
        this.eulerIncr = eulerIncr;
	this.useGridSearch = useGridApproach;
        this.useOptSearch = useOptApproach;
        this.realCD = realCDBlow;
        this.realDD = realDDBlow;
        this.bobConf = new GDMBOBYQAConfig(); // make tunable!
        this.fullyRelaxed = fullyRelaxed;
    }
    
    private GraphBasedDirMut(final GraphBasedDirMut orig){
        this.DEBUG = (GlobalConfig.DEBUGLEVEL > 1);
        this.blowBonds = orig.blowBonds;
        this.collDetect = orig.collDetect.clone();
        this.blowColl = orig.blowColl;
        this.back = orig.back.clone();
        this.newton = orig.newton.clone();
        this.gridHalfLength = orig.gridHalfLength;
        this.gridIncr = orig.gridIncr;
        this.eulerIncr = orig.eulerIncr;
	this.useGridSearch = orig.useGridSearch;
        this.useOptSearch = orig.useOptSearch;
        this.realCD = orig.realCD;
        this.realDD = orig.realDD;
        this.doLocOpt = orig.doLocOpt;
        this.bobConf = orig.bobConf;
        this.fullyRelaxed = orig.fullyRelaxed;
        this.doCDFirst = orig.doCDFirst;
    }
    
    @Override
    public GraphBasedDirMut clone(){
        return new GraphBasedDirMut(this);
    }
    
    @Override
    public String getMyID(){
        return "graph based directed mutation: \n\t gridbased? " + useGridSearch 
                + "\n\t blow factor bonds " + blowBonds + "\n\t blow factor collision " + blowColl 
                + "\n\t grid half " + gridHalfLength + "\n\t grid incr " + gridIncr
                + "\n\t Euler incr " + eulerIncr + "\n\t opt based? " + useOptSearch
                + "\n\t CD first? " + doCDFirst + "\n\t opt first? " + doLocOpt;
    }
    
    @Override
    public Geometry mutate(final Geometry geom){
        
        if(DEBUG){System.out.println("DEBUG: Running our fancy graph based directed mutation.");}
                
        Geometry work = new Geometry(geom);
        
        if(doLocOpt){
            if(DEBUG){System.out.println("DEBUG: doing a local optimization first.");}
            work = newton.localOptimization(work);
        }
        
        if(doCDFirst){
            if(DEBUG){System.out.println("DEBUG: doing a CD now.");}
            final boolean hasColl = collDetect.checkOnlyForCollision(work.getCartesians(),
                    blowColl, work.getBondInfo());
            if(hasColl){
                if(DEBUG){System.out.println("DEBUG: found a collision, returning null.");}
                return null;
            }
        }
        
        final CartesianCoordinates c = work.getCartesians();
        final BondInfo bonds = work.getBondInfo();
        
        // get the energy before we manipulate
        final double[] eparts = new double[c.getNoOfMolecules()];
        final double eBefore = back.energyCalculation(geom.getID(), 42, c.getAll1DCartes(),
                c.getAllAtomTypes(), c.getAllAtomNumbers(), c.getAllAtomsPerMol(), eparts,
                c.getNoOfAtoms(), c.getAllCharges(), c.getAllSpins(), bonds);
        
        // create the connectivity matrix for the geometry (so in between molecules)
        // 1) get the full one
        final SimpleBondInfo fullBonds = CoordTranslation.checkForBonds(c, blowBonds);
        final boolean[][] full = fullBonds.getFullBondMatrix();
        
        if(DEBUG){
            System.out.println("Bonds with blow fator " + blowBonds);
            for (final boolean[] fullRow : full) {
                System.out.println(Arrays.toString(fullRow));
            }
        }
        
        // figure out the number of connections from each molecule to another one (do not count internal bonds)
        final int[] connections = new int[c.getNoOfMolecules()];
        final int[] atsPerMol = c.getAllAtomsPerMol();
        int[][] connCounter = null;
        if(DEBUG){
            connCounter = new int[c.getNoOfMolecules()][c.getNoOfMolecules()];
        }
        int atCounter = 0;
        for(int i = 0; i < c.getNoOfMolecules(); i++){
            final int ats = atsPerMol[i];
            for(int at = atCounter; at < atCounter+ats; at++){
                // go over the bonds (remove self interactions!)
                // we DO count if multiple atoms of the same molecule connect to another molecule
                // not only is this easier, but it makes more sense for the application at hand
                // which is to find the least connected molecule
                for(int x = 0; x < c.getNoOfAtoms(); x++){
                    if(full[at][x]){
                        if(x < atCounter || x >= atCounter+ats){
                            if(DEBUG){
                                
                                // figure out the molecule this atom belongs to
                                int otherMol = -1;
                                int count = 0;
                                for(int y = 0; y < c.getNoOfAtoms(); y++){
                                    count += atsPerMol[y];
                                    if(count > x){
                                        otherMol = y;
                                        break;
                                    }
                                }
                                
                                System.out.println("Connection found for molecule " + i + " "
                                    + "to molecule " + otherMol);
                                connCounter[i][otherMol]++;
                                connCounter[otherMol][i]++;
                            }
                            connections[i]++;
                        }
                    }
                }
            }
            atCounter += ats;
        }
        
        // figure out the least connected molecule in the cluster (easy task...)
        int leastConnID = -1;
        int noLeast = Integer.MAX_VALUE;
        for(int i = 0; i < c.getNoOfMolecules(); i++){
            if(noLeast > connections[i]){
                leastConnID = i;
                noLeast = connections[i];
            }
        }
        
        if(DEBUG){
            
            int secLeastConnID = -1;
            int secNoLeast = Integer.MAX_VALUE;
            for(int i = 0; i < c.getNoOfMolecules(); i++){
                if(secNoLeast > connections[i]){
                    if(i == leastConnID){continue;}
                    secLeastConnID = i;
                    secNoLeast = connections[i];
                }
            }
            
            final float[] colEnOne = new float[]{0.0f,0.0f,1.0f};
            final float[] colEnTwo = new float[]{1.0f,0.0f,0.0f};
            final float[] colTrOne = new float[]{0.0f,0.0f,1.0f};
            final float[] colTrTwo = new float[]{1.0f,0.0f,0.0f};
            final float transparency = 1.0f;
            final String backGround = "white";
            final String circleRing = "black";
            final String vertexShape = "circle";
            final String vertexStyle = "filled";
            final String edgeStyle = "filled";
            final int penWidth = 3;
            
            // dot file output :-)
            System.out.println("graph G {");
            System.out.println("   bgcolor = \"" + backGround + "\";");
            
            // first the vertices
            int maxConns = 0;
            for(int i = 0; i < connections.length; i++){maxConns = Math.max(maxConns,connections[i]/2);}
            final float[] colDef = new float[3];
            for(int i = 0; i < connections.length; i++){
                // calculate our color in RGB space
                /*final double p = connections[i]/(2*maxConns);
                colDef[0] = (float) (colEnTwo[0] * p + colEnOne[0] * (1 - p));
                colDef[1] = (float) (colEnTwo[1] * p + colEnOne[1] * (1 - p));
                colDef[2] = (float) (colEnTwo[2] * p + colEnOne[2] * (1 - p));
                if(colDef[0] > 1.0f){colDef[0] = 1.0f;}
                if(colDef[1] > 1.0f){colDef[1] = 1.0f;}
                if(colDef[2] > 1.0f){colDef[2] = 1.0f;}
                
                final String hex = colorToHex(colDef, transparency);*/
                String hex;
                if(i == leastConnID){
                    hex = "#FF0000";
                } else if(i == secLeastConnID){
                    hex = "#003366";
                } else {
                    hex = "#FFFFFF";
                }
                
                final String labelStr = "" + connections[i];
                System.out.println("   " + (i+1) + " [label=\"" + labelStr + "\", shape = \"" + vertexShape + "\", style = \"" + vertexStyle
                    + "\", color = \"" + circleRing + "\", fillcolor = \"" + hex + "\"];");
            }
            
            // now the edges
            for(int i = 0; i < connections.length-1; i++){
                for(int j = i+1; j < connections.length; j++){
                    if(connCounter[i][j] <= 0) {continue;}
                    
                    /*final double p = connCounter[i][j]/(2*maxConns);
                    colDef[0] = (float) (colTrTwo[0] * p + colTrOne[0] * (1 - p));
                    colDef[1] = (float) (colTrTwo[1] * p + colTrOne[1] * (1 - p));
                    colDef[2] = (float) (colTrTwo[2] * p + colTrOne[2] * (1 - p));
                    if(colDef[0] > 1.0f){colDef[0] = 1.0f;}
                    if(colDef[1] > 1.0f){colDef[1] = 1.0f;}
                    if(colDef[2] > 1.0f){colDef[2] = 1.0f;}
                    
                    final String hex = colorToHex(colDef, transparency);*/
                    final String hex = "#000000";
            
                    System.out.println("   " + (i+1) + " -- " + (j+1) + " [style = \""
                            + edgeStyle + "\", color = \"" + hex + "\", penwidth = \""
                            + penWidth + "\"];");
                }
            }
            
            System.out.println("}");
        }
        
        if(DEBUG){
            System.out.println("DEBUG: Following molecular connection information obtained:");
            for(int i = 0; i < c.getNoOfMolecules(); i++){
                System.out.println(" " + i + "\t" + connections[i]);
            }
        }
        
        if(DEBUG){
            System.out.println("DEBUG: Identified molecule " + leastConnID
                    + " to be the least connected molecule in the cluster with only "
                    + noLeast + " connections.");
        }
        
        // now the hard task: figure out where this molecule should go
        // as this is currently developed code, it is not quite as clean as it
        // could be
        
        /*
         * ATTEMPT ONE: FIND THE SECOND WORST MOLECULE AND MOVE THE MOLECULES IN A
         * GRID AROUND THAT COM
         */
        int secLeastConnID = -1;
        int secNoLeast = Integer.MAX_VALUE;
        for(int i = 0; i < c.getNoOfMolecules(); i++){
            if(secNoLeast > connections[i]){
                if(i == leastConnID){continue;}
                secLeastConnID = i;
                secNoLeast = connections[i];
            }
        }
        
        if(DEBUG){
            System.out.println("DEBUG: Identified molecule " + secLeastConnID
                    + " with " + secNoLeast + " connections to be the second least favourite.");
        }
        
        if(useOptSearch){
            
            if(VERBOSEDEBUG){System.out.println("DEBUG: Starting with energy " + work.getFitness());}
            
            // use bobyqa for a quick optimization of that one least connected molecule
            
            final double[] comFirst = work.getCOM(leastConnID); // this we will tinker with
            final double[] eulerFirst = work.getEulers(leastConnID); // this we will also tinker with
            final double[] comSec = work.getCOM(secLeastConnID);
            final double[] guess = new double[6];
            
            final double[][] boundaries = new double[2][6];
            boundaries[0][0] = -gridHalfLength + comSec[0];
            boundaries[1][0] = gridHalfLength + comSec[0];
            
            boundaries[0][1] = -gridHalfLength + comSec[1];
            boundaries[1][1] = gridHalfLength + comSec[1];
            
            boundaries[0][2] = -gridHalfLength + comSec[2];
            boundaries[1][2] = gridHalfLength + comSec[2];
            // phi
            boundaries[0][3] = -Math.PI;
            boundaries[1][3] = Math.PI;
            // omega
            boundaries[0][4] = -0.5*Math.PI;
            boundaries[1][4] = 0.5*Math.PI;
            // psi
            boundaries[0][5] = -Math.PI;
            boundaries[1][5] = Math.PI;
            
            final BOBYQAOptimizer bobyqa = new BOBYQAOptimizer(boundaries, bobConf.interPoints,
                    bobConf.initialTrust, bobConf.stoppingTrust, bobConf.maxIter,
                    true); // we do normalize
            
            // do something...
            final Lottery random = Lottery.getInstance();
            guess[0] = (random.nextBoolean()) ? comSec[0] - 2.0*random.nextDouble() :
                    comSec[0] + 2.0*random.nextDouble();
            guess[1] = (random.nextBoolean()) ? comSec[1] - 2.0*random.nextDouble() :
                    comSec[1] + 2.0*random.nextDouble();
            guess[2] = (random.nextBoolean()) ? comSec[2] - 2.0*random.nextDouble() :
                    comSec[2] + 2.0*random.nextDouble();
            final double[] eulers = new double[3];
            RandomUtils.randomEulers(eulers);
            guess[3] = eulers[0];
            guess[4] = eulers[1];
            guess[5] = eulers[2];
            
            final GDMAdaptor adap = new GDMAdaptor(work, leastConnID,
                newton, collDetect, guess, boundaries, blowColl, fullyRelaxed);
            
            final double[] init = new double[6];
            adap.normalizeFromBounds(boundaries, guess, init);
            
            final Tuple<Double, double[]> opt = bobyqa.doOptimize(init, adap);

            final double currBestE = opt.getObject1();
            
            if (currBestE >= FixedValues.NONCONVERGEDENERGY) {
                System.out.println("INFO: Could not find any place that was collision free and had a reasonable energy. Complaining.");
                return new Geometry(geom);
            }

            // work is already setup properly, except for the energy :-)
            work = adap.getCurrentGeom();

            if (currBestE >= eBefore) {
                System.out.println("INFO: Could not find any place with a better energy than we had BEFORE spending CPU cycles. Continuing but complaining.");
                // however, this comparison MAY BE unfair since the energy before MAY BE locopted, now it is not!
                // additionally, even if the new energy is somewhat worse than the previous one, this is 
                // NOT necessarily bad but it may still be usable _and_ in line with the intentions of
                // this routine: What we aim at is moving badly positioned molecules (that are "stuck" where
                // they are during normal xover/mutation events) to non-stupid new locations. Of course
                // ideally these new locations could be instantly better, but we still profit if this effect
                // manifests itself only after several further xover/mutation events.
            }
            
            work.setFitness(currBestE);
            
            return work;
        }

        /*
         * DEFAULT: USE GRID SEARCH. ALTHOUGH BEING RATHER EXPENSIVE, THIS IS THE
         * GENERIC WAY TO MUTATE.
         */
        if(useGridSearch){
            
            final double[] comFirst = work.getCOM(leastConnID); // this we will tinker with
            final double[] eulerFirst = work.getEulers(leastConnID); // this we will also tinker with
            final double[] comSec = work.getCOM(secLeastConnID);
            double currBestE = Double.MAX_VALUE;
            final double[] currBest = new double[6];

            /*
             * Here is the issue: the euler angles. it is fairly simple to find the
             * best new spot for an atomic cluster, as we only need to evaluate,
             * which spots are collision free and (from those) which one has the
             * lowest energy. Already figuring out which spot is collision free and
             * which one is not is more difficult with a molecule where the Euler
             * orientation plays a role. I guess one COULD discretize the Euler
             * angles with a rather crude grid themselves (say 30 degree increments)
             * and then go over these as well? That DOES remove quite a bit of
             * problems at a rather large cost of numbers of backend evaluations
             * (well, and obviously coll detect, but that is cheap)
             * lets try that...
             * 
             * another option would be to run a BOBYQA in that grid and find a set
             * of COM and Euler coordinates. Problem: will only find some minimum,
             * not the best one which may matter for non-trivial molecules
             */
            final int noAts = work.getMolecules().get(leastConnID).getNumberOfAtoms();
            final int noPoints = (int) Math.ceil(2*gridHalfLength / gridIncr);
            final int noPointsPhi = (noAts == 1) ? 1 : (int) Math.floor(2 * Math.PI / eulerIncr);
            final int noPointsOmega = (noAts == 1) ? 1 : (int) Math.floor(Math.PI / eulerIncr);
            final int noPointsPsi = (noAts == 1) ? 1 : (int) Math.floor(2 * Math.PI / eulerIncr);
            int evalCounter = 0;
            int collCounter = 0;
            for (int x = 0; x < noPoints; x++) {
                final double currX = comSec[0] - gridHalfLength + x * gridIncr;
                comFirst[0] = currX;
                for (int y = 0; y < noPoints; y++) {
                    final double currY = comSec[1] - gridHalfLength + y * gridIncr;
                    comFirst[1] = currY;
                    for (int z = 0; z < noPoints; z++) {
                        final double currZ = comSec[2] - gridHalfLength + z * gridIncr;
                        comFirst[2] = currZ;
                        for (int phi = 0; phi < noPointsPhi; phi++) {
                            final double currPhi = -Math.PI + phi * eulerIncr;
                            eulerFirst[0] = currPhi;
                            for (int omega = 0; omega < noPointsOmega; omega++) {
                                final double currOmega = -Math.PI / 2 + omega * eulerIncr;
                                eulerFirst[1] = currOmega;
                                for (int psi = 0; psi < noPointsPsi; psi++) {
                                    final double currPsi = -Math.PI + psi * eulerIncr;
                                    eulerFirst[2] = currPsi;

                                    // ok, check for collisions and energy
                                    final CartesianCoordinates currC = work.getCartesians();
                                    final boolean hasColl = collDetect.checkOnlyForCollision(currC, blowColl, bonds);
                                    if (DEBUG) {
                                        System.out.println("DEBUG: For " + currX + "\t"
                                                + currY + "\t" + currZ + "\t" + currPhi + "\t"
                                                + currOmega + "\t" + currPsi + ": has collision? "
                                                + hasColl);
                                    }
                                    if (!hasColl) {
                                        if (DEBUG) {
                                            System.out.println("DEBUG: Evaluating energy of this set...");
                                        }
                                        final double eCurr = back.energyCalculation(geom.getID(), evalCounter, currC.getAll1DCartes(),
                                                c.getAllAtomTypes(), c.getAllAtomNumbers(), c.getAllAtomsPerMol(), eparts,
                                                c.getNoOfAtoms(), c.getAllCharges(), c.getAllSpins(), bonds);
                                        evalCounter++;
                                        if (DEBUG) {
                                            System.out.println("DEBUG: Evaluated energy is " + eCurr);
                                        }
                                        if (eCurr < currBestE) {
                                            // store this state
                                            currBestE = eCurr;
                                            currBest[0] = currX;
                                            currBest[1] = currY;
                                            currBest[2] = currZ;
                                            currBest[3] = currPhi;
                                            currBest[4] = currOmega;
                                            currBest[5] = currPsi;
                                        }
                                    } else {
                                        collCounter++;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (DEBUG) {
                System.out.println("DEBUG: Number of collisions found " + collCounter);
                System.out.println("DEBUG: Number of energy evaluations done " + evalCounter);
            }

            if (currBestE >= FixedValues.NONCONVERGEDENERGY) {
                System.out.println("INFO: Could not find any place that was collision free and had a reasonable energy. Complaining.");
                return new Geometry(geom);
            }

            comFirst[0] = currBest[0];
            comFirst[1] = currBest[1];
            comFirst[2] = currBest[2];
            eulerFirst[0] = currBest[3];
            eulerFirst[1] = currBest[4];
            eulerFirst[2] = currBest[5];

            if (currBestE >= eBefore) {
                System.out.println("INFO: Could not find any place with a better energy than we had BEFORE spending CPU cycles. Continuing but complaining.");
                // however, this comparison MAY BE unfair since the energy before MAY BE locopted, now it is not!
                // additionally, even if the new energy is somewhat worse than the previous one, this is 
                // NOT necessarily bad but it may still be usable _and_ in line with the intentions of
                // this routine: What we aim at is moving badly positioned molecules (that are "stuck" where
                // they are during normal xover/mutation events) to non-stupid new locations. Of course
                // ideally these new locations could be instantly better, but we still profit if this effect
                // manifests itself only after several further xover/mutation events.
            }
            
            work.setFitness(currBestE);

            return work;
        }
        
        /*
         * make sure that we are dealing with water molecules in the expected atom
         * order...
         */
        for(final Molecule m : work.getMolecules()){
            final String[] atoms = m.getAtomTypes();
            if(!(atoms.length == 3 &&
                    atoms[0].equalsIgnoreCase("O") &&
                    atoms[1].equalsIgnoreCase("H") &&
                    atoms[2].equalsIgnoreCase("H"))){
                throw new RuntimeException("Water specialized graph-based directed mutation should better be"
                        + " used for a WATER system in O/H/H order. This system does not qualify.");
            }
        }
        
	/*
         * THIS IS THE OTHER OPTION: BERND'S WATER SPECIFIC SOLUTION.
         * first find s.th. like a representative COM-COM distance,
         * by scanning all these distances
         */
        final int numMolecules = c.getNoOfMolecules();
        final double[] allCOMdists = new double[(int)(numMolecules*(numMolecules-1)*0.5)];
        int counter = 0;
        for(int i = 0; i < numMolecules; i++){
            for(int j = 0; j < i; j++){
                allCOMdists[counter] = distance(work.getCOM(i),work.getCOM(j));
                counter++;
            }
        }
        Arrays.sort(allCOMdists);
        
        if(DEBUG){
            try (final PrintWriter comOutput = new PrintWriter("com.distances")) {
                System.out.println("DEBUG: COM-COM distances: written to file com.distances.");
                for(int i = 0; i < counter-1; i++){
                    comOutput.println(" " + i + "\t" + allCOMdists[i] 
                            + "\t" + (allCOMdists[i+1]-allCOMdists[i]));
                }
            } catch (FileNotFoundException e){
                System.out.println("FileNotFoundException: probably this file could not be created.");
                System.out.println(e.getMessage());
                e.printStackTrace(System.err);
            }
        }
        
        // now nail down the first "jump" in these distances
        double maxDiff = 0.0;
        double maxDiffOld = 0.0;
        int iJump = -1;
        for(int i = 0; i < counter-1; i++){
            double thisDiff = Math.abs(allCOMdists[i+1]-allCOMdists[i]);
            if(thisDiff>maxDiff){
                maxDiffOld = maxDiff;
                maxDiff = thisDiff;
            }
            // this factor is NOT universal, sadly...
            if(maxDiffOld>0 && maxDiff> JUMPDETECTFACTOR*maxDiffOld){
                iJump = i;
                break;
            }
        }
        // calculate the mean nearest-neighbor distance as the average up to the jump:
        // (the next line tries to make the above foolproof, despite the non-universal scaling factor)
        double distCOMaverage = allCOMdists[0]*DISTSCALEFAC;
        if(iJump>-1){
            distCOMaverage = 0.0;
            for(int i=0; i < iJump; i++){
                distCOMaverage += allCOMdists[i];
            }
            distCOMaverage /= (double)iJump;
        }
        if(DEBUGCDDD){
            System.out.println("DEBUG: average nearest-neighbor COM-COM distance: " + distCOMaverage);
        }
                
        // now find those molecules that are within a * distCOMaverage from
        // the second-least-connected molecule (note that a>2 makes no sense for
        // the following, and neither does a<1, but within 1<a<2 my choice of
        // a=1.5 is fairly arbitrary...
        final double factor = 1.5;
        final double distNeigh = factor * distCOMaverage;
        final List<Integer> neighborHood = new ArrayList<>();
        for(int i = 0; i < secLeastConnID; i++){
            if(distance(work.getCOM(secLeastConnID),work.getCOM(i)) < distNeigh){
                neighborHood.add(i);
            }
        }
         for(int i = secLeastConnID+1; i < numMolecules; i++){
            if(distance(work.getCOM(secLeastConnID),work.getCOM(i)) < distNeigh){
                neighborHood.add(i);
            }
        }
        if(DEBUG){
            System.out.println("DEBUG: " + neighborHood.size() + 
                    " molecules found near to 2nd-least-connected:" + neighborHood.toString());
        }
        
        // now loop over unique pairs of these, and try to place the worst one over or below
        // each triangle of molecules (triangle = the 2nd-least-connected and the other two)
        // JMD NOTE: THIS IS INDEED RATHER UGLY. LET ME THINK HOW WE CAN DO THIS PROPERLY...
        final double[] comSec = work.getCOM(secLeastConnID);
        double currBestE = Double.MAX_VALUE;
        int evalCounter = 0;
        int collCounter = 0;
        int tryCounter = 0;
        final double[][] comTriple = new double[3][3];
        final double[][] newPos = new double[2][3];
        final double[] xyzNewO = new double[3];
        final double[][] xyzNewHH = new double[2][3];
        final double[][] xyzBest = new double[3][3];
        // from HERE to THERE this is a probably stupid way of generating
        // a sub-loop over index pairs (1,2), (1,3) and (2,3), used below
        final List<int[]> indexPair = new ArrayList<>();
        final int[] iPair1 = {0,1};
        indexPair.add(iPair1);
        final int[] iPair2 = {0,2};
        indexPair.add(iPair2);
        final int[] iPair3 = {1,2};
        indexPair.add(iPair3);
        // THERE
        final CartesianCoordinates currC = work.getCartesians();
        for(final int i : neighborHood){
            for (final int j : neighborHood){
                if(i<j){ // I suspect this can be done more gracefully with foreach(list) but I do not know how...
                    if(DEBUG){
                        System.out.println("DEBUG: " + i + "\t"+ j);
                    }
                    comTriple[0] = comSec;
                    comTriple[1] = work.getCOM(i);
                    comTriple[2] = work.getCOM(j);
                    
                    final short ret = giveNewPos(comTriple[0],comTriple[1],comTriple[2],distCOMaverage,newPos);
                    if(ret < 0){continue;}
                    // put new points in work/Coords and calc E...
                    // do this for both positions returned in newPos
                    for(int pair = 0; pair < 2; pair++){
                        // and for all pairs 12,13,23 in the comTriple
                        for(final int[] nmPair : indexPair){
                            final int n = nmPair[0];
                            final int m = nmPair[1];
                            // here smuggle in the new coords for molecule leastConnID;
                            // this is done in 2 VERY crude and H2O-specific(!!) steps at the moment:
                            // 1) assume that the O-atom is at(!) COM:
                            System.arraycopy(newPos[pair],0,xyzNewO,0,3);
                            currC.setXYZCoordinatesOfAtom(xyzNewO,3*leastConnID);
                            assert(!Double.isNaN(xyzNewO[0]));
                            assert(!Double.isNaN(xyzNewO[1]));
                            assert(!Double.isNaN(xyzNewO[2]));
                            // 2) place the two H-atoms
                            giveNewHH(xyzNewO,comTriple[n],comTriple[m], xyzNewHH);
                            currC.setXYZCoordinatesOfAtom(xyzNewHH[0],3*leastConnID+1);
                            assert(!Double.isNaN(xyzNewHH[0][0]));
                            assert(!Double.isNaN(xyzNewHH[0][1]));
                            assert(!Double.isNaN(xyzNewHH[0][2]));
                            currC.setXYZCoordinatesOfAtom(xyzNewHH[1],3*leastConnID+2);
                            assert(!Double.isNaN(xyzNewHH[1][0]));
                            assert(!Double.isNaN(xyzNewHH[1][1]));
                            assert(!Double.isNaN(xyzNewHH[1][2]));
                            
                            // testing:
                            //final double distOH1 = distance(xyzNewO,xyzNewH1)*Constants.BOHRTOANG;
                            //final double distOH2 = distance(xyzNewO,xyzNewH2)*Constants.BOHRTOANG;
                            //final double halfDistHH = 0.5*distance(xyzNewHH[0],xyzNewHH[1])*Constants.BOHRTOANG;
                            tryCounter++;
                            final boolean hasColl = collDetect.checkOnlyForCollision(currC, blowColl, bonds);
                            if(DEBUG){
                                System.out.println("DEBUG: triple (" + secLeastConnID + "," + i + "," + j + "), pair:" 
                                                + pair + ", idxPair:" + Arrays.toString(nmPair)+ ": has collision? "
                                                + hasColl);
                            }
                            if(!hasColl){
                                if(DEBUG){
                                    System.out.println("DEBUG: Evaluating energy of this set...");
                                }
                                final double eCurr = back.energyCalculation(geom.getID(), evalCounter, currC.getAll1DCartes(),
                                    c.getAllAtomTypes(), c.getAllAtomNumbers(), c.getAllAtomsPerMol(), eparts,
                                    c.getNoOfAtoms(), c.getAllCharges(), c.getAllSpins(), bonds);
                                evalCounter++;
                                if(DEBUG){
                                    System.out.println("DEBUG: Evaluated energy is " + eCurr);
                                }
                                if(eCurr < currBestE){
                                    currBestE = eCurr;
                                    System.arraycopy(xyzNewO,0,xyzBest[0],0,3);
                                    System.arraycopy(xyzNewHH[0],0,xyzBest[1],0,3);
                                    System.arraycopy(xyzNewHH[1],0,xyzBest[2],0,3);
                                }
                            } else{
                                collCounter++;
                            }
                        }
                    }
                }
            }
        }

        currC.setXYZCoordinatesOfAtom(xyzBest[0],3*leastConnID);
        currC.setXYZCoordinatesOfAtom(xyzBest[1],3*leastConnID+1);
        currC.setXYZCoordinatesOfAtom(xyzBest[2],3*leastConnID+2);

        // more testing: print out final coords, directly from the Cartesian object,
        // simply using the createPrintableCartesians method!
        if(DEBUG){
            System.out.println("DEBUG inside GraphDirMut: modified geometry coming...");
            final String[] saXYZfile = currC.createPrintableCartesians();
            for(final String s : saXYZfile) System.out.println(s);            
        }
        
        // now translate these modified coordinates to the geometry object work
        // which is then returned to the caller:
        work = CoordTranslation.cartesianToGeometry(currC, currC.getNoOfMolecules(), currC.getAllAtomsPerMol(), 
                work.getAllFlexies(), work.getExplicitDoFs(), work.getAllConstraints(true), work.getAllConstraintsXYZ(true), 
                work.getSIDs(), bonds);
        
        if(DEBUG){
            System.out.println("DEBUG: SUMMARY (FOR H2O GRAPH DIR MUT):");
            System.out.println("molecule no. " + (leastConnID+1) + " was moved, "
                    + "with atoms " + (leastConnID*3+1) + ", " + (leastConnID*3+2)
                    + ", " + (leastConnID*3+3));
            System.out.println("total new positions tested: " + (collCounter+evalCounter)
                    + ". Of these, " + collCounter + " had collisions,"
                    + evalCounter + " energy evaluations.");
            System.out.println("old energy (potentially after locopt) : " + 
                    (eBefore*Constants.HARTREETOKJ) + " kJ/mol");
            System.out.println("new energy (before locopt): "+ 
                    (currBestE*Constants.HARTREETOKJ) + "kJ/mol");
        }

        if(currBestE >= FixedValues.NONCONVERGEDENERGY){
            System.out.println("INFO: Could not find any place that was collision free and had a reasonable energy. Complaining.");
            return new Geometry(geom);
        }
        
        if(currBestE >= eBefore){
            System.out.println("INFO: Could not find any place with a better energy than we had BEFORE spending CPU cycles. Continuing but complaining.");
            // however, this comparison MAY BE unfair since the energy before MAY BE locopted, now it is not!
            // additionally, even if the new energy is somewhat worse than the previous one, this is 
            // NOT necessarily bad but it may still be usable _and_ in line with the intentions of
            // this routine: What we aim at is moving badly positioned molecules (that are "stuck" where
            // they are during normal xover/mutation events) to non-stupid new locations. Of course
            // ideally these new locations could be instantly better, but we still profit if this effect
            // manifests itself only after several further xover/mutation events.
        }
        
        if(currBestE < eBefore && true){
            System.out.println("DEBUG: GraphDirMut ran successfully and obtained a lower energy " + currBestE + " vs " + eBefore + " improvement " + (eBefore-currBestE));
        }
        
        work.setFitness(currBestE);
        
        if(DEBUGCDDD){
            
            System.out.println("INFO: We found " + collCounter + " collisions in  " + tryCounter  + " tries.");
            // just for kicks, we are running CD/DD on the current work geometry
            final GeometrySanityCheck sanity = new GeometrySanityCheck(realCD, realDD, true, true);
            
            final boolean isSane = sanity.isSane(work);
            if(!isSane){
                System.err.println("WARNING: Geometry " + work.getID() + " is deemed "+
                        "insane by the GDM-internal sanity check and will be removed "+
                        "in the algorithm.");
            } else {
                System.out.println("INFO: Geometry " + work.getID() + " is deemed sane. "+
                        "Everything fine.");
            }
        }
        
        return work;
    }
    
    /**
     * Returns two new 3D-positions, given three initial positions and a distance dist.
     * Both of the new positions have distance dist from the three given positions.
     * The three given positions may have pairwise distances different from dist
     * but not exceeding 2*dist (in which case we signal an error)
     * @param aCoord cartesian 3D-coords of point1
     * @param bCoord cartesian 3D-coords of point2
     * @param cCoord cartesian 3D-coords of point3
     * @param dist the desired distance between the new points and the 3 given ones
     * @param newPos the two new points, the first one in [0][0:2], the second one in [1][0:2]
     * @return 0 if everything ok, -1 otherwhise.
     */
    private static short giveNewPos(final double[] aCoord, 
            final double[] bCoord, final double[] cCoord, final double dist,
            final double[][] newPos){
        
        assert(newPos.length == 2);
        assert(newPos[0].length == 3);
        
        // check input coords for their pair distances, see comments below
        final double dLimit = 2.0 * dist;
        if(distance(aCoord,bCoord)>=dLimit){
//            System.err.println("WARNING: GraphBasedDirMut distance a-b too large");
            return -1;
        }
        if(distance(aCoord,cCoord)>=dLimit){
//            System.err.println("WARNING: GraphBasedDirMut distance a-c too large");
            return -1;
        }
        if(distance(bCoord,cCoord)>=dLimit){
//            System.err.println("WARNING: GraphBasedDirMut distance b-c too large");
            return -1;
        }
        
        // now we are sure that a,b,c are close enough to each other,
        // so now calculate possible new positions: 
        // first find the center p of the circumscribed circle of a,b,c
        // (rCircCirc is the radius of this circle)
        final double[] ab = new double[3];
        final double[] ac = new double[3];
        //final double[] bc = new double[3];
        final double abDist = distance(aCoord,bCoord);
        ab[0] = (bCoord[0]-aCoord[0])/abDist;
        ab[1] = (bCoord[1]-aCoord[1])/abDist;
        ab[2] = (bCoord[2]-aCoord[2])/abDist;
        final double acDist = distance(aCoord,cCoord);
        ac[0] = (cCoord[0]-aCoord[0])/acDist;
        ac[1] = (cCoord[1]-aCoord[1])/acDist;
        ac[2] = (cCoord[2]-aCoord[2])/acDist;
        final double bcDist = distance(bCoord,cCoord);
        //bc[0] = (cCoord[0]-bCoord[0])/bcDist; not needed
        //bc[1] = (cCoord[1]-bCoord[1])/bcDist;
        //bc[2] = (cCoord[2]-bCoord[2])/bcDist;
        final double cosa = TrivialLinearAlgebra.dotProduct(ab,ac);
        // final double sina = Math.sin(Math.acos(cosa)); or alternatively:
        // (note that by definition of sin and cos, cosa*cosa<=1.0 is ALWAYS true,
        // therefore this sqrt argument can never become negative)
        final double sina = Math.sqrt(1.0-cosa*cosa);
        final double rCircCirc = 0.5 * bcDist / sina;
        // for the statement h=Math.sqrt(e*e-b*b) below to NOT fail,
        // we need to ensure that the INPUT dist is larger than the
        // radius of the circumscribed circle. Since the latter can become
        // arbitrarily large even if all triangle distance remain small,
        // we have to supply this additional failure return:
        if(rCircCirc>dist){
//            System.err.println("WARNING: circumscribed circle too large in GDM");
            return -1;
        }
        // (the distance check at the beginning of this routine ensures
        // that the argument of this sqrt cannot be negative either)
        final double b = Math.sqrt(rCircCirc*rCircCirc-abDist*abDist*0.25);
        final double factor1 = acDist*cosa;
        final double[] nCoord = new double[3];
        nCoord[0] = aCoord[0] + factor1 * ab[0];
        nCoord[1] = aCoord[1] + factor1 * ab[1];
        nCoord[2] = aCoord[2] + factor1 * ab[2];
        final double[] nc = new double[3];
        final double ncDist = distance(cCoord,nCoord);
        nc[0] = (cCoord[0]-nCoord[0])/ncDist;
        nc[1] = (cCoord[1]-nCoord[1])/ncDist;
        nc[2] = (cCoord[2]-nCoord[2])/ncDist;
        final double factor2 = 0.5 * abDist;
        final double[] pCoord = new double[3];
        pCoord[0] = aCoord[0] + factor2 * ab[0] + b * nc[0];
        pCoord[1] = aCoord[1] + factor2 * ab[1] + b * nc[1];
        pCoord[2] = aCoord[2] + factor2 * ab[2] + b * nc[2];
        
        // now find the distance h of the new positions from the a,b,c-plane
        // (again, for this sqrt, the argument cannot be negative due to the
        // dLimit distance check at the beginning of this routine)
        final double e = Math.sqrt(dist*dist-abDist*abDist*0.25);
        // (to ensure non-negativity of this sqrt argument, we do the
        // second test and return(-1) above, with rCircCirc...)
        final double h = Math.sqrt(e*e-b*b);
        
        // construct a unit normal vector on the a,b,c,-plane
        final double[] normal = new double[3];
        TrivialLinearAlgebra.crossProduct(ab,ac,normal);
        final double normalDist = Math.sqrt(TrivialLinearAlgebra.dotProduct(normal, normal));
        normal[0] /= normalDist;
        normal[1] /= normalDist;
        normal[2] /= normalDist;
        
        // finally place the new positions atop and below the triangle a,b,c
        newPos[0][0] = pCoord[0] + h * normal[0];
        newPos[0][1] = pCoord[1] + h * normal[1];
        newPos[0][2] = pCoord[2] + h * normal[2];
        newPos[1][0] = pCoord[0] - h * normal[0];
        newPos[1][1] = pCoord[1] - h * normal[1];
        newPos[1][2] = pCoord[2] - h * normal[2];
        
        return 0;
    }
    
    private static final double OHDIST = 0.5898*Constants.ANGTOBOHR;
    private static final double DHDIST = 0.7827*Constants.ANGTOBOHR;
    
    /**
     * Given three 3D-points a,b,c, this returns two new points that are 
     * (1) in the a,b,c-plane, and
     * (2) as close as possible to the a-b,a-c connection lines while still
     * (3) fulfilling hard-wired(!!) requirements on distances between these two
     * new points and between them and point a.
     * Intended usage: find two H-atom positions for a given O-atom position...
     * (the hard-wired distances OD and DH are taken from what appears to be a
     * fairly TTM3F-relaxed water monomer, with OH-distance=0.98\AA and
     * bond angle 106 degrees; "D" is the midpoint between the two H-atoms)
     * @param aCoord point1, the central new O-atom
     * @param bCoord point2, a neighboring O-atom
     * @param cCoord point3, another neighboring O-atom
     * @param newHH two sets of 3D cartesian coordinates for the H-atoms
     */
    private static void giveNewHH(final double[] aCoord, final double[] bCoord,
            final double[] cCoord, final double[][] newHH){
        
        assert(newHH.length == 2);
        assert(newHH[0].length == 3);
        
        final double abDist = distance(aCoord,bCoord);
        assert(abDist > 0.0);
        
        final double[] ab = new double[3];
        ab[0] = (bCoord[0]-aCoord[0])/abDist;
        ab[1] = (bCoord[1]-aCoord[1])/abDist;
        ab[2] = (bCoord[2]-aCoord[2])/abDist;
        final double acDist = distance(aCoord,cCoord);
        final double[] ac = new double[3];
        ac[0] = (cCoord[0]-aCoord[0])/acDist;
        ac[1] = (cCoord[1]-aCoord[1])/acDist;
        ac[2] = (cCoord[2]-aCoord[2])/acDist;
        
        // construct unit vector w bisecting the ab,ac angle
        final double[] w = new double[3];
        w[0] = 0.5*(ab[0]+ac[0]);
        w[1] = 0.5*(ab[1]+ac[1]);
        w[2] = 0.5*(ab[2]+ac[2]);
        final double wNorm = Math.sqrt(TrivialLinearAlgebra.dotProduct(w,w));
        w[0] /= wNorm;
        w[1] /= wNorm;
        w[2] /= wNorm;
        
        // construct auxiliary point D on the angle bisector
        final double[] dCoord = new double[3];
        dCoord[0] = aCoord[0] + OHDIST*w[0];
        dCoord[1] = aCoord[1] + OHDIST*w[1];
        dCoord[2] = aCoord[2] + OHDIST*w[2];
        
        // construct unit vector mn perpendicular to w but in a,b,c-plane
        final double[] mn = new double[3];
        mn[0] = ac[0] - ab[0];
        mn[1] = ac[1] - ab[1];
        mn[2] = ac[2] - ab[2];
        final double mnNorm = Math.sqrt(TrivialLinearAlgebra.dotProduct(mn,mn));
        mn[0] /= mnNorm;
        mn[1] /= mnNorm;
        mn[2] /= mnNorm;
        
        // obtain new H-atom positions from w and mn
        newHH[0][0] = dCoord[0] + DHDIST * mn[0];
        newHH[0][1] = dCoord[1] + DHDIST * mn[1];
        newHH[0][2] = dCoord[2] + DHDIST * mn[2];
        newHH[1][0] = dCoord[0] - DHDIST * mn[0];
        newHH[1][1] = dCoord[1] - DHDIST * mn[1];
        newHH[1][2] = dCoord[2] - DHDIST * mn[2];
    }
    
    private static String colorToHex(final float[] colDef, final float trans){
        assert(colDef.length == 3);
        final Color col = new Color(colDef[0], colDef[1], colDef[2], trans);
        return "#" + toBrowserHexValue(col.getRed()) + toBrowserHexValue(col.getGreen()) + toBrowserHexValue(col.getBlue()) + toBrowserHexValue(col.getAlpha());
    }
    
    private static String toBrowserHexValue(int number) {
        final StringBuilder builder = new StringBuilder(Integer.toHexString(number & 0xff));
        while (builder.length() < 2) {
            builder.append("0");
        }
        return builder.toString().toUpperCase();
    }
    
    private static class GDMAdaptor extends AbstractBOBYQAMethod {
        
        private Geometry work;
        private final int leastConn;
        private final Newton newton;
        private final CartesianFullBackend back;
        private final CollisionDetectionEngine coll;
        private final double[] point;
        private final double[][] bounds;
        private final double blowCD;
        private final boolean fullyRelaxed;
        private int cIter = 0;
        
        GDMAdaptor(final Geometry work, final int leastConn,
                final Newton newton, final CollisionDetectionEngine coll,
                final double[] point, final double[][] bounds, final double blowCD,
                final boolean fullyRelaxed){
            super();
            this.coll = coll;
            this.newton = newton;
            this.back = newton.getBackend();
            assert(back != null);
            this.leastConn = leastConn;
            this.work = work;
            this.bounds = bounds;
            this.point = point;
            this.blowCD = blowCD;
            this.fullyRelaxed = fullyRelaxed;
        }

        @Override
        public double computeObjectiveValue(final double[] normalized) {
            
            cIter++;
            
            // denormalize
            denormalizeToBounds(bounds, normalized, point);
            
            // put in
            work.setExtCoordMolecule(point, leastConn);
            
            // run a quick CD
            final CartesianCoordinates cartes = work.getCartesians();
            final boolean hasColl = coll.checkOnlyForCollision(cartes, blowCD, work.getBondInfo());
            if(hasColl){
                if(VERBOSEDEBUG) {System.out.println("DEBUG: Collision found.");}
                return FixedValues.NONCONVERGEDENERGY;
            }
            
            // now try the backend
            final String[] atoms = cartes.getAllAtomTypes();
            final double[] xyz1D = cartes.getAll1DCartes();
            final int[] atsPerMol = cartes.getAllAtomsPerMol();
            final long lID = work.getID();
            final short[] atomNos = cartes.getAllAtomNumbers();
            final double[] energyParts = new double[cartes.getNoOfMolecules()];
            final int noAtoms = cartes.getNoOfAtoms();
            final short[] spins = cartes.getAllSpins();
            final float[] charges = cartes.getAllCharges();
            
            if(fullyRelaxed){
                
                work = newton.localOptimization(work);
                final double e = work.getFitness();
                if(VERBOSEDEBUG){System.out.println("DEBUG: new relaxed energy is " + e);}
                
                return e;
            } else {
                final double e =  back.energyCalculation(lID, cIter, xyz1D, atoms, atomNos, atsPerMol, energyParts, noAtoms, charges, spins, work.getBondInfo());
                if(VERBOSEDEBUG){System.out.println("DEBUG: new energy is " + e);}
                
                return e;
            }
        }

        @Override
        public boolean doesNormalize() {
            return true;
        }
        
        Geometry getCurrentGeom (){
            return work;
        }
    }
    
    private static class GDMBOBYQAConfig implements Serializable {
        
        private static final long serialVersionUID = (long) 20140525;
        
        // XXX this should ideally be configurable
        final double initialTrust = 1E-1;
        final double stoppingTrust = 1E-5;
        final int maxIter = 500;
        final int interPoints = 13; // 2N+1 w/ 6 degrees of freedom
    }
}