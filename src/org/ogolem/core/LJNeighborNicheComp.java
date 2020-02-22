/**
Copyright (c) 2013, J. M. Dieterich
              2015-2019, J. M. Dieterich and B. Hartke
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

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;
import static org.ogolem.core.Constants.BOHRTOANG;
import org.ogolem.helpers.IndexSort;

/**
 * An attempt to fully separate all 4 major LJ classes,
 * icosahedral, decahedral, tetrahedral, and fcc,
 * based on the neighborhoods of selected or all atoms
 * @author Bernd Hartke
 * @version 2019-12-29
 */
class LJNeighborNicheComp implements NicheComputer<Molecule,Geometry> {
    
    private static final long serialVersionUID = (long) 20140107;
    public static final int defWidth = 20;
    public static final String defMode = "StandardTypes";
    
    private final int width;
    private final String mode;
    
    private final boolean DEBUG = false;

    LJNeighborNicheComp(final int width, final String mode){
        this.width = width;
        this.mode = mode;
    }
    
    private LJNeighborNicheComp(final LJNeighborNicheComp orig){
        this.width = orig.width;
        this.mode = orig.mode;
    }
    
    @Override
    public LJNeighborNicheComp clone(){
        return new LJNeighborNicheComp(this);
    }
    
    /**
     * checks geometrical configurations of neighbors of several or all interior atoms,
     * and tries to assign the typical Lennard-Jones structural characters
     * "icosahedral", "decahedral", "tetrahedral" or "fcc" to the whole structure,
     * based on the occurrence patterns of 4 local neighborhood types:
     * - the actual icosahedron (occurs only in the icosahedral character!),
     * - the eclipsed version of it (the icosahedron is staggered)
     * - and the typical fcc neighborhood, which can be ABA or ABC
     *   (i.e., when viewed as layers, the 3rd layer is the same as the first or different)
     * These local neighborhoods themselves are categorized via their numbers
     * of right angles and their numbers of longest pair distances. And here things
     * start to get ugly, since we also have to deal with distortions from these 
     * ideal patterns. Therefore, the rules implemented below for extracting these
     * neighborhood configuration types from right-angle and longest-distance counts
     * do look messy and non-intuitive; they simply result from examining a lot of
     * structures occurring in many ogolem runs for LJ clusters...
     * Also, the final rules for deciding on one of the overall labels, given certain
     * counts or percentages of neighborhood configuration types, appear to be somewhat
     * less clear than they could be. This results from attempts to also include special
     * cases (like the curious ico LJ38) and a "mixed" case which apparently occurs
     * frequently (consisting of a mainly decahedral-looking cluster with an icosahedral cap)
     * @param g
     * @return the niche this geometry can be associated to 
     */
    @Override
    public Niche computeNiche(final Geometry g){
        
        final CartesianCoordinates c = g.getCartesians();
        c.moveCoordsToCOM(); // now the COM is at (d,d,d), with d <= 5.e-16 or so (checked that)
//        final double[] comCoords = c.calculateTheCOM();
//        System.out.println("center of mass coords: " + Arrays.toString(comCoords));
        final int noAtoms = c.getNoOfAtoms();
        double[][] xyz = c.getAllXYZCoordsCopy();
        final double[][] xyzSave = new double[3][noAtoms];
// find all distances from the COM
        final double[] dist2COM = new double[noAtoms];
        final int[] dist2COMIndex = new int[noAtoms];
        for(int i=0; i < noAtoms; i++){
            dist2COM[i] = Math.sqrt(xyz[0][i] * xyz[0][i] + xyz[1][i] * xyz[1][i] + xyz[2][i] * xyz[2][i]);
            dist2COMIndex[i] = i;
        }
// index sort, see remarks below...
        IndexSort.quicksort(dist2COM.clone(), dist2COMIndex);
        /*Arrays.sort(dist2COMIndex, new Comparator<Integer>() {
            @Override public int compare(final Integer i1, final Integer i2) {
            return Double.compare(dist2COM[i1], dist2COM[i2]);
            }
        });*/
        if(DEBUG){
            System.out.println("closest to COM is atom " + (dist2COMIndex[0]+1) + " with dist = " + dist2COM[dist2COMIndex[0]]*Constants.BOHRTOANG);
            System.out.println("next closest to COM is atom " + (dist2COMIndex[1]+1) + " with dist = " + dist2COM[dist2COMIndex[1]]*Constants.BOHRTOANG);
            System.out.println("farthest from COM is atom " + (dist2COMIndex[noAtoms-1]+1) + " with dist = " + dist2COM[dist2COMIndex[noAtoms-1]]*Constants.BOHRTOANG);
        }
// to arrive at a resonable estimate for a typical pair distance, calculate all distances between
// the atom closest to the COM and all others, and find the smallest one
        double guessedPairDist = Double.MAX_VALUE;
        for(int i=0; i < noAtoms; i++){
            if(i==dist2COMIndex[0]) continue;
            final double delX = xyz[0][dist2COMIndex[0]] - xyz[0][i];
            final double delY = xyz[1][dist2COMIndex[0]] - xyz[1][i];
            final double delZ = xyz[2][dist2COMIndex[0]] - xyz[2][i]; 
            double nextPairDist = Math.sqrt(delX*delX + delY*delY + delZ*delZ);            
            if(nextPairDist < guessedPairDist) guessedPairDist = nextPairDist;
        }
        if(DEBUG) System.out.println("guessed pair distance = " + guessedPairDist*Constants.BOHRTOANG);
// now focus on the atoms in the "interior", defined as the largest COM-distance
// minus fudge-factor times guessed-pair-distance
        final double fudgeFactor = 0.8; // was 1.1; this is really nasty :-(; and 1.3 does not work anymore for LJ75ico!!
        final double cutCOMdist = dist2COM[dist2COMIndex[noAtoms-1]] - fudgeFactor*guessedPairDist;
        int noCoreAtoms = -1;
        for(int i=0; i < noAtoms; i++){
            if(dist2COM[dist2COMIndex[i]] > cutCOMdist){
                noCoreAtoms = i;
                break;
            }
        }
        if(noCoreAtoms < 0){
            System.err.println("no core atoms in LJNeighborNicheComp!");
            return new Niche("ljneighborniche: no core atoms!");
        }
        if(DEBUG) System.out.println(noCoreAtoms + " atoms are in the interior, namely the following:");
// loop over the following, for all interior atoms...
        final List<Integer> bla90= new ArrayList<>(noCoreAtoms);
        final List<Integer> longestDist = new ArrayList<>(noCoreAtoms);
        int numberSkippedAtomsTooFewNeighbors = 0;
        boolean skippedAll = true;
        for(int n=0; n < noCoreAtoms; n++){
            int centralAtomIndex = dist2COMIndex[n];
            if(DEBUG) System.out.println("core atom " + (centralAtomIndex+1));
//        }
        
// calculate distances of all atoms w.r.t. this interior atom
        final double[] dist2Central = new double[noAtoms];
        int[] distIndexOrdered = new int[noAtoms];
        for(int i=0; i < noAtoms; i++){
            final double dx = xyz[0][i] - xyz[0][centralAtomIndex];
            final double dy = xyz[1][i] - xyz[1][centralAtomIndex];
            final double dz = xyz[2][i] - xyz[2][centralAtomIndex];
            dist2Central[i] = dx*dx + dy*dy + dz*dz;
            distIndexOrdered[i] = i;
        }
// direct sort is not what we want:        
        //Arrays.sort(dist2Central);
        IndexSort.quicksort(dist2Central.clone(), distIndexOrdered);
// instead, we want an index sort. Surprisingly, this is NOT already available in java!
// copied this code piece from 
// http://stackoverflow.com/questions/951848/java-array-sort-quick-way-to-get-a-sorted-list-of-indices-of-an-array
// Tested it for a simple abstract example and for a LJ_75 cluster: works!
        /*Arrays.sort(distIndexOrdered, new Comparator<Integer>() {
            @Override public int compare(final Integer i1, final Integer i2) {
            return Double.compare(dist2Central[i1], dist2Central[i2]);
            }
        });*/
// for testing purposes, print out the ordered distances and their indices
//        for(int i=0; i < noAtoms; i++){
//            System.out.println(i + " " + distIndexOrdered[i] + " " + dist2Central[distIndexOrdered[i]]);
//        }
// now find the set of nearest neighbors, by locating the "first jump" in the distances;
// fortunately, the threshold factor for that can be set within a wide range:
// between 1.02 and 3.0, for the ideal LJ75_ico geometry and for central atoms
        final double nnThreshFactor = 1.5;
        int lastNNIndex = -1;
        double dist_old = dist2Central[distIndexOrdered[1]]; // 2nd, b/c the first is the self-distance = 0
        for(int i=2; i < noAtoms; i++){
            if(dist2Central[distIndexOrdered[i]] > nnThreshFactor * dist_old){
                lastNNIndex = i-1;
                break;
            }
            dist_old = dist2Central[distIndexOrdered[i]];
        }
        if(lastNNIndex < 0){
// this may happen with fairly "wild" geometries occuring directly after a random initialization;
// for now, just try to skip this item in the core atom loop...
            if(DEBUG) System.err.println("FAILURE in LJNeighborNicheComp: lastNNIndex = " + lastNNIndex + " is negative!");
            continue;
//            return new Niche("ljneighborniche: " + 1);
//            System.exit(17);
        }
        if(DEBUG) System.out.println("Core Atom " + (centralAtomIndex+1) + " has " + lastNNIndex + " nearest neighbors ****CoreNN***");
        if(lastNNIndex != 12){ // <: incomplete shell = this is not really an "interior" atom; >: weird distorted geo, skip it!
                               // (otherwise, this may ruin the type assignment rules...)
            if(DEBUG) System.out.println("skipped NON-interior atom with " + lastNNIndex + " nearest neighbors!");
            numberSkippedAtomsTooFewNeighbors++;
            continue;
        }
        skippedAll = false; // if _every_ loop pass ended up in the "continue" statements above, 
// this will remain true, and we will test for that below...
        
// testing: print out lastNNIndex
//        System.out.println("ordered distance entry number " + lastNNIndex + " is the last nearest neighbor.");
//        System.out.println("It corresponds to atom number " + distIndexOrdered[lastNNIndex]);
// for testing, write the nearest neighbors into a file:
        if(DEBUG){
            try{
                PrintWriter outputNN = new PrintWriter("nearest_neighbors" + (centralAtomIndex+1) + ".xyz");
                final String[] sa1 = new String[lastNNIndex+1]; // weird counting, hence +1
                for(int i=0; i <= lastNNIndex; i++){ // weird counting, hence <= instead of <
                    sa1[i] = " K \t" + xyz[0][distIndexOrdered[i]]*BOHRTOANG + "\t"
                                     + xyz[1][distIndexOrdered[i]]*BOHRTOANG + "\t"
                                     + xyz[2][distIndexOrdered[i]]*BOHRTOANG;
                }
                outputNN.println(lastNNIndex+1); // again, weird counting, hence +1
                outputNN.println("nearest neighbors of atom " + (centralAtomIndex+1));
                for(final String s : sa1) {outputNN.println(s);}
                outputNN.flush();
            }catch (FileNotFoundException e){
                System.out.println("FileNotFoundException when writing nearest_neighbors.xyz");
                System.out.println(e.getMessage());
                e.printStackTrace(System.err);
            }
        }
// now, for all(!) surface atoms (i.e., all atoms other than the 1st), count the number of right angles
// in its nearest neighbors
// But first, generate connectivity info for this subset of surface atoms
// (otherwise, we would re-compute several distances several times)
        boolean[][] bonded = new boolean[lastNNIndex][lastNNIndex];
        final double cutoff = dist2Central[distIndexOrdered[lastNNIndex]]*1.5; // yes, this is a nasty fudge factor!
// use this double loop to also calculate all pair distances; problem: they then need 2b ordered,
// w/o loosing their relation to the i,j-indices... I am trying to solve that with a minimum of thought ;-)
        final int halfSize = (lastNNIndex*(lastNNIndex-1))/2;
        final double[] peripheryDists = new double[halfSize];
        final int[] peripheryIndex1 = new int[halfSize];
        final int[] peripheryIndex2 = new int[halfSize];
        int[] peripheryIndexOrdered = new int[halfSize];
        int kk = 0;
        for(int i=0; i < lastNNIndex; i++){
            bonded[i][i] = false;
            for(int j=i+1; j < lastNNIndex; j++){
                final double dx = xyz[0][distIndexOrdered[i+1]] - xyz[0][distIndexOrdered[j+1]];
                final double dy = xyz[1][distIndexOrdered[i+1]] - xyz[1][distIndexOrdered[j+1]];
                final double dz = xyz[2][distIndexOrdered[i+1]] - xyz[2][distIndexOrdered[j+1]];
                final double ddist = dx*dx+dy*dy+dz*dz;
                peripheryDists[kk] = ddist;
                peripheryIndex1[kk] = i;
                peripheryIndex2[kk] = j;
                peripheryIndexOrdered[kk] = kk;
                kk++;
                bonded[i][j] = ddist <= cutoff;
                bonded[j][i] = bonded[i][j];
            }
        }
        IndexSort.quicksort(peripheryDists.clone(), peripheryIndexOrdered);
        int countLongest = 0;
        if(DEBUG) System.out.println((halfSize-1) + ": "
                + peripheryIndex1[peripheryIndexOrdered[halfSize-1]] + " , " + peripheryIndex2[peripheryIndexOrdered[halfSize-1]] + " :  " 
                + peripheryDists[peripheryIndexOrdered[halfSize-1]]);
        for(int i=halfSize-2; i >= 0; i--){
// no abs() needed in gap, since peripheryDists are ordered ;-)            
            double gap=(peripheryDists[peripheryIndexOrdered[i+1]]-peripheryDists[peripheryIndexOrdered[i]])/peripheryDists[peripheryIndexOrdered[i]];
            if(DEBUG) System.out.println(i + ": "
                    + peripheryIndex1[peripheryIndexOrdered[i]] + " , " + peripheryIndex2[peripheryIndexOrdered[i]] + " :  " 
                    + peripheryDists[peripheryIndexOrdered[i]] + ", gap: " + gap);
            countLongest++;
            if(gap>0.05){
                break;
            }
        }
        longestDist.add(countLongest);
// now finally really loop through all surface atoms
        int count90 = 0;
        final double radtodeg = 180.0/Math.PI;
// All these loops should start at 1, not at 0, since we want to skip the 1st = central atom!!        
        for(int i=1; i < lastNNIndex; i++){
            for(int j=1; j < lastNNIndex; j++){
                if(bonded[i][j]){
                    for(int k=1; k < lastNNIndex; k++){
                        if(bonded[i][k] && j != k){ // if we also test if j and k are bonded, we _eliminate_ all right angles by that!
// we want to have "i" "in the middle" (since ij and ik are bound, but possibly not jk)                            
//                            System.out.println((j+2) + " " + (i+2) + " " + (k+2));
                            final double angle = radtodeg * CoordTranslation.calcAngle(xyz, 
                                    distIndexOrdered[j+1], distIndexOrdered[i+1], distIndexOrdered[k+1]);
//                            System.out.println("angle: " + angle);
// again, also this tolerance is somewhat arbitrary; non-central neighborhoods are somewhat distorted,
// even in the perfect lowest-energy minima, to within about 3 degrees...                            
                            if(Math.abs(angle - 90.0) < 3.0) count90++;
                        }
                    }
                }
            }
        }
        bla90.add(count90);
        if(DEBUG) System.out.println("count90 = " + count90);
        }    
// we need to trap the possibility that noCoreAtoms was zero (unusual in practice, but de facto
// it does depend on fudgeFactor, which has a somewhat wildly guessed value; so, this CAN go wrong!)
// The other way it can go wrong is that all(!) entries in the above loop get skipped,
// so we test for that, too...        
        if(noCoreAtoms < 1 || skippedAll){
// here, we always collapse weird cases into the ico niche,
// since this is the default LJ structure anyway...            
            return new Niche("ljneighborniche: " + "ico");
        }
        final int trueCoreAtoms = noCoreAtoms - numberSkippedAtomsTooFewNeighbors;
        if(DEBUG) System.out.println("true number of interior atoms: " + trueCoreAtoms + " (" + numberSkippedAtomsTooFewNeighbors + " expunged)");
        if(DEBUG) System.out.println("numbers of longest distances in all neighborhoods: " + longestDist.toString());
        if(DEBUG) System.out.println("numbers of 90deg angles in all neighborhoods: " + bla90.toString());
        final int blaMax = Collections.max(bla90);
        final int blaMin = Collections.min(bla90);
        if(DEBUG) System.out.println("min and max value of 90count: " + blaMin + " , " + blaMax);
        int[] finalIndex = new int[bla90.size()];
        double[] angle90 = new double [bla90.size()];
        for(int i=0; i < bla90.size(); i++){
            finalIndex[i] = i;
            angle90[i] = bla90.get(i);
        }
        IndexSort.quicksort(angle90, finalIndex);
        if(DEBUG){
            System.out.println("counter    #90deg-angles   #longest-dists");
            for(int i=0; i < bla90.size(); i++){
                System.out.println(i + ": " + bla90.get(finalIndex[i]) + ", " + longestDist.get(finalIndex[i]));
            }
        }
// from the info collected in bla90 and longestDist,
// deduce for each neighborhood to which type it belongs, and count these types
        int countABCfcc = 0;
        int countABAfcc = 0;
        int countEclipsed = 0;
        int countStaggered = 0;
        for(int i=0; i < bla90.size(); i++){
            final int num90 = bla90.get(finalIndex[i]);
            if(num90 == 0){
                countStaggered++;
                continue;
            }
            final int numLong = longestDist.get(finalIndex[i]);
            if(num90 == 36 || num90 == 32){
                if(num90 == 36 && numLong == 6){
                    countABCfcc++;
                }else if(numLong < 10){
                    countABAfcc++;
                }else{
                    countEclipsed++;
                }
                continue;
            }
            countEclipsed++;
        }
        if(DEBUG){
            System.out.println("ABC-fcc:   " + countABCfcc);
            System.out.println("ABA-fcc:   " + countABAfcc);
            System.out.println("eclipsed:  " + countEclipsed);
            System.out.println("staggered: " + countStaggered);
        }
// from these counts, deduce the geo-type-label for the whole structure        
        final String ljType;
        if(mode.equalsIgnoreCase("StandardTypes")){
// this convoluted logic attempts to deduce one of the standard LJ structural types (plus a "mixed" type);
// (the niching parameter "width" is simply ignored in this case)            
            if(countStaggered > 0){ // then ico or mixed
                if (countABCfcc > 0){
                    ljType = "mixed";
                }else{
                    ljType = "ico";
                }
            }else if(countABCfcc == 0 && countABAfcc == 0){ // catch the TOTALLY WEIRD 38ico: only eclipsed
                if(countEclipsed == 0){ // this should NOT happen...
                    ljType = "mixed";
                }else{
                    ljType = "ico";
                }
            }else{ // then deca, fcc, t_d or mixed
                if(countEclipsed > 0){ // then deca, t_d or mixed
                    if (countABCfcc > 0){ // then deca or mixed
                        final double percABCfcc = (double) countABCfcc / (double) trueCoreAtoms;
                        final double percEclipsed = (double) countEclipsed / (double) trueCoreAtoms;
                        if(percABCfcc > percEclipsed){
                            ljType = "deca";
                        }else{
                            ljType = "mixed";
                        }
                    }else{
                        ljType = "t_d";
                    }
                }else{
                    ljType = "fcc";
                }
            }
        }else if(mode.equalsIgnoreCase("Any")){
// this mode assumes NO prior knowledge about LJ structural types,
// instead, a niche string is built from the above 4 numbers (as percentages),
// each of them binned into bins of width "width", between 0 and 100
            final int percABCfcc = (int) Math.round((double) countABCfcc / (double) trueCoreAtoms * 100.0);
            final int percABAfcc = (int) Math.round((double) countABAfcc / (double) trueCoreAtoms * 100.0);
            final int percEclipsed = (int) Math.round((double) countEclipsed / (double) trueCoreAtoms * 100.0);
            final int percStaggered = (int) Math.round((double) countStaggered / (double) trueCoreAtoms * 100.0);
            final String partABCfcc = Integer.toString(percABCfcc/width);
            final String partABAfcc = Integer.toString(percABAfcc/width);
            final String partEclipsed = Integer.toString(percEclipsed/width);
            final String partStaggered = Integer.toString(percStaggered/width);
            ljType = partABCfcc + "-" + partABAfcc + "-" + partEclipsed + "-" + partStaggered;
        }else{
// catching other modes should have been done in GlobalConfig, but just in case...            
            System.err.println("Wrong mode " + mode + " in ljneighbor nicher; this should NOT happen!");
            ljType = "WRONG";
        }
        return new Niche("ljneighborniche: " + ljType);
    }
}