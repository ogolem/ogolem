/**
Copyright (c) 2012-2015, J. M. Dieterich and B. Hartke             
              2017, J. M. Dieterich and B. Hartke
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
import java.util.Arrays;
import java.util.List;
import static org.ogolem.core.Constants.BOHRTOANG;
import org.ogolem.generic.GenericMutation;
import org.ogolem.helpers.Tuple;
/**
 * Directed mutation.
 * @author Bernd Hartke
 * @author Johannes Dieterich
 * @version 2017-03-03
 */
public class DirectedMutation implements GenericMutation<Molecule,Geometry> {
    
    private static final long serialVersionUID = (long) 20140401;
    private final static boolean DEBUG = false;
    private final CollisionDetectionEngine collDetect;
    private final double blowColl;
    private final CartesianFullBackend back;
    private final double scaleFactor;
    private final double rangeModFactor;

    DirectedMutation(final CartesianFullBackend back, final CollisionDetection.CDTYPE whichCollDetect, final double blowColl,
            final double scaleFactor, final double rangeModFactor){
        this.back = back;
        this.blowColl = blowColl;
        this.scaleFactor = scaleFactor;
        this.rangeModFactor = rangeModFactor;
        this.collDetect = new CollisionDetection(whichCollDetect);
    }
    
    DirectedMutation(final DirectedMutation orig){
        this.back = orig.back.clone();
        this.blowColl = orig.blowColl;
        this.collDetect = orig.collDetect.clone();
        this.rangeModFactor = orig.rangeModFactor;
        this.scaleFactor = orig.scaleFactor;
    }
    
    @Override
    public DirectedMutation clone() {
        return new DirectedMutation(this);
    }

    @Override
    public String getMyID() {
        return "DIRECTED GEOMETRY MUTATION\n\tscale factor: " + scaleFactor
                + "\n\trange mod factor: " + rangeModFactor;
    }

    @Override
    public Geometry mutate(final Geometry geom){
            
        //please note that the backend probably will NOT support to evaluate anything but the *same* system as always.
        // so no mixing, erasing of molecules, adding of molecules,...
        
        // nonsense message just to see if we get here at all...
        if(DEBUG) {System.out.println("DEBUG: hi there, this is DirMut!");}
        
        final Geometry work = new Geometry(geom);        
        
        // OPTION: try to scale coords to reduce collision probability for inner points
        // scaleFactor = 1.2 is a reasonable default that should be passed to DirMut
        // (of course, 1.0 does no scaling, it merely is some testing nonsense;
        // 1.2 is sufficient to make some room inside my simple test cluster)
        //
        // allow for some more room in the search box, to make placement of molecules at
        // the cluster surface more likely <-> if this is not wanted, this may be skipped.
        // This simply tacks on additional room outside of the cluster edges. The amount
        // of this room is a multiple (given by rangeModFactor) of the estimated pair distance.
        // rangeModFactor = 0.5 is a reasonable default that should be passed to DirMut
        final CoordRanges cRanges = setupGrid(work, scaleFactor, rangeModFactor);
        
        // energy-per-molecule needed for finding worst molecule(s)
        // (this could be hidden inside worstIndex if is it not really needed as input for
        // the backend (see below), which I suspect... Then, worstIndex could also be
        // modularly replaced by other methods that use other criteria for "worst" ! )
        final CartesianCoordinates c = work.getCartesians();
        final double[] energyParts = new double[c.getNoOfMolecules()];
        final double e = back.energyCalculation(work.getID(), -1, c.getAll1DCartes(), c.getAllAtomTypes(), c.getAllAtomNumbers(), 
                     c.getAllAtomsPerMol(), energyParts, c.getNoOfAtoms(), c.getAllCharges(), c.getAllSpins(), work.getBondInfo(),
                     c.containedEnvType() == CartesianCoordinates.ENVTYPE.RIGID);
        if(DEBUG){System.out.println("DEBUG: e is " + e);}
        //TODO expand DirectedMutation to do its stuff for n>1 molecules.
        // for the moment, just pretend that we only want to move one molecule...
        // (to move more, worstIndex could simply return an integer array of the N worst molecules,
        // and then the remainder could be looped over N -- with only some slight adjustments
        // in a few locations, hopefully ;-) ;-) 
        // The more tricky considerations are: If the 1st molecule is moved, this should result
        // in a changed geometry, which would then effectively block the found vacancy for further
        // placement attempts (hopefully). But does it really make sense to proceed further w/o
        // re-calculating the energy contributions? And w/o doing a locopt? But if all that is done,
        // this definitely becomes too expensive for something that is done at random for rather
        // arbitrary pool members.)
        final int indexMolMove = worstIndex(work, energyParts);
        
        if(DEBUG) {System.out.println("DEBUG: index of molecule with worst Econtrib: " + indexMolMove);}
        
        // this is the molecule that will be moved to a new position below
        final Molecule mol = work.getMoleculeAtPosition(indexMolMove);
        
        // here we save the "old" position of this molecule; this is used later to avoid
        // rediscovering it in the grid search for a new position -- which would render the
        // whole DirMut rather pointless... (on the other hand, one could think about NOT doing
        // this check, which would be equivalent to assuming that there _are_ better vacancies
        // than the current/old point, with the fallback option of rediscovering the old point
        // if there are _no_ better vacancies... Not sure if this really makes sense, though :-( )
        final double[] oldPoint = mol.getExternalCenterOfMass();
        
        // it now needs to get assigned new coords&EulerAngles. 
        // so we proceed to find the best vacancy by attempting to re-place the molecule on a grid.
        
        // this vacancy search needs to cover lots of empty space around and within
        // the geometry, this is why we attempted to restrict it to as few test points as possible...
        // (if lots of space outside of the cluster is included in the above scan, it may make
        // sense to also run everything through DD (which then should be done within the above loop);
        // otherwise this may be a waste of time...)
        //
        final List<double[]> alGoodPoints = getCollFreePoints(cRanges, cRanges.avPairDist, oldPoint, mol,
                work, collDetect, blowColl, work.getBondInfo());
        
        // at this point, we have a list of potentially good points in alGoodPoints,
        // so now try to find the best one of them, by calculating their energies,
        // i.e.: loop over this point list, calling some energy calculation (and locopt?!),
        // and remember which one is best

        if(alGoodPoints.isEmpty()){
            if(DEBUG) {System.out.println("DEBUG: No fitting new position found, returning previous geometry");}
            return new Geometry(geom);
        }
        
        // now we have the most promising vacancy => move the worst molecule there:
        // are the following two lines...
        final Tuple<Double,double[]> tup = findBestPoint(back, alGoodPoints, mol, work, energyParts, work.getBondInfo());
        final double eBest = tup.getObject1();
        if(eBest > e){
            System.out.println("INFO: Best point on dirmut grid actually worse than initial. Continuing but complaining...");
            if(DEBUG){System.out.println("DEBUG: eBest " + eBest + " vs initial " + e);}
        }
        final double[] bestPoint = tup.getObject2();
        mol.setExternalCenterOfMass(bestPoint);
        work.setMoleculeAtPosition(indexMolMove, mol); // needed?
        // equivalent to the next few lines? (except for the issue of the Eulers, of course)
//        final double[] best6Point = new double[6];
//        best6Point[0] = bestPoint[0];
//        best6Point[1] = bestPoint[1];
//        best6Point[2] = bestPoint[2];
//        best6Point[3] = 0.0; // here, some random Euler needed
//        best6Point[4] = 0.0; // here, some random Euler needed
//        best6Point[5] = 0.0; // here, some random Euler needed
//        work.setExtCoordMolecule(best6Point,indexMolMove);        
        //TODO: of course, this Euler stuff needs a _real_ solution: naively, one should try placements
        // with many different orientations for each vacancy test point, but obviously this would be
        // a big effort in computer time which would not really pay off.
        // The cheapest "solution" is to ignore the Eulers (what I am probably doing here),
        // the next less cheap "solution" would be to try a _very_ limited number of different Eulers
        // (which could be a standard set or a random set).

        return work;
        // REMARK please note that this *will* be afterwards optimized (if it passes CD/DD)
        // so no previous extra optimization if possible
    }
    
    private static int worstIndex(final Geometry geom, final double[] energyParts){
        
        // returns index of molecule with worst (largest = "most positive") energy contribution to a geometry
        // (actually, this could be a specific implementation of a far more general type which would return
        // the element that is best/worst by any criterion, from any collection (or also the best/worst N entries)
        int index = -1;
        double maxE = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < energyParts.length; i++) {
            if(DEBUG) {System.out.println("DEBUG:" + i + ": " + energyParts[i]);}
            if (energyParts[i] > maxE) {
                index = i;
                maxE = energyParts[i];
            }
        }
        
        if(DEBUG) {System.out.println("DEBUG: Worst index is " + index);}
        
        return index;
    }
    
    private static void scaleMe(final double[][] c, final double factor){
        
        for (final double[] c1 : c) {
            for (int i = 0; i < c[0].length; i++) {
                c1[i] *= factor;
            }
        }
    }

    private static class CoordRanges {
        
        double xBegin, xEnd, yBegin, yEnd, zBegin, zEnd;
        double avPairDist;
        
        void modify(final double offset){
            xBegin -= offset;
            xEnd += offset;
            yBegin -= offset;
            yEnd += offset;
            zBegin -= offset;
            zEnd += offset;
        }
    }

    private static CoordRanges getMinMax(final double[][] c){
        
        final CoordRanges cRanges = new CoordRanges();
        // find ranges in x
        cRanges.xBegin = Double.POSITIVE_INFINITY;
        cRanges.xEnd = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < c[0].length; i++){
            cRanges.xBegin = Math.min(cRanges.xBegin, c[0][i]);
            cRanges.xEnd = Math.max(cRanges.xEnd, c[0][i]);
        }
        // ... in y
        cRanges.yBegin = Double.POSITIVE_INFINITY;
        cRanges.yEnd = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < c[1].length; i++){
            cRanges.yBegin = Math.min(cRanges.yBegin, c[1][i]);
            cRanges.yEnd = Math.max(cRanges.yEnd, c[1][i]);
        }
        // ... and in z
        cRanges.zBegin = Double.POSITIVE_INFINITY;
        cRanges.zEnd = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < c[2].length; i++){
            cRanges.zBegin = Math.min(cRanges.zBegin, c[2][i]);
            cRanges.zEnd = Math.max(cRanges.zEnd, c[2][i]);
        }
        if(DEBUG) {System.out.println("DEBUG: x coord range: " + cRanges.xBegin * BOHRTOANG + " , " + cRanges.xEnd * BOHRTOANG);}
        if(DEBUG) {System.out.println("DEBUG: y coord range: " + cRanges.yBegin * BOHRTOANG + " , " + cRanges.yEnd * BOHRTOANG);}
        if(DEBUG) {System.out.println("DEBUG: z coord range: " + cRanges.zBegin * BOHRTOANG + " , " + cRanges.zEnd * BOHRTOANG);}
        
        return cRanges;
    }
    
    private static void estimatePairDist(double[][] xyzCoords, CoordRanges cRanges){
        
        final double cubeRootNoCOM = Math.pow(xyzCoords[0].length, 1.0/3.0);
        final double dx = (cRanges.xEnd - cRanges.xBegin); // would be pairDist estimate if divided by cubeRootNoCOM
        final double dy = (cRanges.yEnd - cRanges.yBegin);
        final double dz = (cRanges.zEnd - cRanges.zBegin);
        cRanges.avPairDist = (dx + dy + dz) / (3.0*cubeRootNoCOM);
        if(DEBUG) {System.out.println("DEBUG: avPairDist: " + cRanges.avPairDist * BOHRTOANG);}
    }
    
    private static List<double[]> getCollFreePoints(final CoordRanges cRanges, final double pointDist, 
            final double[] oldPoint, Molecule mol, Geometry work, final CollisionDetectionEngine collDetect, 
            final double blowColl, final BondInfo bonds){
        
        final List<double[]> alGoodPoints = new ArrayList<>();
        // number of points to cover in each direction (enlargement factor at the end is a wild guess...)
        final int iPointsX = (int) Math.round(Math.abs(cRanges.xEnd - cRanges.xBegin) / pointDist * 2.0);
        final int iPointsY = (int) Math.round(Math.abs(cRanges.yEnd - cRanges.yBegin) / pointDist * 2.0);
        final int iPointsZ = (int) Math.round(Math.abs(cRanges.zEnd - cRanges.zBegin) / pointDist * 2.0);
        
        if(DEBUG){
            System.out.println("DEBUG: points in x: " + iPointsX);
            System.out.println("DEBUG: points in y: " + iPointsY);
            System.out.println("DEBUG: points in z: " + iPointsZ);
            System.out.println("DEBUG: total number of positions to test: " + iPointsX * iPointsY * iPointsZ);
        }
        
        final double xIncr = Math.abs(cRanges.xEnd - cRanges.xBegin) / (double) iPointsX;
        final double yIncr = Math.abs(cRanges.yEnd - cRanges.yBegin) / (double) iPointsY;
        final double zIncr = Math.abs(cRanges.zEnd - cRanges.zBegin) / (double) iPointsZ;
        double testX, testY, testZ;
        double checkX, checkY, checkZ, checkSquared;
        double pointDistSquared = pointDist * pointDist;
        testX = cRanges.xBegin - xIncr;
        int countSelfColl = 0;
        for (int ix = 0; ix < iPointsX; ix++){
            testX += xIncr;
            testY = cRanges.yBegin - yIncr;
            for (int iy = 0; iy < iPointsY; iy++){
                testY += yIncr;
                testZ = cRanges.zBegin - zIncr;
                for (int iz = 0; iz < iPointsZ; iz++){
                    testZ += zIncr;
                    // test for distance w.r.t. old point
                    checkX = testX - oldPoint[0];
                    checkY = testY - oldPoint[1];
                    checkZ = testZ - oldPoint[2];
                    checkSquared = checkX * checkX + checkY * checkY + checkZ * checkZ;
                    // do following stuff only if distance w.r.t. old point is large enough,
                    // to avoid re-discovery of old point...
                    if (checkSquared > pointDistSquared){
                        // choose point w/in interval and test for collision
                        double[] point = new double[3];
                        point[0] = testX;
                        point[1] = testY;
                        point[2] = testZ;
                        if(DEBUG) {System.out.println(point[0] * BOHRTOANG + " , " + point[1] * BOHRTOANG + " , " + point[2] * BOHRTOANG);}
                        // Try to weed out the worst positioning attempts by CD:
                        mol.setExternalCenterOfMass(point); // see below: what are the Eulers here?!
                        final CartesianCoordinates cart = work.getCartesians();
                        
                        final CollisionInfo collinfo = collDetect.checkForCollision(cart, blowColl, bonds);

                        if(!collinfo.hasCollision()){
                            if (DEBUG) {System.out.println("DEBUG: added nonColl point " + Arrays.toString(point));}
                            alGoodPoints.add(point);
                        }
                    } else {
                        countSelfColl += 1;
                        if (DEBUG) {System.out.println("DEBUG: Self-collision #" + countSelfColl + ": " 
                                           + Math.sqrt(checkSquared) * BOHRTOANG + " <= "
                                           + Math.sqrt(pointDistSquared) * BOHRTOANG);}
                    }
                }
            }
        }
        
        if(DEBUG) {System.out.println("DEBUG: number of collision-free positions: " + alGoodPoints.size());}
        
        return alGoodPoints;
    }
    
    private static Tuple<Double,double[]> findBestPoint(final CartesianFullBackend back, final List<double[]> alGoodPoints, Molecule mol,
            Geometry work, final double[] energyParts, final BondInfo bonds){
        
        if(DEBUG){
            String methodID = back.getMethodID();
            System.out.println("DEBUG: Method in DirectedMutation " + methodID);
        }// just for debugging, to verify backend ID
        
        final int numberOfAtoms = work.getNumberOfAtoms();
        double[] bestPoint = new double[3];
        double eBest = Double.POSITIVE_INFINITY;
        int iter = 0;
        for (final double[] aPoint : alGoodPoints){
            if(DEBUG) {
                System.out.println("DEBUG: Testing point " + java.util.Arrays.toString(aPoint));
                System.out.println("  XX " + aPoint[0]*Constants.BOHRTOANG + "  " + aPoint[1]*Constants.BOHRTOANG + "  " + aPoint[2]*Constants.BOHRTOANG);
            }
            mol.setExternalCenterOfMass(aPoint);
            final CartesianCoordinates c = work.getCartesians();
            final double[] xyz1D = c.getAll1DCartes();
            final double e = back.energyCalculation(work.getID(), iter, xyz1D, c.getAllAtomTypes(), c.getAllAtomNumbers(), 
                     c.getAllAtomsPerMol(), energyParts, numberOfAtoms, c.getAllCharges(), c.getAllSpins(), bonds, c.containedEnvType() == CartesianCoordinates.ENVTYPE.RIGID);
            if(DEBUG) {System.out.println("DEBUG:    has energy: " + e);}
            if (e < eBest){
                eBest = e;
                bestPoint = aPoint;
            }
            iter++;
        }
        if(DEBUG) {System.out.println("DEBUG: best point was: " + Arrays.toString(bestPoint));}
        if(DEBUG) {System.out.println("DEBUG:    with energy: " + eBest);}
        
        return new Tuple<>(eBest,bestPoint);
    }
    
    /** the stuff below has not really anything to do with DirectedMutation;
        it is all about a cartesian coordinate grid "around" the cluster,
        so maybe it is more logical to place this in connection with the
        CartesianCoordinates or with Geometry or ... ? */
    private static CoordRanges setupGrid(final Geometry geom, final double scaleFactor,
            final double rangeModFactor){
        
        // New place(s) for molecules are found on a grid,
        // which in turn is based on the COM cartesians, so get them:
        final double[][] xyzCoords = geom.getAllCOMs();
        // OPTION: try to scale coords to reduce collision probability for inner points
        if (Math.abs(scaleFactor-1.0) > 1.0e-5){
            scaleMe(xyzCoords, scaleFactor);            
            geom.setAllCOMs(xyzCoords); 
        }
        // trying to pack xBegin,xEnd,... into an inner class cRanges
        // (which is a lame excuse for not using a more general and/or prettier construction)
        final CoordRanges cRanges = getMinMax(xyzCoords);
        
        // try to estimate COM-COM distances w/o recourse to info on molecules
        // (could be replaced by direct info on average pair distances, if available in some way)
        estimatePairDist(xyzCoords, cRanges);
        
        // allow for some more room in the search box, to make placement of molecules at
        // the cluster surface more likely <-> if this is not wanted, this may be skipped
        if(rangeModFactor > 0.0) {
            final double rangeMod = rangeModFactor * cRanges.avPairDist; // should be passed to DirMut
            cRanges.modify(rangeMod);
        }
        
        return cRanges;
    }
}

