/**
Copyright (c) 2014-2016, J. M. Dieterich and B. Hartke
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
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import static org.ogolem.core.CoordTranslation.distance;
import org.ogolem.helpers.Tuple;
import org.ogolem.math.TrivialLinearAlgebra;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;
import static org.ogolem.core.CoordTranslation.distance;

/**
 * A collection of point providers.
 * @author Johannes Dieterich
 * @version 2014-06-18
 */
public class DirMutPointProviders implements Serializable {
    
    private static final long serialVersionUID = (long) 20140528;
    private static final boolean DEBUG = false;
    
    static interface PointProvider extends Cloneable, Serializable {
        
        public PointProvider clone();
        
        String getMyID();
        
        boolean hasNextPoint();
        
        int putNextPointIn(final Geometry geom, final List<Tuple<Integer,Integer>> sortedMols, final int moveMol, final int movePartner);
        
        void reset();
    }
    
    static PointProvider parseProvider (final String configString) throws Exception {
        
        PointProvider provider = null;
        
        if(configString.startsWith("comonlygrid")){
            final String[] tokens = configString.substring(11).trim().split("/");
            for(int i = 0; i < tokens.length; i++){tokens[i] = tokens[i].trim();}
            
            double gridHalf = COMOnlyCoordGrid.DEFAULTGRIDSIZE;
            double gridIncr = COMOnlyCoordGrid.DEFAULTGRIDINCR;
            
            for(final String token : tokens){
                if(token.isEmpty()){continue;}
                if(token.startsWith("halfcomgrid=")){
                    gridHalf = Double.parseDouble(token.substring(12).trim());
                } else if(token.startsWith("comincrement=")){
                    gridIncr = Double.parseDouble(token.substring(13).trim());
                } else {
                    throw new RuntimeException("Unknown COM only option " + token);
                }
            }
            
            provider = new COMOnlyCoordGrid(gridHalf,gridIncr);
        } else if(configString.startsWith("fullgrid")){
            final String[] tokens = configString.substring(8).trim().split("/");
            for(int i = 0; i < tokens.length; i++){tokens[i] = tokens[i].trim();}
            
            double gridHalf = FullExternalCoordGrid.DEFAULTGRIDSIZE;
            double gridIncr = FullExternalCoordGrid.DEFAULTGRIDINCR;
            double eulerIncr = FullExternalCoordGrid.DEFAULTEULERINCR;
            
            for(final String token : tokens){
                if(token.isEmpty()){continue;}
                if(token.startsWith("halfcomgrid=")){
                    gridHalf = Double.parseDouble(token.substring(12).trim());
                } else if(token.startsWith("comincrement=")){
                    gridIncr = Double.parseDouble(token.substring(13).trim());
                } else if(token.startsWith("eulerincrement=")){
                    eulerIncr = Double.parseDouble(token.substring(15).trim());
                } else {
                    throw new RuntimeException("Unknown fullgrid option " + token);
                }
            }
            
            provider = new FullExternalCoordGrid(gridHalf,gridIncr,eulerIncr);
        } else if(configString.startsWith("randompoint")){
            
            final String[] tokens = configString.substring(11).trim().split("/");
            for(int i = 0; i < tokens.length; i++){tokens[i] = tokens[i].trim();}
            
            int noPoints = 1;
            double maxCOMDiff = 2.0*Constants.ANGTOBOHR;
            
            for(final String token : tokens){
                if(token.isEmpty()){continue;}
                if(token.startsWith("nopoints=")){
                    noPoints = Integer.parseInt(token.substring(9).trim());
                } else if(token.startsWith("maxcomdiff=")){
                    maxCOMDiff = Double.parseDouble(token.substring(11).trim());
                } else {
                    throw new RuntimeException("Unknown option in randompoint: " + token);
                }
            }
            
            assert(noPoints >= 0);
            
            provider = new LotteryPointProvider(noPoints,maxCOMDiff);
        } else if(configString.equalsIgnoreCase("waterspecific")){
            
            provider = new WaterSpecificProvider();
        } else if(configString.equalsIgnoreCase("comdistaverager")){
            
            provider = new AverageCOMDistProvider();
        } else if(configString.startsWith("eulerchanger")){
            
            final String[] tokens = configString.substring(11).trim().split("/");
            for(int i = 0; i < tokens.length; i++){tokens[i] = tokens[i].trim();}
            
            short mode = 2; // some is default
            for(final String token : tokens){
                if(token.isEmpty()){continue;}
                if(token.startsWith("mode=")){
                    final String modeStr = token.substring(5).trim();
                    if(modeStr.equalsIgnoreCase("one")){
                        mode = 0;
                    } else if(modeStr.equalsIgnoreCase("some")){
                        mode = 2;
                    } else if(modeStr.equalsIgnoreCase("all")){
                        mode = 1;
                    } else {
                        throw new RuntimeException("Unknown mode option " + modeStr);
                    }
                } else {
                    throw new RuntimeException("Unknown option in eulerchanger: " + token);
                }
            }
            
            provider = new EulerChanger(mode);
        } else {
            throw new RuntimeException("Unknown GDM provider: " + configString);
        }
        
        return provider;
    }
    
    static class FullExternalCoordGrid implements PointProvider {
        
        public static final double DEFAULTEULERINCR = Math.toRadians(30);
        public static final double DEFAULTGRIDINCR = 0.5*Constants.ANGTOBOHR;
        public static final double DEFAULTGRIDSIZE = 5.0*Constants.ANGTOBOHR;
        
        private static final long serialVersionUID = (long) 20140528;
        private final double gridHalfLength;
        private final double gridIncr;
        private final double eulerIncr;

        private List<double[]> trialPoints = null;
        private Iterator<double[]> trialPIter = null;
        
        FullExternalCoordGrid(final double gridHalf, final double gridIncr,
                final double eulerIncr){
            this.gridHalfLength = gridHalf;
            this.gridIncr = gridIncr;
            this.eulerIncr = eulerIncr;
        }
        
        @Override
        public PointProvider clone() {
            return new FullExternalCoordGrid(this.gridHalfLength, this.gridIncr, this.eulerIncr);
        }

        @Override
        public String getMyID() {
            return "all external coordinate grid with half length: " + gridHalfLength
                    + "\n\t and grid increment: " + gridIncr
                    + "\n\t and Euler increment: " + eulerIncr;
        }

        private List<double[]> trialPoints(final Geometry geom, final int moveMol, final int movePartner) {
            
            final double[] comSec = geom.getCOM(movePartner);
            
            final int noAts = geom.getMolecules().get(moveMol).getNumberOfAtoms();
            final int noPoints = (int) Math.ceil(2*gridHalfLength / gridIncr);
            final int noPointsPhi = (noAts == 1) ? 1 : (int) Math.floor(2 * Math.PI / eulerIncr);
            final int noPointsOmega = (noAts == 1) ? 1 : (int) Math.floor(Math.PI / eulerIncr);
            final int noPointsPsi = (noAts == 1) ? 1 : (int) Math.floor(2 * Math.PI / eulerIncr);
            final List<double[]> pointsList = new ArrayList<>(noPoints*noPointsPhi*noPointsOmega*noPointsPsi);
            for (int x = 0; x < noPoints; x++) {
                final double currX = comSec[0] - gridHalfLength + x * gridIncr;
                for (int y = 0; y < noPoints; y++) {
                    final double currY = comSec[1] - gridHalfLength + y * gridIncr;
                    for (int z = 0; z < noPoints; z++) {
                        final double currZ = comSec[2] - gridHalfLength + z * gridIncr;
                        for (int phi = 0; phi < noPointsPhi; phi++) {
                            final double currPhi = -Math.PI + phi * eulerIncr;
                            for (int omega = 0; omega < noPointsOmega; omega++) {
                                final double currOmega = -Math.PI / 2 + omega * eulerIncr;
                                for (int psi = 0; psi < noPointsPsi; psi++) {
                                    final double currPsi = -Math.PI + psi * eulerIncr;

                                    final double[] thisPoint = new double[]{currX,currY,currZ,currPhi,currOmega,currPsi};
                                    pointsList.add(thisPoint);

                                }
                            }
                        }
                    }
                }
            }
            
            return pointsList;
        }
        
        @Override
        public int putNextPointIn(final Geometry geom, final List<Tuple<Integer, Integer>> sortedMols,
                final int moveMol, final int movePartner) {
            
            if(moveMol < 0 || movePartner < 0){throw new RuntimeException("this point provider must be correctly used w.r.t. to move mol and move partner!");}
            
            if(trialPoints == null){
                // XXX very inefficient
                trialPoints = trialPoints(geom, moveMol, movePartner);
                trialPIter = trialPoints.iterator();
                if(!trialPIter.hasNext()){
                    if(DEBUG) System.err.println("DEBUG: Found NOTHING!!!");
                    return moveMol;
                }
            }
            
            final double[] point = trialPIter.next();
            geom.setExtCoordMolecule(point, moveMol);
            
            return moveMol;
        }

        @Override
        public void reset() {
            this.trialPIter = null;
            this.trialPoints = null;
        }

        @Override
        public boolean hasNextPoint() {
            if(trialPIter == null  && trialPoints == null){return true;}
            return trialPIter.hasNext();
        }
    }
    
    static class EulerChanger implements PointProvider {
        
        private static final long serialVersionUID = (long) 20140614;
        
        private final Lottery random = Lottery.getInstance();
        private final short mode;
        private List<Integer> changeMols = null;
        
        EulerChanger(final short mode){
            this.mode = mode;
        }

        EulerChanger(final EulerChanger orig){
            this.mode = orig.mode;
        }
        
        @Override
        public PointProvider clone() {
            return new EulerChanger(this);
        }

        @Override
        public String getMyID() {
            return "Euler angle changer, mode: " + mode;
        }

        @Override
        public int putNextPointIn(final Geometry geom, final List<Tuple<Integer, Integer>> sortedMols,
                final int moveMol, final int movePartner) {
            
            if(changeMols == null){
                // figure out which euler angles to change...
                changeMols = new ArrayList<>();
                if(mode == 0){
                    // one
                    final int which = random.nextInt(geom.getNumberOfIndieParticles());
                    changeMols.add(which);
                } else if(mode == 1) {
                    // all
                    for(int i = 0; i < geom.getNumberOfIndieParticles(); i++){
                        changeMols.add(i);
                    }
                } else if(mode == 2) {
                    // some
                    final int moveable = random.nextInt(geom.getNumberOfIndieParticles());
                    RandomUtils.listOfPoints(moveable, 0, geom.getNumberOfIndieParticles(), changeMols);
                } else {
                    throw new RuntimeException("Mode " + mode + " is not supported for Euler.");
                }
            }
            
            // pop next molcule off...
            final int which = changeMols.remove(0);
            // change Euler angles
            RandomUtils.randomEulers(geom.getMoleculeAtPosition(which).getOrientation());
            
            return which;
        }

        @Override
        public void reset() {
            changeMols = null;
        }

        @Override
        public boolean hasNextPoint() {
            if(changeMols == null){return true;}
            return !changeMols.isEmpty();
        }
    }
    
    static class COMOnlyCoordGrid implements PointProvider {

        public static final double DEFAULTGRIDINCR = 0.5*Constants.ANGTOBOHR;
        public static final double DEFAULTGRIDSIZE = 5.0*Constants.ANGTOBOHR;
    
        private static final long serialVersionUID = (long) 20140528;
        private final double gridHalfLength;
        private final double gridIncr;

        private List<double[]> trialPoints = null;
        private Iterator<double[]> trialPIter = null;
        
        COMOnlyCoordGrid(final double gridHalf, final double gridIncr){
            this.gridHalfLength = gridHalf;
            this.gridIncr = gridIncr;
        }
        
        @Override
        public PointProvider clone() {
            return new COMOnlyCoordGrid(this.gridHalfLength, this.gridIncr);
        }

        @Override
        public String getMyID() {
            return "COM only coordinate grid with half length: " + gridHalfLength
                    + "\n\t and grid increment: " + gridIncr;
        }

        private List<double[]> trialPoints(final Geometry geom, final int movePartner) {
            
            final double[] comSec = geom.getCOM(movePartner);
            
            final int noPoints = (int) Math.ceil(2*gridHalfLength / gridIncr);
            final List<double[]> pointsList = new ArrayList<>(noPoints*noPoints*noPoints);
            for (int x = 0; x < noPoints; x++) {
                final double currX = comSec[0] - gridHalfLength + x * gridIncr;
                for (int y = 0; y < noPoints; y++) {
                    final double currY = comSec[1] - gridHalfLength + y * gridIncr;
                    for (int z = 0; z < noPoints; z++) {
                        final double currZ = comSec[2] - gridHalfLength + z * gridIncr;
                        
                        final double[] thisPoint = new double[]{currX,currY,currZ};
                        pointsList.add(thisPoint);
                    }
                }
            }
            
            return pointsList;
        }

        @Override
        public int putNextPointIn(final Geometry geom, final List<Tuple<Integer, Integer>> sortedMols,
                final int moveMol, final int movePartner) {
            
            if(moveMol < 0 || movePartner < 0){throw new RuntimeException("this point provider must be correctly used w.r.t. to move mol and move partner!");}
            
            if(trialPoints == null){
                // XXX very inefficient
                trialPoints = trialPoints(geom, movePartner);
                trialPIter = trialPoints.iterator();
                if(!trialPIter.hasNext()){
                    if(DEBUG) System.err.println("DEBUG: Found NOTHING!!!");
                    return moveMol;
                }
            }
            
            final double[] point = trialPIter.next();
            geom.getMoleculeAtPosition(moveMol).setExternalCenterOfMass(point);
            
            return moveMol;
        }
        
        @Override
        public void reset() {
            this.trialPIter = null;
            this.trialPoints = null;
        }

        @Override
        public boolean hasNextPoint() {
            if(trialPIter == null && trialPoints == null){return true;}
            return trialPIter.hasNext();
        }
    }
    
    static class WaterSpecificProvider implements PointProvider {
        
        private static final long serialVersionUID = (long) 20140528;
        public static final double JUMPDETECTFACTOR = 3.0;
        public static final double DISTSCALEFAC = 1.042;
        
        private List<double[][]> trialPoints = null;
        private Iterator<double[][]> trialPIter = null;
        
        WaterSpecificProvider(){
        }

        @Override
        public PointProvider clone() {
            return new WaterSpecificProvider();
        }

        @Override
        public String getMyID() {
            return "Hartke's finest water specific point provider";
        }

        
        @Override
        public int putNextPointIn(Geometry geom, final List<Tuple<Integer, Integer>> sortedMols,
                final int moveMol, final int movePartner) {
            
            if(moveMol < 0 || movePartner < 0){throw new RuntimeException("this point provider must be correctly used w.r.t. to move mol and move partner!");}
            
            if(trialPoints == null){
                this.trialPoints = trialPoints(geom,movePartner);
                this.trialPIter = trialPoints.iterator();
                if(!trialPIter.hasNext()){
                    if(DEBUG) System.err.println("DEBUG: Found NOTHING!!!");
                    return moveMol;
                }
            }
            
            final double[][] thisSet = trialPIter.next();
            
            // put in
            final Molecule mol = geom.getMoleculeAtPosition(moveMol);
            final double[] com = mol.getExternalCenterOfMass();
            com[0] = 0.0; com[1] = 0.0; com[2] = 0.0;
            
            final CartesianCoordinates currC = geom.getCartesians();
            
            currC.setXYZCoordinatesOfAtom(thisSet[0],3*moveMol);
            currC.setXYZCoordinatesOfAtom(thisSet[1],3*moveMol+1);
            currC.setXYZCoordinatesOfAtom(thisSet[2],3*moveMol+2);
            
            currC.calculateMoleculeCOM(moveMol, com); // and directly it is in :-)
            
            final CartesianCoordinates cMol = mol.getCartesians();
            final double[][] xyzMol = cMol.getAllXYZCoord();
            // now put the rest in
            // the O
            xyzMol[0][0] = thisSet[0][0] - com[0];
            xyzMol[1][0] = thisSet[0][1] - com[1];
            xyzMol[2][0] = thisSet[0][2] - com[2];
            // H1
            xyzMol[0][1] = thisSet[1][0] - com[0];
            xyzMol[1][1] = thisSet[1][1] - com[1];
            xyzMol[2][1] = thisSet[1][2] - com[2];
            
            // H2
            xyzMol[0][0] = thisSet[2][0] - com[0];
            xyzMol[1][0] = thisSet[2][1] - com[1];
            xyzMol[2][0] = thisSet[2][2] - com[2];
            
            //XXX kinda nasty and error prone :-(
            
            return moveMol;
        }
        
        @Override
        public void reset() {
            this.trialPIter = null;
            this.trialPoints = null;
        }
        
        @Override
        public boolean hasNextPoint(){
            if(trialPIter == null && trialPoints == null){return true;}
            return trialPIter.hasNext();
        }
        
        private List<double[][]> trialPoints(final Geometry geom, int movePartner) {
            
            final CartesianCoordinates c = geom.getCartesians();
            
            /*
             * make sure that we are dealing with water molecules in the expected atom
             * order...
             */
            for(final Molecule m : geom.getMolecules()){
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
             * first find s.th. like a representative COM-COM distance,
             * by scanning all these distances
             */
            final int numMolecules = c.getNoOfMolecules();
            final double[] allCOMdists = new double[(int)(numMolecules*(numMolecules-1)*0.5)];
            int counter = 0;
            for(int i = 0; i < numMolecules; i++){
                for(int j = 0; j < i; j++){
                    allCOMdists[counter] = distance(geom.getCOM(i),geom.getCOM(j));
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
            if(DEBUG){
                System.out.println("DEBUG: average nearest-neighbor COM-COM distance: " + distCOMaverage);
            }
                
            // now find those molecules that are within a * distCOMaverage from
            // the second-least-connected molecule (note that a>2 makes no sense for
            // the following, and neither does a<1, but within 1<a<2 my choice of
            // a=1.5 is fairly arbitrary...
            final double factor = 1.5;
            final double distNeigh = factor * distCOMaverage;
            final List<Integer> neighborHood = new ArrayList<>();
            assert(movePartner >= 0);
            for(int i = 0; i < movePartner; i++){
                if(distance(geom.getCOM(movePartner),geom.getCOM(i)) < distNeigh){
                    neighborHood.add(i);
                }
            }
             for(int i = movePartner+1; i < numMolecules; i++){
                if(distance(geom.getCOM(movePartner),geom.getCOM(i)) < distNeigh){
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
            final double[] comSec = geom.getCOM(movePartner);
            final double[][] comTriple = new double[3][3];
            final double[][] newPos = new double[2][3];
            final double[] xyzNewO = new double[3];
            final double[][] xyzNewHH = new double[2][3];
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
            final List<double[][]> pointsList = new ArrayList<>();
            for(final int i : neighborHood){
                for (final int j : neighborHood){
                    if(i<j){ // I suspect this can be done more gracefully with foreach(list) but I do not know how...
                        if(DEBUG){
                            System.out.println("DEBUG: " + i + "\t"+ j);
                        }
                        comTriple[0] = comSec;
                        comTriple[1] = geom.getCOM(i);
                        comTriple[2] = geom.getCOM(j);
                    
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
                                
                                final double[][] thisWaterBlock = new double[3][];
                                // first: O
                                System.arraycopy(newPos[pair],0,xyzNewO,0,3);
                                thisWaterBlock[0] = xyzNewO.clone();
                                assert(!Double.isNaN(xyzNewO[0]));
                                assert(!Double.isNaN(xyzNewO[1]));
                                assert(!Double.isNaN(xyzNewO[2]));
                                // 2) place the two H-atoms
                                giveNewHH(xyzNewO,comTriple[n],comTriple[m], xyzNewHH);
                                assert(!Double.isNaN(xyzNewHH[0][0]));
                                assert(!Double.isNaN(xyzNewHH[0][1]));
                                assert(!Double.isNaN(xyzNewHH[0][2]));
                                assert(!Double.isNaN(xyzNewHH[1][0]));
                                assert(!Double.isNaN(xyzNewHH[1][1]));
                                assert(!Double.isNaN(xyzNewHH[1][2]));
                                thisWaterBlock[1] = xyzNewHH[0].clone();
                                thisWaterBlock[2] = xyzNewHH[1].clone();
                                
                                // add to list of points
                                pointsList.add(thisWaterBlock);
                            }
                        }
                    }
                }
            }
            
            return pointsList;
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
         * @return 0 if everything ok, -1 otherwise.
         */
        private static short giveNewPos(final double[] aCoord, 
            final double[] bCoord, final double[] cCoord, final double dist,
            final double[][] newPos){
        
            assert(newPos.length == 2);
            assert(newPos[0].length == 3);
        
            // check input coords for their pair distances, see comments below
            final double dLimit = 2.0 * dist;
            if(distance(aCoord,bCoord)>=dLimit){
//                 System.err.println("WARNING: GraphBasedDirMut distance a-b too large");
                return -1;
            }
            if(distance(aCoord,cCoord)>=dLimit){
//                System.err.println("WARNING: GraphBasedDirMut distance a-c too large");
                return -1;
            }
            if(distance(bCoord,cCoord)>=dLimit){
//                System.err.println("WARNING: GraphBasedDirMut distance b-c too large");
                return -1;
            }
        
            // now we are sure that a,b,c are close enough to each other,
            // so now calculate possible new positions: 
            // first find the center p of the circumscribed circle of a,b,c
            // (rCircCirc is the radius of this circle)
            final double[] ab = new double[3];
            final double[] ac = new double[3];
            //final double[] bc = new double[3];
            final double abDist = distance(aCoord, bCoord);
            ab[0] = (bCoord[0] - aCoord[0]) / abDist;
            ab[1] = (bCoord[1] - aCoord[1]) / abDist;
            ab[2] = (bCoord[2] - aCoord[2]) / abDist;
            final double acDist = distance(aCoord, cCoord);
            ac[0] = (cCoord[0] - aCoord[0]) / acDist;
            ac[1] = (cCoord[1] - aCoord[1]) / acDist;
            ac[2] = (cCoord[2] - aCoord[2]) / acDist;
            final double bcDist = distance(bCoord, cCoord);
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
//                System.err.println("WARNING: circumscribed circle too large in GDM");
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
    }
    
    static class AverageCOMDistProvider implements PointProvider {
        
        private static final long serialVersionUID = (long) 20140528;
        public static final double JUMPDETECTFACTOR = 3.0;
        public static final double DISTSCALEFAC = 1.042;
        
        private List<double[]> trialPoints = null;
        private Iterator<double[]> trialPIter = null;
        
        AverageCOMDistProvider(){
        }

        @Override
        public PointProvider clone() {
            return new AverageCOMDistProvider();
        }

        @Override
        public String getMyID() {
            return "Hartke's finest average COM distance point provider";
        }

        
        @Override
        public int putNextPointIn(final Geometry geom, final List<Tuple<Integer, Integer>> sortedMols,
                final int moveMol, final int movePartner) {
            
            if(moveMol < 0 || movePartner < 0){throw new RuntimeException("this point provider must be correctly used w.r.t. to move mol and move partner!");}
            
            if(trialPoints == null){
                this.trialPoints = trialPoints(geom,movePartner);
                this.trialPIter = trialPoints.iterator();
                if(!trialPIter.hasNext()){
                    if(DEBUG) System.err.println("DEBUG: Found NOTHING!!!");
                    return moveMol;
                }
            }
            
            final double[] thisCOM = trialPIter.next();
            
            // put in
            final Molecule mol = geom.getMoleculeAtPosition(moveMol);
            mol.setExternalCenterOfMass(thisCOM);
            
            return moveMol;
        }
        
        @Override
        public void reset() {
            this.trialPIter = null;
            this.trialPoints = null;
        }
        
        @Override
        public boolean hasNextPoint(){
            if(trialPIter == null && trialPoints == null){return true;}
            return trialPIter.hasNext();
        }
        
        private List<double[]> trialPoints(final Geometry geom, int movePartner) {
            
            final CartesianCoordinates c = geom.getCartesians();
                    
            /*
             * first find s.th. like a representative COM-COM distance,
             * by scanning all these distances
             */
            final int numMolecules = c.getNoOfMolecules();
            final double[] allCOMdists = new double[(int)(numMolecules*(numMolecules-1)*0.5)];
            int counter = 0;
            for(int i = 0; i < numMolecules; i++){
                for(int j = 0; j < i; j++){
                    allCOMdists[counter] = distance(geom.getCOM(i),geom.getCOM(j));
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
            if(DEBUG){
                System.out.println("DEBUG: average nearest-neighbor COM-COM distance: " + distCOMaverage);
            }
                
            // now find those molecules that are within a * distCOMaverage from
            // the second-least-connected molecule (note that a>2 makes no sense for
            // the following, and neither does a<1, but within 1<a<2 my choice of
            // a=1.5 is fairly arbitrary...
            final double factor = 1.5;
            final double distNeigh = factor * distCOMaverage;
            final List<Integer> neighborHood = new ArrayList<>();
            assert(movePartner >= 0);
            for(int i = 0; i < movePartner; i++){
                if(distance(geom.getCOM(movePartner),geom.getCOM(i)) < distNeigh){
                    neighborHood.add(i);
                }
            }
             for(int i = movePartner+1; i < numMolecules; i++){
                if(distance(geom.getCOM(movePartner),geom.getCOM(i)) < distNeigh){
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
            final double[] comSec = geom.getCOM(movePartner);
            final double[][] comTriple = new double[3][3];
            final double[][] newPos = new double[2][3];
            final double[] xyzNewO = new double[3];
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
            final List<double[]> pointsList = new ArrayList<>();
            for(final int i : neighborHood){
                for (final int j : neighborHood){
                    if(i<j){ // I suspect this can be done more gracefully with foreach(list) but I do not know how...
                        if(DEBUG){
                            System.out.println("DEBUG: " + i + "\t"+ j);
                        }
                        comTriple[0] = comSec;
                        comTriple[1] = geom.getCOM(i);
                        comTriple[2] = geom.getCOM(j);
                    
                        final short ret = giveNewPos(comTriple[0],comTriple[1],comTriple[2],distCOMaverage,newPos);
                        if(ret < 0){continue;}
                        // put new points in work/Coords and calc E...
                        // do this for both positions returned in newPos
                        for(int pair = 0; pair < 2; pair++){
                            // and for all pairs 12,13,23 in the comTriple
                            for(final int[] nmPair : indexPair){
                                // here smuggle in the new coords for molecule leastConnID;                                
                                System.arraycopy(newPos[pair],0,xyzNewO,0,3);
                                pointsList.add(newPos[pair].clone());
                                
                            }
                        }
                    }
                }
            }
            
            return pointsList;
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
         * @return 0 if everything ok, -1 otherwise.
         */
        private static short giveNewPos(final double[] aCoord, 
            final double[] bCoord, final double[] cCoord, final double dist,
            final double[][] newPos){
        
            assert(newPos.length == 2);
            assert(newPos[0].length == 3);
        
            // check input coords for their pair distances, see comments below
            final double dLimit = 2.0 * dist;
            if(distance(aCoord,bCoord)>=dLimit){
//                 System.err.println("WARNING: GraphBasedDirMut distance a-b too large");
                return -1;
            }
            if(distance(aCoord,cCoord)>=dLimit){
//                System.err.println("WARNING: GraphBasedDirMut distance a-c too large");
                return -1;
            }
            if(distance(bCoord,cCoord)>=dLimit){
//                System.err.println("WARNING: GraphBasedDirMut distance b-c too large");
                return -1;
            }
        
            // now we are sure that a,b,c are close enough to each other,
            // so now calculate possible new positions: 
            // first find the center p of the circumscribed circle of a,b,c
            // (rCircCirc is the radius of this circle)
            final double[] ab = new double[3];
            final double[] ac = new double[3];
            //final double[] bc = new double[3];
            final double abDist = distance(aCoord, bCoord);
            ab[0] = (bCoord[0] - aCoord[0]) / abDist;
            ab[1] = (bCoord[1] - aCoord[1]) / abDist;
            ab[2] = (bCoord[2] - aCoord[2]) / abDist;
            final double acDist = distance(aCoord, cCoord);
            ac[0] = (cCoord[0] - aCoord[0]) / acDist;
            ac[1] = (cCoord[1] - aCoord[1]) / acDist;
            ac[2] = (cCoord[2] - aCoord[2]) / acDist;
            final double bcDist = distance(bCoord, cCoord);
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
//                System.err.println("WARNING: circumscribed circle too large in GDM");
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
    }
    
    static class LotteryPointProvider implements PointProvider {
        
        private static final long serialVersionUID = (long) 20140528;
        private final int noPoints;
        private final double maxMove;
        private final Lottery r = Lottery.getInstance();
        private int counter = 0;
        
        LotteryPointProvider(final int noPoints, final double maxCOMDiff){
            this.maxMove = maxCOMDiff;
            this.noPoints = noPoints;
        }

        @Override
        public PointProvider clone() {
            return new LotteryPointProvider(noPoints,maxMove);
        }

        @Override
        public String getMyID() {
            return "random provider: no points " + noPoints + " max move " + maxMove;
        }

        @Override
        public int putNextPointIn(final Geometry geom, final List<Tuple<Integer, Integer>> sortedMols, 
                final int moveMol, final int movePartner) {
            
            if(moveMol < 0 || movePartner < 0){throw new RuntimeException("this point provider must be correctly used w.r.t. to move mol and move partner!");}
            
            final double[] com = geom.getMoleculeAtPosition(movePartner).getExternalCenterOfMass().clone();
            com[0] = (r.nextBoolean()) ? com[0] + r.nextDouble()*maxMove : com[0] - r.nextDouble()*maxMove;
            com[1] = (r.nextBoolean()) ? com[1] + r.nextDouble()*maxMove : com[1] - r.nextDouble()*maxMove;
            com[2] = (r.nextBoolean()) ? com[2] + r.nextDouble()*maxMove : com[2] - r.nextDouble()*maxMove;
            
            geom.getMoleculeAtPosition(moveMol).setExternalCenterOfMass(com);
            
            counter++;
            
            return moveMol;
        }
        
        @Override
        public boolean hasNextPoint(){
            return (counter < noPoints);
        }
        
        @Override
        public void reset() {
            this.counter = 0;
        }
    }
}
