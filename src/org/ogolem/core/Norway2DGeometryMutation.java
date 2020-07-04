/**
Copyright (c) 2014, J. M. Dieterich
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

import java.util.ArrayList;
import java.util.List;
import static org.ogolem.core.GlobOptAtomics.calculateMolecularSizes;
import org.ogolem.generic.GenericMutation;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * This is a an "experienced guess" algorithm. Since we
 * know from the 3D knappsack problem that one should arrange "big" molecules
 * first to have an optimal packing, we do exactly this.
 * @author Johannes Dieterich
 * @version 2016-12-18
 */
public class Norway2DGeometryMutation implements GenericMutation<Molecule,Geometry> {

    private static final long serialVersionUID = (long) 20160402;
    private static final boolean DEBUG = false;
    
    private static final int FAILEDATTEMPTSTOINCR = 500;
    private static final double INCRBOHR = 7;
    private static final int TRIESBEFORERESET = 10000;
    private static final int MAXRESETS = 5;
    
    public static enum MUTMODE{ASCENDING,RANDOM,BYSIZE};
    
    private final Lottery random = Lottery.getInstance();
    private final CollisionDetection colldetect;
    private final double blowColl;
    private final double blowDiss;
    private final DissociationDetection.DDTYPE whichDissDetect;
    private final MUTMODE mode;
    
    //TODO explicit DoF are not initialized
    
    Norway2DGeometryMutation(final CollisionDetection.CDTYPE whichCollDetect, final double blowColl, final double blowDiss,
            final DissociationDetection.DDTYPE whichDissDetect, final MUTMODE mode){
        this.colldetect = new CollisionDetection(whichCollDetect);
        this.blowColl = blowColl;
        this.blowDiss = blowDiss;
        this.whichDissDetect = whichDissDetect;
        this.mode = mode;
    }
    
    Norway2DGeometryMutation(final Norway2DGeometryMutation orig){
        this.blowColl = orig.blowColl;
        this.blowDiss = orig.blowDiss;
        this.colldetect = orig.colldetect.clone();
        this.whichDissDetect = orig.whichDissDetect;
        this.mode = orig.mode;
    }
    
    @Override
    public Norway2DGeometryMutation clone() {
        return new Norway2DGeometryMutation(this);
    }

    @Override
    public String getMyID() {
        return "NORWAY 2D PACKING MUTATION";
    }

    @Override
    public Geometry mutate(final Geometry orig) {
        
        final List<Integer> order;
        switch(mode){
            case ASCENDING:
                // ascending order
                order = new ArrayList<>();
                for(int i = 0; i < orig.getNumberOfIndieParticles(); i++){order.add(i);}
                break;
            case RANDOM:
                // random order
                order = new ArrayList<>();
                RandomUtils.randomList(orig.getNumberOfIndieParticles(), order);
                break;
            case BYSIZE:
                // by size
                order = calculateMolecularSizes(orig);
                break;
            default:
                throw new RuntimeException("Invalid mode " + mode + " in norway!");
        }

        /*
         * now we start one by one with the molecular IDs as given and put the
         * clusters together.
         */
        Geometry mutated = new Geometry(orig);
        final GeometryConfig gc = mutated.returnMyConfig();
        
        // copy the molecular configs (we will need them later)
        final List<MoleculeConfig> allMCs = new ArrayList<>();
        for(int i = 0; i < mutated.getNumberOfIndieParticles(); i++){
            allMCs.add(gc.geomMCs.get(i).clone());
        }
        
        // now empty all
        gc.geomMCs.clear();
        
        // for housekeeping...
        final List<Integer> whereIsWhat = new ArrayList<>();
        
        // put the first molecule into the center of the new cluster
        whereIsWhat.add(order.get(0));
        gc.geomMCs.add(allMCs.get(order.get(0)));
        final double[] com = {0.0,0.0,0.0};
        gc.geomMCs.get(0).externalCOM = com;
        RandomUtils.randomEulers(gc.geomMCs.get(0).externalOrient);
        
        // make the geometry config one molecule long
        gc.noOfParticles = 1;
        
        assert(gc.geomMCs.size() == 1); // just to make sure...

        // boolean for emergency exit
        boolean hasEmergency = false;

        // check the box size
        double[] cellSize = returnInitialBoxSize(mutated.getCartesians());

        // get the collision info object
        //XXX although of course only a single evaluation makes sense, this is somewhat unsafe...
        final CollisionInfo collInfo = new SingleCollisionInfo();
        
        /*
         * loop over all molecules
         */
        final double[] rndCOM = new double[2];
        int noAtoms = gc.geomMCs.get(0).noOfAtoms;
        int[] atsPerMol = new int[]{noAtoms};
        for(int molCounter = 1; molCounter < orig.getNumberOfIndieParticles(); molCounter++){
            
            gc.noOfParticles++;
            final int wherePlaced = place(whereIsWhat, allMCs, order.get(molCounter), gc);
            
            assert(gc.noOfParticles == gc.geomMCs.size());

            final int atsThisMol = gc.geomMCs.get(wherePlaced).noOfAtoms;
            noAtoms += atsThisMol;
            
            int jumpIndex = 0;
            final int[] newAtsPerMol = new int[molCounter+1];
            for(int x = 0; x < wherePlaced; x++){
                newAtsPerMol[x] = atsPerMol[x];
                jumpIndex += atsPerMol[x];
            }
            newAtsPerMol[wherePlaced] = atsThisMol;
            for(int x = wherePlaced+1; x < molCounter; x++){
                newAtsPerMol[x] = atsPerMol[x-1];
            }
            atsPerMol = newAtsPerMol;
            
            // put it somewhere in space with some random orientation
            int countFailedAtt = 0;
            int countTotalFailed = 0;
            int countResets = 0;
            int countIters = 0;
            
            // prepare a cartesian as cache, add the molecular data that does not change by altering the Euler/COM
            final CartesianCoordinates cartesCache =  new CartesianCoordinates(noAtoms, molCounter+1, atsPerMol);
            final double[][] xyz = cartesCache.getAllXYZCoord();
            final String[] atoms = cartesCache.getAllAtomTypes();
            int off = 0;
            for(int x = 0; x < molCounter+1; x++){
                if(x == wherePlaced){
                    final MoleculeConfig mc = gc.geomMCs.get(x);
                    final int noMolAts = mc.noOfAtoms;
                    System.arraycopy(mc.atomTypes, 0, atoms, off, noMolAts);
                    off += noMolAts;
                    continue;
                }
                final Molecule mol = new Molecule(gc.geomMCs.get(x));
                final double[][] molXYZ = mol.giveRotTransCartesians();
                final String[] molAtoms = mol.getAtomTypes();
                final int noMolAts = mol.getNumberOfAtoms();
                System.arraycopy(molXYZ[0], 0, xyz[0], off, noMolAts);
                System.arraycopy(molXYZ[1], 0, xyz[1], off, noMolAts);
                System.arraycopy(molXYZ[2], 0, xyz[2], off, noMolAts);
                
                System.arraycopy(molAtoms, 0, atoms, off, noMolAts);
                
                off += noMolAts;
            }
            
            // to make sure as we cannot be sure that the MCs contain the atom numbers yet
            cartesCache.recalcAtomNumbersForced();
            
            collInfo.resizeDistsAndClearState(noAtoms);
            
            if(DEBUG){
                final List<String> bondInfo = mutated.getBondInfo().translateToInput();
                System.out.println("BondInfo coming");
                for(final String s : bondInfo){
                    System.out.println(s);
                }
            }
            
            final double[][] rotResult = new double[3][atsThisMol];
            while(!hasEmergency){
                countIters++;
                
                if(countTotalFailed >= TRIESBEFORERESET){
                    System.err.println("WARNING: Needing too many iterations to pack this geometry " + orig.getID() + ". Resetting cell size.");
                    
                    countResets++;
                    if(countResets >= MAXRESETS){
                        System.err.println("WARNING: Reset of cell size has taken place " + MAXRESETS +  " times now. Returning geometry " + orig.getID() + " as is.");
                        break;
                    }
                    
                    // task: figure out the cell size w/o considering the location of the last added molecule
                    // also, we KNOW that there is one molecule close to or in 0/0/0
                    // so, set the coordinates of the last mol in our cache cartesian to that value
                    for(int coord = 0; coord < 3; coord++){
                        for(int at = jumpIndex; at < jumpIndex+atsThisMol; at++){
                            xyz[coord][at] = 0.0;
                        }
                    }
                    
                    // we also know, that below routine does not actually do anything but looping over the coordinates
                    cellSize = returnInitialBoxSize(cartesCache);
                    
                    // reset counters
                    countFailedAtt = 0;
                    countTotalFailed = 0;
                }
                
                // set a random orientation
                RandomUtils.randomEulers(gc.geomMCs.get(wherePlaced).externalOrient);
                
                for(int j = 0; j < 2; j++){
                    rndCOM[j] = random.nextDouble()*cellSize[j];
                    // get the sign
                    if(random.nextBoolean()){
                        // minus
                        rndCOM[j] = -rndCOM[j];
                    }
                }
                
                if(DEBUG){
                    System.out.println("DEBUG: x-y COM for " + wherePlaced + " is " + rndCOM[0] + ", " + rndCOM[1]);
                }
                    

                // set the COM in
                final double[] extCOM = gc.geomMCs.get(wherePlaced).externalCOM;
                extCOM[0] = rndCOM[0];
                extCOM[1] = rndCOM[1];
                extCOM[2] = 0.0; /* KEEP IT IN THE x-y-PLANE*/
                
                // get the cartesian of this mol and put it into the cache object
                final Molecule mol = new Molecule(gc.geomMCs.get(wherePlaced));
                mol.giveRotTransCartesians(rotResult);
                final int noMolAts = mol.getNumberOfAtoms();
                System.arraycopy(rotResult[0], 0, xyz[0], jumpIndex, noMolAts);
                System.arraycopy(rotResult[1], 0, xyz[1], jumpIndex, noMolAts);
                System.arraycopy(rotResult[2], 0, xyz[2], jumpIndex, noMolAts);

                // collision and dissociation detection
                /*XXX we DO NOT NEED FULL ONES! WE ONLY NEED INCREMENTAL ONES!!!*/
                collInfo.resizeDistsAndClearState(noAtoms);
                colldetect.checkForCollision(cartesCache, blowColl, mutated.getBondInfo(),collInfo);
                if(!collInfo.hasCollision()){
                    if(DEBUG){System.out.println("DEBUG: NO collision detected.");}
                    // do a dissociation detection
                    final boolean hasDissociation = DissociationDetection.checkForDissociation(collInfo.getPairWiseDistances(),
                            cartesCache.getAllAtomTypes(), cartesCache.getAllAtomNumbers(), blowDiss, whichDissDetect);

                    if(!hasDissociation){
                        // we found a geometry
                        if(DEBUG){System.out.println("DEBUG: NO Dissociation detected.");}
                        countFailedAtt = 0;
                        countIters = 0;
                        break;
                    } else {
                        countTotalFailed++;
                        if(DEBUG){System.out.println("DEBUG: Dissociation detected.");}
                    }
                } else{
                    // for incrementing the cell size (if needed)
                    countFailedAtt++;
                    countTotalFailed++;
                    if(DEBUG){
                        System.out.println("DEBUG: At least one collision detected.");
                        final List<CollisionInfo.Collision> colls = collInfo.getCollisions();
                        for(final CollisionInfo.Collision coll : colls){
                            System.out.println("DEBUG: Collision detected between " + coll.getAtomOne() + " and " + coll.getAtomTwo());
                        }
                    }
                }

                
                if(countFailedAtt >= FAILEDATTEMPTSTOINCR){
                    // inflate cell
                    for (int j = 0; j < 2; j++) {
                        final double incr = random.nextDouble()*INCRBOHR;
                        cellSize[j] += incr;
                        if(DEBUG){System.out.println("DEBUG: Incrementing cell size: " + j + " " + incr + " " + FAILEDATTEMPTSTOINCR + " " + countFailedAtt);}
                    }
                    countFailedAtt = 0;
                }

                /*if(countIters >= FixedValues.MAXTOEMERGENCY){
                    System.err.println("WARNING: Emergency caused in norway.");
                    hasEmergency = true;
                    break;
                }*/
            }

            /*if(hasEmergency){
                break;
            }*/

            // we should really hand the packed geometries back...
            mutated = new Geometry(gc);
        }

        /*if(hasEmergency){
            // return a null'd set of geometries
            System.out.println("WARNING: Noticing emergency in Norway packing. Bailing out. No big problem though.");
            return orig;
        }*/
        
        if(DEBUG){
            final String[] xyz = mutated.makePrintableAbsoluteCoord(false);
            System.out.println("DEBUG: packed geometry before fitting the environment.");
            for(String s : xyz){
                System.out.println(s);
            }
        }

        if(orig.containsEnvironment()){
            // the environment
            int emerCounter = 0;
            boolean fitsToEnv = false;
            while (!fitsToEnv && emerCounter < FixedValues.MAXTOEMERGENCY) {

                // init that one
                mutated.initializeEnvironment();

                // check whether they fit together
                fitsToEnv = mutated.doesFitWithEnvironment();

                emerCounter++;
                // in case of an "emergency", the last geom will be returned and will fail to optimize, too bad ;-)
            }
        }
        
        if(DEBUG){System.out.println("DEBUG: ONE 2D NORWAY PACKING SUCCESSFULLY DONE.");}

        return mutated;
    }
    
    private double[] returnInitialBoxSize(final CartesianCoordinates cartes){

        final double[][] coords = cartes.getAllXYZCoord();
        final double[] box = new double[2];

        // find the max/min values and round them up
        for(int i = 0; i < 2; i++){
            // loop over the coordinates
            for(int j = 0; j < cartes.getNoOfAtoms(); j++){
                // loop over the atoms
                final double d = Math.abs(coords[i][j]);
                if(d > box[i]){
                    box[i] = d;
                }
            }
        }

        // add something (up to INCRBOHR bohr) on top
        for(int i = 0; i < 2; i++){
            box[i] += random.nextDouble()*INCRBOHR;
        }

        return box;
    }
    
    private static int place(final List<Integer> whereIsWhat, final List<MoleculeConfig> allMCs,
            final int mol, final GeometryConfig gc){
        
        assert(whereIsWhat.size() == gc.geomMCs.size());
        for(int i = 0; i < whereIsWhat.size(); i++){
            if(mol < whereIsWhat.get(i)){
                // place here
                gc.geomMCs.add(i, allMCs.get(mol));
                whereIsWhat.add(i, mol);
                return i;
            }
        }
        
        // if we end up here: append
        gc.geomMCs.add(allMCs.get(mol));
        whereIsWhat.add(mol);
        
        return whereIsWhat.size()-1;
    }
}
