/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.properties.Property;

/**
 * This is (most likely) the central thing of the whole program. The geometry.
 * It needs to provide a lot of states, getters and setters
 * as well as methods.
 * @author Johannes Dieterich
 * @version 2016-08-29
 */
public class Geometry extends ContinuousProblem<Molecule> {
	
    /*
     * Warning concerning declaring serialVersionUID could be safely ignored since no
     * further control on the serialization process is probably needed.
     * BUT: Since a CLEAN coding is wanted, we implement it!
     * Whenever something changes in this class, increment the serialVersionUID field.
     * AND: Nasty errors when sucking in old serialized objects may follow otherwise.
     */
    private static final long serialVersionUID = (long) 20140402;

    // a boolean saying whether or not it is already locally optimized
    private boolean isLocOptimized = false;

    // number of independent particles (molecules)
    private int noOfIndieParticles;

    // Vector of Molecules
    private List<Molecule> molecules;

    private boolean coordRepresSynced = false;

    // the environment
    private Environment env;
    
    // the bonds
    private final BondInfo bonds;
	
    // properties
    private final List<Property> properties;
    
    /*
     * Constructor madness
     */
    public Geometry(final GeometryConfig gc) {
        assert(gc.noOfParticles > 0);
        this.noOfIndieParticles = gc.noOfParticles;
        this.myID = gc.lID;
        this.fatherID = gc.fatherID;
        this.motherID = gc.motherID;
        this.fitness = gc.fitness;
        this.molecules = new ArrayList<>(gc.noOfParticles);
        assert(gc.geomMCs != null);
        assert(gc.noOfParticles == gc.geomMCs.size());
        for (int i = 0; i < noOfIndieParticles; i++) {
            final MoleculeConfig mc = gc.geomMCs.get(i);
            final Molecule molecule = new Molecule(mc);
            molecules.add(molecule);
        }
        this.env = (gc.env == null) ? null : gc.env.clone();
        this.bonds = gc.bonds;
        assert(this.bonds != null);
        this.properties = new ArrayList<>();
    }

    /**
     * Extracting a geometry from a set of Cartesian coordinates.
     * @param cartes the Cartesian coordinates of this geometry
     * @param id the id
     * @param noOfMols number of molecules in it
     * @param atsPerMol atoms per molecules
     * @param molFlexies are the molecules flexible
     * @param allFlexies the explicitly optimized coordinates as a list of boolean arrays, each array representing the allowed states in that z-Matrix. List must be of right length, but entries can be null if there is no flexible coordinate in that molecule.
     * @param molConstraints are there constraints on the molecules
     * @param constrXYZ constraints as a Cartesian set
     * @param sids the string IDs of the molecules
     * @param bonds the bond information of this geometry
     */
    public Geometry(final CartesianCoordinates cartes, final long id, final int noOfMols, final int[] atsPerMol,
            final boolean[] molFlexies, final List<boolean[][]> allFlexies, final boolean[] molConstraints,
            final boolean[][] constrXYZ, final String[] sids, final BondInfo bonds) {
        assert(noOfMols > 0);
        this.noOfIndieParticles = noOfMols;
        this.myID = id;
        final GeometryConfig gc = CoordTranslation.cartesianToGeomConfig(cartes, noOfMols,
                atsPerMol, molFlexies, allFlexies, molConstraints, constrXYZ, sids, bonds);
        assert(gc.noOfParticles == noOfMols);
        this.fatherID = gc.fatherID;
        this.motherID = gc.motherID;
        this.fitness = gc.fitness;
        this.molecules = new ArrayList<>(gc.noOfParticles);
        assert(gc.geomMCs != null);
        assert(gc.noOfParticles == gc.geomMCs.size());
        assert(gc.noOfParticles == noOfIndieParticles);
        for (int i = 0; i < noOfIndieParticles; i++) {
            final MoleculeConfig mc = gc.geomMCs.get(i);
            final Molecule molecule = new Molecule(mc);
            molecules.add(molecule);
        }
        this.env = (gc.env == null) ? null : gc.env.clone();
        this.bonds = bonds;
        assert(this.bonds != null);
        this.properties = new ArrayList<>();
    }

    /**
     * A copy constructor.
     * @param orig the original geometry
     */
    public Geometry(final Geometry orig) {
        final int indies = orig.noOfIndieParticles;
        this.noOfIndieParticles = indies;
        this.myID = orig.myID;
        this.fatherID = orig.fatherID;
        this.motherID = orig.motherID;
        this.fitness = orig.fitness;
        this.molecules = new ArrayList<>(indies);
        assert(indies == orig.molecules.size());
        for (int i = 0; i < indies; i++) {
            final Molecule mol = orig.molecules.get(i);
            final Molecule molecule = new Molecule(mol);
            molecules.add(molecule);
        }
        this.env = (orig.env == null) ? null : orig.env.clone();
        assert(orig.bonds != null);
        this.bonds = orig.bonds.clone();
        this.properties = new ArrayList<>();
        for(int i = 0; i < orig.properties.size(); i++){
            final Property prop = orig.properties.get(i);
            properties.add(prop.clone());
        }
    }
    
    @Override
    public Geometry clone(){
        return new Geometry(this);
    }

    /*
     * Getters and Setters
     */
    public int getNumberOfIndieParticles() {
        return noOfIndieParticles;
    }
    
    public Molecule getMoleculeAtPosition(final int mol) {
        assert(mol >= 0);
        assert(mol < noOfIndieParticles);
        assert(molecules.size() == noOfIndieParticles);
        return molecules.get(mol);
    }
    
    @Override
    public long getFatherID(){
        return fatherID;
    }
    
    @Override
    public long getMotherID(){
        return motherID;
    }

    @Override
    public long getID() {
        return myID;
    }
    
    public boolean isLocalOptimized() {
        return isLocOptimized;
    }

    boolean areCoordsSynced() {
        return coordRepresSynced;
    }

    public boolean containsEnvironment(){
        return !(env == null);
    }

    @Override
    public void setFitness(final double fitn) {
        fitness = fitn;
    }

    public void setLocalOptimized(final boolean isOpt) {
        isLocOptimized = isOpt;
    }

    void setMoleculeAtPosition(final int whichMol, final Molecule molecule) {
        assert(whichMol >= 0);
        assert(whichMol < noOfIndieParticles);
        assert(noOfIndieParticles == molecules.size());
        molecules.set(whichMol, molecule);
    }

    @Override
    public double getFitness() {
        return fitness;
    }

    List<Molecule> getMolecules() {
        return molecules;
    }

    public CartesianCoordinates getCartesians() {

        final CartesianCoordinates cartesians = CoordTranslation.geometryToCartesian(this, false);
        cartesians.setEnergy(fitness);

        return cartesians;
    }

    CartesianCoordinates getCartesiansWithEnvironment(){

        final CartesianCoordinates cartesians = CoordTranslation.geometryToCartesian(this, true);
        cartesians.setEnergy(fitness);

        return cartesians;
    }

    Environment getEnvironmentCopy(){
        return env.clone();
    }

    Environment getEnvironment(){
        return env;
    }
    
    public BondInfo getBondInfo(){
        assert(bonds != null);
        return this.bonds;
    }

    void setEnvironment(final Environment environment){
        env = environment;
    }

    void setMolecules(final List<Molecule> mols) {
        assert(mols != null);
        assert(mols.size() == noOfIndieParticles);
        this.molecules = mols;
    }

    @Override
    public void setID(final long id) {
        myID = id;
    }

    public void setFather(final long id) {
        fatherID = id;
    }

    public void setMother(final long id) {
        motherID = id;
    }

    void setCoordsSynced(final boolean b) {
        coordRepresSynced = b;
    }

    /*
     * Methods
     */
    public int getNumberOfAtoms(){
        
        int noOfAtoms = 0;
        for(int part = 0; part < noOfIndieParticles; part++){
            final Molecule mol = molecules.get(part);
            noOfAtoms += mol.getNumberOfAtoms();
        }

        return noOfAtoms;
    }
    
    Molecule removeMolecule(final int id){
        assert(id >= 0);
        assert(id < noOfIndieParticles);
        assert(molecules.size() == noOfIndieParticles);
        this.noOfIndieParticles--;
        
        return molecules.remove(id);
    }
    
    void addMolecule(final Molecule mol, final int pos){
        assert(pos >= 0);
        this.noOfIndieParticles++;
        molecules.add(pos,mol);
        assert(molecules.size() == noOfIndieParticles);
    }

    /**
     * Sets random COMs for all molecules.
     * @param cellSize The cell in which the molecules are allowed to be.
     * @param whichCollisionDetection which CD
     * @param whichDissDetect which DD
     * @param blowColl blow factor for CD
     * @param blowDiss blow factor for DD
     */
    void setRandomCOMOfMolecules(final double[] cellSize, final CollisionDetection.CDTYPE whichCollisionDetection,
            final DissociationDetection.DDTYPE whichDissDetect, final double blowColl,
            final double blowDiss) {
        
        int count = 0;
        while (true) {
            for (final Molecule mol : molecules) {
                mol.setRandomCOM(cellSize);
            }

            final CartesianCoordinates cartes = getCartesians();
            final CollisionDetection collDetect = new CollisionDetection(whichCollisionDetection);
            final CollisionInfo collInfo = collDetect.checkForCollision(cartes, blowColl, bonds);
            
            if (!collInfo.hasCollision()) {
                // collision detection without results

                // now checking for dissociation
                final boolean diss = DissociationDetection.checkForDissociation(collInfo.getPairWiseDistances(),
                        cartes.getAllAtomTypes(), cartes.getAllAtomNumbers(), blowDiss, whichDissDetect);
                if (!diss) {
                    // no dissociation, we break the loop here
                    coordRepresSynced = true;
                    break;
                } // else: dissociation detected, on we go, stays unchanged
            } // else: stays unchanged, collision detected. -> print a line?
            
            count++;
            if(count >= FixedValues.MAXTOEMERGENCY){
                System.err.println("WARNING: MAXTOEMERGENCY reached in setRandomCOMOfMolecules. Breaking out.");
                break;
            }
        }
    }

    /**
     * Initializes the Molecules (orientation) in a single attempt.
     * Does not check whether or not the chosen orientations cause clashes and/or
     * dissociation.
     */
    void setRandomOrientationOfMolecules() {
        
        for (final Molecule mol : molecules) {
            mol.setRandomOrient();
        }
    }

    /**
     * Initializes the Molecules orientation and considers
     * a maximum amount of trials.
     * @param whichCollDetect
     * @param whichDissDetect
     * @param blowDissDetect
     * @param emergency
     * @throws Exception an exception in case that no solution could be found within the number of steps
     */
    void setRandomOrientationOfMolecules(final CollisionDetection.CDTYPE whichCollDetect, final DissociationDetection.DDTYPE whichDissDetect,
            final double blowDissDetect, final double blowBonds, final int emergency) throws Exception{

        boolean repeat = true;
        int counter = 0;
        while (repeat) {


            for(final Molecule mol : molecules){
                mol.setRandomOrient();
            }
            
            // check sanity
            final CartesianCoordinates cartes = getCartesians();
            boolean sanity = GeometrySanityCheck.checkSanity(cartes, bonds, blowBonds);
            if(!sanity){
                continue;
            }

            final CollisionDetection collDetect = new CollisionDetection(whichCollDetect);
            final CollisionInfo collInfo = collDetect.checkForCollision(cartes, blowBonds, bonds);
            if (collInfo.hasCollision() == false) {
                // collision detection without results

                // checking for dissociation
                final boolean diss = DissociationDetection.checkForDissociation(collInfo.getPairWiseDistances(),
                        cartes.getAllAtomTypes(), cartes.getAllAtomNumbers(), blowDissDetect, whichDissDetect);
                if (!diss){
                    // no dissociation, we break the loop here
                    coordRepresSynced = true;
                    break;
                } // else: dissociation detected, on we go, bRepeat stays unchanged
            } // repeat stays unchanged, collision or dissociation detected
            counter++;
            if(counter >= emergency) {throw new Exception("No possible solution found w/o hitting emergency.");}
        }
    }

    public void initializeEnvironment(){
        if(env != null){
            env.initializeConnections(getCartesians());
        }
    }

    public boolean doesFitWithEnvironment(){
        if(env != null){
            return env.doesItFit(getCartesians());
        } else {
            return true;
        }
    }

    public void setRandomExtCoordMolecule(final CollisionDetection.CDTYPE whichCollisionDetection, 
            final DissociationDetection.DDTYPE whichDissociationDetection,
            double[] cellSize, double blowDissDetect, double blowCollDetect,
            boolean bAnyDissDetect, float explicitDoFRandRatio, boolean molecularCD) {

        // TODO MAKE TUNABLE!!
        boolean beAggressive = true;

        final CollisionDetection collDetect = new CollisionDetection(whichCollisionDetection);
        final Random random = new Random();
        final ArrayList<Molecule> alBackupMoles = new ArrayList<>(noOfIndieParticles);

        // initialize it
        for(final Molecule mol : molecules) {alBackupMoles.add(new Molecule(mol));}

        boolean repeat = true;
        while (repeat) {
            int i = 0;
            for(final Molecule mol : molecules){
                /*
                 * due to final fields in the Molecule class, one needs to create a fresh
                 * Molecule every time.
                 */
                final Molecule molecule = new Molecule(mol);
                if (molecule.getFlexy()) {
                    // try to mutate
                    ZMatrix zmat = molecule.getZMatrix();

                    final boolean[][] baDoFs = molecule.getDegreesOfFreedom();
                    final BondInfo tempBonds = CoordTranslation.checkForBonds(molecule.getCartesians(), blowCollDetect);

                    int iHowManyDoFs = 0;
                    for (int k = 0; k < baDoFs.length; k++) {
                        for (int j = 0; j < baDoFs[0].length; j++) {
                            if (baDoFs[k][j]) {iHowManyDoFs++;}
                        }
                    }

                    if (iHowManyDoFs == 0) {
                        /*
                         * we do nothing since this case can happen. the user
                         * might, e.g., specify a molecule as being flexible
                         * but forget to say which degree of freedom is used.
                         */
                        System.err.println("WARNING: The total number of DoF is "
                                + "zero but the molecule is flexy. This makes no "
                                + "sense.");
                    } else {

                        // first drag a random, how many DoFs should be manipulated
                        final int howMany = random.nextInt((int)explicitDoFRandRatio*iHowManyDoFs)+1;

                        // now drag a series of randoms which coordinates these are
                        final ArrayList<Integer> whichDoFs = new ArrayList<>(howMany);
                        for(int j = 0; j < howMany; j++){
                            int em = 0;
                            do{
                                final int r = random.nextInt(iHowManyDoFs)+1;
                                if(!whichDoFs.contains(r)){
                                    whichDoFs.add(r);
                                    break;
                                }
                                em++;
                            } while(em < FixedValues.MAXTOEMERGENCY);
                            if(em >= FixedValues.MAXTOEMERGENCY) {System.err.println("WARNING: Emergency detected in setRandomExtCoordMolecule()!");}
                        }

                        // now translate back to what these mean
                        for(int iWhichCoord : whichDoFs){

                            // now translate back to what this means
                            int temp = 0;
                            int iWhichKind = 2;
                            int iWhichValue = 0;

                            OuterLoop:
                            for (int k = 0; k < baDoFs.length; k++) {
                                for (int j = 0; j < baDoFs[0].length; j++) {
                                    if (baDoFs[k][j]) {temp++;}

                                    if (temp >= iWhichCoord) {
                                        iWhichKind = j;
                                        iWhichValue = k;
                                        break OuterLoop;
                                    }
                                }
                            }

                            int emerg = 0;
                            double backup = 0.0;
                            switch(iWhichKind){
                                case 0: backup = zmat.getABondLength(iWhichValue); break;
                                case 1: backup = zmat.getABondAngle(iWhichValue); break;
                                case 2: backup = zmat.getADihedral(iWhichValue); break;
                                default:
                                    System.err.println("ERROR: This is a bug in geom init. Contact author(s)!");
                                    break;
                            }

                            SingleDOFLoop:
                            while(emerg < FixedValues.MAXTOEMERGENCY){
                                emerg++;

                                // reset the "previous state"
                                switch (iWhichKind) {
                                    case 0: zmat.setABondLength(iWhichValue, backup); break;
                                    case 1: zmat.setABondAngle(iWhichValue, backup); break;
                                    case 2: zmat.setADihedral(iWhichValue, backup); break;
                                    default:
                                        System.err.println("ERROR: This is a bug in geom init (2). Contact author(s)!");
                                        break;
                                }

                                double dFactor = random.nextDouble();

                                if(beAggressive){
                                    if (iWhichKind == 0) {zmat.setABondLength(iWhichValue, 2.5 * dFactor);}
                                    else if (iWhichKind == 0) {zmat.setABondAngle(iWhichValue, Math.PI * dFactor);}
                                    else if (iWhichKind == 2){
                                        final boolean b = random.nextBoolean();
                                        if(b) {zmat.setADihedral(iWhichValue, Math.PI*dFactor);}
                                        else {zmat.setADihedral(iWhichValue, -Math.PI*dFactor);}
                                    }
                                } else{
                                    final double d = random.nextDouble();
                                    double val = 0.0;
                                    if(iWhichKind == 0){
                                        //XXX this might need some more thinking and tinkering...
                                        val = zmat.getABondLength(iWhichValue) * d;
                                    } else if (iWhichKind == 1) {
                                        val = Math.PI * d;
                                    } else if (iWhichKind == 2) {
                                        final boolean b = random.nextBoolean();
                                        val = (b) ? Math.PI * d : -Math.PI * d;
                                    }

                                    zmat.setAValue(iWhichKind, iWhichValue, val);
                                }

                                if (molecularCD) {
                                    final CartesianCoordinates ct = CoordTranslation.zMatToCartesians(zmat);
                                    final boolean cd = collDetect.checkForCollision(ct, blowCollDetect, tempBonds).hasCollision();
                                    if(cd){
                                        continue SingleDOFLoop;
                                    } else{
                                        // wonderful! :-)
                                        break SingleDOFLoop;
                                    }
                                } else {
                                    // we just break the loop
                                    break SingleDOFLoop;
                                }
                            }
                        }

                        molecule.setZMatrix(zmat);
                        CartesianCoordinates refCartes = CoordTranslation.zMatToCartesians(zmat);
                        molecule.setReferenceCartesian(refCartes);

                    }

                }
                molecule.setRandomCOM(cellSize);
                molecule.setRandomOrient();
                molecules.set(i, molecule);
                i++;
            }
            CartesianCoordinates cartes = getCartesians();
            assert(cartes.getAllXYZCoord() != null);
            CollisionInfo collInfo = collDetect.checkForCollision(cartes, blowCollDetect, bonds);

            if (!collInfo.hasCollision()) {

                // collision detection without results

                if(bAnyDissDetect){
                // checking for dissociation
                boolean bDiss = DissociationDetection.checkForDissociation(collInfo.getPairWiseDistances(),
                        cartes.getAllAtomTypes(), cartes.getAllAtomNumbers(), blowDissDetect, whichDissociationDetection);
                    if (bDiss) {
                        // dissociation detected, on we go, bRepeat stays unchanged
                        molecules = alBackupMoles;
                    } else {
                        // no dissociation, we break the loop here
                        coordRepresSynced = true;
                        break;
                    }
                } else {
                    // we just break the loop w/o DD
                    coordRepresSynced = true;
                    break;
                }
            } else {
                molecules = alBackupMoles;
                // XXX if(DEBUG) bRepeat stays unchanged, collision detected. -> print a line?
            }
        }
    }

    /**
     * Sets the external coordinates of a molecule.
     * @param extCoords The external coordinates first the COM and then the euler angles. Needs to be 6 doubles long.
     * @param whichMol Which molecule to use.
     */
    public void setExtCoordMolecule(final double[] extCoords, final int whichMol){

        final Molecule mol = molecules.get(whichMol);
        mol.setExternalCenterOfMass(extCoords[0],extCoords[1],extCoords[2]);
        mol.setOrientation(extCoords[3],extCoords[4],extCoords[5]);
    }
    
    public double[] getCOM(final int i){
        
        assert(i < noOfIndieParticles);
        return molecules.get(i).getExternalCenterOfMass();
    }
    
    public double[] getEulers(final int i){
        
        assert(i < noOfIndieParticles);
        return molecules.get(i).getOrientation();
    }

    public double[][] getAllCOMs(){
        final double[][] all = new double[3][noOfIndieParticles];

        int counter = 0;
        double[] com;
        for(final Molecule mol : molecules){
            com = mol.getExternalCenterOfMass();
            all[0][counter] = com[0];
            all[1][counter] = com[1];
            all[2][counter] = com[2];
            counter++;
        }

        return all;
    }

    void setAllCOMs(final double[][] all){

        int counter = 0;
        for(final Molecule mol : molecules){
            final double[] com = mol.getExternalCenterOfMass();
            com[0] = all[0][counter];
            com[1] = all[1][counter];
            com[2] = all[2][counter];
            counter++;
        }
    }

    double[] getAllExtCoords(){

        final double[] all = new double[noOfIndieParticles*6];

        int counter = 0;
        double[] com;
        double[] orient;
        for(final Molecule mol : molecules){
            com = mol.getExternalCenterOfMass();
            orient = mol.getOrientation();
            all[counter    ] = com[0];
            all[counter + 1] = com[1];
            all[counter + 2] = com[2];
            all[counter + 3] = orient[0];
            all[counter + 4] = orient[1];
            all[counter + 5] = orient[2];

            counter += 6;
        }

        return all;
    }

    void setAllExtCoords(final double[] all){

        int counter = 0;
        for(final Molecule mol : molecules){
            final double[] com    = mol.getExternalCenterOfMass();
            final double[] orient = mol.getOrientation();
            com[0] = all[counter    ];
            com[1] = all[counter + 1];
            com[2] = all[counter + 2];
            orient[0] = all[counter + 3];
            orient[1] = all[counter + 4];
            orient[2] = all[counter + 5];

            counter += 6;
        }
    }
    
    /**
     * Determines a cell enclosing a randomized version of this geometry in a
     * rather crude and approximative way.
     * @param cellSize The cell size. Must be a 3D vector, overwritten at exit.
     */
    public void determineCellSize(final double[] cellSize, final double blowBonds) {
        
        assert(cellSize != null);
        assert(cellSize.length == 3);
        
        int tries = 0;
        int successes = 0;

        final Geometry gTempGeom = new Geometry(this);

        boolean repeat = true;
        while (repeat) {
            for (final Molecule molecule : gTempGeom.molecules) {
                molecule.setRandomCOM(cellSize);
                molecule.setRandomOrient();
            }

            final CartesianCoordinates cartes = getCartesians();

            // ALWAYS take the simple pairwise check for this!
            final CollisionDetection collDetect = new CollisionDetection(CollisionDetection.CDTYPE.SIMPLEPAIRWISE);
            final CollisionInfo collInfo = collDetect.checkForCollision(cartes, blowBonds, bonds);
            if (!collInfo.hasCollision()) {
                // collision detection without problem
                successes++;
                if (successes >= 20) {
                    // it should at least complete 20 times successfully before we accept the cell
                    repeat = false;
                    // found the cell!
                }// else: not yet having the proper cell but we are one step closer to it.
            } else {
                if (tries > 1000) {
                    // increase cell size by two in each direction
                    cellSize[0] += 2;
                    cellSize[1] += 2;
                    cellSize[2] += 2;
                    //out.PrintMisc(sOutputPath, sToBePrinted);
                    //System.OUT.println("New cell size is: " + iaCellSize[0] +"," + iaCellSize[1] +"," + iaCellSize[2]);
                    tries = 0;
                } // else: not yet reached the maximum number of attempts
            }
            tries++;
        }
    }

    public String[] makePrintableAbsoluteCoord(final boolean withEnv) {

        if(withEnv){
            return getCartesiansWithEnvironment().createPrintableCartesians();
        } else {
            return getCartesians().createPrintableCartesians();
        }
    }

    GeometryConfig returnMyConfig() {
        GeometryConfig gc = new GeometryConfig();
        gc.fitness = fitness;
        int iNoOfIndies = molecules.size();
        gc.noOfParticles = iNoOfIndies;
        gc.fatherID = fatherID;
        gc.motherID = motherID;
        gc.lID = myID;
        ArrayList<MoleculeConfig> alMC = new ArrayList<>(iNoOfIndies);
        for (int i = 0; i < iNoOfIndies; i++) {
            MoleculeConfig mc = molecules.get(i).returnMyConfig();
            alMC.add(i, mc);
        }
        alMC.trimToSize();
        gc.geomMCs = alMC;
        if(env != null){
            gc.env = env.clone();
        } else{
            gc.env = null;
        }
        gc.bonds = bonds.clone();
        
        return gc;
    }

    void mutateACoord(final int whichMol, final int whichCoord, final double mutFactorCOM, 
            final double mutFactorEuler, final int whichWay) {
        /*
         * which coordinate should be mutated?
         * 0: x
         * 1: y
         * 2: z
         * 3: Euler1
         * 4: Euler2
         * 5: Euler3
         */
        final Molecule mTemp = molecules.get(whichMol);
        mTemp.mutateACoord(whichCoord, mutFactorCOM, mutFactorEuler, whichWay);
    }

    void mutateEnvironment(){
        env.mutateConnections(getCartesians());
    }

    public boolean[] getAllFlexies() {
        boolean[] baFlexies = new boolean[noOfIndieParticles];
        for (int i = 0; i < noOfIndieParticles; i++) {
            baFlexies[i] = molecules.get(i).getFlexy();
        }
        return baFlexies;
    }

    /**
     * Checks all molecules in this geometry for constraints.
     * @return An array of booleans which are true if the molecule has a contstraint and false otherwise.
     */
    @Deprecated
    public boolean[] getAllConstraints() {
        boolean[] baConstraints = new boolean[noOfIndieParticles];
        for (int i = 0; i < noOfIndieParticles; i++) {
            baConstraints[i] = molecules.get(i).isConstricted();
        }
        return baConstraints;
    }

    /**
     * Checks all molecules in this geometry for constraints.
     * @return An array of booleans which are true if the molecule has a constraint and false otherwise.
     */
    public boolean[] getAllConstraints(final boolean withEnv) {

        final boolean b = (containsEnvironment() && withEnv) ? true : false;
        
        boolean[] baConstraints = new boolean[noOfIndieParticles + ((b)?1:0)];
        for (int i = 0; i < noOfIndieParticles; i++) {
            baConstraints[i] = molecules.get(i).isConstricted();
        }

        if(b) {
            baConstraints[baConstraints.length-1] = env.isEnvironmentRigid();
        }

        return baConstraints;
    }

    /**
     * Returns all constraints corresponding to a complete set of cartesian coordinates of this geometry.
     * @return An array field with the constraints.
     */
    @Deprecated
    public boolean[][] getAllConstraintsXYZ(){
        boolean[][] baConstraints = new boolean[3][getNumberOfAtoms()];

        int iCounter = 0;
        for(Molecule mol : molecules){
            final int iNoOfAtoms = mol.getNumberOfAtoms();
            if(!mol.isConstricted()){
                for(int i = 0; i < iNoOfAtoms; i++){
                    baConstraints[0][iCounter] = false;
                    baConstraints[1][iCounter] = false;
                    baConstraints[2][iCounter] = false;
                    iCounter++;
                }
            } else{
                // ok, slightly different
                final boolean[][] baMoleculeConstraints =  mol.getConstraints();
                for(int i = 0; i < iNoOfAtoms; i++){
                    baConstraints[0][iCounter] = baMoleculeConstraints[0][i];
                    baConstraints[1][iCounter] = baMoleculeConstraints[1][i];
                    baConstraints[2][iCounter] = baMoleculeConstraints[2][i];
                    iCounter++;
                }
            }
        }

        return baConstraints;
    }

    /**
     * Returns all constraints corresponding to a complete set of cartesian coordinates of this geometry.
     * @return An array field with the constraints.
     */
    public boolean[][] getAllConstraintsXYZ(final boolean withEnv){
        int iEntries;
        if(withEnv && containsEnvironment()){
            iEntries = getNumberOfAtoms() + env.atomsInEnv();
        } else{
            iEntries = getNumberOfAtoms();
        }
        boolean[][] baConstraints = new boolean[3][iEntries];

        int iCounter = 0;
        for(Molecule mol : molecules){
            final int iNoOfAtoms = mol.getNumberOfAtoms();
            if(!mol.isConstricted()){
                for(int i = 0; i < iNoOfAtoms; i++){
                    baConstraints[0][iCounter] = false;
                    baConstraints[1][iCounter] = false;
                    baConstraints[2][iCounter] = false;
                    iCounter++;
                }
            } else{
                // ok, slightly different
                final boolean[][] baMoleculeConstraints =  mol.getConstraints();
                for(int i = 0; i < iNoOfAtoms; i++){
                    baConstraints[0][iCounter] = baMoleculeConstraints[0][i];
                    baConstraints[1][iCounter] = baMoleculeConstraints[1][i];
                    baConstraints[2][iCounter] = baMoleculeConstraints[2][i];
                    iCounter++;
                }
            }
        }

        if(withEnv && containsEnvironment()){
            // default allocation is false
            if(env.isEnvironmentRigid()){
                for(int i = 0; i < env.atomsInEnv(); i++){
                    baConstraints[0][iCounter] = true;
                    baConstraints[1][iCounter] = true;
                    baConstraints[2][iCounter] = true;
                    iCounter++;
                }
            } else{
                // we need it more fine-grained
                //TODO implement it
            }
        }

        return baConstraints;
    }
    
    public String[] getSIDs(){
        
        final String[] sids = new String[noOfIndieParticles];
        for(int i = 0; i < noOfIndieParticles; i++) {sids[i] = molecules.get(i).getSID();}
        
        return sids;
    }

    public boolean isThereAFlexy() {
        boolean bAnyFlexy = false;
        for (int i = 0; i < noOfIndieParticles; i++) {
            if (molecules.get(i).getFlexy()) {
                bAnyFlexy = true;
                break;
            }
        }
        return bAnyFlexy;
    }
    
    public List<boolean[][]> getExplicitDoFs() {
        
        final List<boolean[][]> allFlexies = new ArrayList<>(getNumberOfIndieParticles());
        for (int i = 0; i < noOfIndieParticles; i++) {
            if (!molecules.get(i).getFlexy()) {
                allFlexies.add(null);
            } else {
                final boolean[][] flexiesForMol = molecules.get(i).getDegreesOfFreedom();
                allFlexies.add(flexiesForMol);
            }
        }
        
        return allFlexies;
    }

    public boolean isThereAConstraint(){
        boolean bAnyConstraint = false;
        for (int i = 0; i < noOfIndieParticles; i++) {
            if (molecules.get(i).isConstricted()) {
                bAnyConstraint = true;
                break;
            }
        }
        return bAnyConstraint;
    }
    
    public double[] energyParts(){
        
        //final int no = (env == null) ? iNumberOfIndieParticles : iNumberOfIndieParticles+1;
        // XXX do we need the env?
        final double[] parts = new double[noOfIndieParticles];
        
        for(int i = 0; i < noOfIndieParticles; i++){
            parts[i] = molecules.get(i).getEnergy();
        }
        
        return parts;
    }
    
    /**
     * Gets the FULL genom.
     * @return The 1D cartesian coordinates, order x/y/z.
     */
    @Override
    public double[] getGenomeAsDouble(){
        return getCartesians().getAll1DCartes();
    }
    
    @Override
    public Molecule[] getGenomeCopy(){
        
        final Molecule[] mols = new Molecule[noOfIndieParticles];
        for(int i = 0; i < molecules.size(); i++){
            mols[i] = molecules.get(i).clone();
        }
        
        return mols;
    }
    
    @Override
    public void setGenome(final Molecule[] mols){
        
        assert(mols != null);
        assert(mols.length == noOfIndieParticles);
        for(int i = 0; i < mols.length; i++){
            molecules.set(i, mols[i]);
        }
    }

    @Override
    public void setFatherID(final long fatherID) {
        this.fatherID = fatherID;
    }

    @Override
    public void setMotherID(final long motherID) {
        this.motherID = motherID;
    }
    
    /**
     * Adds a property value object to this geometry
     * @param prop the property to be added
     * @return true if an old instance of this property got replaced, false otherwise
     */
    public boolean addProperty(final Property prop){
        
        for(int i = 0; i < properties.size(); i++){
            final Property p = properties.get(i);
            if(prop.getClass() == p.getClass()){
                // replace
                properties.set(i, prop);
                return true;
            }
        }
        
        properties.add(prop);
        return false;
    }
    
    public void resetProperties(){
        properties.clear();
    }
    
    public Iterator<Property> getPropertyIterator(){
        return properties.iterator();
    }
    
    /**
     * THIS IS A VERY DANGEROUS ROUTINE. USE WITH CARE AFTER STUDYING THE CODE W.R.T. WHAT IT DOES (AND IMPLIES!)
     * @param xyz  the xyz coordinates IN THE CORRECT ORDER AND LENGTH
     */
    void updateAllCoordinates(final double[][] xyz){
        
        int offset = 0;
        for(int particle = 0; particle < noOfIndieParticles; particle++){
            
            final Molecule mol = molecules.get(particle);
            mol.updateAllCoordinates(xyz, offset);
            
            offset += mol.getNumberOfAtoms();
        }
        
    }
}
