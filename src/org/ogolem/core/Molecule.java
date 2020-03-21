/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2015, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Random;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.helpers.RandomUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This describes a full molecules consisting of atoms and further information.
 * @author Johannes Dieterich
 * @version 2020-02-09
 */
public class Molecule extends ContinuousProblem<Double> {

    // serialVersion so that serialization works flawlessly
    private static final long serialVersionUID = (long) 20200209;
    private static final Logger l = LoggerFactory.getLogger(Molecule.class);

    // String ID of the molecule
    private String sID;

    // integer ID of the molecule
    private int iID;

    // a single boolean describing whether or not the *whole* molecule is flexible or not
    private boolean isFlexy;

    // a single boolean describing whether or not the *whole* molecule is constraint or not
    private boolean isConstrained;

    // double array for the external COM
    private double[] extCOM;
    // double array describing the orientation of the molecule through Euler angles
    private double[] extOrient;

    private double[][] refXYZ;
    private String[] atomTypes;
    private short[] atomNumbers;

    // the degrees of freedom
    private boolean[][] dofs;

    // the constraints
    private boolean[][] constraints;

    // z matrix, just populated for a flexible molecule, otherwise null
    private ZMatrix zmat;

    // the number of atoms
    private int noOfAtoms;

    // the charges on the atoms
    private float[] charges;

    // the spins on the atoms
    private short[] spins;

    // kind of a private portion of the energy
    private double energy;

    /**
     * Simple constructor
     * @param mc the molecule configuration describing this molecule
     */
    Molecule(final MoleculeConfig mc) {
        assert(mc.charges != null);
        assert(mc.spins != null);
        assert(mc.charges.length == mc.spins.length);
        assert(mc.charges.length == mc.noOfAtoms);
        assert(mc.externalOrient != null);
        assert(mc.atomTypes != null);
        assert(mc.refXYZ != null);
        assert(mc.refXYZ.length == 3);
        assert(mc.refXYZ[0].length == mc.noOfAtoms);
        this.myID = -1l;
        this.fatherID = -1l;
        this.motherID = -1l;
        this.sID = mc.sID;
        this.iID = mc.iID;
        this.isFlexy = mc.flexy;
        this.isConstrained = mc.constricted;
        this.extCOM = mc.externalCOM;
        this.extOrient = mc.externalOrient;
        this.atomTypes = mc.atomTypes;
        this.refXYZ = mc.refXYZ;
        if(mc.atomNumbers == null){
            // recalc always
            this.atomNumbers = null;
        } else{
            assert(mc.atomNumbers != null);
            assert(mc.atomNumbers.length == mc.noOfAtoms);
            this.atomNumbers = mc.atomNumbers;
        }
        this.noOfAtoms = mc.noOfAtoms;
        this.charges = mc.charges;
        this.spins = mc.spins;
        this.energy = mc.molecularEnergy;
        if(isFlexy){
            assert(mc.zmat != null);
            this.zmat = mc.zmat;
            if(mc.degreesOfFreedom == null){
                System.out.println("found oneee...");
            }
            assert(mc.degreesOfFreedom != null);
            this.dofs = mc.degreesOfFreedom;
        }
        if(isConstrained){
            assert(mc.constraints != null);
            this.constraints = mc.constraints;
        }
    }

    /**
     * A copy constructor.
     * @param original the original molecule
     */
    Molecule(final Molecule original) {
        
        this.myID = original.myID;
        this.fatherID = original.fatherID;
        this.motherID = original.motherID;
        this.noOfAtoms = original.noOfAtoms;
        this.sID = original.sID;
        this.iID = original.iID;
        this.isFlexy = original.isFlexy;
        this.isConstrained = original.isConstrained;
        this.extCOM = original.extCOM.clone();
        this.extOrient = original.extOrient.clone();
        this.atomTypes = original.atomTypes.clone();
        this.refXYZ = new double[3][this.noOfAtoms];
        System.arraycopy(original.refXYZ[0], 0, this.refXYZ[0], 0, this.noOfAtoms);
        System.arraycopy(original.refXYZ[1], 0, this.refXYZ[1], 0, this.noOfAtoms);
        System.arraycopy(original.refXYZ[2], 0, this.refXYZ[2], 0, this.noOfAtoms);
        if(original.atomNumbers == null){
            // recalc always
            this.atomNumbers = null;
        } else{
            assert(original.atomNumbers.length == original.noOfAtoms);
            this.atomNumbers = original.atomNumbers.clone();
        }
        this.charges = original.charges.clone();
        this.spins = original.spins.clone();
        this.energy = original.energy;
        if (isFlexy) {
            this.zmat = new ZMatrix(original.zmat);
            if(original.dofs == null){
                System.out.println("FOUND ONE!!");
            }
            this.dofs = original.dofs.clone();
        }
        if(isConstrained){
            this.constraints = original.constraints;
        }
    }
    
    /**
     * A constructor working with a set of cartesian coordinates.
     * @param cartes The cartesian coordinates for this molecule.
     * @param iID the ID of this molecule as an integer
     * @param sID the ID of this molecule as a String
     */
    public Molecule(final CartesianCoordinates cartes, final int iID, final String sID) {
        this(cartes,iID,sID,false,null,false,null);
    }
    
    /**
     * A constructor working with a set of cartesian coordinates.
     * @param cartes The cartesian coordinates for this molecule.
     * @param iID the ID of this molecule as an integer
     * @param sID the ID of this molecule as a String
     * @param isFlexyMol if we deal with a flexible molecule or not
     * @param zmatDoFs the degrees of freedom of this molecule
     * @param isConstrained if this molecule is constrained or not
     * @param molConstraints the molecular constraints
     */
    public Molecule(final CartesianCoordinates cartes, final int iID, final String sID,
            final boolean isFlexyMol, final boolean[][] zmatDoFs, final boolean isConstrained,
            final boolean[][] molConstraints) {

        final MoleculeConfig mc = CoordTranslation.cartesianToBareMolConfig(cartes, iID, sID, isFlexyMol,
                zmatDoFs, isConstrained, molConstraints);
        
        this.myID = -1l;
        this.fatherID = -1l;
        this.motherID = -1l;
        
        this.sID = mc.sID;
        this.iID = iID;
        this.isFlexy = isFlexyMol;
        this.isConstrained = isConstrained;
        this.extCOM = mc.externalCOM;
        this.energy = mc.molecularEnergy;
        this.extOrient = mc.externalOrient;
        this.charges = mc.charges;
        this.spins = mc.spins;
        this.noOfAtoms = mc.noOfAtoms;
        this.atomNumbers = mc.atomNumbers;
        this.atomTypes = mc.atomTypes;
        this.refXYZ = mc.refXYZ;
        if (isFlexy) {
            this.zmat = mc.zmat;
            this.dofs = mc.degreesOfFreedom;
        }
        if(isConstrained){
            this.constraints = mc.constraints;
        }
    }
    
    @Override
    public Molecule clone(){
        return new Molecule(this);
    }

    /*
     * Getters and Setters
     */
    String[] getAtomTypes(){
        return atomTypes;
    }
    
    /**
     * Returns the reference set of cartesian coordinates for this molecule as
     * an xyz[3][noAtoms] field. Be careful!
     * @return the reference cartesians
     */
    double[][] getReferenceCartesians(){
        return refXYZ;
    }
    
    public CartesianCoordinates getCartesians() {
        return CoordTranslation.moleculeToCartesian(this, true);
    }

    double[] getExternalCenterOfMass() {
        return extCOM;
    }

    int getNumberOfAtoms() {
        return noOfAtoms;
    }

    double getEnergy() {
        return energy;
    }

    short[] getAtomNumbers(){
        if(atomNumbers == null){
            System.err.println("EXPENSIVE: Needing to reculate atomic numbers.");
            atomNumbers = new short[noOfAtoms];
            for(int x = 0; x < noOfAtoms; x++){
                atomNumbers[x] = AtomicProperties.giveAtomicNumber(atomTypes[x]);
            }
        }
        
        return atomNumbers;
    }
    
    /**
     * Returns a boolean[noAtoms][3] array of the internal degrees of freedom w.r.t. the z-Matrix.
     * @return degrees of freedom
     */
    public boolean[][] getDegreesOfFreedom() {
        return dofs;
    }

    void setExternalCenterOfMass(final double[] extCenterOfMass) {
        assert(extCenterOfMass != null);
        assert(extCenterOfMass.length == 3);
        extCOM = extCenterOfMass;
    }
    
    void setExternalCenterOfMass(final double comX, final double comY, final double comZ) {
        extCOM[0] = comX;
        extCOM[1] = comY;
        extCOM[2] = comZ;
    }
    
    void setReferenceCartesian(final double[][] xyz) {
        assert(xyz != null);
        assert(xyz.length == 3);
        assert(xyz[0].length == noOfAtoms);
        assert(xyz[1].length == noOfAtoms);
        assert(xyz[2].length == noOfAtoms);
        
        System.arraycopy(xyz[0], 0, refXYZ[0], 0, noOfAtoms);
        System.arraycopy(xyz[1], 0, refXYZ[1], 0, noOfAtoms);
        System.arraycopy(xyz[2], 0, refXYZ[2], 0, noOfAtoms);
    }

    void setReferenceCartesian(final CartesianCoordinates cartes) {
        assert(cartes != null);
        assert(cartes.getNoOfMolecules() == 1);
        assert(cartes.getNoOfAtoms()== this.noOfAtoms);
        final double[][] xyz = cartes.getAllXYZCoord();
        setReferenceCartesian(xyz);
    }

    double[] getOrientation() {
        
        extOrient[0] = CoordTranslation.sanitizePhi(extOrient[0]);
        extOrient[1] = CoordTranslation.sanitizeOmega(extOrient[1]);
        extOrient[2] = CoordTranslation.sanitizePsi(extOrient[2]);

        assert(extOrient[0] <= Math.PI && extOrient[0] >= -Math.PI);
        assert(extOrient[1] <= 0.5*Math.PI && extOrient[1] >= -0.5*Math.PI);
        assert(extOrient[2] <= Math.PI && extOrient[2] >= -Math.PI);
        
        return extOrient;
    }

    float[] getAllCharges() {
        return charges;
    }

    short[] getAllSpins() {
        return spins;
    }

    int getTotalCharge() {
        float fTotalCharge = 0;
        for (int i = 0; i < charges.length; i++) {
            fTotalCharge += charges[i];
        }
        return Math.round(fTotalCharge);
    }

    public void setEnergy(double energy) {
        this.energy = energy;
    }

    /**
     * This not only sets the new orientation of the molecule but also checks whether the new values
     * are within the valid range of the yaw-pitch-roll notation. If not, it sanitizes the value.
     * @param orientation The new orientation that should be set. Order: phi, omega, psi.
     */
    void setOrientation(final double[] orientation) {
        assert(orientation != null);
        assert(orientation.length == 3);
        /*
         * test whether we are within the valid range which is as follows for the yaw-pitch-roll
         * notation:
         * phi: -pi < phi <= pi
         * omega: -pi/2 < omega <= pi/2
         * psi: -pi < psi <= pi
         */
        extOrient[0] = orientation[0];
        extOrient[1] = orientation[1];
        extOrient[2] = orientation[2];
        sanitizeEulers();
    }
    
    void setOrientation(final double phi, final double omega, final double psi) {
        extOrient[0] = phi;
        extOrient[1] = omega;
        extOrient[2] = psi;
        sanitizeEulers();
    }
    
    void sanitizeEulers() {
        /*
         * test whether we are within the valid range which is as follows for the yaw-pitch-roll
         * notation:
         * phi: -pi < phi <= pi
         * omega: -pi/2 < omega <= pi/2
         * psi: -pi < psi <= pi
         */
        
        extOrient[0] = CoordTranslation.sanitizePhi(extOrient[0]);
        extOrient[1] = CoordTranslation.sanitizeOmega(extOrient[1]);
        extOrient[2] = CoordTranslation.sanitizePsi(extOrient[2]);

        assert(extOrient[0] <= Math.PI && extOrient[0] >= -Math.PI);
        assert(extOrient[1] <= 0.5*Math.PI && extOrient[1] >= -0.5*Math.PI);
        assert(extOrient[2] <= Math.PI && extOrient[2] >= -Math.PI);
    }

    void setZMatrix(ZMatrix zMatrix) {
        zmat = zMatrix;
    }

    public int getMolPosition() {
        return iID;
    }

    public String getSID() {
        assert(sID != null);
        assert(!sID.equalsIgnoreCase("N/A"));
        assert(!sID.isEmpty());
        return sID;
    }

    public void setID(int ID) {
        this.iID = ID;
    }

    public void setSID(String sID) {
        this.sID = sID;
    }

    public boolean getFlexy() {
        return isFlexy;
    }

    public boolean isConstricted(){
        return isConstrained;
    }

    public boolean[][] getConstraints(){
        return constraints;
    }

    ZMatrix getZMatrix() {
        return new ZMatrix(zmat);
    }

    /*
     * Methods
     */
    void setRandomCOM(final double[] maxCellSize) {
        
        final Random random = new Random();
        final double[] cell = new double[3];
        cell[0] = (maxCellSize[0]-1)/2;
        cell[1] = (maxCellSize[1]-1)/2;
        cell[2] = (maxCellSize[2]-1)/2;
        
        for (int i = 0; i < extCOM.length; i++) {
            final double d = (random.nextBoolean()) ? -random.nextDouble()*cell[i] : random.nextDouble()*cell[i];            
            extCOM[i] = d;
        }
    }

    void setRandomOrient() {
        if (getNumberOfAtoms() == 0) {
            extOrient[0] = 0.0;
            extOrient[1] = 0.0;
            extOrient[2] = 0.0;
        } else {
            RandomUtils.randomEulers(extOrient);
        }
    }

    MoleculeConfig returnMyConfig() {
        final MoleculeConfig mc = new MoleculeConfig(true);
        mc.flexy = isFlexy;
        mc.constricted = isConstrained;
        mc.externalCOM = extCOM.clone();
        mc.externalOrient = extOrient.clone();
        mc.iID = iID;
        mc.sID = sID;
        int iNoOfAtoms = noOfAtoms;
        //System.out.println("Molecule: vAtoms.size reports: " + iNoOfAtoms);
        mc.noOfAtoms = iNoOfAtoms;
        if (isFlexy) {
            mc.zmat = new ZMatrix(zmat);
            mc.degreesOfFreedom = dofs.clone();
        }
        if(isConstrained){
            mc.constraints = constraints.clone();
        }

        mc.charges = charges.clone();
        mc.spins = spins.clone();
        mc.atomTypes = this.atomTypes.clone();
        mc.atomNumbers = (this.atomNumbers == null) ? null : this.atomNumbers.clone();
        mc.refXYZ = new double[3][noOfAtoms];
        System.arraycopy(refXYZ[0], 0, mc.refXYZ[0], 0, noOfAtoms);
        System.arraycopy(refXYZ[1], 0, mc.refXYZ[1], 0, noOfAtoms);
        System.arraycopy(refXYZ[2], 0, mc.refXYZ[2], 0, noOfAtoms);
        if(isConstrained){
            mc.constraints = constraints.clone();
        }

        return mc;
    }

    /**
     * This works like a "mini-decorator" hiding different ways of mutating behind it.
     * @param whichCoord The coordinate to be mutated (or where the mutating should start).
     * @param mutFacCOM How strong the mutation should be for COMs.
     * @param mutFacEuler How strong the mutation should be for Euler angles.
     * @param whichWay Which algorithm to use for the mutation.
     */
    void mutateACoord(final int whichCoord, final double mutFacCOM, final double mutFacEuler, final int whichWay) {
        switch (whichWay) {
            case 0:
                mutateSomeCoords(whichCoord, mutFacCOM, mutFacEuler);
                break;
            case 1:
                mutateOneCoord(whichCoord, mutFacCOM, mutFacEuler);
                break;
            case 2:
                mutateSomeCoordsOfType(whichCoord, mutFacCOM, mutFacEuler);
                break;
            default:
                l.error("Unknown way " + whichWay + " for mutation!");
                break;// no mutation, something went wrong.
        }
        sanitizeEulers();
    }
    
    /**
     * Actually, this method does NOT just mutate a single gene but the gene
     * that is specified by the choice and ALL the genes of same type (COM or Euler)
     * after that. This has proven to provide a substantial speedup in comparison 
     * to the "cleaner" version.
     * Therefore the fallthrough warning can be safely suppressed. It's not a bug,
     * it's a feature! ;-)
     * @param whichCoord Where to begin with the mutation.
     * @param mutFacCOM How much the COM genes should be mutated as a factor.
     * @param mutFacEuler How much the COM genes should be mutated as a factor.
     */
    @SuppressWarnings(value = "fallthrough")
    private void mutateSomeCoordsOfType(final int whichCoord, final double mutFacCOM, final double mutFacEuler) {
        /*
         * which coordinate should be mutated?
         * 0: x
         * 1: y
         * 2: z
         * 3: Euler1: Phi
         * 4: Euler2: Omega
         * 5: Euler3: Psi
         */
        switch (whichCoord) {
            case 0:
                extCOM[0] *= mutFacCOM;
            case 1:
                extCOM[1] *= mutFacCOM;
            case 2:
                extCOM[2] *= mutFacCOM;
                break;
            case 3:
                mutateEuler(0, mutFacEuler);
            case 4:
                mutateEuler(1, mutFacEuler);
            case 5:
                mutateEuler(2, mutFacEuler);
                break;
            default:
                l.error("Unknown coordinate " + whichCoord + " for mutation!");
                //nothing
        }
    }

    /**
     * Actually, this method does NOT just mutate a single gene but the gene
     * that is specified by the choice and ALL the genes after that. This has
     * proven to provide a substantial speedup in comparison to the "cleaner"
     * version.
     * Therefore the fallthrough warning can be safely suppressed. It's not a bug,
     * it's a feature! ;-)
     * @param whichCoord Where to begin with the mutation.
     * @param mutFacCOM How much the COM genes should be mutated as a factor.
     * @param mutFacEuler How much the Euler genes should be mutated as a factor.
     */
    @SuppressWarnings(value = "fallthrough")
    private void mutateSomeCoords(final int whichCoord, double mutFacCOM, final double mutFacEuler) {
        /*
         * which coordinate should be mutated?
         * 0: x
         * 1: y
         * 2: z
         * 3: Euler1: Phi
         * 4: Euler2: Omega
         * 5: Euler3: Psi
         */
        switch (whichCoord) {
            case 0:
                extCOM[0] *= mutFacCOM;
            case 1:
                extCOM[1] *= mutFacCOM;
            case 2:
                extCOM[2] *= mutFacCOM;
            case 3:
                mutateEuler(0, mutFacEuler);
            case 4:
                mutateEuler(1, mutFacEuler);
            case 5:
                mutateEuler(2, mutFacEuler);
                break;
            default:
                l.error("Unknown coordinate " + whichCoord + " for mutation!");
                //nothing
        }
    }

    /**
     * This really just mutates a single coordinate.
     * @param whichCoord Where to begin with the mutation.
     * @param mutFacCOM How much the COM genes should be mutated as a factor.
     * @param mutFacEuler How much the Euler genes should be mutated as a factor.
     */
    private void mutateOneCoord(final int whichCoord, final double mutFacCOM, final double mutFacEuler) {
        /*
         * which coordinate should be mutated?
         * 0: x
         * 1: y
         * 2: z
         * 3: Euler1: Phi
         * 4: Euler2: Omega
         * 5: Euler3: Psi
         */
        switch (whichCoord) {
            case 0:
                extCOM[0] *= mutFacCOM;
                break;
            case 1:
                extCOM[1] *= mutFacCOM;
                break;
            case 2:
                extCOM[2] *= mutFacCOM;
                break;
            case 3:
                mutateEuler(0, mutFacEuler);
                break;
            case 4:
                mutateEuler(1, mutFacEuler);
                break;
            case 5:
                mutateEuler(2, mutFacEuler);
                break;
            default:
                l.error("Unknown coordinate " + whichCoord + " for mutation!");
                break; //nothing
        }
    }

    private void mutateEuler(final int whichEuler, final double mutFacEuler) {
        
        assert(whichEuler < 3);
        
        /*if(extOrient[whichEuler] == 0.0){
            final double[] eulers = new double[3];
            RandomUtils.randomEulers(eulers);
            extOrient[whichEuler] = eulers[whichEuler];
            return;
        }
        at the current moment the above seems too intrusive...
        */
        
        final double[] eulerCopy = extOrient.clone();
        RandomUtils.randomEulerIncrements(eulerCopy, mutFacEuler);
        extOrient[whichEuler] = eulerCopy[whichEuler];
        
        switch(whichEuler){
            case 0:
                extOrient[0] = CoordTranslation.sanitizePhi(extOrient[0]); break;
            case 1:
                extOrient[1] = CoordTranslation.sanitizeOmega(extOrient[1]); break;
            case 2:
                extOrient[2] = CoordTranslation.sanitizePsi(extOrient[2]); break;
            default:
                System.err.println("ERROR: Impossible condition in mutate euler.");
        }

        assert(extOrient[0] <= Math.PI && extOrient[0] >= -Math.PI);
        assert(extOrient[1] <= 0.5*Math.PI && extOrient[1] >= -0.5*Math.PI);
        assert(extOrient[2] <= Math.PI && extOrient[2] >= -Math.PI);
    }

    double getMyTotalWeight() {

        double totalWeight = 0.0;
        if(this.atomNumbers == null){
            for (final String atom : atomTypes) {
                totalWeight += AtomicProperties.giveWeight(atom);
            }
        } else {
            for (final short atom : atomNumbers) {
                totalWeight += AtomicProperties.giveWeight(atom);
            }
        }
        
        return totalWeight;
    }

    private void writeObject(final ObjectOutputStream oos) throws IOException {
        
        oos.writeInt(iID);
        oos.writeUTF(sID);        
        oos.writeInt(noOfAtoms);
        if(noOfAtoms == 1){
            // just save type, no, position, spin, charge, energy, constraint info
            oos.writeUTF(atomTypes[0]);
            oos.writeShort(atomNumbers[0]);
            oos.writeDouble(extCOM[0]);
            oos.writeDouble(extCOM[1]);
            oos.writeDouble(extCOM[2]);
            oos.writeShort(spins[0]);
            oos.writeFloat(charges[0]);
            oos.writeDouble(energy);
            oos.writeBoolean(isConstrained);
            if(isConstrained){
                oos.writeBoolean(constraints[0][0]);
                oos.writeBoolean(constraints[1][0]);
                oos.writeBoolean(constraints[2][0]);
            }
            
            return;
        }
        
        oos.writeBoolean(isFlexy);
        oos.writeBoolean(isConstrained);
        if(isFlexy){
            for(int i = 0; i < noOfAtoms; i++){
                for(int j = 0; j < 3; j++){
                    oos.writeBoolean(dofs[i][j]);
                }
            }
        }
        if(isConstrained){
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < noOfAtoms; j++){
                    oos.writeBoolean(constraints[i][j]);
                }
            }
        }
        oos.writeDouble(energy);
        oos.writeDouble(extCOM[0]);
        oos.writeDouble(extCOM[1]);
        oos.writeDouble(extCOM[2]);
        oos.writeDouble(extOrient[0]);
        oos.writeDouble(extOrient[1]);
        oos.writeDouble(extOrient[2]);
        for(int i = 0; i < noOfAtoms; i++){
            oos.writeShort(spins[i]);
        }
        for(int i = 0; i < noOfAtoms; i++){
            oos.writeFloat(charges[i]);
        }
        for(int i = 0; i < noOfAtoms; i++){
            oos.writeUTF(atomTypes[i]);
        }
        for(int i = 0; i < noOfAtoms; i++){
            oos.writeShort(atomNumbers[i]);
        }
        for(int j = 0; j < 3; j++){
            for(int i = 0; i < noOfAtoms; i++){
                oos.writeDouble(refXYZ[j][i]);
            }
        }
        if(isFlexy){
            oos.writeObject(zmat);
        }
    }

    private void readObject(final ObjectInputStream ois) throws IOException,
            ClassNotFoundException, ClassCastException,
            IllegalAccessException, NoSuchFieldException {

        final Class<? extends Molecule> cl = this.getClass();
        iID = ois.readInt();
        sID = ois.readUTF();
        
        noOfAtoms = ois.readInt();
        
        if(noOfAtoms == 1){
            // read the data actually needed
            final String at = ois.readUTF();
            atomTypes = new String[]{at};
        
            final short no = ois.readShort();
            atomNumbers = new short[]{no};
            
            extCOM = new double[3];
            extCOM[0] = ois.readDouble();
            extCOM[1] = ois.readDouble();
            extCOM[2] = ois.readDouble();
            
            final short sp = ois.readShort();
            spins = new short[]{sp};
            final float ch = ois.readFloat();
            charges = new float[]{ch};
            
            energy = ois.readDouble();
            
            isConstrained = ois.readBoolean();
            
            constraints = null;
            if(isConstrained){
                constraints = new boolean[3][1];
                constraints[0][0] = ois.readBoolean();
                constraints[1][0] = ois.readBoolean();
                constraints[2][0] = ois.readBoolean();
            }
            
            // initialize what is not read
            isFlexy = false;
            dofs = null;
            zmat = null;
            
            extOrient = new double[3];            
            refXYZ = new double[3][1];
            
            return;
        }
        
        isFlexy = ois.readBoolean();

        isConstrained= ois.readBoolean();

        dofs = null;
        if(isFlexy){
            dofs = new boolean[3][noOfAtoms];
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < noOfAtoms; j++){
                    dofs[i][j] = ois.readBoolean();
                }
            }
        }
        constraints = null;
        if(isConstrained){
            constraints = new boolean[noOfAtoms][3];
            for(int i = 0; i < noOfAtoms; i++){
                for(int j = 0; j < 3; j++){
                    constraints[i][j] = ois.readBoolean();
                }
            }
        }
        
        energy = ois.readDouble();
        extCOM = new double[3];
        extCOM[0] = ois.readDouble();
        extCOM[1] = ois.readDouble();
        extCOM[2] = ois.readDouble();
        extOrient = new double[3];
        extOrient[0] = ois.readDouble();
        extOrient[1] = ois.readDouble();
        extOrient[2] = ois.readDouble();

        charges = new float[noOfAtoms];
        spins = new short[noOfAtoms];
        for(int i = 0; i < noOfAtoms; i++){
            spins[i] = ois.readShort();
        }
        for(int i = 0; i < noOfAtoms; i++){
            charges[i] = ois.readFloat();
        }
        
        atomTypes = new String[noOfAtoms];
        for(int i = 0; i < noOfAtoms; i++){
            atomTypes[i] = ois.readUTF();
        }
        atomNumbers = new short[noOfAtoms];
        for(int i = 0; i < noOfAtoms; i++){
            atomNumbers[i] = ois.readShort();
        }
        refXYZ = new double[3][noOfAtoms];
        for(int j = 0; j < 3; j++){
            for(int i = 0; i < noOfAtoms; i++){
                refXYZ[j][i] = ois.readDouble();
            }
        }
        
        zmat = null;
        if(isFlexy){
            zmat = (ZMatrix) ois.readObject();
        }
    }
    
    public void moveReferenceToCOM(){
        
        final double[] com = new double[3];
        assert(atomNumbers != null);
        assert(refXYZ != null);
        assert(refXYZ.length == 3);
        assert(refXYZ[0].length == noOfAtoms);
        
        double denominator = 0.0;
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < noOfAtoms; i++) {
                // get weight
                final double weight = AtomicProperties.giveWeight(atomNumbers[i]);

                // denominator: simply add up all weights for all atoms
                if(j == 0) {denominator += weight;}
                // numerator: multiply position with weight and add up
                com[j] += weight * refXYZ[j][i];
            }
        }
        
        // divide numerator by denominator
        com[0] /= denominator;
        com[1] /= denominator;
        com[2] /= denominator;
    
        // move coordinates
        for (int i = 0; i < 3; i++) {
            final double d = com[i];
            for (int j = 0; j < noOfAtoms; j++) {
                refXYZ[i][j] -= d;
            }
        }
    }
        
    @Override
    public long getID(){
        return myID;
    }

    @Override
    public long getFatherID(){
        return this.fatherID;
    }

    @Override
    public long getMotherID(){
        return this.motherID;
    }

    @Override
    public double getFitness(){
        return this.energy;
    }
    
    @Override
    public double[] getGenomeAsDouble(){
        
        if(!this.isFlexy) {return null;}
        else{
            throw new UnsupportedOperationException("Not supported yet.");
        }
    }

    @Override
    public Double[] getGenomeCopy(){
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public void setGenome(final Double[] genome){
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setID(final long id) {
        this.myID = id;
    }

    @Override
    public void setFatherID(final long fatherID) {
        this.fatherID = fatherID;
    }

    @Override
    public void setMotherID(long motherID) {
        this.motherID = motherID;
    }

    @Override
    public void setFitness(final double fitness) {
        this.energy = fitness;
    }
    
    /**
     * THIS IS A VERY DANGEROUS ROUTINE. USE WITH CARE AFTER STUDYING THE CODE W.R.T. WHAT IT DOES (AND IMPLIES!)
     * @param xyz  the xyz coordinates IN THE CORRECT ORDER AND LENGTH
     */
    void updateAllCoordinates(final double[][] xyz, final int offset){
        
        extCOM[0] = 0.0;
        extCOM[1] = 0.0;
        extCOM[2] = 0.0;
        extOrient[0] = 0.0;
        extOrient[1] = 0.0;
        extOrient[2] = 0.0;
            
        // first copy the non-COM'd coordinates in
        System.arraycopy(xyz[0], offset, refXYZ[0], 0, noOfAtoms);
        System.arraycopy(xyz[1], offset, refXYZ[1], 0, noOfAtoms);
        System.arraycopy(xyz[2], offset, refXYZ[2], 0, noOfAtoms);
            
        // then calculate the COM
        double denominator = 0.0;
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < noOfAtoms; i++) {
                // get weight
                final double weight = AtomicProperties.giveWeight(atomNumbers[i]);

                // denominator: simply add up all weights for all atoms
                if(j == 0) {denominator += weight;}
                // numerator: multiply position with weight and add up
                extCOM[j] += weight * refXYZ[j][i];
            }
        }
        
        // divide numerator by denominator
        extCOM[0] /= denominator;
        extCOM[1] /= denominator;
        extCOM[2] /= denominator;
            
        // then substract the now calculated COM
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < noOfAtoms; i++) {
                refXYZ[j][i] -= extCOM[j];
            }
        }
        
        if(isFlexy){
            // update the zmatrix
            CoordTranslation.updateZMat(refXYZ, zmat);
        }
    }
    
    /**
     * Give the fully rotated and translated cartesian back.
     * @return a fully rotated and translated set of cartesian coordinates
     */
    double[][] giveRotTransCartesians(){
        
        final double[][] res = new double[3][noOfAtoms];
        giveRotTransCartesians(res);
        
        return res;
    }
    
    /**
     * Give the fully rotated and translated cartesian back.
     * @return a fully rotated and translated set of cartesian coordinates
     */
    void giveRotTransCartesians(final double[][] result){
        
        assert(result.length == 3);
        assert(result[0].length == noOfAtoms);
        assert(result[1].length == noOfAtoms);
        assert(result[2].length == noOfAtoms);
        
        if(noOfAtoms == 1){
            result[0][0] = refXYZ[0][0] + extCOM[0];
            result[1][0] = refXYZ[1][0] + extCOM[1];
            result[2][0] = refXYZ[2][0] + extCOM[2];
            
            return;
        }
        
        // first rotate
        CoordTranslation.rotateXYZ(refXYZ, extOrient, result);
        
        // then translate
        for(int coord = 0; coord < 3; coord++){
            for(int at = 0; at < noOfAtoms; at++){
                result[coord][at] += extCOM[coord];
            }
        }
    }
}
