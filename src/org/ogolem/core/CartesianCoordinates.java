/**
Copyright (c) 2009-2016, J. M. Dieterich and B. Hartke
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

import static java.lang.Math.cos;
import static java.lang.Math.sin;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Locale;
import static org.ogolem.core.Constants.BOHRTOANG;
import static org.ogolem.core.Constants.HARTREETOKJ;

/**
 * This provides a complete set of cartesian coordinates.
 * @author Johannes Dieterich
 * @version 2016-07-16
 */
public class CartesianCoordinates implements StructuralData, Cloneable {

    // serialVersion so that serialization works flawlessly
    private static final long serialVersionUID = (long) 2010917;

    private final int noOfMols;
    private final int noOfAtoms;
    private final int[] atomsPerMol;
    private double[][] xyzCoords;
    private String[] atomNames;
    private String method;
    private double energy;
    private float[] charges;
    private short[] spins;
    private short[] atomNumbers;
    private ZMatrix[] refZMats;
    private double[] partialEnergies;
    private Environment refEnv;

    /*
     * Constructors
     */

    /**
     * A specific constructor to initialize all sizes correct.
     * @param noOfAtoms The total number of atoms in the geometry.
     * @param noOfMols The number of molecules in the geometry.
     * @param iaAtPerMol The number of atoms per molecule.
     */
    public CartesianCoordinates(final int noOfAtoms,final int noOfMols, final int[] iaAtPerMol){
	this.atomsPerMol = iaAtPerMol;
	this.noOfAtoms = noOfAtoms;
	this.noOfMols = noOfMols;
	this.xyzCoords = new double[3][noOfAtoms];
	this.atomNames = new String[noOfAtoms];
        this.charges = new float[noOfAtoms];
        this.spins = new short[noOfAtoms];
        this.refZMats = new ZMatrix[noOfMols];
        this.partialEnergies = new double[noOfMols];
    }
    
    /**
     * A specific constructor which does NOT initialize any array if wanted.
     * @param noOfAtoms The total number of atoms in the geometry.
     * @param noOfMols The number of molecules in the geometry.
     * @param iaAtPerMol The number of atoms per molecule.
     * @param noArrayAlloc if the arrays should be allocated or not. be careful if you specify false!
     */
    public CartesianCoordinates(final int noOfAtoms, final int noOfMols,
            final int[] iaAtPerMol, final boolean noArrayAlloc){
        if(noArrayAlloc){
            this.atomsPerMol = iaAtPerMol;
            this.noOfAtoms = noOfAtoms;
            this.noOfMols = noOfMols;
        } else{
            this.atomsPerMol = iaAtPerMol;
            this.noOfAtoms = noOfAtoms;
            this.noOfMols = noOfMols;
            this.xyzCoords = new double[3][noOfAtoms];
            this.atomNames = new String[noOfAtoms];
            this.charges = new float[noOfAtoms];
            this.spins = new short[noOfAtoms];
            this.refZMats = new ZMatrix[noOfMols];
            this.partialEnergies = new double[noOfMols];
        }
    }

    /**
     * A copy constructor.
     * @param original the original Cartesian coordinates object to be copied.
     */
    public CartesianCoordinates(final CartesianCoordinates original){
        this.noOfAtoms = original.noOfAtoms;
        this.noOfMols = original.noOfMols;
        this.atomsPerMol = original.atomsPerMol.clone();
        this.atomNames = original.atomNames.clone();
        this.xyzCoords = original.getAllXYZCoordsCopy();
        this.energy = original.energy;
        this.method = original.method;
        this.charges = original.charges.clone();
        if(original.atomNumbers != null) {this.atomNumbers = original.atomNumbers.clone();}
        else {recalcAtomNumbers();} // slow
        
        this.refZMats = new ZMatrix[noOfMols];
        if(original.refZMats != null){
            for(int i = 0; i < noOfMols; i++){
                refZMats[i] = (original.refZMats[i] == null) ? null : new ZMatrix(original.refZMats[i]);
            }
        }
        
        this.partialEnergies = (original.partialEnergies == null) ? new double[noOfMols] : original.partialEnergies.clone();
        
        this.spins = original.spins.clone();
    }
    
    @Override
    public CartesianCoordinates clone(){
        return new CartesianCoordinates(this);
    }
    
    /*
     * Getters and Setters
     */
    @Override
    public int getNoOfAtoms() {
        return noOfAtoms;
    }

    @Override
    public int getNoOfMolecules() {
        return noOfMols;
    }

    int getAtomsPerMol(int i) {
        return atomsPerMol[i];
    }

    public int[] getAllAtomsPerMol() {
        return atomsPerMol;
    }

    public int getTotalCharge() {
        float fTotalCharge = 0;
        for (int i = 0; i < noOfAtoms; i++) {
            fTotalCharge += charges[i];
        }

        return Math.round(fTotalCharge);
    }

    public float getExactTotalCharge() {
        float fTotalCharge = 0;
        for (int i = 0; i < noOfAtoms; i++) {
            fTotalCharge += charges[i];
        }

        return fTotalCharge;
    }

    public int getTotalSpin() {
        int iTotalSpin = 0;
        for (int i = 0; i < noOfAtoms; i++) {
            iTotalSpin += spins[i];
        }
        return iTotalSpin;
    }

    @Override
    public float[] getAllCharges() {
        return charges;
    }

    @Override
    public short[] getAllSpins() {
        return spins;
    }

    float getChargeOfAtom(final int atom) {
        assert(atom >= 0);
        assert(atom < noOfAtoms);
        return charges[atom];
    }

    short getSpinOfAtom(final int atom) {
        assert(atom >= 0);
        assert(atom < noOfAtoms);
        return spins[atom];
    }
    
    @Override
    public short[] getAllAtomNumbers(){
        if(atomNumbers == null){
            System.err.println("INFO: Please fill the atom numbers in explicitly in the calling code. Calculating now...");
            recalcAtomNumbersForced();
        }
        return atomNumbers;
    
    }
    
    short getAtomNumberOfAtom(final int which){
        return atomNumbers[which];
    }
    
    void allocateNumbers(){
        atomNumbers = new short[noOfAtoms];
    }
    
    void setAtomNumbers(final short[] numbers){
        assert(numbers.length == noOfAtoms);
        atomNumbers = numbers;
    }

    public double[] getXYZCoordinatesOfAtom(final int atom) {
        assert(atom >= 0);
        assert(atom < noOfAtoms);
        return new double[]{xyzCoords[0][atom],xyzCoords[1][atom],xyzCoords[2][atom]};
    }

    boolean containsEnvironment(){
        return (refEnv != null);
    }

    boolean containedEnvIsFlexy(){
        return ((refEnv == null) ? false : refEnv.isEnvironmentRigid());
    }

    Atom getAtomAtPos(final int at){
        assert(at >= 0);
        assert(at < noOfAtoms);
        final AtomConfig ac = new AtomConfig();
        ac.sID = atomNames[at];
        ac.iID = at;
        ac.atomNo = atomNumbers[at];
        ac.energypart = 0.0;
        ac.position[0] = xyzCoords[0][at];
        ac.position[1] = xyzCoords[1][at];
        ac.position[2] = xyzCoords[2][at];
        final Atom atom = new Atom(ac);
        
        return atom;
    }

    /**
     * Returns a copy of the reference environment and null in case that there
     * is none.
     * @return Reference environment.
     */
    Environment getReferenceEnvironmentCopy(){
        return (refEnv == null) ? null : refEnv.clone();
    }

    public double[][] getAllXYZCoord() {
        return xyzCoords;
    }
    
    public void getAllXYZCoord(final double[][] xyz) {
        assert(xyz.length == 3);
        assert(xyz[0].length == xyzCoords[0].length);
        System.arraycopy(xyzCoords[0], 0, xyz[0], 0, noOfAtoms);
        System.arraycopy(xyzCoords[1], 0, xyz[1], 0, noOfAtoms);
        System.arraycopy(xyzCoords[2], 0, xyz[2], 0, noOfAtoms);
    }

    public double[][] getAllXYZCoordsCopy() {

        final double[][] copy = new double[3][noOfAtoms];
        System.arraycopy(xyzCoords[0], 0, copy[0], 0, noOfAtoms);
        System.arraycopy(xyzCoords[1], 0, copy[1], 0, noOfAtoms);
        System.arraycopy(xyzCoords[2], 0, copy[2], 0, noOfAtoms);
        
        return copy;
    }

    String getAtomType(final int atom) {
        assert(atom >= 0);
        assert(atom < noOfAtoms);
        return atomNames[atom];
    }

    public String[] getAllAtomTypes() {
        return atomNames;
    }

    String getMethod() {
        return method;
    }

    public double getEnergy() {
        return energy;
    }

    
    public ZMatrix[] getZMatrices(){
        return refZMats;
    }

    public ZMatrix getMolecularRefZMatrix(final int which) {
        return refZMats[which];
    }

    public double[] getAll1DCartes() {
        final double[] xyz1D = new double[noOfAtoms * 3];
        getAll1DCartes(xyz1D);
        return xyz1D;
    }
    
    public void getAll1DCartes(final double[] xyz1D) {
        assert(xyz1D != null);
        assert(xyz1D.length >= noOfAtoms*3);
        System.arraycopy(xyzCoords[0], 0, xyz1D, 0, noOfAtoms);
        System.arraycopy(xyzCoords[1], 0, xyz1D, noOfAtoms, noOfAtoms);
        System.arraycopy(xyzCoords[2], 0, xyz1D, noOfAtoms * 2, noOfAtoms);
    }

    public void setAll1DCartes(final double[] da1DXYZ, int iNoOfAts) {
        assert(xyzCoords != null);
        assert(xyzCoords[0].length == iNoOfAts);
        assert(da1DXYZ.length == 3*iNoOfAts);
        System.arraycopy(da1DXYZ, 0, xyzCoords[0], 0, iNoOfAts);
        System.arraycopy(da1DXYZ, iNoOfAts, xyzCoords[1], 0, iNoOfAts);
        System.arraycopy(da1DXYZ, iNoOfAts * 2, xyzCoords[2], 0, iNoOfAts);
    }

    public void setXYZCoordinatesOfAtom(final double[] xyz, final int i) {
        assert(xyz.length == 3);
        assert(i >= 0);
        assert(!Double.isNaN(xyz[0]));
        assert(!Double.isNaN(xyz[1]));
        assert(!Double.isNaN(xyz[2]));
        assert(!Double.isInfinite(xyz[0]));
        assert(!Double.isInfinite(xyz[1]));
        assert(!Double.isInfinite(xyz[2]));
        this.xyzCoords[0][i] = xyz[0];
        this.xyzCoords[1][i] = xyz[1];
        this.xyzCoords[2][i] = xyz[2];
    }

    public void setAllXYZ(final double[][] xyz) {
        assert(xyz.length == 3);
        assert(xyz[0].length == noOfAtoms);
        this.xyzCoords = xyz;
    }
    
    public void setAllXYZAsCopy(final double[][] xyz) {
        assert(xyz.length == 3);
        assert(xyz[0].length == noOfAtoms);
        System.arraycopy(xyz[0], 0, xyzCoords[0], 0, noOfAtoms);
        System.arraycopy(xyz[1], 0, xyzCoords[1], 0, noOfAtoms);
        System.arraycopy(xyz[2], 0, xyzCoords[2], 0, noOfAtoms);
    }

    public void setAtom(final String atom, final int pos) {
        assert(pos >= 0);
        assert(atom != null);
        this.atomNames[pos] = atom;
    }

    public void setAllAtomTypes(final String[] atomTypes) {
        assert(atomTypes != null);
        assert(atomTypes.length == noOfAtoms);
        this.atomNames = atomTypes;
        recalcAtomNumbers();
    }
    
    void setRefEnvironment(final Environment env){
        this.refEnv = env;
    }

    void setMethod(final String method) {
        this.method = method;
    }

    public void setEnergy(final double e) {
        assert(!Double.isNaN(e));
        assert(!Double.isInfinite(e));
        this.energy = e;
    }

    @Override
    public void setAllCharges(final float[] newCharges) {
        assert(newCharges != null);
        assert(newCharges.length == noOfAtoms);
        this.charges = newCharges;
    }

    @Override
    public void setAllSpins(final short[] newSpins) {
        assert(newSpins != null);
        assert(newSpins.length == noOfAtoms);
        this.spins = newSpins;
    }

    void setChargeAtAtom(final float charge, final int atom) {
        assert(atom >= 0);
        charges[atom] = charge;
    }

    void setSpinAtAtom(final short spin, final int atom) {
        assert(atom >= 0);
        spins[atom] = spin;
    }
    
    public void setZMatrices(final ZMatrix[] refZMats){
        this.refZMats = refZMats;
    }

    void setMolecularRefZMatrix(final ZMatrix zmat) {
        if(refZMats == null) {refZMats = new ZMatrix[1];}
        refZMats[0] = zmat;
    }

    /*
     * Methods
     */
    // print cartesians to StringArray
    public String[] createPrintableCartesians() {

        final String[] saToBePrinted = new String[noOfAtoms + 2];

        String sMethodID = "";
        if(method != null){
            sMethodID = method + " ";
        }

        saToBePrinted[0] = Integer.toString(noOfAtoms);
        saToBePrinted[1] = sMethodID + "Energy: " + String.format(Locale.US, "%20.4f", energy * HARTREETOKJ);

        for (int i = 2; i < noOfAtoms + 2; i++) {
            final double dXValue = xyzCoords[0][i - 2] * BOHRTOANG;
            final double dYValue = xyzCoords[1][i - 2] * BOHRTOANG;
            final double dZValue = xyzCoords[2][i - 2] * BOHRTOANG;
            saToBePrinted[i] = atomNames[i - 2] + "   " + String.format(Locale.US, "%20.7f", dXValue) + "   " + String.format(Locale.US, "%20.7f", dYValue) + "   " + String.format(Locale.US, "%20.7f", dZValue);
        }

        return saToBePrinted;
    }

    /**
     * returns a molecular Cartesian for a molecule at position i and checks for
     * problems
     * @param i the i-th molecular cartesian to return
     * @param flexy if the molecular cartesian is flexy, if it should contain its reference zmatrix
     * @return the i-th molecular cartesian, if i is wrong (too big or negative) returns a clone of this cartesian instead
     */
    public CartesianCoordinates giveMolecularCartes(final int i, final boolean flexy) {
        // if else clause to consider strange cases!
        if (i > noOfMols) {
            // wrong i, return the original
            return this.clone();
        } else if (i < 0) {
            // another wrong i, return the original
            return this.clone();
        } else if (noOfMols == 1) {
            // anyway just one molecule, return it
            return this.clone();
        } else {
            return giveMolecularCartesForced(i,flexy);
        }
    }
    
    /**
     * The non error-checking version of the molecular cartesian creating routine.
     * @param i the i-th molecular cartesian to return
     * @param flexy if the molecular cartesian is flexy, if it should contain its reference zmatrix
     * @return CartesianCoordinate object of specified molecule. NON ERROR CHECKING!
     */
    CartesianCoordinates giveMolecularCartesForced(final int i, final boolean flexy) {
        
        /*
         *  since the molecularCartes contains just one Molecule, the number of molecules is
         *  1 and the number of atoms per molecule equals the total amount of atoms
         */
        final int atoms = atomsPerMol[i];
        final int[] atomPerMol = {atoms};
        final CartesianCoordinates molecularCartes = new CartesianCoordinates(atoms, 1, atomPerMol);

        // calculate the offset of the molecule at position i
        int offset = 0;
        for (int j = 0; j < i; j++) {
            offset += atomsPerMol[j];
        }

        System.arraycopy(xyzCoords[0], offset, molecularCartes.xyzCoords[0], 0, atoms);
        System.arraycopy(xyzCoords[1], offset, molecularCartes.xyzCoords[1], 0, atoms);
        System.arraycopy(xyzCoords[2], offset, molecularCartes.xyzCoords[2], 0, atoms);
        System.arraycopy(atomNames, offset, molecularCartes.atomNames, 0, atoms);
        System.arraycopy(charges, offset, molecularCartes.charges, 0, atoms);
        System.arraycopy(spins, offset, molecularCartes.spins, 0, atoms);
        molecularCartes.allocateNumbers();
        System.arraycopy(atomNumbers, offset, molecularCartes.atomNumbers, 0, atoms);

        /*
         *  set the rest of information, method and energy
         *  (no idea whether this is actually needed...)
         */
        molecularCartes.setMethod(method);
        molecularCartes.setEnergy(partialEnergies[i]);
        
        /*
         * set the reference zmatrix
         */
        if (flexy && refZMats[i] != null) {
            final ZMatrix[] molZMat = new ZMatrix[1];
            molZMat[0] = new ZMatrix(refZMats[i]);
            molecularCartes.setZMatrices(molZMat);
        }
        
        return molecularCartes;
    }
    
    public void updateCartesians(){
        if(this.noOfMols != 1){throw new RuntimeException("ERROR: CartesianCoordinate update of Cartesians only allowed for molecular cartesians! This object holds " + noOfMols + " molecules.");}
        CoordTranslation.updateCartesians(refZMats[0], this);
    }

    public ZMatrix calculateZMatrix() {
        ZMatrix zmat = CoordTranslation.cartesToZMat(this);
        return zmat;
    }
    
    public void updateRefZMatrix(){
        if(this.noOfMols != 1){throw new RuntimeException("ERROR: CartesianCoordinate update of ZMatrix only allowed for molecular cartesians! This object holds " + noOfMols + " molecules.");}
        CoordTranslation.updateZMat(this, refZMats[0]);
    }

    String[] returnPrintablePDB(final String[] moleculeNames) {
        
        if (moleculeNames.length != noOfMols) {
            System.err.println("Too little molecular names specified in PDB creation.");
            return null;
        }

        final ArrayList<String> pdb = new ArrayList<>(noOfAtoms * 2);

        // a hash value of the coordinates used as a name
        pdb.add("COMPND    " + java.util.Arrays.hashCode(xyzCoords));
        pdb.add("AUTHOR    GENERATED BY OGOLEM");

        DecimalFormat format = new DecimalFormat("0.000");

        int iOffSet = 0;
        int iEndSet = 0;
        int iAtomID = 0;

        for (int i = 0; i < noOfMols; i++) {
            iEndSet += atomsPerMol[i];
            for (int j = iOffSet; j < iEndSet; j++) {
                /*
                 * first ATOM, then running atom number,
                 * then atom-String, then domain, then molecule running, then
                 * xyz
                 * */

                iAtomID++;
                String sTemp = "ATOM";

                String sTemp2 = "" + iAtomID;
                while (sTemp2.length() < 7) {
                    sTemp2 = " " + sTemp2;
                }

                sTemp = sTemp + sTemp2;

                // atom string: whilst pdb starts with 1, we start with 0, so -1
                sTemp2 = atomNames[iAtomID - 1] + " ";
                while (sTemp2.length() < 6) {
                    sTemp2 = " " + sTemp2;
                }

                sTemp = sTemp + sTemp2;

                // the domain
                sTemp2 = moleculeNames[i];

                sTemp = sTemp + sTemp2;

                // the molecular number
                sTemp2 = "" + i;
                while (sTemp2.length() < 6) {
                    sTemp2 = " " + sTemp2;
                }

                sTemp = sTemp + sTemp2;

                // the cartesian coordinates
                sTemp2 = format.format(xyzCoords[0][iAtomID - 1] * BOHRTOANG);
                while (sTemp2.length() < 12) {
                    sTemp2 = " " + sTemp2;
                }

                sTemp = sTemp + sTemp2;

                sTemp2 = format.format(xyzCoords[1][iAtomID - 1] * BOHRTOANG);
                while (sTemp2.length() < 8) {
                    sTemp2 = " " + sTemp2;
                }

                sTemp = sTemp + sTemp2;

                sTemp2 = format.format(xyzCoords[2][iAtomID - 1] * BOHRTOANG);
                while (sTemp2.length() < 8) {
                    sTemp2 = " " + sTemp2;
                }

                sTemp = sTemp + sTemp2;

                pdb.add(sTemp);
            }

            // after every molecule a TER
            pdb.add("TER");
            iOffSet = iEndSet;
        }

        final String[] saPDB = new String[pdb.size()];
        for (int i = 0; i < pdb.size(); i++) {
            //copy it
            saPDB[i] = pdb.get(i);
        }

        return saPDB;
    }
    
    /**
     * Rotates a single molecule in this set of Cartesians by the specified Euler angles.
     * @param molID the molecule to be rotated
     * @param eulers the Euler angles for rotation
     * @param scr a scratch field of size [3][numberOfAtoms]
     * @param rot a scratch field for the rotation matrix, size [3][3]
     * @param scr1 a scratch array of size [3]
     * @param scr2 a scratch array of size [3]
     */
    public void rotateMolecule(final int molID, final double[] eulers, final double[][] scr,
            final double[][] rot, final double[] scr1, final double[] scr2){
        
        assert(molID < noOfMols);
        assert(scr.length == 3);
        assert(scr1.length >= 3);
        assert(scr2.length >= 3);
        
        final double[] com = scr1;
        calculateMoleculeCOM(molID,com);
        moveMolCoordsFromPoint(molID, com);
        
        final int noAts = atomsPerMol[molID];
        assert(scr[0].length >= noAts);
        
        // rotate the points w/ our custom routine
        int offXYZ = 0;
        for(int i = 0; i < molID; i++){offXYZ += atomsPerMol[i];}
        rotateXYZ(xyzCoords,eulers,offXYZ,noAts,scr,rot,scr2);        
        
        // copy the rotated data back
        System.arraycopy(scr[0], 0, xyzCoords[0], offXYZ, noAts);
        System.arraycopy(scr[1], 0, xyzCoords[1], offXYZ, noAts);
        System.arraycopy(scr[2], 0, xyzCoords[2], offXYZ, noAts);
        
        // translate back
        moveMolCoordsToPoint(molID, com);
    }
    
    private static void rotateXYZ(double[][] xyz, final double[] eulers, final int offXYZ, final int lengthXYZ, final double[][] scr,
            final double[][] rot, final double[] cache) {
        
        assert(rot.length >= 3);
        assert(rot[0].length >= 3);
        assert(cache.length >= 3);
            
        // rotate using yaw-pitch-roll notation
        final double phi = eulers[0];
        final double omega = eulers[1];
        final double psi = eulers[2];
        
        if(phi == 0.0 && omega == 0.0 && psi == 0.0) {return;}
        
        final double pc = cos(phi);
        final double oc = cos(omega);
        final double sc = cos(psi);
        
        final double ps = sin(phi);
        final double os = sin(omega);
        final double ss = sin(psi);
        
        rot[0][0] = oc*sc;
        rot[1][0] = ps*os*sc-pc*ss;
        rot[2][0] = pc*os*sc+ps*ss;
        rot[0][1] = oc*ss;
        rot[1][1] = ps*os*ss+pc*sc;
        rot[2][1] = pc*os*ss-ps*sc;
        rot[0][2] = -os;
        rot[1][2] = ps*oc;
        rot[2][2] = pc*oc;
 
        // call the matrix multiplication
        org.ogolem.math.TrivialLinearAlgebra.matMult(rot, xyz, scr, 3, 3, lengthXYZ, 0, 0, 0, offXYZ);        
    }
    
    public double[] calculateMoleculeCOM(final int molID){
        
        final double[] com = new double[3];
        calculateMoleculeCOM(molID, com);
        
        return com;
    }
    
    public void calculateMoleculeCOM(final int molID, final double[] com){
        
        assert(molID < noOfMols);
        assert(atomNumbers != null);
        assert(com.length == 3);
        
        int off = 0;
        for(int i = 0; i < molID; i++){
            off += atomsPerMol[i];
        }
        
        com[0] = 0.0; com[1] = 0.0; com[2] = 0.0;
        double denominator = 0.0;
        for(int at = off; at < off+atomsPerMol[molID]; at++){
            // get weight
            final double weight = AtomicProperties.giveWeight(atomNumbers[at]);

            // denominator: simply add up all weights for all atoms
            denominator += weight;
            // numerator: multiply position with weight and add up
            for (int j = 0; j < 3; j++) {
                com[j] += weight * xyzCoords[j][at];
            }
        }
        for (int k = 0; k < 3; k++) {
            // divide numerator by denominator
            com[k] /= denominator;
        }
    }

    /**
     * This method returns the center of mass as a double array with three fields,
     * the coordinates in cartesian form.
     * @return the calculated COM of this Cartesian set as an array of size [3]
     */
    public double[] calculateTheCOM(){
        
        final double[] com = new double[3];
        calculateTheCOM(com);
        
        return com;
    }
    
    public void calculateTheCOM(final double[] com) {
        
        // formula: \mathbf{R}=\dfrac{\sum_i mi_i\cdot \mathbf{r}_i}{\sum_i m_i}
        if(atomNumbers == null){recalcAtomNumbersForced();} // inefficient, yes. XXX print warning?
        assert(com.length == 3);
        assert(xyzCoords != null);
        assert(xyzCoords.length == 3);
        assert(xyzCoords[0].length == noOfAtoms);
        
        double denominator = 0.0;
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < noOfAtoms; i++) {
                // get weight
                final double weight = AtomicProperties.giveWeight(atomNumbers[i]);

                // denominator: simply add up all weights for all atoms
                if(j == 0) {denominator += weight;}
                // numerator: multiply position with weight and add up
                com[j] += weight * xyzCoords[j][i];
            }
        }
        
        // divide numerator by denominator
        com[0] /= denominator;
        com[1] /= denominator;
        com[2] /= denominator;        
    }

    public double[] calculateCOMOfMol(final int molecule){

        assert(molecule >= 0);
        assert(molecule < noOfMols);
        assert(atomNumbers != null);
        assert(xyzCoords != null);
        assert(xyzCoords.length == 3);
        assert(xyzCoords[0].length == noOfAtoms);
        
        // find offset and end
        int offset = 0;
        int end = 0;
        for(int i = 0; i <= molecule; i++){
            end += atomsPerMol[i];
            if(i != 0) {offset += atomsPerMol[i-1];}
        }

        final double[] com = new double[3];
        // formula: \mathbf{R}=\dfrac{\sum_i mi_i\cdot \mathbf{r}_i}{\sum_i m_i}
        double denominator = 0.0;
        for (int j = 0; j < 3; j++) {
            for (int i = offset; i < end; i++) {

                // get weight for all atoms
                final double weight = AtomicProperties.giveWeight(atomNames[i]);

                // denominator: simply add up all weights for all atoms
                if(j == 0) {denominator += weight;}

                // numerator: multiply position with weight and add up
                com[j] += weight *  xyzCoords[j][i];
            }
        }
        
        // divide all numerators by denominator
        com[0] /= denominator;
        com[1] /= denominator;
        com[2] /= denominator;

        return com;
    }
    
    public void moveMolCoordsFromPoint(final int molID, final double[] point){
        
        assert(point != null);
        assert(point.length == 3);
        assert(xyzCoords != null);
        assert(xyzCoords.length == 3);
        assert(noOfAtoms == xyzCoords[0].length);
        assert(molID >= 0);
        assert(molID < noOfMols);
        
        int off = 0;
        for(int i = 0; i < molID; i++){
            off += atomsPerMol[i];
        }
        
        // this order is way more performant.
        for (int i = 0; i < 3; i++) {
            final double d = point[i];
            for (int j = off; j < off+atomsPerMol[molID]; j++) {
                xyzCoords[i][j] -= d;
            }
        }
    }

    public void moveCoordsFromPoint(final double[] point) {
        
        assert(point != null);
        assert(point.length == 3);
        assert(xyzCoords != null);
        assert(xyzCoords.length == 3);
        assert(noOfAtoms == xyzCoords[0].length);
        // this order is way more performant.
        for (int i = 0; i < 3; i++) {
            final double d = point[i];
            for (int j = 0; j < noOfAtoms; j++) {
                xyzCoords[i][j] -= d;
            }
        }
    }
    
    public void moveMolCoordsToPoint(final int molID, final double[] point){
        
        assert(point != null);
        assert(point.length == 3);
        assert(xyzCoords != null);
        assert(xyzCoords.length == 3);
        assert(noOfAtoms == xyzCoords[0].length);
        assert(molID >= 0);
        assert(molID < noOfMols);
        
        int off = 0;
        for(int i = 0; i < molID; i++){
            off += atomsPerMol[i];
        }
        
        // this order is way more performant.
        for (int i = 0; i < 3; i++) {
            final double d = point[i];
            for (int j = off; j < off+atomsPerMol[molID]; j++) {
                xyzCoords[i][j] += d;
            }
        }
    }

    public void moveCoordsToPoint(final double[] point){
        
        assert(point != null);
        assert(point.length == 3);
        assert(xyzCoords != null);
        assert(xyzCoords.length == 3);
        assert(noOfAtoms == xyzCoords[0].length);
        // this order is way more performant.
        for (int i = 0; i < 3; i++) {
            final double d = point[i];
            for (int j = 0; j < noOfAtoms; j++) {
                xyzCoords[i][j] += d;
            }
        }
    }

    /**
     * Moves this set of Cartesian coordinates to the center of mass. Meaning:
     * after calling this routine, the coordinates will have their COM in 0/0/0.
     */
    public void moveCoordsToCOM() {
        // first calculate the COM
        final double[] com = calculateTheCOM();

        // then move them
        moveCoordsFromPoint(com);
    }
    
    public final void recalcAtomNumbers(){
        if(atomNumbers == null) {recalcAtomNumbersForced();}
    }
    
    public final void recalcAtomNumbersForced(){
        if(atomNumbers == null) {atomNumbers = new short[noOfAtoms];}
        for(int i = 0; i < noOfAtoms; i++){
            atomNumbers[i] = AtomicProperties.giveAtomicNumber(atomNames[i]);
        }
    }

    /**
     * Uses the same layout as the createPrintableCartesians() method. Just that
     * this is stored in a single String with newline characters.
     * @return The cartesian set as a single String.
     */
    @Override
    public String toString(){

        // we use the String[] provided for this task
        final String[] saCartes = createPrintableCartesians();

        // now we turn this into a String

        final StringBuffer sBuff = new StringBuffer();

        for(final String sTemp : saCartes){
            sBuff.append(sTemp);
            sBuff.append(System.getProperty("line.separator"));
        }

        return sBuff.toString();
    }
    
    
    String createDefaultMolecularType(){
        
        assert(atomNames != null);
        
        String type = "";
        for (final String at : atomNames) {type += at;}

        return type;
    }
    
    void setPartialEnergy(final int i, final double energy){
        
        assert(partialEnergies != null);
        partialEnergies[i] = energy;
    }
    
    void setPartialEnergies(final double[] partialEnergies){
        assert(partialEnergies != null);
        assert(partialEnergies.length == noOfMols);
        this.partialEnergies = partialEnergies;
    }
        
    public double[] getPartialEnergies(){
        return partialEnergies;
    }
    
    /**
     * Mirror the coordinates.
     * @param mode 2: on the xy plane, 1: on the xz plane, 0: on the yz plane.
     */
    void mirrorCoordinates(final int mode){
        
        assert(mode < 3 && mode >= 0);
        
        for(int i = 0; i < 3; i++){
            if(i == mode){
                for(int j = 0; j < noOfAtoms; j++){
                    xyzCoords[i][j] = -xyzCoords[i][j];
                }
            }
        }
    }
        
    public CartesianCoordinates shuffleCartesians(final int[] mapping){
        
        final CartesianCoordinates shuffled = this.clone();
        shuffleCartesians(mapping,shuffled);
        
        return shuffled;
    }
    
    /**
     * USE WITH CAUTION! Prefer to use the object generating method! ONLY USE IN
     * TIGHT LOOPS!
     * @param mapping The mapping for the shuffling
     * @param shuffled The matrix in which the shuffling will be put into. Needs to be almost identical to this Cartesian, best use an old clone.
     */
    public void shuffleCartesians(final int[] mapping, final CartesianCoordinates shuffled){
        
        assert(mapping != null);
        assert(mapping.length == noOfMols);
        assert(shuffled.noOfAtoms == this.noOfAtoms);
        assert(shuffled.noOfMols == this.noOfMols);
        
        shuffled.method = this.method;
        shuffled.energy = this.energy;
        shuffled.refEnv = (this.refEnv == null) ? null : this.refEnv.clone();
        shuffled.refZMats = new ZMatrix[this.refZMats.length];
        
        final double[][] xyz = shuffled.getAllXYZCoord();
        final String[] ats = shuffled.getAllAtomTypes();
        final short[] spinsSh = shuffled.getAllSpins();
        final float[] chargesSh = shuffled.getAllCharges();
        final int[] atsPerMol = shuffled.getAllAtomsPerMol();
        final short[] atomNos = shuffled.getAllAtomNumbers();
        
        final int[] offsOrig = new int[noOfMols];
        final int[] offsShuffle = new int[noOfMols];
        // first shuffle the atom numbers and z matrices
        for(int i = 0; i < noOfMols; i++){
            atsPerMol[i] = atomsPerMol[mapping[i]];
            shuffled.refZMats[i] = (refZMats[mapping[i]] == null) ? null : new ZMatrix(refZMats[mapping[i]]);
        }
        
        // now compute the offsets
        int offOrig = 0;
        int offShuffle = 0;
        for(int i = 0; i < noOfMols; i++){
            offsOrig[i] = offOrig;
            offsShuffle[i] = offShuffle;
            offOrig += atomsPerMol[i];
            offShuffle += atsPerMol[i];
        }
        
        // now do the actual shuffling...
        for(int i = 0; i < noOfMols; i++){
            assert(atsPerMol[i] == atomsPerMol[mapping[i]]);
            final int offOr = offsOrig[mapping[i]];
            final int offSh = offsShuffle[i];
            final int atoms = atsPerMol[i];
            // this ought to be the fastest way to shuffle
            System.arraycopy(atomNames, offOr, ats, offSh, atoms);
            System.arraycopy(xyzCoords[0], offOr, xyz[0], offSh, atoms);
            System.arraycopy(xyzCoords[1], offOr, xyz[1], offSh, atoms);
            System.arraycopy(xyzCoords[2], offOr, xyz[2], offSh, atoms);
            System.arraycopy(spins, offOr, spinsSh, offSh, atoms);
            System.arraycopy(charges, offOr, chargesSh, offSh, atoms);
            System.arraycopy(atomNumbers, offOr, atomNos, offSh, atoms);
        }
    }

    @Override
    public CartesianCoordinates getCartesianCoordinates() {
        return this;
    }
}
