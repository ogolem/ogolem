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

import static java.lang.Math.*;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.helpers.Tuple;
import org.ogolem.math.TrivialLinearAlgebra;
import static org.ogolem.math.TrivialLinearAlgebra.crossProduct;

/**
 * This class transforms coordinates back and forth from different coordinate
 * systems and object representations.
 * @author Johannes Dieterich
 * @version 2016-08-29
 */
public final class CoordTranslation {

    private static final boolean DEBUG = false;
    private CoordTranslation(){};

    /**
     * Translates a Geometry into a CartesianCoordinates object.
     * @param geom Our geometry object that we want to transform.
     * @return The returned complete set of cartesian coordinates.
     */
    static CartesianCoordinates geometryToCartesian(final Geometry geom,
            final boolean hasEnvironment) {

        final int noOfMols = geom.getNumberOfIndieParticles();
        final int[] atomsPerMol = new int[noOfMols];

        // loop over all molecules to get the total number of atoms per geometry
        // and put the number of atoms per molecule into the array
        int totNoAtoms = 0;
        for (int i = 0; i < noOfMols; i++) {
            // sometimes Java code can be ugly too... ;-)
            final int iAtomsPerMol = geom.getMoleculeAtPosition(i).getNumberOfAtoms();
            totNoAtoms += iAtomsPerMol;
            atomsPerMol[i] = iAtomsPerMol;
        }

        final boolean flexyGeom = geom.isThereAFlexy();
        final boolean[] allFlexies = geom.getAllFlexies();
        final ZMatrix[] allZMatrices = new ZMatrix[totNoAtoms];
        
        CartesianCoordinates cartesian = new CartesianCoordinates(totNoAtoms, noOfMols, atomsPerMol);

        final float[] charges = cartesian.getAllCharges();
        final short[] spins = cartesian.getAllSpins();
        final double[][] xyz = cartesian.getAllXYZCoord();
        final String[] atoms = cartesian.getAllAtomTypes();
        final short[] nos = new short[totNoAtoms];

        int offset = 0;
        for (int i = 0; i < noOfMols; i++) {

            final Molecule molecule = geom.getMoleculeAtPosition(i);

            if(atomsPerMol[i] == 1){
                // shortcut for performance, removing a lot of object allocations
                charges[offset] = molecule.getAllCharges()[0];
                spins[offset] = molecule.getAllSpins()[0];
                atoms[offset] = molecule.getAtomTypes()[0];
                nos[offset] = molecule.getAtomNumbers()[0];
                final double[] coords = molecule.getExternalCenterOfMass();
                xyz[0][offset] = coords[0];
                xyz[1][offset] = coords[1];
                xyz[2][offset] = coords[2];
            } else{
                final int noAtThisMol = molecule.getNumberOfAtoms();
                System.arraycopy(molecule.getAllCharges(), 0, charges, offset, noAtThisMol);
                System.arraycopy(molecule.getAllSpins(), 0, spins, offset, noAtThisMol);
                System.arraycopy(molecule.getAtomTypes(), 0, atoms, offset, noAtThisMol);
                System.arraycopy(molecule.getAtomNumbers(), 0, nos, offset, noAtThisMol);

                //XXX ideally, I'd like shortcuts here as well, but I admit they are not going to be even close as pretty.
                final CartesianCoordinates molecularCartes = moleculeToCartesian(molecule, false);
                assert(molecularCartes.getNoOfAtoms() == atomsPerMol[i]);
                final double[][] xyzMol = molecularCartes.getAllXYZCoord();
                System.arraycopy(xyzMol[0], 0, xyz[0], offset, noAtThisMol);
                System.arraycopy(xyzMol[1], 0, xyz[1], offset, noAtThisMol);
                System.arraycopy(xyzMol[2], 0, xyz[2], offset, noAtThisMol);

                if (flexyGeom && allFlexies[i]) {
                    // flexy geom and flexy molecule, set refZMat
                    allZMatrices[i] = molecule.getZMatrix().clone();
                } else if (flexyGeom && !allFlexies[i]) {
                    // in principle a flexy geom, set null'd zmat for this molecule
                    allZMatrices[i] = null;
                }// else: nothing flexy, nothing to be done
            }
            cartesian.setPartialEnergy(i, molecule.getEnergy());
            offset += atomsPerMol[i];
        }

        cartesian.setAtomNumbers(nos);
        cartesian.setZMatrices(allZMatrices);
        cartesian.setEnergy(geom.getFitness());

        if(geom.containsEnvironment() && hasEnvironment){
            // we need to act on this
            final Environment env = geom.getEnvironmentCopy();
            cartesian = env.marryThem(cartesian);
        }

        
        return cartesian;
    }

    static CartesianCoordinates moleculeToCartesian(final Molecule molecule, final boolean coordsAsCopy) {

        final int noOfAtoms = molecule.getNumberOfAtoms();
        final int[] atsPerMol = {noOfAtoms};

        // since it is a molecular cartesian, the number of molecules is 1
        CartesianCoordinates cartesian;
        if (noOfAtoms == 1) {
            cartesian = new CartesianCoordinates(noOfAtoms, 1, atsPerMol, true);
            cartesian.setAllAtomTypes(molecule.getAtomTypes()); // a recalc of the atom types will automagically take place here
            final double[] coords = molecule.getExternalCenterOfMass();
            final double[][] xyz = new double[3][1];
            xyz[0][0] = coords[0];
            xyz[1][0] = coords[1];
            xyz[2][0] = coords[2];
            cartesian.setAllXYZ(xyz);
        } else {
            if (!molecule.getFlexy()){
                cartesian = new CartesianCoordinates(noOfAtoms, 1, atsPerMol, true);
                cartesian.setAtomNumbers(molecule.getAtomNumbers());
                cartesian.setAllAtomTypes(molecule.getAtomTypes());
                if(coordsAsCopy){
                    final double[][] xyz = molecule.getReferenceCartesians();
                    final double[][] xyzC = new double[3][noOfAtoms];
                    System.arraycopy(xyz[0], 0, xyzC[0], 0, noOfAtoms);
                    System.arraycopy(xyz[1], 0, xyzC[1], 0, noOfAtoms);
                    System.arraycopy(xyz[2], 0, xyzC[2], 0, noOfAtoms);
                    cartesian.setAllXYZ(xyzC);
                } else {
                    cartesian.setAllXYZ(molecule.getReferenceCartesians());
                }
                // OK, this should not be necessary. But god knows what has happened to the coordinates in the meanwhile: better be safe than sorry
                cartesian.moveCoordsToCOM();
            } else {
                final ZMatrix zmat = molecule.getZMatrix();
                molecule.moveReferenceToCOM();

                // an initial, arbitrary cartesian set
                cartesian = zMatToCartesians(zmat);

                final double[][] xyz = molecule.getReferenceCartesians();
                final CartesianCoordinates refCa = new CartesianCoordinates(noOfAtoms, 1, new int[]{noOfAtoms},true);
                refCa.setAtomNumbers(molecule.getAtomNumbers());
                refCa.setAllXYZ(xyz);
                refCa.setAllAtomTypes(molecule.getAtomTypes());

                cartesian = alignTwoCartes(refCa, cartesian);

            }

            // rotate the geometry with the euler angles
            final double[][] daXYZ = rotateXYZ(cartesian.getAllXYZCoord(), molecule.getOrientation());

            // move it to the center of mass
            moveToCOM(daXYZ, molecule.getExternalCenterOfMass());

            // put the new coordinates back into a cartesian
            cartesian.setAllXYZ(daXYZ);
            cartesian.recalcAtomNumbers();
        }
        cartesian.setAllCharges(molecule.getAllCharges());
        cartesian.setAllSpins(molecule.getAllSpins());
        cartesian.setEnergy(molecule.getEnergy());

        return cartesian;
    }
    
    /**
     * Updates a geometry object using the data in a CartesianCoordinates one. Typically a good
     * way to save some memory traffic for quickly recurring translations. Care must be taken
     * by the calling code that the two objects do indeed belong to each other as only limited
     * sanity checks are made in that regard. In particular, we do NOT currently check atom types, charges, spins.
     * Hence, this must NOT be exposed as a public API! Will NOT update energies, constraints, ...
     * @param cartes Cartesian coordinates, must not be null.
     * @param geom geometry object to be updated, must not be null, must fit the Cartesian coordinates
     */
    static void updateGeometryFromCartesian(final CartesianCoordinates cartes, final Geometry geom){
        
        assert(cartes != null);
        assert(geom != null);
        
        
        // environment (if applicable) must happen first
        CartesianCoordinates cartesian = cartes;
        if(geom.containsEnvironment()){
            // separate the environment off
            final Environment env = cartes.getReferenceEnvironmentCopy();
            cartesian = env.divorceThem(cartes);
            // put the environment in
            geom.setEnvironment(env);
        }
        
        final int noMols = geom.getNumberOfIndieParticles();
        final int noMolsCartes = cartesian.getNoOfMolecules();
        if(noMols != noMolsCartes){
            throw new RuntimeException("Geometry and CartesianCoordinates have different number of molecules in them. " + noMols + " vs " + noMolsCartes);
        }
        
        final int[] atsPerMol = cartesian.getAllAtomsPerMol();
        final double[][] allXYZ = cartesian.getAllXYZCoord();
        int xyzIndex = 0;
        for(int mol = 0; mol < noMols; mol++){
            
            final int atsMol = atsPerMol[mol];
            final int atsMolGeom = geom.getMoleculeAtPosition(mol).getNumberOfAtoms();
            if(atsMol != atsMolGeom){
                throw new RuntimeException("Geometry and CartesianCoordinates have different atom numbers at position " + mol + ": " + atsMol + " vs " + atsMolGeom);
            }
            
            if(atsPerMol[mol] == 1){
                // short cut: simply put in the external COM
                final double[] com = geom.getMoleculeAtPosition(mol).getExternalCenterOfMass();
                com[0] = allXYZ[0][xyzIndex];
                com[1] = allXYZ[1][xyzIndex];
                com[2] = allXYZ[2][xyzIndex];
            } else {
                // not so simple case
                final Molecule m = geom.getMoleculeAtPosition(mol);
                
                // set the external quantities
                final double[] comMol = cartesian.calculateCOMOfMol(mol);
                m.setExternalCenterOfMass(comMol);
                final double[] eulers = m.getOrientation();
                eulers[0] = 0.0;
                eulers[1] = 0.0;
                eulers[2] = 0.0;
                
                // set the internal quantities
                final double[][] molCartes = m.getReferenceCartesians();
                for(int coord = 0; coord < 3; coord++){
                    for(int at = 0; at < atsMol; at++){
                        molCartes[coord][at] = allXYZ[coord][xyzIndex+at] - comMol[coord];
                    }
                }
                
                // z-Matrix (if applicable)
                if(m.getFlexy()){
                    final ZMatrix zmat = m.getZMatrix();
                    updateZMat(molCartes, zmat);
                }
                
            }
            
            xyzIndex += atsMol;
        }
    }

    /**
     * This fully translates a set of cartesian coordinates into a geometry featuring
     * our internal coordinates and almost all the stuff we need.
     * @param cartesian The set of cartesian coordinates.
     * @param noOfMols number of molecules
     * @param atsPerMol atoms per molecule
     * @param molFlexies are the molecules flexible?
     * @param degreesOfFreedom a list of the degrees of freedom per molecule. Must not be null. Must be of length #molecules.
     * @param molConstraints are the molecules constrained?
     * @param constraintsXYZ the Cartesian constraints of this geometry
     * @param bonds the bond information
     * @param sids the molecular IDs as strings
     * @return The returned geometry object.
     */
    public static Geometry cartesianToGeometry(final CartesianCoordinates cartesian, final int noOfMols,
            final int[] atsPerMol, final boolean[] molFlexies,
            final List<boolean[][]> degreesOfFreedom, final boolean[] molConstraints,
            final boolean[][] constraintsXYZ, final String[] sids, final BondInfo bonds) {

        final GeometryConfig gc = cartesianToGeomConfig(cartesian, noOfMols,
                atsPerMol, molFlexies, degreesOfFreedom, molConstraints, constraintsXYZ, sids, bonds);
        final Geometry geom = new Geometry(gc);
        geom.setCoordsSynced(true);

        return geom;
    }

    static GeometryConfig cartesianToGeomConfig(final CartesianCoordinates cartes,
            final int noOfMols, final int[] atsPerMol, final boolean[] flexies,
            final List<boolean[][]> degreesOfFreedom,
            final boolean[] constraints, final boolean[][] constraintsXYZ,
            final String[] sids, final BondInfo bonds) {
        
        assert(flexies != null);
        assert(sids != null);
        assert(atsPerMol.length == sids.length);
        assert(flexies.length == sids.length);
        assert(degreesOfFreedom != null);
        assert(degreesOfFreedom.size() == cartes.getNoOfMolecules());
        
        final GeometryConfig gc = new GeometryConfig();

        CartesianCoordinates cartesian = cartes;
        if(cartes.containsEnvironment()){
            // separate the environment off
            final Environment env = cartes.getReferenceEnvironmentCopy();
            cartesian = env.divorceThem(cartes);
            // put the environment into the gc
            gc.env = env;
        }

        gc.fitness = cartesian.getEnergy();
        gc.noOfParticles = noOfMols;
        gc.bonds = bonds;
        
        /*
         * offsets for molecules are important!
         * how to deal with external ones? difficult, non trivial!
         * approach: separate the molecules and calculate for each one the COM and
         * set the orientation. Then one can move on to the internal coordinates.
         * translation and rotation:
         * set the internal COM to the origin of the cartesian coordinate system
         */
        
        final String[] allTypes = cartesian.getAllAtomTypes();
        final short[] allAtomNos = cartesian.getAllAtomNumbers();
        final short[] allSpins = cartesian.getAllSpins();
        final float[] allCharges = cartesian.getAllCharges();
        final double[] allPartialE = cartesian.getPartialEnergies();

        final ArrayList<MoleculeConfig> alMCs = new ArrayList<>(noOfMols);
        int offset = 0;
        for (int i = 0; i < noOfMols; i++) {
            
            final int molAts = atsPerMol[i];
            // keep track of the constraints
            boolean[][] molConstraints = null;
            if(constraints[i]){
                molConstraints = new boolean[3][molAts];
                for(int j = 0; j < molAts; j++){
                    molConstraints[0][j] = constraintsXYZ[0][offset+j];
                    molConstraints[1][j] = constraintsXYZ[1][offset+j];
                    molConstraints[2][j] = constraintsXYZ[2][offset+j];
                }
            }

            if(molAts == 1){
                // shortcut for performance
                final MoleculeConfig mc = new MoleculeConfig(true);
                mc.noOfAtoms = molAts;
                mc.iID = i;
                mc.sID = sids[i];
                mc.externalCOM = cartesian.getXYZCoordinatesOfAtom(offset);
                mc.externalOrient = new double[3];
                mc.refXYZ = new double[3][1];
                mc.charges = new float[]{allCharges[offset]};
                mc.spins = new short[]{allSpins[offset]};
                mc.atomTypes = new String[]{allTypes[offset]};
                mc.atomNumbers = new short[]{allAtomNos[offset]};
                mc.molecularEnergy = allPartialE[i];
                // No flexibility! (not needed to set it)
                // No constraint! (not needed to set it)
                alMCs.add(mc);
            } else{
                //XXX a shortcut here would also be nice?
                /*
                 * First step: separation
                 */
                final CartesianCoordinates molecularCartesian = cartesian.giveMolecularCartesForced(i, flexies[i]);
                /*
                 * Second step: translation
                 */
                final MoleculeConfig mc = cartesianToBareMolConfig(molecularCartesian, i, sids[i],
                        flexies[i], degreesOfFreedom.get(i), constraints[i], molConstraints);
                if(DEBUG){
                    final double[] com2 = molecularCartesian.calculateTheCOM();
                    System.out.println("DEBUG: COM in cartes2gc " + com2[0] + " " + com2[1] + " " + com2[2]);
                }
                alMCs.add(mc);
            }

            offset += molAts;
        }
        gc.geomMCs = alMCs;
        
        return gc;
    }

    static Molecule cartesianToMolecule(final CartesianCoordinates cartes, final int moleculeID,
            final String sID, final boolean flexyMolecule, final boolean[][] zmatDoFs,
            final boolean constrictedMol, final boolean[][] molConstraints) {

        final MoleculeConfig mc = cartesianToBareMolConfig(cartes, moleculeID, sID,
                flexyMolecule, zmatDoFs, constrictedMol, molConstraints);
        final Molecule molecule = new Molecule(mc);
        
        return molecule;
    }

    /**
     * This returns a MoleculeConfig when handing over cartesian coordinates and some needed
     * parameters. It will NOT set any degrees of freedom that is NOT it's purpose!
     * @param cartesian the cartesian coordinates
     * @param moleculeID this molecules ID
     * @param sID this molecules String-based ID
     * @param constrictedMol is it a constraint molecule?
     * @param flexyMolecule any explicit degrees of freedom?
     * @param zmatDoFs the flexible points of the molecule in its zmatrix format
     * @param molConstraints the cartesian-based molecular constraints (may be null if not constraint)
     * @return mc A MoleculeConfig object which can be turned into a molecule.
     */
    public static MoleculeConfig cartesianToBareMolConfig(final CartesianCoordinates cartesian,
            final int moleculeID, final String sID, final boolean flexyMolecule, final boolean[][] zmatDoFs,
            final boolean constrictedMol, final boolean[][] molConstraints) {
        
        // 1. STEP: CREATE THE NEEDED CONFIG OBJECTS
        final MoleculeConfig mc = new MoleculeConfig(true);
        final int noAtoms = cartesian.getNoOfAtoms();
        mc.noOfAtoms = noAtoms;
        mc.iID = moleculeID;
        mc.sID = sID;
        mc.charges = cartesian.getAllCharges();
        mc.spins = cartesian.getAllSpins();

        // 1a. STEP: ONE ATOM OR MORE?
        if (noAtoms == 1) {
            // No Orientation
            mc.externalOrient = new double[3];
            // No flexibility! (not needed to set it)
            // COM equals the coordinates
            mc.externalCOM = cartesian.getXYZCoordinatesOfAtom(0);
            // we need to set the reference cartesians and atom types
            mc.atomNumbers = cartesian.getAllAtomNumbers();
            mc.atomTypes = cartesian.getAllAtomTypes();
            mc.refXYZ = cartesian.getAllXYZCoord();
            // this energy part
            mc.molecularEnergy = cartesian.getEnergy();
            
            return mc;
        } else {
            // more than one atom
            // 2. STEP: COM
            final double[] com = cartesian.calculateTheCOM();
            mc.externalCOM = com;
            // move coordinates to COM
            cartesian.moveCoordsFromPoint(com);
            
            /*
             * 3. STEP: Use the CartesianCoordinates object for the molecule as a reference
             * therefore we can leave our euler angles at 0.0/0.0/0.0 . Easy, ain't it? :-)
             * we do not need to explicitly set them
             */
            // we need to set the reference cartesians and atom types
            mc.atomNumbers = cartesian.getAllAtomNumbers();
            mc.atomTypes = cartesian.getAllAtomTypes();
            mc.refXYZ = cartesian.getAllXYZCoord();

            /*
             * we just need to actually create internal coordinates for molecules that we also want
             * to optimize internally. Otherwise we can well live with the reference cartesian and
             * benefit from the reduced amount of cpu cycles needed for coordinate translation.
             */
            if (flexyMolecule) {
                final ZMatrix zmat = cartesToZMat(cartesian);
                // now set the zmatrix
                mc.flexy = true;
                mc.zmat = zmat;
                assert(zmatDoFs != null);
                assert(zmatDoFs.length == 3);
                assert(zmatDoFs[0].length == cartesian.getNoOfMolecules());
                mc.degreesOfFreedom = zmatDoFs;
            }

            if(constrictedMol){
                assert(molConstraints != null);
                assert(molConstraints.length == 3);
                assert(molConstraints[0].length == cartesian.getNoOfMolecules());
                mc.constricted = true;
                mc.constraints = molConstraints;
            }
            
            // this energy part
            mc.molecularEnergy = cartesian.getEnergy();

            return mc;
        }
    }
    
    static ZMatrix cartesToZMat(final CartesianCoordinates cartesian) {
        // get the reference ZMatrix: 0; since we are already at the molecule
        final ZMatrix zmat = cartesian.getMolecularRefZMatrix(0);
        return cartesToZMat(cartesian, zmat);
    }

    static ZMatrix cartesToZMat(final CartesianCoordinates cartesian, final ZMatrix zmat) {
        
        final ZMatrix zmat2 = new ZMatrix(zmat);
        updateZMat(cartesian,zmat2);
        return zmat2;
    }
    
    static void updateZMat(final CartesianCoordinates cartesian, final ZMatrix zmat) {

        final double[][] xyz = cartesian.getAllXYZCoord();
        updateZMat(xyz,zmat);
    }
    
    static void updateZMat(final double[][] xyz, final ZMatrix zmat) {

        final int[] bps = zmat.getAllBondConnects();
        final int[] aps = zmat.getAllAnglesConnects();
        final int[] dps = zmat.getAllDihedralConnects();
        
        final double[] scr1 = new double[3];
        final double[] scr2 = new double[3];
        final double[] scr3 = new double[3];
        final double[] scr4 = new double[3];
        final double[] scr5 = new double[3];
        for (int i = 1; i < zmat.getNoOfAtoms(); i++) {
            
            // bond length
            final double bo = distance(xyz,i,bps[i]);
            zmat.setABondLength(i, bo);
            
            if(i == 1) {continue;}
            
            // bond angle
            final double ang = calcAngle(xyz, i, bps[i], aps[i], scr1, scr2);
            zmat.setABondAngle(i, ang);
            
            if(i == 2) {continue;}
            
            // dihedral
            final double di = calcDihedral(xyz, i, bps[i], aps[i], dps[i], scr1, scr2, scr3, scr4, scr5);
            zmat.setADihedral(i, di);
        }
    }
    
    /**
     * Rotate a set of Cartesian coordinates around the x-axis.
     * @param xyz the set of Cartesian coordinates to be rotated. Must be of dimensions[3][N].
     * @param rotation the rotation angle in rad.
     * @return the rotated set of Cartesian coordinates. Will be of dimensions [3][N].
     */
    public static double[][] rotateXYZAroundX(final double[][] xyz, final double rotation){
        
        final double[][] rotated = new double[3][xyz[0].length];
        final double[][] rot = new double[3][3];
        
        rotateXYZAroundX(xyz,rotation,rotated,rot);
        
        return rotated;
    }
    
    /**
     * Rotate a set of Cartesian coordinates around the x-axis. Memory-efficient version.
     * @param xyz the set of Cartesian coordinates to be rotated. Must be of dimensions[3][N].
     * @param rotation the rotation angle in rad.
     * @param rotated a matrix of dimension [3][N] for the rotated Cartesian coordinates. Will be overwritten on return.
     * @param scr1 a scratch matrix of dimension [3][3].
     */
    public static void rotateXYZAroundX(final double[][] xyz, final double rotation,
            final double[][] rotated, final double[][] scr1){
        
        assert(xyz.length == 3);
        assert(scr1.length == 3);
        assert(scr1[0].length == 3);
        assert(scr1[1].length == 3);
        assert(scr1[2].length == 3);
        assert(rotated.length == 3);
        assert(rotated[0].length == xyz[0].length);
        assert(rotated[1].length == xyz[1].length);
        assert(rotated[2].length == xyz[2].length);
        
        final double sinR = sin(rotation);
        final double cosR = cos(rotation);
        
        final double[][] rot = scr1;
        rot[0][0] = 1.0;
        rot[0][1] = 0.0;
        rot[0][2] = 0.0;
        rot[1][0] = 0.0;
        rot[1][1] = cosR;
        rot[1][2] = -sinR;
        rot[2][0] = 0.0;
        rot[2][1] = sinR;
        rot[2][2] = cosR;
        
        // call the matrix multiplication
        org.ogolem.math.TrivialLinearAlgebra.matMult(rot, xyz, rotated);
    }
    
    /**
     * Rotate a set of Cartesian coordinates around the y-axis.
     * @param xyz the set of Cartesian coordinates to be rotated. Must be of dimensions[3][N].
     * @param rotation the rotation angle in rad.
     * @return the rotated set of Cartesian coordinates. Will be of dimensions [3][N].
     */
    public static double[][] rotateXYZAroundY(final double[][] xyz, final double rotation){
        
        final double[][] rotated = new double[3][xyz[0].length];
        final double[][] rot = new double[3][3];
        
        rotateXYZAroundY(xyz,rotation,rotated,rot);
        
        return rotated;
    }
    
    /**
     * Rotate a set of Cartesian coordinates around the y-axis. Memory-efficient version.
     * @param xyz the set of Cartesian coordinates to be rotated. Must be of dimensions[3][N].
     * @param rotation the rotation angle in rad.
     * @param rotated a matrix of dimension [3][N] for the rotated Cartesian coordinates. Will be overwritten on return.
     * @param scr1 a scratch matrix of dimension [3][3].
     */
    public static void rotateXYZAroundY(final double[][] xyz, final double rotation,
            final double[][] rotated, final double[][] scr1){
        
        assert(xyz.length == 3);
        assert(scr1.length == 3);
        assert(scr1[0].length == 3);
        assert(scr1[1].length == 3);
        assert(scr1[2].length == 3);
        assert(rotated.length == 3);
        assert(rotated[0].length == xyz[0].length);
        assert(rotated[1].length == xyz[1].length);
        assert(rotated[2].length == xyz[2].length);
        
        final double sinR = sin(rotation);
        final double cosR = cos(rotation);
        
        final double[][] rot = scr1;
        rot[0][0] = cosR;
        rot[0][1] = 0.0;
        rot[0][2] = sinR;
        rot[1][0] = 0.0;
        rot[1][1] = 1.0;
        rot[1][2] = 0.0;
        rot[2][0] = -sinR;
        rot[2][1] = 0.0;
        rot[2][2] = cosR;
        
        // call the matrix multiplication
        org.ogolem.math.TrivialLinearAlgebra.matMult(rot, xyz, rotated);
    }

    /**
     * Rotate a set of Cartesian coordinates around the z-axis.
     * @param xyz the set of Cartesian coordinates to be rotated. Must be of dimensions[3][N].
     * @param rotation the rotation angle in rad.
     * @return the rotated set of Cartesian coordinates. Will be of dimensions [3][N].
     */
    public static double[][] rotateXYZAroundZ(final double[][] xyz, final double rotation){
        
        final double[][] rotated = new double[3][xyz[0].length];
        final double[][] rot = new double[3][3];
        
        rotateXYZAroundZ(xyz,rotation,rotated,rot);
        
        return rotated;
    }
    
    /**
     * Rotate a set of Cartesian coordinates around the z-axis. Memory-efficient version.
     * @param xyz the set of Cartesian coordinates to be rotated. Must be of dimensions[3][N].
     * @param rotation the rotation angle in rad.
     * @param rotated a matrix of dimension [3][N] for the rotated Cartesian coordinates. Will be overwritten on return.
     * @param scr1 a scratch matrix of dimension [3][3].
     */
    public static void rotateXYZAroundZ(final double[][] xyz, final double rotation,
            final double[][] rotated, final double[][] scr1){
        
        assert(xyz.length == 3);
        assert(scr1.length == 3);
        assert(scr1[0].length == 3);
        assert(scr1[1].length == 3);
        assert(scr1[2].length == 3);
        assert(rotated.length == 3);
        assert(rotated[0].length == xyz[0].length);
        assert(rotated[1].length == xyz[1].length);
        assert(rotated[2].length == xyz[2].length);
        
        final double sinR = sin(rotation);
        final double cosR = cos(rotation);
        
        final double[][] rot = scr1;
        rot[0][0] = cosR;
        rot[0][1] = -sinR;
        rot[0][2] = 0.0;
        rot[1][0] = sinR;
        rot[1][1] = cosR;
        rot[1][2] = 0.0;
        rot[2][0] = 0.0;
        rot[2][1] = 0.0;
        rot[2][2] = 1.0;
        
        // call the matrix multiplication
        org.ogolem.math.TrivialLinearAlgebra.matMult(rot, xyz, rotated);
    }
    
    static void rotateXYZAroundZ(final double[][] xyz, final double rotation,
            final double[][] rotated, final double[][] rotMat, final double[] cache){
        
        assert(xyz.length == 3);
        assert(cache.length >= 3);
        assert(rotated.length == 3);
        assert(rotated[0].length == xyz[0].length);
        assert(rotated[1].length == xyz[1].length);
        assert(rotated[2].length == xyz[2].length);
        
        final double sinR = sin(rotation);
        final double cosR = cos(rotation);
        
        final double[][] rot = rotMat;
        assert(rot.length == 3);
        assert(rot[0].length == 3);
        assert(rot[1].length == 3);
        assert(rot[2].length == 3);
        rot[0][0] = cosR;
        rot[0][1] = -sinR;
        rot[0][2] = 0.0;
        rot[1][0] = sinR;
        rot[1][1] = cosR;
        rot[1][2] = 0.0;
        rot[2][0] = 0.0;
        rot[2][1] = 0.0;
        rot[2][2] = 1.0;
        
        // call the matrix multiplication
        org.ogolem.math.TrivialLinearAlgebra.matMult(rot, xyz, rotated, 3, 3, xyz[0].length, cache);
    }

    /**
     * Rotates a matrix of xyz-coordinates using the yaw-pitch-roll definition.
     * Sanitizes the Euler angles if need be.
     * @param xyz The matrix to be rotated in XYZ coordinates (double[3][n] though!).
     * @param eulers The Euler angles. (phi, omega, psi)
     * @return The rotated coordinates. Is the same object if psi=omega=phi=0.0 .
     */
    public static double[][] rotateXYZ(final double[][] xyz, final double[] eulers) {
        
        final double[][] rot = new double[3][3];
        return rotateXYZ(xyz, eulers, rot);
    }
    
    /**
     * Rotates a matrix of xyz-coordinates using the yaw-pitch-roll definition.
     * Sanitizes the Euler angles if need be.
     * @param xyz The matrix to be rotated in XYZ coordinates (double[3][n] though!).
     * @param eulers The Euler angles. (phi, omega, psi)
     * @param rot the rotation matrix as a scratch object, size [3][3]
     * @return The rotated coordinates. Is the same object if psi=omega=phi=0.0 .
     */
    public static double[][] rotateXYZ(final double[][] xyz, final double[] eulers, final double[][] rot) {
        
        final double phi = eulers[0];
        final double omega = eulers[1];
        final double psi = eulers[2];
        if(phi == 0.0 && omega == 0.0 && psi == 0.0) {return xyz;}
        
        final double[][] rotated = new double[3][xyz[0].length];
        final double[] cache = new double[3];
        rotateXYZ(xyz,eulers,rot,rotated, cache);
        
        return rotated;
    }
        
        
    /**
     * Rotates a matrix of xyz-coordinates using the yaw-pitch-roll definition.
     * Sanitizes the Euler angles if need be.
     * @param xyz The matrix to be rotated in XYZ coordinates (double[3][n] though!).
     * @param eulers The Euler angles. (phi, omega, psi)
     * @param rot the rotation matrix as a scratch object, size [3][3]
     * @param rotated The rotated coordinates.
     * @param cache a scratch array of length n or longer
     */
    public static void rotateXYZ(final double[][] xyz, final double[] eulers, final double[][] rot,
            final double[][] rotated, final double[] cache) {
        
        assert(rot.length == 3);
        assert(rot[0].length == 3);
        assert(xyz.length == 3);
        assert(rotated.length == 3);
        assert(rotated[0].length == xyz[0].length);
        assert(cache.length >= 3);
        
        // rotate using yaw-pitch-roll notation
        double phi = eulers[0];
        double omega = eulers[1];
        double psi = eulers[2];
        
        if(phi == 0.0 && omega == 0.0 && psi == 0.0) {
            System.arraycopy(xyz[0], 0, rotated[0], 0, xyz[0].length);
            System.arraycopy(xyz[1], 0, rotated[1], 0, xyz[0].length);
            System.arraycopy(xyz[2], 0, rotated[2], 0, xyz[0].length);
        }
        
        // first step: sanitize the eulers
        phi = sanitizePhi(phi);
        omega = sanitizeOmega(omega);
        psi = sanitizePsi(psi);

        assert(phi <= Math.PI && phi >= -Math.PI);
        assert(omega <= 0.5*Math.PI && omega >= -0.5*Math.PI);
        assert(psi <= Math.PI && psi >= -Math.PI);
                
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
        org.ogolem.math.TrivialLinearAlgebra.matMult(rot, xyz, rotated, 3, 3, xyz[0].length, cache);
    }
    
    /**
     * Derives a matrix of xyz-coordinates using the yaw-pitch-roll definition
     * with respect to phi (yaw).
     * Sanitizes the Euler angles if need be.
     * @param xyz The matrix to be rotated in XYZ coordinates (double[3][n] though!).
     * @param eulers The Euler angles. (phi, omega, psi)
     * @param rot the rotation matrix as a scratch object, size [3][3]
     * @param rotated The derivative of the rotated coordinates w.r.t. phi.
     * @param cache a scratch array of length n or longer
     */
    public static void rotateXYZ_dphi(final double[][] xyz, final double[] eulers, final double[][] rot,
            final double[][] rotated, final double[] cache) {
        
        assert(rot.length == 3);
        assert(rot[0].length == 3);
        assert(xyz.length == 3);
        assert(rotated.length == 3);
        assert(rotated[0].length == xyz[0].length);
        assert(cache.length >= xyz[0].length);
        
        // rotate using yaw-pitch-roll notation
        double phi = eulers[0];
        double omega = eulers[1];
        double psi = eulers[2];
        
        // first step: sanitize the eulers
        phi = sanitizePhi(phi);
        omega = sanitizeOmega(omega);
        psi = sanitizePsi(psi);

        assert(phi <= Math.PI && phi >= -Math.PI);
        assert(omega <= 0.5*Math.PI && omega >= -0.5*Math.PI);
        assert(psi <= Math.PI && psi >= -Math.PI);
                
        final double pc = cos(phi);
        final double oc = cos(omega);
        final double sc = cos(psi);
        
        final double ps = sin(phi);
        final double os = sin(omega);
        final double ss = sin(psi);
        
        rot[0][0] = 0.0;
        rot[1][0] = pc*os*sc+ps*ss;
        rot[2][0] = -ps*os*sc+pc*ss;
        rot[0][1] = 0.0;
        rot[1][1] = pc*os*ss-ps*sc;
        rot[2][1] = -ps*os*ss-pc*sc;
        rot[0][2] = 0.0;
        rot[1][2] = pc*oc;
        rot[2][2] = -ps*oc;
 
        // call the matrix multiplication
        org.ogolem.math.TrivialLinearAlgebra.matMult(rot, xyz, rotated, 3, 3, xyz[0].length, cache);
    }
    
    /**
     * Derives a matrix of xyz-coordinates using the yaw-pitch-roll definition
     * with respect to omega (pitch).
     * Sanitizes the Euler angles if need be.
     * @param xyz The matrix to be rotated in XYZ coordinates (double[3][n] though!).
     * @param eulers The Euler angles. (phi, omega, psi)
     * @param rot the rotation matrix as a scratch object, size [3][3]
     * @param rotated The derivative of the rotated coordinates w.r.t. omega.
     * @param cache a scratch array of length n or longer
     */
    public static void rotateXYZ_domega(final double[][] xyz, final double[] eulers, final double[][] rot,
            final double[][] rotated, final double[] cache) {
        
        assert(rot.length == 3);
        assert(rot[0].length == 3);
        assert(xyz.length == 3);
        assert(rotated.length == 3);
        assert(rotated[0].length == xyz[0].length);
        assert(cache.length >= xyz[0].length);
        
        // rotate using yaw-pitch-roll notation
        double phi = eulers[0];
        double omega = eulers[1];
        double psi = eulers[2];
        
        // first step: sanitize the eulers
        phi = sanitizePhi(phi);
        omega = sanitizeOmega(omega);
        psi = sanitizePsi(psi);

        assert(phi <= Math.PI && phi >= -Math.PI);
        assert(omega <= 0.5*Math.PI && omega >= -0.5*Math.PI);
        assert(psi <= Math.PI && psi >= -Math.PI);
                
        final double pc = cos(phi);
        final double oc = cos(omega);
        final double sc = cos(psi);
        
        final double ps = sin(phi);
        final double os = sin(omega);
        final double ss = sin(psi);
        
        rot[0][0] = -os*sc;
        rot[1][0] = ps*oc*sc;
        rot[2][0] = pc*oc*sc;
        rot[0][1] = -os*ss;
        rot[1][1] = ps*oc*ss;
        rot[2][1] = pc*oc*ss;
        rot[0][2] = -oc;
        rot[1][2] = -ps*os;
        rot[2][2] = -pc*os;
        
        // call the matrix multiplication
        org.ogolem.math.TrivialLinearAlgebra.matMult(rot, xyz, rotated, 3, 3, xyz[0].length, cache);
    }
    
    /**
     * Derives a matrix of xyz-coordinates using the yaw-pitch-roll definition
     * with respect to psi (roll).
     * Sanitizes the Euler angles if need be.
     * @param xyz The matrix to be rotated in XYZ coordinates (double[3][n] though!).
     * @param eulers The Euler angles. (phi, omega, psi)
     * @param rot the rotation matrix as a scratch object, size [3][3]
     * @param rotated The derivative of the rotated coordinates w.r.t. psi.
     * @param cache a scratch array of length n or longer
     */
    public static void rotateXYZ_dpsi(final double[][] xyz, final double[] eulers, final double[][] rot,
            final double[][] rotated, final double[] cache) {
        
        assert(rot.length == 3);
        assert(rot[0].length == 3);
        assert(xyz.length == 3);
        assert(rotated.length == 3);
        assert(rotated[0].length == xyz[0].length);
        assert(cache.length >= xyz[0].length);
        
        // rotate using yaw-pitch-roll notation
        double phi = eulers[0];
        double omega = eulers[1];
        double psi = eulers[2];
        
        // first step: sanitize the eulers
        phi = sanitizePhi(phi);
        omega = sanitizeOmega(omega);
        psi = sanitizePsi(psi);

        assert(phi <= Math.PI && phi >= -Math.PI);
        assert(omega <= 0.5*Math.PI && omega >= -0.5*Math.PI);
        assert(psi <= Math.PI && psi >= -Math.PI);
                
        final double pc = cos(phi);
        final double oc = cos(omega);
        final double sc = cos(psi);
        
        final double ps = sin(phi);
        final double os = sin(omega);
        final double ss = sin(psi);
        
        rot[0][0] = -oc*ss;
        rot[1][0] = -ps*os*ss-pc*sc;
        rot[2][0] = -pc*os*ss+ps*sc;
        rot[0][1] = oc*sc;
        rot[1][1] = ps*os*sc-pc*ss;
        rot[2][1] = pc*os*sc+ps*ss;
        rot[0][2] = 0.0;
        rot[1][2] = 0.0;
        rot[2][2] = 0.0;
         
        // call the matrix multiplication
        org.ogolem.math.TrivialLinearAlgebra.matMult(rot, xyz, rotated, 3, 3, xyz[0].length, cache);
    }
    
    public static final double sanitizePhi(final double phi){
        return sanitizePeriodic(phi,-Math.PI,Math.PI,2*Math.PI);
    }
    
    public static final double sanitizeOmega(final double omega){
        return sanitizePeriodic(omega,-0.5*Math.PI,0.5*Math.PI,Math.PI);
    }
    
    public static final double sanitizePsi(final double psi){
        return sanitizePeriodic(psi,-Math.PI,Math.PI,2*Math.PI);
    }
    
    private static double sanitizePeriodic(final double val, final double lower, final double upper, final double period){
        
        if(val > upper){return lower+Math.abs((val-lower)%period);}
        else if(val < lower){return upper-Math.abs((val+lower)%period);}
        
        return val;
    }
    
    /**
     * Moves a set of cartesian coordinates to the center of mass.
     * @param daXYZ The cartesian coordinates to be moved.
     * @param daCOM The center of mass.
     * @return The moved xyz coordinates.
     */
    private static void moveToCOM(final double[][] daXYZ, final double[] daCOM) {
        
        // add the COM values to all XYZ coordinates
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < daXYZ[0].length; j++) {
                daXYZ[i][j] += daCOM[i];
            }
        }
    }

    /**
     * Checks a geometry for bonds and returns the informations as a boolean[][] array.
     * @param cartes The geometry.
     * @param blowFacBonds factor for bond detection.
     * @return The bond information as a simple bond info object
     */
    static SimpleBondInfo checkForBonds(final Geometry geom, final double blowFacBonds) {

        final int noMols = geom.getNumberOfIndieParticles();
        final int noOfAtoms = geom.getNumberOfAtoms();
        final SimpleBondInfo bonds = new SimpleBondInfo(noOfAtoms);
        
        for (int i = 0; i < noOfAtoms; i++) {
            // of course between two identical atoms there is a bond. kind of...
            bonds.setBond(i, i, BondInfo.UNCERTAIN);
        }
        
        final List<double[][]> rotTransXYZ = new ArrayList<>();
        for(int part = 0; part < noMols; part++){
            final Molecule mol = geom.getMoleculeAtPosition(part);
            rotTransXYZ.add(mol.giveRotTransCartesians());
        }
        
        int off1 = 0;
        for(int part1 = 0; part1 < noMols; part1++){
            
            final Molecule mol1 = geom.getMoleculeAtPosition(part1);
            final double[][] xyz1 = rotTransXYZ.get(part1);
            final short[] atomNos1 = mol1.getAtomNumbers();
            final int noAts1 = mol1.getNumberOfAtoms();
            
            // inside each molecule:
            for (int i = 0; i < noAts1 - 1; i++) {

                final double rad1 = AtomicProperties.giveRadius(atomNos1[i]);

                for (int j = i + 1; j < noAts1; j++) {
                    final double radii = blowFacBonds * (rad1
                        + AtomicProperties.giveRadius(atomNos1[j]));
                    final double radiiSq = radii*radii;

                    final double dx = xyz1[0][i] - xyz1[0][j];
                    final double dy = xyz1[1][i] - xyz1[1][j];
                    final double dz = xyz1[2][i] - xyz1[2][j];
                
                    final double distanceSq = dx*dx+dy*dy+dz*dz;
                
                    if(distanceSq <= radiiSq){
                        bonds.setBond(off1+i,off1+j, BondInfo.UNCERTAIN);
                    }
                }
            }
            
            
            // in between molecules
            int off2 = off1+noAts1;
            for(int part2 = part1+1; part2 < noMols; part2++){
                
                final Molecule mol2 = geom.getMoleculeAtPosition(part2);
                final double[][] xyz2 = rotTransXYZ.get(part2);
                final short[] atomNos2 = mol2.getAtomNumbers();
                final int noAts2 = mol2.getNumberOfAtoms();
                
                for (int i = 0; i < noAts1; i++) {

                    final double rad1 = AtomicProperties.giveRadius(atomNos1[i]);

                    for (int j = 0; j < noAts2; j++) {
                        final double radii = blowFacBonds * (rad1
                            + AtomicProperties.giveRadius(atomNos2[j]));
                        final double radiiSq = radii*radii;

                        final double dx = xyz1[0][i] - xyz2[0][j];
                        final double dy = xyz1[1][i] - xyz2[1][j];
                        final double dz = xyz1[2][i] - xyz2[2][j];
                
                        final double distanceSq = dx*dx+dy*dy+dz*dz;
                
                        if(distanceSq <= radiiSq){
                            bonds.setBond(off1+i,off2+j, BondInfo.UNCERTAIN);
                        }
                    }
                }
                
                off2 += noAts2;
            }
            
            off1 += noAts1;
        }
        
        return bonds;
    }
    
    /**
     * Checks a set of CartesianCoordinates for bonds and returns the informations as a boolean[][] array.
     * @param cartes The cartesian coordinates.
     * @param blowFacBonds factor for bond detection.
     * @return The bond information as a simple bond info object
     */
    public static SimpleBondInfo checkForBonds(final CartesianCoordinates cartes, final double blowFacBonds) {

        final int noOfAtoms = cartes.getNoOfAtoms();
        final SimpleBondInfo bonds = new SimpleBondInfo(noOfAtoms);

        for (int i = 0; i < noOfAtoms; i++) {
            // of course between two identical atoms there is a bond. kind of...
            bonds.setBond(i, i, BondInfo.UNCERTAIN);
        }

        final double[][] xyz = cartes.getAllXYZCoord();
        final short[] atomNos = cartes.getAllAtomNumbers();
        
        for (int i = 0; i < noOfAtoms - 1; i++) {

            final double rad1 = AtomicProperties.giveRadius(atomNos[i]);

            for (int j = i + 1; j < noOfAtoms; j++) {
                final double radii = blowFacBonds * (rad1
                        + AtomicProperties.giveRadius(atomNos[j]));
                final double radiiSq = radii*radii;

                final double dx = xyz[0][i] - xyz[0][j];
                final double dy = xyz[1][i] - xyz[1][j];
                final double dz = xyz[2][i] - xyz[2][j];
                
                final double distanceSq = dx*dx+dy*dy+dz*dz;
                
                if(distanceSq <= radiiSq){
                    bonds.setBond(i,j, BondInfo.UNCERTAIN);
                }
            }
        }
        
        return bonds;
    }

    /**
     * Returns an arbitrary cartesian for the given z matrix, with 0,0,0 being the COM though. Therefore, aligning needs to take place after this.
     * @param zmat The to be transformed z matrix.
     * @return cartes The set of cartesian coordinates.
     */
    public static CartesianCoordinates zMatToCartesians(final ZMatrix zmat) {
        
        final int noOfAtoms = zmat.getNoOfAtoms();
        final int[] atsPerMol = {noOfAtoms};
        final CartesianCoordinates cartes = new CartesianCoordinates(noOfAtoms, 1, atsPerMol);
        final String[] orig = zmat.getAllAtomNames();
        final String[] next = cartes.getAllAtomTypes();
        System.arraycopy(orig, 0, next, 0, noOfAtoms);
        cartes.recalcAtomNumbers();
        
        updateCartesians(zmat,cartes);
        
        return cartes;
    }
    
    public static void updateCartesians(final ZMatrix zmat, final CartesianCoordinates cartes) {
        
        assert(cartes.getNoOfMolecules() == 1);
        assert(cartes.getNoOfAtoms() == zmat.getNoOfAtoms());

        final int noOfAtoms = zmat.getNoOfAtoms();
        final double[][] xyz = cartes.getAllXYZCoord();

        final double[] bondLengths = zmat.getAllBondLengths();
        final double[] bondAngles = zmat.getAllBondAngles();
        final double[] dihedrals = zmat.getAllDihedrals();
        final int[] bondConnects = zmat.getAllBondConnects();
        final int[] angleConnects = zmat.getAllAnglesConnects();
        final int[] dihedralConnects = zmat.getAllDihedralConnects();
        
        // first atom: in 0/0/0 -> DO NOT IGNORE! WE MAY HAVE OLD STUFF IN THERE.
        xyz[0][0] = 0.0;
        xyz[1][0] = 0.0;
        xyz[2][0] = 0.0;
        if(noOfAtoms < 2){
            cartes.moveCoordsToCOM();
            cartes.setMolecularRefZMatrix(zmat);
        
            return;
        }
        // second atom: 0/0/z (bondlength to atom 1)
        assert(bondConnects[1] == 0);
        xyz[0][1] = 0.0;
        xyz[1][1] = 0.0;
        xyz[2][1] = bondLengths[1];
        if(noOfAtoms < 3){
            cartes.moveCoordsToCOM();
            cartes.setMolecularRefZMatrix(zmat);
        
            return;
        }

        final double[] aVec = new double[3];
        final double[] vec1 = new double[3];
        final double[] vec2 = new double[3];
        final double[] vec3 = new double[3];
        final double[] nVec = new double[3];
        final double[] nnVec = new double[3];
        for (int i = 2; i < noOfAtoms; i++) {
            
            if(DEBUG){System.out.println("DEBUG: Atom " + i);}
            
            final int bConn = bondConnects[i];
            aVec[0] = xyz[0][bConn];
            aVec[1] = xyz[1][bConn];
            aVec[2] = xyz[2][bConn];
            final int aConn = angleConnects[i];
            vec2[0] = aVec[0] - xyz[0][aConn];
            vec2[1] = aVec[1] - xyz[1][aConn];
            vec2[2] = aVec[2] - xyz[2][aConn];
            final double bond = bondLengths[i];
            final double angle = bondAngles[i];
            if(DEBUG){System.out.println("DEBUG: bond and connect " + bConn + " " + bond*Constants.BOHRTOANG);}
            if(DEBUG){System.out.println("DEBUG: angle and connect " + aConn + " " + Math.toDegrees(angle));}
            
            final double dihedral;
            if (i == 2) {
                //the third atom:
                vec3[0] = aVec[0];
                vec3[1] = aVec[1] - 1.0;
                vec3[2] = aVec[2];
                dihedral = Math.PI / 2;
            } else {
                // any other atom
                final int dConn = dihedralConnects[i];
                vec3[0] = aVec[0] - xyz[0][dConn];
                vec3[1] = aVec[1] - xyz[1][dConn];
                vec3[2] = aVec[2] - xyz[2][dConn];
                if(DEBUG){System.out.println("DEBUG: dihedral and connect " + dConn + " " + Math.toDegrees(dihedrals[i]));}
                dihedral = dihedrals[i];
            }

            TrivialLinearAlgebra.crossProduct(vec2, vec3, nVec);
            TrivialLinearAlgebra.crossProduct(vec2, nVec, nnVec);
            
            // normalize the new vectors and manipulate them with the dihedral
            final double normN = (nVec[0]*nVec[0] + nVec[1]*nVec[1] + nVec[2]*nVec[2]);
            final double normNN = (nnVec[0]*nnVec[0] + nnVec[1]*nnVec[1] + nnVec[2]*nnVec[2]);
            final double tempN = -Math.sin(dihedral)/normN;
            final double tempNN = Math.cos(dihedral)/normNN;
            for (int j = 0; j < 3; j++) {
                // create the next vector
                vec1[j] = nVec[j]*tempN + nnVec[j]*tempNN;
            }

            // normalize and throw bond length and angle on it
            final double norm1 = Math.sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
            final double norm2 = Math.sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2]);
            final double tmp1 = bond * Math.sin(angle) / norm1;
            final double tmp2 = bond * Math.cos(angle) / norm2;
            for (int j = 0; j < 3; j++) {
                vec1[j] *= tmp1;
                vec2[j] *= tmp2;
            }
            
            if(DEBUG){
                System.out.println("aVec " + aVec[0]*Constants.BOHRTOANG + " " + aVec[1]*Constants.BOHRTOANG + " " + aVec[2]*Constants.BOHRTOANG);
                System.out.println("vec1 " + vec1[0]*Constants.BOHRTOANG + " " + vec1[1]*Constants.BOHRTOANG + " " + vec1[2]*Constants.BOHRTOANG);
                System.out.println("vec2 " + vec2[0]*Constants.BOHRTOANG + " " + vec2[1]*Constants.BOHRTOANG + " " + vec2[2]*Constants.BOHRTOANG);
            }
            
            // calculate the actual postions
            for (int j = 0; j < 3; j++) {
                final double d = aVec[j] + vec1[j] - vec2[j];
                if(Double.isInfinite(d) || Double.isNaN(d)){
                    throw new RuntimeException("Coordinate is NaN/Infinite: " + i + " " + j + " " + d);
                }
                xyz[j][i] = d;
            }
            if(DEBUG){System.out.println(cartes.getAllAtomTypes()[i] + "\t" + xyz[0][i] + "\t" + xyz[1][i] + "\t" + xyz[2][i]);}

        }
        
        if(DEBUG){
            System.out.println("DEBUG: Before moving to COM!");
            for(final String s : cartes.createPrintableCartesians()){
                System.out.println(s);
            }
        }
        
        cartes.recalcAtomNumbers();
        cartes.moveCoordsToCOM();
        cartes.setMolecularRefZMatrix(zmat);
        
        if(DEBUG){
            System.out.println("DEBUG: After moving to COM!");
            for(final String s : cartes.createPrintableCartesians()){
                System.out.println(s);
            }
        }
    }
    
    /**
     * Determines the 3D-Euclidean distance between cartesian coordinate points a and b
     * @param a point a
     * @param b point b
     * @return distance
     */
    public static double distance(final double[] a, final double[] b){
        
        final double dx = a[0]-b[0];
        final double dy = a[1]-b[1];
        final double dz = a[2]-b[2];
        
        return Math.sqrt(dx*dx+dy*dy+dz*dz);
    }
    
    /**
     * Determines the 3D-Euclidean distance between cartesian coordinates in points i and j
     * @param xyz set of cartesian coordinates
     * @param i point i
     * @param j point j
     * @return distance
     */
    public static double distance(final double[][] xyz, final int i, final int j){
        
        assert(i >= 0);
        assert(j >= 0);
        assert(xyz != null);
        assert(xyz.length == 3);
        assert(xyz[0].length > i);
        assert(xyz[0].length > j);
        final double dx = xyz[0][i]-xyz[0][j];
        final double dy = xyz[1][i]-xyz[1][j];
        final double dz = xyz[2][i]-xyz[2][j];
        
        return Math.sqrt(dx*dx+dy*dy+dz*dz);
    }

    public static double calcAngle(final double[][] xyz, final int i, final int j, final int k){

        // the angle is given by acos of the dot product of the two (normalised) direction vectors: v1 v2 = |v1||v2| cos(angle)
        final double[] daVectorOne = new double[3];
        final double[] daVectorTwo = new double[3];
        // calculate the direction vectors
        for (int c = 0; c < 3; c++) {
            daVectorOne[c] = xyz[c][i] - xyz[c][j];
            daVectorTwo[c] = xyz[c][k] - xyz[c][j];
        }

        final double dAngle = angle(daVectorOne, daVectorTwo);

        return dAngle;

    }
    
    public static double calcAngle(final double[][] xyz, final int i, final int j,
            final int k, final double[] scr1, final double[] scr2){

        // the angle is given by acos of the dot product of the two (normalised) direction vectors: v1 v2 = |v1||v2| cos(angle)
        final double[] daVectorOne = scr1;
        final double[] daVectorTwo = scr2;
        // calculate the direction vectors
        for (int c = 0; c < 3; c++) {
            daVectorOne[c] = xyz[c][i] - xyz[c][j];
            daVectorTwo[c] = xyz[c][k] - xyz[c][j];
        }

        final double dAngle = angle(daVectorOne, daVectorTwo);

        return dAngle;

    }
    
    public static double calcAngle(final double[] pos1, final double[] pos2, final double[] pos3){

        // the angle is given by acos of the dot product of the two (normalised) direction vectors: v1 v2 = |v1||v2| cos(angle)
        final double[] daVectorOne = new double[3];
        final double[] daVectorTwo = new double[3];
        // calculate the direction vectors
        for (int c = 0; c < 3; c++) {
            daVectorOne[c] = pos1[c] - pos2[c];
            daVectorTwo[c] = pos3[c] - pos2[c];
        }

        final double dAngle = angle(daVectorOne, daVectorTwo);

        return dAngle;

    }
    
    public static double calcAngle(final double[] pos1, final double[] pos2,
            final double[] pos3, final double[] scr1, final double[] scr2){

        // the angle is given by acos of the dot product of the two (normalised) direction vectors: v1 v2 = |v1||v2| cos(angle)
        final double[] daVectorOne = scr1;
        final double[] daVectorTwo = scr2;
        // calculate the direction vectors
        for (int c = 0; c < 3; c++) {
            daVectorOne[c] = pos1[c] - pos2[c];
            daVectorTwo[c] = pos3[c] - pos2[c];
        }

        final double dAngle = angle(daVectorOne, daVectorTwo);

        return dAngle;

    }

    private static double angle(final double[] v1, final double[] v2) {
        
        // the angle is given by acos of the dot product of the two (normalized) direction vectors: v1 v2 = |v1||v2| cos(angle)

        // normalize them
        final double n1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
        final double n2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);

        // calculate the dot product
        final double dp = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])/(n1*n2);

        // calculate the bond angle and return it
        final double ang = Math.acos(dp);
        
        return ang;
    }

    /**
     * Calculates the dihedral angle out of four cartesian coordinates in the order present in a z matrix.
     * see e.g.: http://structbio.biochem.dal.ca/jrainey/dihedralcalc.html for an explanation.
     * @param xyz cartesian set
     * @param k atom 1 for dihedral calculation
     * @param l atom 2 for dihedral calculation
     * @param m atom 3 for dihedral calculation
     * @param n atom 4 for dihedral calculation
     * @param scr1 scratch space, 3 element array
     * @param scr2 scratch space, 3 element array
     * @param scr3 scratch space, 3 element array
     * @param scr4 scratch space, 3 element array
     * @param scr5 scratch space, 3 element array
     * @return The dihedral angle.
     */    
    public static double calcDihedral(final double[][] xyz, final int k, final int l,
            final int m, final int n,
            final double[] scr1, final double[] scr2, final double[] scr3,
            final double[] scr4, final double[] scr5) {
        
        for (int i = 0; i < 3; i++) {
            scr1[i] = xyz[i][k] - xyz[i][l];
            scr2[i] = xyz[i][l] - xyz[i][m];
            scr3[i] = xyz[i][m] - xyz[i][n];
        }

        crossProduct(scr1, scr2, scr4);
        crossProduct(scr2, scr3, scr1);
        crossProduct(scr4, scr1, scr5);

        // calculate the length of the daCThree vector
        final double le = sqrt(scr5[0]*scr5[0] + scr5[1]*scr5[1] + scr5[2]*scr5[2]);

        if (le < 0.001) {
            // take care of numerical inaccuracies
            return 0.0;
        }

        // calculate the dihedral as the vector angle
        final double dihedral = angle(scr4, scr1);

        return (TrivialLinearAlgebra.dotProduct(scr2, scr5) > 0.0) ? -dihedral : dihedral;
    }
    
    public static double calcDihedral(final double[][] xyz, final int k, final int l,
            final int m, final int n) {

        final double[] scr1 = new double[3];
        final double[] scr2 = new double[3];
        final double[] scr3 = new double[3];
        final double[] scr4 = new double[3];
        final double[] scr5 = new double[3];
        
        return calcDihedral(xyz,k,l,m,n,scr1,scr2,scr3,scr4,scr5);
    }
    
    public static double calcDihedral(final double[] pos1, final double[] pos2,
                final double[] pos3, final double[] pos4,
            final double[] scr1, final double[] scr2, final double[] scr3,
            final double[] scr4, final double[] scr5) {
        
        for (int i = 0; i < 3; i++) {
            scr1[i] = pos1[i] - pos2[i];
            scr2[i] = pos2[i] - pos3[i];
            scr3[i] = pos3[i] - pos4[i];
        }

        crossProduct(scr1, scr2, scr4);
        crossProduct(scr2, scr3, scr1);
        crossProduct(scr4, scr1, scr5);

        // calculate the length of the daCThree vector
        final double le = sqrt(scr5[0]*scr5[0] + scr5[1]*scr5[1] + scr5[2]*scr5[2]);

        if (le < 0.001) {
            // take care of numerical inaccuracies
            return 0.0;
        }

        // calculate the dihedral as the vector angle
        final double dihedral = angle(scr4, scr1);

        return (TrivialLinearAlgebra.dotProduct(scr2, scr5) > 0.0) ? -dihedral : dihedral;
    }
    
    public static double calcDihedral(final double[] pos1, final double[] pos2,
            final double[] pos3, final double[] pos4) {

        final double[] scr1 = new double[3];
        final double[] scr2 = new double[3];
        final double[] scr3 = new double[3];
        final double[] scr4 = new double[3];
        final double[] scr5 = new double[3];
        
        return calcDihedral(pos1,pos2,pos3,pos4,scr1,scr2,scr3,scr4,scr5);
    }
    
    /**
     * Aligns a set of cartesian coordinates to a set of reference cartesian coordinates.
     * The aligning uses a least-squares fit algorithm described in
     * S. K. Kearsley, Acta. Cryst., 1989, A45c, 208-210
     * and stems from reference implementations by
     * B. Rupp and S. Parikin (www.structure.illn.gov) in Fortran and
     * M. Eriksson in C
     * @param refCartes The reference set. Will be assumed to be in the COM.
     * @param cartes The to be aligned set. Will not be touched but a copy be moved to the COM.
     * @return A tuple of the aligned cartesian and the RMSD of the aligned to the reference.
     * @throws java.lang.Exception e.g. if the diagonalization fails
     */
    public static Tuple<CartesianCoordinates,Double> alignTwoCartesians(final CartesianCoordinates refCartes,
            final CartesianCoordinates cartes) throws Exception {
        return alignTwoCartesians(refCartes.getAllXYZCoord(), cartes);
    }
        
    /**
     * Aligns a set of cartesian coordinates to a set of reference cartesian coordinates.
     * The aligning uses a least-squares fit algorithm described in
     * S. K. Kearsley, Acta. Cryst., 1989, A45c, 208-210
     * and stems from reference implementations by
     * B. Rupp and S. Parikin (www.structure.illn.gov) in Fortran and
     * M. Eriksson in C
     * @param refCartes The reference set. Will be assumed to be in the COM.
     * @param cartes The to be aligned set. Will not be touched but a copy be moved to the COM.
     * @return A tuple of the aligned cartesian and the RMSD of the aligned to the reference.
     * @throws java.lang.Exception e.g. if the diagonalization fails
     */
    public static Tuple<CartesianCoordinates,Double> alignTwoCartesians(final double[][] refCartes,
            final CartesianCoordinates cartes) throws Exception {
        
        assert(refCartes != null);
        assert(refCartes.length == 3);
        assert(refCartes[0].length > 0);
        assert(cartes != null);
        assert(cartes.getNoOfAtoms() > 0);
        
        final CartesianCoordinates aligned = new CartesianCoordinates(cartes);
        aligned.moveCoordsToCOM();
        
        // set up the kearsley matrix
        final contrib.jama.Matrix matKearsley = new contrib.jama.Matrix(4,4);
        setupKearsleyMatrix(refCartes, aligned.getAllXYZCoord(), matKearsley.getArray());

        final contrib.jama.EigenvalueDecomposition eigen = new contrib.jama.EigenvalueDecomposition(matKearsley, true);
        final contrib.jama.Matrix matEigenVec = eigen.getV();
        final double[] evals = eigen.getRealEigenvalues();
        /*
        for(final double e : evals){
            System.err.println("Eigenval " + e);
        }*/

        // calculate the rotation matrix
        final double[][] rot = new double[3][3];
        calculateRotationMatrix(matEigenVec.getArray(), rot);

        // rotate the cartesians
        final double[][] newCartes = new double[3][cartes.getNoOfAtoms()];
        TrivialLinearAlgebra.matMult(rot, aligned.getAllXYZCoord(), newCartes);

        aligned.setAllXYZ(newCartes);

        final double rmsd = Math.sqrt(evals[0]/Math.min(refCartes[0].length,cartes.getNoOfAtoms()));
        
        return new Tuple<>(aligned, rmsd);
    }

    /**
     * Aligns a set of cartesian coordinates to a set of reference cartesian coordinates.
     * The aligning uses a least-squares fit algorithm described in
     * S. K. Kearsley, Acta. Cryst., 1989, A45c, 208-210
     * and stems from reference implementations by
     * B. Rupp and S. Parikin (www.structure.illn.gov) in Fortran and
     * M. Eriksson in C
     * @param refCartes The reference set.
     * @param cartes The to be aligned set.
     * @return cartes The set of cartesian coordinates aligned to the reference set.
     */
    public static CartesianCoordinates alignTwoCartes(final CartesianCoordinates refCartes,
            final CartesianCoordinates cartes) {

        /*
         * first move the cartesian coordinate set to the COM. this might not be needed in all
         * cases but we do it out of safety
         */
        refCartes.moveCoordsToCOM();
        
        // solve the thing with the other routine
        try {
            final Tuple<CartesianCoordinates, Double> tup = alignTwoCartesians(refCartes,cartes);
            return tup.getObject1();
        } catch (Exception e) {
            // return the non aligned molecule, if the EigenValueDecomposition fails (which might happen)
            System.err.println("WARNING: Failed to align molecules. " + e.toString() + "\nReturning non-aligned cartes!");
            return cartes;
        }

        /*
         * FROM HERE ON IT IS OLD (not-tested!!!) code which attempts to align using moments of
         * intertia.
         * It stayed here since it might even be useful at some later point.

        AtomicProperties props = new AtomicProperties();
        double[] daMasses = new double[iNoOfAtoms];
        double[] daDistances1 = new double[iNoOfAtoms];
        double[] daDistances2 = new double[iNoOfAtoms];

        // get all masses and all distances to the COM
        for(int i = 0; i < iNoOfAtoms; i++){
        daMasses[i] = props.giveWeight(saAtoms[i]);
        daDistances1[i] = Math.sqrt(Math.pow(daRefCoords[0][i], 2) + Math.pow(daRefCoords[1][i], 2) + Math.pow(daRefCoords[2][i], 2));
        daDistances2[i] = Math.sqrt(Math.pow(daCoords[0][i], 2) + Math.pow(daCoords[1][i], 2) + Math.pow(daCoords[2][i], 2));
        }

        // calculate the moment of inertia tensor for both sets
        Matrix tensor1 = new Matrix(3,3);
        Matrix tensor2 = new Matrix(3,3);
        double[][] daTensor1 = tensor1.getArray();
        double[][] daTensor2 = tensor2.getArray();

        int iKroneckerDelta;

        for(int i = 0; i < 3; i++){
        for(int j = i; j < 3; j++){
        if(i == j){
        iKroneckerDelta = 1;
        } else{
        iKroneckerDelta = 0;
        }

        for(int k = 0; k < iNoOfAtoms; k++){
        daTensor1[i][j] += daMasses[k] * (Math.pow(daDistances1[k],2)*iKroneckerDelta - daRefCoords[i][k] * daRefCoords[j][k]);
        daTensor2[i][j] += daMasses[k] * (Math.pow(daDistances2[k],2)*iKroneckerDelta - daCoords[i][k] * daCoords[j][k]);
        }

        // the tensor is actually symmetric
        daTensor1[j][i] = daTensor1[i][j];
        daTensor2[j][i] = daTensor2[i][j];
        }
        }

        // calculate the eigenvectors of the tensors
        EigenvalueDecomposition eigen1 = new EigenvalueDecomposition(tensor1);
        EigenvalueDecomposition eigen2 = new EigenvalueDecomposition(tensor2);

        Matrix eigenvecs1 = eigen1.getV();
        Matrix eigenvecs2 = eigen2.getV();

        //calculate the eigenvalues
        double[] daEigenvals1 = eigen1.getRealEigenvalues();
        double[] daEigenvals2 = eigen2.getRealEigenvalues();

        double[][] daEigenVecs1 = eigenvecs1.getArray();
        double[][] daEigenVecs2 = eigenvecs2.getArray();
        double[][] daOrderedVecs1 = new double[3][3];
        double[][] daOrderedVecs2 = new double[3][3];

        // order them by their eigenvalue
        if(daEigenvals1[0] >= daEigenvals1[1] && daEigenvals1[0] >= daEigenvals1[2]){
        // the first eigenvector is the one having the biggest eigenvalue
        daOrderedVecs1[0] = daEigenVecs1[0];

        if(daEigenvals1[1] >= daEigenvals1[2]){
        daOrderedVecs1[1] = daEigenVecs1[1];
        daOrderedVecs1[2] = daEigenVecs1[2];
        } else{
        daOrderedVecs1[1] = daEigenVecs1[2];
        daOrderedVecs1[2] = daEigenVecs1[1];
        }


        } else if(daEigenvals1[1] >= daEigenvals1[0] && daEigenvals1[1] >= daEigenvals1[2]){
        // the second should be in slot one
        daOrderedVecs1[0] = daEigenVecs1[1];

        if(daEigenvals1[0] >= daEigenvals1[2]){
        daOrderedVecs1[1] = daEigenVecs1[0];
        daOrderedVecs1[2] = daEigenVecs1[2];
        } else{
        daOrderedVecs1[1] = daEigenVecs1[2];
        daOrderedVecs1[2] = daEigenVecs1[0];
        }


        } else if(daEigenvals1[2] >= daEigenvals1[0] && daEigenvals1[2] >= daEigenvals1[1]){
        // the third should be in slot one
        daOrderedVecs1[0] = daEigenVecs1[2];

        if(daEigenvals1[0] >= daEigenvals1[1]){
        daOrderedVecs1[1] = daEigenVecs1[0];
        daOrderedVecs1[2] = daEigenVecs1[1];
        } else{
        daOrderedVecs1[1] = daEigenVecs1[1];
        daOrderedVecs1[2] = daEigenVecs1[0];
        }


        } else {
        System.err.println("Coordinate translation reports problems in the aligning step. First vector.");
        return cartes;
        }


        if(daEigenvals2[0] >= daEigenvals2[1] && daEigenvals2[0] >= daEigenvals2[2]){
        // the first eigenvector is the one having the biggest eigenvalue
        daOrderedVecs2[0] = daEigenVecs2[0];

        if(daEigenvals1[1] >= daEigenvals1[2]){
        daOrderedVecs2[1] = daEigenVecs2[1];
        daOrderedVecs2[2] = daEigenVecs2[2];
        } else{
        daOrderedVecs2[1] = daEigenVecs2[2];
        daOrderedVecs2[2] = daEigenVecs2[1];
        }


        } else if(daEigenvals2[1] >= daEigenvals2[0] && daEigenvals2[1] >= daEigenvals2[2]){
        // the second should be in slot one
        daOrderedVecs2[0] = daEigenVecs2[1];

        if(daEigenvals2[0] >= daEigenvals2[2]){
        daOrderedVecs2[1] = daEigenVecs2[0];
        daOrderedVecs2[2] = daEigenVecs2[2];
        } else{
        daOrderedVecs2[1] = daEigenVecs2[2];
        daOrderedVecs2[2] = daEigenVecs2[0];
        }


        } else if(daEigenvals2[2] >= daEigenvals2[0] && daEigenvals2[2] >= daEigenvals2[1]){
        // the third should be in slot one
        daOrderedVecs2[0] = daEigenVecs2[2];

        if(daEigenvals2[0] >= daEigenvals2[1]){
        daOrderedVecs2[1] = daEigenVecs2[0];
        daOrderedVecs2[2] = daEigenVecs2[1];
        } else{
        daOrderedVecs2[1] = daEigenVecs2[1];
        daOrderedVecs2[2] = daEigenVecs2[0];
        }


        } else {
        System.err.println("Coordinate translation reports problems in the aligning step. Second vector.");
        return cartes;
        }


        /*
         * now the eigenvectors are ordered by there "size", aligning
         * should take place now. For this we need the rotation matrix, that
         * translates Ref = Rot*A (all are matrices).
         * Therefore: Ref*A^{-1} = Rot
         */
        /*        Matrix rotMatrix;
        Matrix refInertias = new Matrix(daOrderedVecs1);
        Matrix newInertias = new Matrix(daOrderedVecs2);

        newInertias = newInertias.inverse();

        rotMatrix = refInertias.times(newInertias);


        /*
         * now that we have our rotation matrix, we can easily rotate the whole
         * set of cartesian coordinates.
         * Aligned = rot*old
         */
        /*        Matrix alignedCartes;
        Matrix oldCartes = new Matrix(daCoords);

        alignedCartes = rotMatrix.times(oldCartes);

        double[][] daNewCartes = alignedCartes.getArrayCopy();

        cartes.setAllXYZ(daNewCartes);
         */        
    }

    /**
     * Translates a set of cartesian coordinates into a set of spherical coordinates.
     * @param xyz cartesian coordinates
     * @return spherical coordinates as a double[][]
     */
    public static double[][] cartesianToSphericalCoord(final double[][] xyz) {

        final int noOfAtoms = xyz[0].length;
        final double[][] sphericalCoords = new double[3][noOfAtoms];

        for (int i = 0; i < noOfAtoms; i++) {

            final double x = xyz[0][i];
            final double y = xyz[1][i];
            final double z = xyz[2][i];

            // first the r value
            sphericalCoords[0][i] = sqrt(x*x + y*y + z*z);
            // now phi
            sphericalCoords[1][i] = atan2(y, x);
            // finally omega
            sphericalCoords[2][i] = PI / 2 - atan(z / sqrt(x*x + y*y));
        }

        return sphericalCoords;
    }
    
    //TODO connect this with the upper one
    static void cartes2Spherical(final double x, final double y, final double z,
            final double[] spher){

        spher[0] = sqrt(x*x+y*y+z*z);

        if(x > 0){
            spher[1] = atan(y/x);
        } else if(x == 0){
            spher[1] = signum(y)*PI/2;
        } else if(x < 0 && y >= 0){
            spher[1] = atan(y/x) + PI;
        } else if(x < 0 && y < 0){
            spher[1] = atan(y/x) - PI;
        }
   
        spher[2] = acos(z/spher[0]);
    }

    /**
     * Translates a set of spherical coordinates into a set of cartesian coordinates.
     * @param spherical spherical coordinates
     * @return cartesian coordinates as a double[][]
     */
    public static double[][] sphericalToCartesianCoord(final double[][] spherical) {

        final int noOfAtoms = spherical[0].length;
        final double[][] xyz = new double[3][noOfAtoms];

        for (int i = 0; i < noOfAtoms; i++) {

            final double r   = spherical[0][i];
            final double phi = spherical[1][i];
            final double om  = spherical[2][i];
            
            final double sinOm = Math.sin(om);
            xyz[0][i] = r * sinOm * cos(phi);
            xyz[1][i] = r * sinOm * sin(phi);
            xyz[2][i] = r * cos(om);
        }

        return xyz;
    }
    
    public static void sphericalToCartesianCoord(final double[] spherical, final double[] xyz) {

        assert(spherical != null);
        assert(spherical.length == 3);
        assert(xyz != null);
        assert(xyz.length == 3);
        
        final double r = spherical[0];
        final double phi = spherical[1];
        final double om = spherical[2];

        final double sinOm = sin(om);
        xyz[0] = r * sinOm * cos(phi);
        xyz[1] = r * sinOm * sin(phi);
        xyz[2] = r * cos(om);
    }
    

    /**
     * Sets the Kearsley matrix up.
     * @param refCoords The reference coordinate set, with the COM being 0.0,0.0,0.0.
     * @param coords The to be moved coordinate set, with the COM being 0.0,0.0,0.0.
     * @param kearsley The kearsley matrix, needs to be initialized with zeros!
     */
    private static void setupKearsleyMatrix(final double[][] refCoords, final double[][] coords, final double[][] kearsley) {
        
        assert(kearsley != null);
        assert(kearsley.length == 4);
        assert(kearsley[0].length == 4);

        final int noOfAtoms = refCoords[0].length;

        // set the actual kearsley matrix up
        final double[] diff = new double[3];
        final double[] summ = new double[3];
        for (int i = 0; i < noOfAtoms; i++) {

            for(int j = 0; j < 3; j++){
                diff[j] = refCoords[j][i] - coords[j][i];
                summ[j] = refCoords[j][i] + coords[j][i];
            }

            kearsley[0][0] += diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
            kearsley[0][1] += summ[1]*diff[2] - diff[1]*summ[2];
            kearsley[0][2] += diff[0]*summ[2] - summ[0]*diff[2];
            kearsley[0][3] += summ[0]*diff[1] - diff[0]*summ[1];

            kearsley[1][1] += summ[1]*summ[1] + summ[2]*summ[2] + diff[0]*diff[0];
            kearsley[1][2] += diff[0]*diff[1] - summ[0]*summ[1];
            kearsley[1][3] += diff[0]*diff[2] - summ[0]*summ[2];

            kearsley[2][2] += summ[0]*summ[0] + summ[2]*summ[2] + diff[1]*diff[1];
            kearsley[2][3] += diff[1]*diff[2] - summ[1]*summ[2];

            kearsley[3][3] += summ[0]*summ[0] + summ[1]*summ[1] + diff[2]*diff[2];
        }

        /*
         * only upper half is filled so far, we transpose this down now since the matrix is
         * symmetric
         */
        for (int i = 0; i < 3; i++) {
            for (int j = i+1; j < 4; j++) {
                kearsley[j][i] = kearsley[i][j];
            }
        }
    }

    /**
     * Calculate the rotation matrix from the quaternion representing the optimal superposition.
     * @param eigen The eigenvalue matrix.
     * @param rot The minimum [3,3] rotation matrix, changed on exit.
     * @return The rotation matrix.
     */
    private static void calculateRotationMatrix(final double[][] eigen, final double[][] rot) {
        
        assert(rot != null);
        assert(rot.length >= 3);
        assert(rot[0].length >= 3);

        final double q1 = eigen[0][0];
        final double q2 = eigen[1][0];
        final double q3 = eigen[2][0];
        final double q4 = eigen[3][0];
        
        rot[0][0] = q1*q1 + q2*q2 - q3*q3 - q4*q4;
        rot[0][1] = 2 * (q2*q3 + q1*q4);
        rot[0][2] = 2 * (q2*q4 - q1*q3);

        rot[1][0] = 2 * (q2*q3 - q1*q4);
        rot[1][1] = q1*q1 + q3*q3 - q2*q2 - q4*q4;
        rot[1][2] = 2 * (q3*q4 + q1*q2);

        rot[2][0] = 2 * (q2*q4 + q1*q3);
        rot[2][1] = 2 * (q3*q4 - q1*q2);
        rot[2][2] = q1*q1 + q4*q4 - q2*q2 - q3*q3;
    }
    
    /**
     * rotates a set of cartesian coordinates such that the point is on
     * the z-axis. Taken from phenix and optimized (sorry, Bernd!).
     * @param xyz set of cartesians
     * @param point the point to be rotated onto the z axis
     * @param ats number of atoms
     * @return the rotated cartesians
     */
    public static double[][] rotatePointToZAxis(final double[][] xyz, final double[] point, final int ats){
        
        assert(xyz.length == 3);
        assert(ats == xyz[0].length);
        
        // strategy: use the Euler angle rotation and determine the necessary Euler angles

        /*
         * Euler angle rotations are done with a given handedness, and a normal to a
         * plane can have two different orientations; in order to remedy the possible
         * error situation arising from this, let alpha=180-alpha and beta=180-beta
         * if the yy is negative; this is most easily achieved by changing the signs
         * of xx and zz
         */
         if (point[1] < 0.0){
            point[0] *= -1;
            point[2] *= -1;
         }
         
         /*
          * determine Euler angle alpha as angle between vector (xx,yy,zz) and 
          * "line of nodes" (terminus technicus of the Euler angle experts..., it is
          * the position of the y-axis after the first of the three rotations, and the
          * axis about which the z-axis is rotated; more simply: the line of nodes is
          * perpendicular to the plane spanned by the z-axis and the vector (xx,yy,zz))
          */
         final double pXSq = point[0]*point[0];
         final double pYSq = point[1]*point[1];
         final double pZSq = point[2]*point[2];
         final double dist2 = sqrt(pXSq+pYSq);
         final double ca = point[0]/dist2;
         final double alpha = acos(ca);
         
         /*
          * determine angle beta as angle between vector (xx,yy,zz) and current z-axis
          * (note: after the first rotation, the values of xx and yy have changed, but
          * not zz and not the length of the vector (xx,yy,zz), and this determines the
          * angle, so the first rotation does not affect this angle.)
          */
         final double dist = sqrt(pXSq+pYSq+pZSq);
         final double cb = point[2]/dist;
         final double beta = acos(cb);
         
         // third Euler angle is irrelevant
         final double gamma = 0.0;
         
         if(DEBUG) {System.out.println("DEBUG: alpha " + alpha + " beta " + beta + " gamma " + gamma);}
         
         // rotation using given Euler angles
         final double sa = sin(alpha);
         //final double ca = cos(alpha);
         final double sb = sin(beta);
         //final double cb = cos(beta);
         //final double sg = 0.0; sin(gamma);
         //final double cg = 1.0; cos(gamma);
         
         assert(abs(ca - cos(alpha)) < 1E-6);
         assert(abs(cb - cos(beta)) < 1E-6);
         assert(abs(sin(gamma)) < 1E-6);
         assert(abs(cos(gamma)-1.0) < 1E-6);
         
         // from rotdet
         final double a11 = ca * cb;// * cg - sa * sg;
         final double a12 = sa * cb;// * cg + ca * sg;
         final double a13 = -sb;// * cg;
         final double a21 = - sa;//-ca * cb * sg - sa * cg;
         final double a22 = ca;//-sa * cb * sg + ca * cg;
         //final double a23 = 0.0;//sb * sg;
         final double a31 = ca * sb;
         final double a32 = sa * sb;
         final double a33 = cb;

         final double[][] rotcoord = new double[3][ats];

         for (int i = 0; i < ats; i++) {
            rotcoord[0][i] = a11 * xyz[0][i] + a12 * xyz[1][i] + a13 * xyz[2][i];
            rotcoord[1][i] = a21 * xyz[0][i] + a22 * xyz[1][i];// + a23 * xyz[2][i];
            rotcoord[2][i] = a31 * xyz[0][i] + a32 * xyz[1][i] + a33 * xyz[2][i];
         }
                  
         /*// from rotdt2
         final double[][] tmpcoord = new double[3][ats];
         final double[][] rotcoord2 = new double[3][ats];

         for(int i = 0; i < ats; i++){
            tmpcoord[0][i]=ca*xyz[0][i]+sa*xyz[1][i];
            tmpcoord[1][i]=-sa*xyz[0][i]+ca*xyz[1][i];
            tmpcoord[2][i]=xyz[2][i];
         }
         for(int i = 0; i < ats; i++){
            rotcoord2[0][i]=cb*tmpcoord[0][i]-sb*tmpcoord[2][i];
            rotcoord2[1][i]=tmpcoord[1][i];
            rotcoord2[2][i]=sb*tmpcoord[0][i]+cb*tmpcoord[2][i];
         }*/
         
         return rotcoord;         
    }
}
