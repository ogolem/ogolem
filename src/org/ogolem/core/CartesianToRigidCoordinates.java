/**
Copyright (c) 2015, J. M. Dieterich and B. Hartke
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
import static org.ogolem.core.RigidBodyCoordinates.MAXMOVECOMCOORD;
import org.ogolem.generic.GenericBackend;

/**
 * Turns a Cartesian backend into a rigid body one.
 * @author Johannes Dieterich
 * @version 2015-03-28
 */
public class CartesianToRigidCoordinates implements CoordinateRepresentation {
    
    private static final long serialVersionUID = (long) 20150218;
    private static final boolean DEBUG = false;
    
    private final CartesianFullBackend back;
    
    private final double[] euler = new double[3];
    private final double[][] rotMat = new double[3][3];

    private Geometry cache;
    private CartesianCoordinates cartes;
    private double[] xyz1D;
    private List<CartesianCoordinates> cartesBackup;
    private List<CartesianCoordinates> cartesians;
    private List<double[][]> rotatedCoordsCache;
    private List<double[][][]> rotatedCache;
    private double[][] coms;
    private double[][] comDisplacementCache;
    private double[] matMultCache;
    
    private Gradient gradObj;

    private double[] energyparts;
    
    CartesianToRigidCoordinates(final CartesianFullBackend cartesback){
        this.back = cartesback;
    }
    
    CartesianToRigidCoordinates(final CartesianToRigidCoordinates orig){
        this.back = orig.back.clone();
    }

    @Override
    public CartesianToRigidCoordinates clone() {
        return new CartesianToRigidCoordinates(this);
    }
    
    @Override
    public String getMyID() {
        return back.getMethodID();
    }
    
    @Override
    public void setUnderlyingBackend() throws CoordRepresentationsNotCompatibleException {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public GenericBackend<Molecule, Geometry> getUnderlyingBackend() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    
    @Override
    public int numberOfActiveCoordinates(final Geometry individual) {
        return 6*individual.getNumberOfIndieParticles();
    }
    
    @Override
    public void resetToStable(final double[] coordinates) {
        
        // basically, we rotate the cartesian backup coordinates so that the Euler parameters can be zero'd
        final int mols = cache.getNumberOfIndieParticles();
        final double[][] rotmat = new double[3][3];
        final double[] eulers = new double[3];
        for(int mol = 0; mol < mols; mol++){
            eulers[0] = coordinates[6*mol+3];
            eulers[1] = coordinates[6*mol+4];
            eulers[2] = coordinates[6*mol+5];
            
            final double[][] xyz = cartesBackup.get(mol).getAllXYZCoord();
            final double[][] rotated = CoordTranslation.rotateXYZ(xyz, eulers, rotmat);
            cartesBackup.get(mol).setAllXYZ(rotated);
            
            coordinates[6*mol+3] = 0.0;
            coordinates[6*mol+4] = 0.0;
            coordinates[6*mol+5] = 0.0;
        }
    }

    @Override
    public double[] getActiveCoordinates(final Geometry individual) {
        
        this.cache = individual;
        this.energyparts = new double[individual.getNumberOfIndieParticles()];
        final int mols = individual.getNumberOfIndieParticles();
        this.coms = new double[mols][3];
        this.comDisplacementCache = new double[mols][];
        this.cartesBackup = new ArrayList<>(individual.getNumberOfIndieParticles());
        this.cartesians = new ArrayList<>(individual.getNumberOfIndieParticles());
        this.rotatedCache = new ArrayList<>(individual.getNumberOfIndieParticles());
        this.rotatedCoordsCache = new ArrayList<>(individual.getNumberOfIndieParticles());
        this.gradObj = new Gradient(3,cartes.getNoOfAtoms());
        final double[] coords = new double[mols*6];
        
        // simply copy COM and Euler coordinates
        int maxAtPerMol = 0;
        for(int mol = 0; mol < mols; mol++){
            final Molecule m = cache.getMoleculeAtPosition(mol);
            
            // no rigidification
            final CartesianCoordinates cartes = new CartesianCoordinates(m.getNumberOfAtoms(),1,new int[]{m.getNumberOfAtoms()});
            cartes.setAllAtomTypes(m.getAtomTypes());
            cartes.recalcAtomNumbers();
            cartes.setAllCharges(m.getAllCharges());
            cartes.setAllSpins(m.getAllSpins());
            cartes.setAllXYZAsCopy(m.getReferenceCartesians());
            
            final double[] displ = cartes.calculateTheCOM(); // ideally, all zeros
            cartes.moveCoordsToCOM(); // we move these to the COM to make rotations easier and then substract the displacement in the end again
            
            maxAtPerMol = Math.max(maxAtPerMol, cartes.getNoOfAtoms());
            
            cartesians.add(cartes);
            cartesBackup.add(cartes.clone());
            
            final double[][][] rotCache = new double[3][3][cartes.getNoOfAtoms()];
            rotatedCache.add(rotCache);
            
            final double[][] rotCoordCache = new double[3][cartes.getNoOfAtoms()];
            rotatedCoordsCache.add(rotCoordCache);
            
            comDisplacementCache[mol] = displ;
            final double[] com = cache.getCOM(mol);
            coms[mol][0] = com[0] + displ[0];
            coms[mol][1] = com[1] + displ[1];
            coms[mol][2] = com[2] + displ[2];
            coords[mol*6] = com[0];
            coords[mol*6+1] = com[1];
            coords[mol*6+2] = com[2];
            final double[] euler = cache.getEulers(mol);
            coords[mol*6+3] = euler[0];
            coords[mol*6+4] = euler[1];
            coords[mol*6+5] = euler[2];
        }
        
        this.matMultCache = new double[maxAtPerMol];
        this.cartes = cache.getCartesiansWithEnvironment();
        this.xyz1D = new double[cache.getNumberOfAtoms()*3];
        
        if(DEBUG){
            System.out.println("DEBUG: Printing out active coordinates for individual " + individual.getID());
            for(final double d : coords){
                System.out.println("DEBUG: " + d);
            }
            System.out.println("DEBUG: Coordinates done");
        }
        
        return coords;
    }
    
    @Override
    public void updateActiveCoordinates(final Geometry individual, final double[] coordinates) {
        
        /*
         * the only problem could be if individual is cache and was somehow touched outside.
         * however, inside here this works since we take the cartesian coordinates but with the
         * very first energy or gradient in the rotation effectively replace the existing ones
         * with a copy.
         */
        final int mols = individual.getNumberOfIndieParticles();
        final double[][] rotmat = new double[3][3];
        for(int mol = 0; mol < mols; mol++){
            final Molecule m = individual.getMoleculeAtPosition(mol);
            final double[] com = m.getExternalCenterOfMass();
            final double[] euler = m.getOrientation();
            
            if(coordinates == null){throw new RuntimeException("Coordinates to be updated are null! Contact author(s)!");}
            if(euler == null){throw new RuntimeException("Eulers to be updated are null! Contact author(s)!");}
            if(com == null){throw new RuntimeException("COM to be updated are null! Contact author(s)!");}
            
            com[0] = coordinates[6*mol];
            com[1] = coordinates[6*mol+1];
            com[2] = coordinates[6*mol+2];
            
            final double[][] xyz = cartesBackup.get(mol).getAllXYZCoord();
            final double[] realEulers = new double[]{coordinates[6*mol+3],coordinates[6*mol+4],coordinates[6*mol+5]};
            final double[][] rotated = CoordTranslation.rotateXYZ(xyz, realEulers, rotmat);
            
            // something interesting can happen here. as the backend may add dummies (or whatever) atoms,
            // we can get an cartesian set LONGER than what the individual has. hence, this will end up as
            // an assertion error. if we use setReferenceCartesian. Ergo: we copy ourselves
            final double[][] finalXYZ = m.getReferenceCartesians();
            final int noAtoms = m.getNumberOfAtoms();
            System.arraycopy(rotated[0], 0, finalXYZ[0], 0, noAtoms);
            System.arraycopy(rotated[1], 0, finalXYZ[1], 0, noAtoms);
            System.arraycopy(rotated[2], 0, finalXYZ[2], 0, noAtoms);
            
            final double[] mEulers = m.getOrientation();
            mEulers[0] = 0.0;
            mEulers[1] = 0.0;
            mEulers[2] = 0.0;
        }
    }
    
    @Override
    public double fitness(final double[] currCoords, final int iteration) {

        if(DEBUG){
            final Geometry scr = cache.clone();
            updateActiveCoordinates(scr,currCoords);
            final String[] cart = scr.makePrintableAbsoluteCoord(true);
            System.out.println("DEBUG: AT ENERGY EVALUTION " + iteration);
            for(final String s : cart){
                System.out.println(s);
            }
        }
        
        updateState(currCoords);
        
        final double e = back.energyCalculation(cache.getID(), iteration, xyz1D, cartes.getAllAtomTypes(),
                cartes.getAllAtomNumbers(), cartes.getAllAtomsPerMol(), energyparts, cartes.getNoOfAtoms(),
                cartes.getAllCharges(), cartes.getAllSpins(), cache.getBondInfo());
        
        assert(!Double.isInfinite(e));
        assert(!Double.isNaN(e));
        
        return e;
    }
    
    @Override
    public Geometry fitness(final Geometry individual, final boolean forceOneEval) {
        
        final double[] coords = getActiveCoordinates(individual.clone());
        final double e = fitness(coords, 42);
        
        assert(!Double.isInfinite(e));
        assert(!Double.isNaN(e));
        
        individual.setFitness(e);
        
        return individual;
    }
    
    @Override
    public double gradient(final double[] currCoords, final double[] gradient, final int iteration) {
        
        /*
         * this basically uses the derivations for the pure rigid body backend and
         * bastardizes them into something that hopefully just works (TM)
         * however, I would expect some buildup of numerical noise throughout the
         * course of the optimization
         */
        
        if(DEBUG){
            final Geometry scr = cache.clone();
            updateActiveCoordinates(scr,currCoords);
            final String[] cart = scr.makePrintableAbsoluteCoord(true);
            System.out.println("DEBUG: AT GRADIENT EVALUTION " + iteration);
            for(final String s : cart){
                System.out.println(s);
            }
        }
        if(DEBUG){
            System.out.println("DEBUG: Printing out current coordinates for individual " + cache.getID() + " at gradient " + iteration);
            for(final double d : currCoords){
                System.out.println("DEBUG: " + d);
            }
            System.out.println("DEBUG: Coordinates done");
        }
        
        updateState(currCoords);
        
        if(DEBUG){
            System.out.println("DEBUG: 1D cartes coming:");
            final String[] atoms = cartes.getAllAtomTypes();
            final int noAtomsTot = cartes.getNoOfAtoms();
            for(int i = 0; i < noAtomsTot; i++){
                System.out.println(atoms[i] + "\t" + xyz1D[i]*Constants.BOHRTOANG + "\t" + xyz1D[noAtomsTot+i]*Constants.BOHRTOANG + "\t" + xyz1D[2*noAtomsTot+i]*Constants.BOHRTOANG);
            }
            System.out.println("DEBUG: 1D cartes done.");
        }
        
        back.gradientCalculation(cache.getID(), iteration, xyz1D, cartes.getAllAtomTypes(),
                cartes.getAllAtomNumbers(), cartes.getAllAtomsPerMol(), energyparts, cartes.getNoOfAtoms(),
                cartes.getAllCharges(), cartes.getAllSpins(), cache.getBondInfo(), gradObj);
        final double energy = gradObj.getTotalEnergy();
        final double[][] gradXYZ = gradObj.getTotalGradient();
        
        // zero the gradient?!
        for(int x = 0; x < gradient.length; x++){
            gradient[x] = 0.0;
        }
        
        final double[] eulers = new double[3];
        final double[][] rotMat = new double[3][3];
        int off = 0;
        for(int mol = 0; mol < cartesians.size(); mol++){
            final double[][] molXYZ = cartesBackup.get(mol).getAllXYZCoord();
            if(DEBUG){
                final int noAtomsInMol = cartesBackup.get(mol).getNoOfAtoms();
                System.out.println("DEBUG: For molecule " + mol + "the molecular gradient is:");
                for(int atom = 0; atom < noAtomsInMol; atom++){
                    System.out.println("" + atom + " " + gradXYZ[0][off+atom] + " \t " + gradXYZ[1][off+atom] + " \t "+ gradXYZ[2][off+atom]);
                    if(Double.isNaN(gradXYZ[0][off+atom]) || Double.isNaN(gradXYZ[0][off+atom]) || Double.isNaN(gradXYZ[0][off+atom])){
                        System.err.println("ERROR: Found a NaN at molecular gradient for atom " + atom);
                    }
                }
                System.out.println("DEBUG: For molecule " + mol + "the basic cartesian is:");
                for(int atom = 0; atom < noAtomsInMol; atom++){
                    System.out.println("" + atom + " " + molXYZ[0][atom] + " \t " + molXYZ[1][atom] + " \t "+ molXYZ[2][atom]);
                }
                
            }
            
            // apply the rotations as per Burnham's idea and Bernd's partial implementation
            eulers[0] = currCoords[6*mol+3];
            eulers[1] = currCoords[6*mol+4];
            eulers[2] = currCoords[6*mol+5];
            
            final int noAtomsInMol = cartesBackup.get(mol).getNoOfAtoms();
            
            final double[][][] thisRotCache = rotatedCache.get(mol);
            final double[][] molXYZ_dPhi = thisRotCache[0];
            CoordTranslation.rotateXYZ_dphi(molXYZ, eulers, rotMat, molXYZ_dPhi, matMultCache);
            final double[][] molXYZ_dOmega = thisRotCache[1];
            CoordTranslation.rotateXYZ_domega(molXYZ, eulers, rotMat,molXYZ_dOmega, matMultCache);
            final double[][] molXYZ_dPsi = thisRotCache[2];
            CoordTranslation.rotateXYZ_dpsi(molXYZ, eulers, rotMat,molXYZ_dPsi, matMultCache);
            
            /* NOTE TO MYSELF:
             * molXYZ: must be correct (energy is right, transl gradient is right)
             * molGrad: likely be right (transl. gradient is right, unrotated Euler is right, phi is right)
             * order of rotateYXZ_dXXX calls does not influence result
             * first or second mol has non-zero Euler: does not matter
             * changing the rotation matrix and its derivative to another yaw pitch roll definition does not change anything
             * all points to the culprit being with the second and third Euler but it simply cannot be
             * or: cartesBackup gets somehow manipulated in such a fashion that it only destroys two Euler gradient elements. nothign else...
             */
                        
            // according to bernd, this is NOT a matrix multiplication but an elementwise multiplication, i.e., a dot product
            double gradCoordPhi = 0.0;
            double gradCoordOmega = 0.0;
            double gradCoordPsi = 0.0;
            for(int coord = 0; coord < 3; coord++){
                for(int atom = 0; atom < noAtomsInMol; atom++){
                    final double gradE = gradXYZ[coord][off+atom];
                    gradCoordPhi += gradE*molXYZ_dPhi[coord][atom];
                    gradCoordOmega += gradE*molXYZ_dOmega[coord][atom];
                    gradCoordPsi += gradE*molXYZ_dPsi[coord][atom];
                }
            }
            gradient[6*mol+3] = gradCoordPhi;
            gradient[6*mol+4] = gradCoordOmega;
            gradient[6*mol+5] = gradCoordPsi;
                        
            // translational part: simple summation of the UNMODIFIED gradient elements
            for(int coord = 0; coord < 3; coord++){
                double gradCoord = 0.0;
                for(int atom = 0; atom < noAtomsInMol; atom++){
                    gradCoord += gradXYZ[coord][off+atom];
                }
                gradient[6*mol+coord] = gradCoord;
            }
            
            off += noAtomsInMol;
        }
                
        return energy;
    }
    
    static double[][] getEulers(final double[] currCoords, final int mols){
        
        final double[][] eulers = new double[mols][3];
        for(int mol = 0; mol < mols; mol++){
            eulers[mol][0] = currCoords[mol*6+3];
            eulers[mol][1] = currCoords[mol*6+4];
            eulers[mol][2] = currCoords[mol*6+5];
        }
        
        return eulers;
    }
    
    static double[][] getCOMs(final double[] currCoords, final int mols){
        
        final double[][] coms = new double[mols][3];
        for(int mol = 0; mol < mols; mol++){
            coms[mol][0] = currCoords[mol*6];
            coms[mol][1] = currCoords[mol*6+1];
            coms[mol][2] = currCoords[mol*6+2];
        }
        
        return coms;
    }
    
    private void updateState(final double[] currCoords){

        final int mols = cache.getNumberOfIndieParticles();
        final int totAtoms = cache.getNumberOfAtoms();
        
        int off = 0;
        for(int mol = 0; mol < mols; mol++){
            // update the coms
            coms[mol][0] = currCoords[mol*6];
            coms[mol][1] = currCoords[mol*6+1];
            coms[mol][2] = currCoords[mol*6+2];
    
            // rotate the coords with respect to the Euler angles
            final CartesianCoordinates cartTmp = cartesians.get(mol);
            euler[0] = currCoords[mol*6+3];
            euler[1] = currCoords[mol*6+4];
            euler[2] = currCoords[mol*6+5];
            final double[][] rotated = rotatedCoordsCache.set(mol, cartTmp.getAllXYZCoord()); // replace the cache object
            CoordTranslation.rotateXYZ(cartesBackup.get(mol).getAllXYZCoord(), euler, rotMat, rotated, matMultCache);
            cartTmp.setAllXYZ(rotated);
            
            // copy cartesian coordinates to xyz1D
            for(int atom = 0; atom < rotated[0].length; atom++){
                xyz1D[           off+atom] = rotated[0][atom] + coms[mol][0];
                xyz1D[  totAtoms+off+atom] = rotated[1][atom] + coms[mol][1];
                xyz1D[2*totAtoms+off+atom] = rotated[2][atom] + coms[mol][2];
            }
            
            off += rotated[0].length;
        }
    }
    
    List<CartesianCoordinates> getCartesCacheCopy(){
        // just to make sure that no funky business with internal state happens
        // OK, this is NOT making a change to the current issue, however, as we
        // only use this functionality for numerical gradients and for unit tests
        // which are anyways both slow: leave it in!
        final List<CartesianCoordinates> cartCopy = new ArrayList<>();
        for(final CartesianCoordinates cartes : cartesians){
            cartCopy.add(cartes.clone());
        }
        
        return cartCopy;
    }
    
    List<CartesianCoordinates> getCartesBackupCopy(){
        // just to make sure that no funky business with internal state happens
        // OK, this is NOT making a change to the current issue, however, as we
        // only use this functionality for numerical gradients and for unit tests
        // which are anyways both slow: leave it in!
        final List<CartesianCoordinates> cartCopy = new ArrayList<>();
        for(final CartesianCoordinates cartes : cartesBackup){
            cartCopy.add(cartes.clone());
        }
        
        return cartCopy;
    }
    
    double[][] getCOMCache(){
        return coms;
    }
    
    CartesianFullBackend getMyBackend(){
        return back;
    }

    @Override
    public BOUNDSTYPE boundariesInRepresentation(final Geometry individual) {
        return BOUNDSTYPE.SOME;
    }

    @Override
    public void bestEstimateBoundaries(final double[] currCoords, final double[] low, final double[] high) {
        
        // first the COM coordinates, then the Euler ones
        assert(low.length == high.length);
        assert(currCoords.length == low.length);
        assert(currCoords.length%6 == 0);
        
        for(int coord = 0; coord < currCoords.length; coord+=6){
            
            low[coord] = currCoords[coord] - MAXMOVECOMCOORD;
            low[coord+1] = currCoords[coord+1] - MAXMOVECOMCOORD;
            low[coord+2] = currCoords[coord+2] - MAXMOVECOMCOORD;
            high[coord] = currCoords[coord] + MAXMOVECOMCOORD;
            high[coord+1] = currCoords[coord+1] + MAXMOVECOMCOORD;
            high[coord+2] = currCoords[coord+2] + MAXMOVECOMCOORD;
            
            // now the Eulers
            low[coord+4] = -Math.PI;
            low[coord+5] = -0.5*Math.PI;
            low[coord+6] = -Math.PI;
            high[coord+4] = Math.PI;
            high[coord+5] = 0.5*Math.PI;
            high[coord+6] = Math.PI;
        }
    }
}
