/**
Copyright (c) 2014, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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
import org.ogolem.generic.GenericBackend;

/**
 * A coordinate set for rigid body optimizations.
 * @author Johannes Dieterich
 * @version 2020-02-09
 */
public class RigidBodyCoordinates implements CoordinateRepresentation {
    
    //XXX be better behaved for single atom molecules!
    
    public static final double MAXMOVECOMCOORD = 10.0; // in bohr
    
    private static final long serialVersionUID = (long) 20200209;
    private static final boolean DEBUG = false;

    private final RigidBodyBackend backend;
    private final double[] euler = new double[3];

    private Geometry cache;
    private List<CartesianCoordinates> cartesBackup;
    private List<CartesianCoordinates> cartesians;
    private List<double[][]> rotatedCoordsCache;
    private List<double[][][]> rotatedCache;
    private double[][] coms;
    private List<double[][]> gradCache;
    private double[][] comDisplacementCache;
    
    public RigidBodyCoordinates(final RigidBodyBackend back){
        this.backend = back;
    }
    
    public RigidBodyCoordinates(final RigidBodyCoordinates orig){
        this.backend = orig.backend.clone();
    }
    
    @Override
    public String getMyID() {
        return "rigid body backend: " + backend.getMethodID();
    }
    
    @Override
    public GenericBackend<Molecule, Geometry> clone() {
        return new RigidBodyCoordinates(this);
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
        if(!backend.suitableForThisBackend(individual)){
            System.err.println("Chosen backend " + backend.getMethodID() + " not suitable for geometry!");
            return 0;
        }
        
        return individual.getNumberOfIndieParticles()*6; // six coordinates per molecule
    }

    @Override
    public double[] getActiveCoordinates(final Geometry individual) {
        
        final int mols = individual.getNumberOfIndieParticles();
        this.cache = individual;//.clone();
        this.coms = new double[mols][3];
        this.comDisplacementCache = new double[mols][];
        this.cartesBackup = new ArrayList<>(individual.getNumberOfIndieParticles());
        this.cartesians = new ArrayList<>(individual.getNumberOfIndieParticles());
        this.gradCache = new ArrayList<>(individual.getNumberOfIndieParticles());
        this.rotatedCache = new ArrayList<>(individual.getNumberOfIndieParticles());
        this.rotatedCoordsCache = new ArrayList<>(individual.getNumberOfIndieParticles());
        final double[] coords = new double[mols*6];
        int maxAtPerMol = 0;
        for(int mol = 0; mol < mols; mol++){
            final Molecule m = cache.getMoleculeAtPosition(mol);
            backend.rigidify(m);
            final CartesianCoordinates cartes = new CartesianCoordinates(m.getNumberOfAtoms(),1,new int[]{m.getNumberOfAtoms()});
            cartes.setAllAtomTypes(m.getAtomTypes());
            cartes.recalcAtomNumbers();
            cartes.setAllCharges(m.getAllCharges());
            cartes.setAllSpins(m.getAllSpins());
            cartes.setAllXYZAsCopy(m.getReferenceCartesians());
            
            final double[] displ = cartes.calculateTheCOM(); // ideally, all zeros
            cartes.moveCoordsToCOM(); // we move these to the COM to make rotations easier and then substract the displacement in the end again
            
            final CartesianCoordinates adjCartes= backend.adjustCartesians(cache, cartes, mol);
            
            maxAtPerMol = Math.max(maxAtPerMol, adjCartes.getNoOfAtoms());
            
            cartesians.add(adjCartes);
            cartesBackup.add(adjCartes.clone());
            
            final double[][][] rotCache = new double[3][3][adjCartes.getNoOfAtoms()];
            rotatedCache.add(rotCache);
            
            final double[][] rotCoordCache = new double[3][adjCartes.getNoOfAtoms()];
            rotatedCoordsCache.add(rotCoordCache);
            
            final double[][] gradPart = new double[3][adjCartes.getNoOfAtoms()];
            gradCache.add(gradPart);
            
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
    public void resetToStable(final double[] coordinates){
        
        // basically, we rotate the cartesian backup coordinates so that the Euler parameters can be zero'd
        final int mols = cache.getNumberOfIndieParticles();
        final double[] eulers = new double[3];
        for(int mol = 0; mol < mols; mol++){
            eulers[0] = coordinates[6*mol+3];
            eulers[1] = coordinates[6*mol+4];
            eulers[2] = coordinates[6*mol+5];
            
            final double[][] xyz = cartesBackup.get(mol).getAllXYZCoord();
            final double[][] rotated = CoordTranslation.rotateXYZ(xyz, eulers);
            cartesBackup.get(mol).setAllXYZ(rotated);
            
            coordinates[6*mol+3] = 0.0;
            coordinates[6*mol+4] = 0.0;
            coordinates[6*mol+5] = 0.0;
        }
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
            final double[][] rotated = CoordTranslation.rotateXYZ(xyz, realEulers);
            
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
        return backend.energy(cache, cartesians, coms, iteration);
    }

    @Override
    public Geometry fitness(final Geometry individual, final boolean forceOneEval) {
        
        final double[] optCoord = getActiveCoordinates(individual.clone());
        final double e = fitness(optCoord, 0);
        
        individual.setFitness(e);
        
        return individual;
    }
    
    @Override
    public double gradient(final double[] currCoords, final double[] gradient, final int iteration) {
        
        if(false){
            return NumericalGradients.gradientForRigidBody(cache, cartesBackup, getCOMs(currCoords,cache.getNumberOfIndieParticles()),
                    getEulers(currCoords,cache.getNumberOfIndieParticles()), gradient, backend, iteration);
        }
        
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
        final double energy = backend.gradient(cache, cartesians, coms, gradCache, iteration);
        
        // zero the gradient?!
        for(int x = 0; x < gradient.length; x++){
            gradient[x] = 0.0;
        }
        
        final double[] eulers = new double[3];
        for(int mol = 0; mol < cartesians.size(); mol++){
            final double[][] molGrad = gradCache.get(mol);
            final double[][] molXYZ = cartesBackup.get(mol).getAllXYZCoord();
            if(DEBUG){
                final int noAtomsInMol = cartesBackup.get(mol).getNoOfAtoms();
                System.out.println("DEBUG: For molecule " + mol + "the molecular gradient is:");
                for(int atom = 0; atom < noAtomsInMol; atom++){
                    System.out.println("" + atom + " " + molGrad[0][atom] + " \t " + molGrad[1][atom] + " \t "+ molGrad[2][atom]);
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
            final double[][] molXYZ_dOmega = thisRotCache[1];
            final double[][] molXYZ_dPsi = thisRotCache[2];
            CoordTranslation.rotateXYZ_dphi_domega_dpsi(molXYZ, eulers, molXYZ_dPhi, molXYZ_dOmega, molXYZ_dPsi);
                        
            // rotational parts:
            // according to bernd, this is NOT a matrix multiplication but an elementwise multiplication, i.e., a dot product
            double gradCoordPhi = 0.0;
            double gradCoordOmega = 0.0;
            double gradCoordPsi = 0.0;
            // translational part: simple summation of the UNMODIFIED gradient elements
            for(int coord = 0; coord < 3; coord++){
                double gradCoord = 0.0;
                for(int atom = 0; atom < noAtomsInMol; atom++){
                    final double gradXYZE = molGrad[coord][atom];
                    gradCoordPhi += gradXYZE*molXYZ_dPhi[coord][atom];
                    gradCoordOmega += gradXYZE*molXYZ_dOmega[coord][atom];
                    gradCoordPsi += gradXYZE*molXYZ_dPsi[coord][atom];
                    gradCoord += gradXYZE;
                }
                gradient[6*mol+coord] = gradCoord;
            }
            gradient[6*mol+3] = gradCoordPhi;
            gradient[6*mol+4] = gradCoordOmega;
            gradient[6*mol+5] = gradCoordPsi;
        }
        
        if(DEBUG){
            final double[] gradNum = new double[gradient.length];
            final double e =  NumericalGradients.gradientForRigidBody(cache, cartesBackup, getCOMs(currCoords,cache.getNumberOfIndieParticles()),
                    getEulers(currCoords,cache.getNumberOfIndieParticles()), gradient, backend, iteration);
            final double NUMPREC = 1.0e-7;
            boolean issue = false;
            if(Math.abs(e-energy) > NUMPREC){
                System.err.println("DEBUG: Difference in energies " + e + " " + energy);
                issue = true;
            }
            for(int i = 0; i < gradient.length; i++){
                if(Math.abs(gradNum[i]-gradient[i]) > NUMPREC){
                    System.err.println("DEBUG: Difference in gradients at " + i + " numerical " + gradNum[i] + " analytical " + gradient[i]);
                    issue = true;
                }
            }
            if(issue){
                //throw new RuntimeException("Failure in RigidBodyCoordinates, see previous output");
            }
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
        
        for(int mol = 0; mol < mols; mol++){
            // update the coms
            coms[mol][0] = currCoords[mol*6];
            coms[mol][1] = currCoords[mol*6+1];
            coms[mol][2] = currCoords[mol*6+2];
    
            // rotate the coords with respect to the Euler angles
            final CartesianCoordinates cartes = cartesians.get(mol);
            euler[0] = currCoords[mol*6+3];
            euler[1] = currCoords[mol*6+4];
            euler[2] = currCoords[mol*6+5];
            final double[][] rotated = rotatedCoordsCache.set(mol, cartes.getAllXYZCoord()); // replace the cache object
            CoordTranslation.rotateXYZ(cartesBackup.get(mol).getAllXYZCoord(), euler, rotated);
            cartes.setAllXYZ(rotated);
        }
    }
    
    List<CartesianCoordinates> getCartesCacheCopy(){
        // just to make sure that no funky business with internal state happens
        // OK, this is NOT making a change to the current issue, however, as we
        // only use this functionality for numerical gradients and for unit tests
        // which are anyways both slow: leave it in!
        final List<CartesianCoordinates> cartCopy = new ArrayList<>();
        cartesians.forEach((cartes) -> {
            cartCopy.add(cartes.clone());
        });
        
        return cartCopy;
    }
    
    List<CartesianCoordinates> getCartesBackupCopy(){
        // just to make sure that no funky business with internal state happens
        // OK, this is NOT making a change to the current issue, however, as we
        // only use this functionality for numerical gradients and for unit tests
        // which are anyways both slow: leave it in!
        final List<CartesianCoordinates> cartCopy = new ArrayList<>();
        cartesBackup.forEach((cartes) -> {
            cartCopy.add(cartes.clone());
        });
        
        return cartCopy;
    }
    
    double[][] getCOMCache(){
        return coms;
    }
    
    RigidBodyBackend getMyBackend(){
        return backend;
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
