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

/**
 * Wraps a force field engine written in native code using the JNI.
 * The force fields come (most likely) in the doomed Fortran language, so we 
 * put another C layer in between (the actual wrapper) which will be the one 
 * and only interface to the OGOLEM java code.
 * @author Johannes Dieterich
 * @version 2016-09-03
 */
class FFEngineWrapper implements CartesianFullBackend, Newton {

    // the ID
    private static final long serialVersionUID = (long) 20101108;
    private static final boolean DEBUG = false;

    private static long noOfGeomLocalOpts = (long) 0;
    private static long noOfMolLocalOpts = (long) 0;
    
    /*
     * The JNI is a two-sided sword! Just primitives are directly handed over,
     * objects are handed over as pointers, which is of course NOT GOOD!
     * Therefore, a conversion into primitives needs to take place here.
     */

    private final int ffID;
    private final double[] grad1D;
    private final double cutE;

    /*
     * Constructors
     * first one for the external local optimization
     * second for gradient calculation
     */
    FFEngineWrapper(final GlobalConfig globConf, final int ffID, final double cutE) {
        this.ffID = ffID;
        this.grad1D = null;
        this.cutE = cutE;
    }

    FFEngineWrapper(int iWhichBackend, int noAtoms, double cutE) {
        this.ffID = iWhichBackend;
        this.grad1D = new double[noAtoms*3];
        this.cutE = cutE;
    }

    /*
     * Methods
     */

    @Override
    public FFEngineWrapper clone(){
        return new FFEngineWrapper(this.ffID, this.grad1D.length/3, this.cutE);
    }

    @Override
    public String getMethodID(){
        return "FF engine wrapper";
    }

    @Override
    public String myIDandMethod(){
        return "FF engine wrapper";
    }

    @Override
    public double energyCalculation(long lID, int iIteration, double[] daXYZ1D,
            String[] saAtomTypes, short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges, short[] iaSpins,
            final BondInfo bonds) {
        
        assert(atomNos.length == grad1D.length/3);
        
        // actually, this might work (despite the name) for molecules as well, depending on the FF backend
        double energy = FixedValues.NONCONVERGEDENERGY;
        try{
            energy = CalculateGeomEnergy(lID, iIteration, ffID, atomNos, daXYZ1D, iNoOfAtoms, energyparts, atsPerMol);
        } catch(Exception e){
            e.printStackTrace(System.err);
        } catch (Error err){
            err.printStackTrace(System.err);
        }
        
        if(Double.isNaN(energy) || Double.isInfinite(energy) || energy >= FixedValues.NONCONVERGEDENERGY || energy < cutE) return FixedValues.NONCONVERGEDENERGY;
        
        return energy;
    }

    @Override
    public void gradientCalculation(long lID, int iIteration, double[] daXYZ1D,
            String[] saAtomTypes, short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges, short[] iaSpins,
            final BondInfo bonds, final Gradient gradient) {
        
        // actually, this might work (despite the name) for molecules as well, depending on the FF backend
        
        if(DEBUG) System.out.println("DEBUG: no of atoms " + atomNos.length);
        
        assert(atomNos.length == grad1D.length/3);
        assert(energyparts.length == atsPerMol.length);
        
        // clean the scrach gradient
        for(int i = 0; i < grad1D.length; i++) grad1D[i] = 0.0;
        
        double dEnergy =  FixedValues.NONCONVERGEDENERGY;
        try{
            dEnergy = CalculateGeomGradient(lID, iIteration, ffID, atomNos, daXYZ1D, grad1D, iNoOfAtoms, energyparts, atsPerMol);
        } catch(Exception e){
            e.printStackTrace(System.err);
        } catch (Error err){
            err.printStackTrace(System.err);
        }
        
        if(Double.isNaN(dEnergy) || Double.isInfinite(dEnergy)){
            if(DEBUG) System.out.println("DEBUG: Correcting NaN or Inf in FFEngineWrapper gradient energy...");
            //dEnergy = FixedValues.NONCONVERGEDENERGY;
            //for(int i = 0; i < grad1D.length; i++) grad1D[i] = 0.0;
            throw new RuntimeException("NaN is bad in energy of gradient in FFEngineWrapper!");
        }
        
        if(dEnergy >= FixedValues.NONCONVERGEDENERGY){
            if(DEBUG) System.out.println("DEBUG: Correcting NONCONV in FFEngineWrapper gradient energy...");
            dEnergy = FixedValues.NONCONVERGEDENERGY;
            for(int i = 0; i < grad1D.length; i++) grad1D[i] = 0.0;
            //throw new RuntimeException("NaN is bad in energy of gradient in FFEngineWrapper!");
        } else if(dEnergy <= cutE){
            System.out.println("WARNING: Unreasonable energy " + dEnergy + " detected. Assuming problem and correcting.");
            dEnergy = FixedValues.NONCONVERGEDENERGY;
        }
        
        for(int i = 0; i < grad1D.length; i++){
            if(Double.isInfinite(grad1D[i]) || Double.isNaN(grad1D[i])){
                if(DEBUG) System.out.println("DEBUG: Correcting NaN or Inf in FFEngineWrapper gradient...");
                grad1D[i] = 0.0;
                throw new RuntimeException("NaN is bad in gradient in FFEngineWrapper!");
            }
        }
        
        if(DEBUG) System.out.println("DEBUG: Energy reported is " + dEnergy);
                
        final double[][] daCartesGrad = gradient.getTotalGradient();
        System.arraycopy(grad1D, 0, daCartesGrad[0], 0, iNoOfAtoms);
        System.arraycopy(grad1D, iNoOfAtoms, daCartesGrad[1], 0, iNoOfAtoms);
        System.arraycopy(grad1D, iNoOfAtoms * 2, daCartesGrad[2], 0, iNoOfAtoms);
        
        gradient.setTotalEnergy(dEnergy);
        
        if(DEBUG){
            // numerical gradient
            final Gradient numGrad = NumericalGradients.numericalGradient(lID, iIteration, daXYZ1D, saAtomTypes, atomNos, atsPerMol, energyparts, iNoOfAtoms, faCharges, iaSpins, bonds, this);
            // compare
            final double numE = numGrad.getTotalEnergy();
            final double[][] g = numGrad.getTotalGradient();
            if(numE != dEnergy) System.err.println("DEBUG: Energy mismatch! " + (numE-dEnergy));
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < iNoOfAtoms; j++){
                    final double diff = daCartesGrad[i][j]-g[i][j];
                    if(Math.abs(diff) > 1E-5) System.out.println("DEBUG: Difference between numerical and analytical found: " + diff);
                    else System.out.println("DEBUG: No difference, showing analytical: " + daCartesGrad[i][j]);
                }
            }
        }
    }

    @Override
    public Geometry localOptimization(Geometry gStartGeom) {
        incrGeom();
        CartesianCoordinates cartes = gStartGeom.getCartesians();
        int iNoOfAtoms = cartes.getNoOfAtoms();
        double[] da1DPos = cartes.getAll1DCartes();
        short[] iaAtomicNum = cartes.getAllAtomNumbers();
        da1DPos = OptimiseMolComplete(ffID, iaAtomicNum, da1DPos);
        cartes.setAll1DCartes(da1DPos, iNoOfAtoms);
        Geometry gEndGeom = new Geometry(cartes, gStartGeom.getID(), gStartGeom.getNumberOfIndieParticles(),
                cartes.getAllAtomsPerMol(), gStartGeom.getAllFlexies(), gStartGeom.getExplicitDoFs(), gStartGeom.getAllConstraints(true),
                gStartGeom.getAllConstraintsXYZ(true), gStartGeom.getSIDs(), gStartGeom.getBondInfo().clone());
        gEndGeom.setFitness(cartes.getEnergy());
        gEndGeom.setFather(gStartGeom.getFatherID());
        gEndGeom.setMother(gStartGeom.getMotherID());

        return gEndGeom;
    }

    @Override
    public Molecule localOptimization(Molecule mStartMolecule) {
        incrMol();
        CartesianCoordinates cartes = mStartMolecule.getCartesians();
        int iNoOfAtoms = cartes.getNoOfAtoms();
        double[] da1DPos = cartes.getAll1DCartes();
        short[] iaAtomicNum = cartes.getAllAtomNumbers();
        da1DPos = OptimiseGeomComplete(ffID, iaAtomicNum, da1DPos);
        cartes.setAll1DCartes(da1DPos, iNoOfAtoms);
        Molecule mEndMolecule = new Molecule(cartes, mStartMolecule.getMolPosition(), mStartMolecule.getSID(),
                mStartMolecule.getFlexy(), mStartMolecule.getDegreesOfFreedom(), mStartMolecule.isConstricted(),
                mStartMolecule.getConstraints());
        mEndMolecule.setID(mStartMolecule.getMolPosition());
        mEndMolecule.setSID(mStartMolecule.getSID());

        return mEndMolecule;
    }

    @Override
    public CartesianCoordinates cartesToCartes(long id, CartesianCoordinates cartes,
            boolean[][] constraints, boolean isConstricted, final BondInfo bonds) throws Exception{
        throw new Exception("Method not implemented in External FF Caller.");
    }
    
    @Override
    public CartesianFullBackend getBackend(){
        return this; //TODO might not be always correct...
    }

    /*
     * The respective native methods
     */
    private native double CalculateGeomEnergy(long lID, int iIteration, int iWhichFF, short[] iaAtomicNumbers, double[] da1DAtomicPos, int noAtoms, double[] energyparts, int[] atsPerMol);

    private native double CalculateGeomGradient(long lID, int iIteration, int iWhichFF, short[] iaAtomicNumbers, double[] da1DAtomicPos, double[] gradient, int noAtoms, double[] energyparts, int[] atsPerMol);

    private native double[] OptimiseGeomComplete(int iWhichFF, short[] iaAtomicNumbers, double[] da1DAtomicPos);

    private native double[] OptimiseMolComplete(int iWhichFF, short[] iaAtomicNumbers, double[] da1DAtomicPos);

    /*
     * Loading of external library libFFEngineWrapper.so (for UNIX) written in C
     */
    static {
        System.loadLibrary("FFEngineWrapper");
    }
    
    private static synchronized void incrGeom(){
        noOfGeomLocalOpts++;
    }
    
    private static synchronized void incrMol(){
        noOfMolLocalOpts++;
    }
    
    @Override
    public long getNumberOfGeomLocalOpts(){
        return noOfGeomLocalOpts;
    }
    
    @Override
    public long getNumberOfMolLocalOpts(){
        return noOfMolLocalOpts;
    }
}
