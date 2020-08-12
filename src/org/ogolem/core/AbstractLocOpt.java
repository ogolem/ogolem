/**
Copyright (c) 2011-2013, J. M. Dieterich
              2016, J. M. Dieterich and B. Hartke
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
 * A blueprint for all local optimization implementations.
 * Should be extended and overidden upon demand.
 * @author Johannes Dieterich
 * @version 2020-07-29
 */
public abstract class AbstractLocOpt implements Newton{
    
    private static final long serialVersionUID = (long) 20200729;
    protected final boolean DEBUG;
    private static long noOfGeomLocalOpts = (long) 0;
    private static long noOfMolLocalOpts = (long) 0;
    protected final boolean doSanityCheck;
    protected final boolean doCollDetect;
    protected final boolean doDissDetect;
    protected final double blowBonds;
    protected final double blowBondsEnv;
    protected final double blowDissoc;
    
    /**
     * The constructor.
     * @param conf The full global config, just to be sure we get everything.
     */
    public AbstractLocOpt(final GlobalConfig conf){
        this.doSanityCheck = conf.doPostSanityCheck;
        this.doCollDetect = conf.doPostCD; 
        this.blowBonds = conf.blowFacBondDetect;
        this.blowBondsEnv = conf.blowFacEnvClusterClashes;
        this.doDissDetect = conf.doPostDD;
        this.blowDissoc = conf.blowFacDissocDetect;
        this.DEBUG = (GlobalConfig.DEBUGLEVEL > 0);
    }
    
    public AbstractLocOpt(final AbstractLocOpt orig){
        this.doSanityCheck = orig.doSanityCheck;
        this.blowBonds = orig.blowBonds;
        this.blowBondsEnv = orig.blowBondsEnv;
        this.doCollDetect = orig.doCollDetect;
        this.doDissDetect = orig.doDissDetect;
        this.blowDissoc = orig.blowDissoc;
        this.DEBUG = orig.DEBUG;
    }
    
    protected abstract String myID();
    
    @Override
    public abstract AbstractLocOpt clone();

    @Override
    public Geometry localOptimization(final Geometry gStart){
        
        incrGeom();
        
        CartesianCoordinates cartes;
        boolean bWithEnv;
        if(gStart.containsEnvironment()){
            cartes = gStart.getCartesiansWithEnvironment();
            bWithEnv = true;
        } else{
            cartes = gStart.getCartesians();
            bWithEnv = false;
        }

        final ZMatrix[] zmatVec = cartes.getZMatrices();
        final Environment refEnv = cartes.getReferenceEnvironmentCopy();
        try{
            cartes = cartesToCartes(gStart.getID(), cartes,
                    gStart.getAllConstraintsXYZ(bWithEnv), gStart.isThereAConstraint(), gStart.getBondInfo());
        } catch(Exception e){
            System.err.println("WARNING: Error in " + myID() + " locopt. Returning non-optimized geometry. " + e.toString());
            if(DEBUG) {
                e.printStackTrace(System.err);
            }
            gStart.setFitness(FixedValues.NONCONVERGEDENERGY);
            
            return gStart;
        }

        // we anyway set the old zmatrices first in, so that we have something in there
        cartes.setZMatrices(zmatVec);

        for(int i = 0; i < zmatVec.length; i++){
            if(zmatVec[i] != null){
                final CartesianCoordinates cartesTemp = cartes.giveMolecularCartes(i, true);
                zmatVec[i] = cartesTemp.calculateZMatrix();
            } // else: nothing needs to happen, unflexible molecule
        }

        cartes.setZMatrices(zmatVec);
        cartes.setRefEnvironment(refEnv);
        cartes.recalcAtomNumbers();

        if (doSanityCheck) {
            // check the cartesian for sanity, also check potentially for collisions
            final boolean sanity = GeometrySanityCheck.checkSanity(cartes, gStart.getBondInfo(), blowBonds, blowBondsEnv, doCollDetect, doDissDetect, blowDissoc);

            if (!sanity) {
                // something's insane ;-)
                if(DEBUG){
                    System.out.println("DEBUG: Something insane in geometry " + gStart.getID());
                }
                gStart.setLocalOptimized(true);
                gStart.setFitness(FixedValues.NONCONVERGEDENERGY);

                return gStart;
            }
        }

        final Geometry gEnd = new Geometry(cartes, gStart.getID(), gStart.getNumberOfIndieParticles(),
                cartes.getAllAtomsPerMol(), gStart.getAllFlexies(), gStart.getExplicitDoFs(),
                gStart.getAllConstraints(false), gStart.getAllConstraintsXYZ(false), gStart.getSIDs(),
                gStart.getBondInfo().clone());
        
        gEnd.setFitness(cartes.getEnergy());
        gEnd.setFather(gStart.getFatherID());
        gEnd.setMother(gStart.getMotherID());
        gEnd.setLocalOptimized(true);
        
        return gEnd;
    }

    @Override
    public Molecule localOptimization(final Molecule mStart) {
        
        incrMol();
        
        CartesianCoordinates cartes = mStart.getCartesians();
        final ZMatrix[] zmatVec = cartes.getZMatrices();
        try{
            // since we can't be certain that there is something really "unique" in the molecules ID, we hash it :-)
            final int hashCode = mStart.hashCode();
            final BondInfo molBonds = null; //TODO this is obviously wrong
            cartes = cartesToCartes((long) hashCode, cartes, mStart.getConstraints(),
                    mStart.isConstricted(), molBonds);
        } catch(Exception e){
            System.err.println("WARNING: Error in " + myID() + " locopt. Returning non-optimized molecule. " + e.toString());
            mStart.setEnergy(FixedValues.NONCONVERGEDENERGY);
            if(DEBUG) {
                e.printStackTrace(System.err);
            }
            
            return mStart;
        }

        // we anyway set the old zmatrix first in, so that we have something in there
        cartes.setZMatrices(zmatVec);

       
        if (zmatVec[0] != null) {
            final CartesianCoordinates cartesTemp = cartes.giveMolecularCartes(0, true);
            zmatVec[0] = cartesTemp.calculateZMatrix();
        } // else: nothing needs to happen, unflexible molecule

        cartes.setZMatrices(zmatVec);
        
        final Molecule mEnd = new Molecule (cartes, mStart.getMolPosition(), mStart.getSID(),
                mStart.getFlexy(), mStart.getDegreesOfFreedom(), mStart.isConstricted(), mStart.getConstraints());
        mEnd.setID(mStart.getMolPosition());
        mEnd.setSID(mStart.getSID());
        
        return mEnd;
    }

    @Override
    public abstract CartesianCoordinates cartesToCartes(long lID, CartesianCoordinates cartes,
            boolean[][] baConstraints, boolean isConstricted, final BondInfo bonds) throws Exception;
    
    @Override
    public CartesianFullBackend getBackend(){
        return null;
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
