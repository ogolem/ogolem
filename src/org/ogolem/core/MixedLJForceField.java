/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

/**
 * This provides a backend for mixed atom type Lennard-Jones calculations.
 * Uses Lorentz-Berthelot combination rules.
 * @author Johannes Dieterich
 * @version 2020-05-25
 */
public class MixedLJForceField implements CartesianFullBackend {
    
    // the ID
    private static final long serialVersionUID = (long) 20200201;

    @Override
    public MixedLJForceField clone(){
        return new MixedLJForceField();
    }

    @Override
    public String getMethodID(){
        return "LJ Force Field for heterogenous cases";
    }

    @Override
    public double energyCalculation(final long lID, final int iIteration,
            final double[] xyz1D, final String[] saAtomTypes, final short[] atomNos,
            final int[] atsPerMol, final double[] energyparts,
            final int iNoOfAtoms, final float[] faCharges, final short[] iaSpins, final BondInfo bonds,
            final boolean hasRigidEnv) {
        
        // zero the partial contributions
        for(int i = 0; i < energyparts.length; i++) energyparts[i] = 0.0;

        // some cutoff constants
        final double t1 = 1.0/0.64;
        final double t1Sq = t1*t1;
        final double t1Hex = t1Sq*t1Sq*t1Sq;
        final double t112 = t1Hex*t1Hex;
        
        // get all LJ parameters in O(N)
        final double[] eps = new double[iNoOfAtoms];
        final double[] sig = new double[iNoOfAtoms];
        for (int i = 0; i < iNoOfAtoms; i++) {
            if(atomNos[i] == 0){continue;} // dummy
            eps[i] = AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[i]);
            sig[i] = AtomicProperties.giveLennardJonesSigma(saAtomTypes[i]);
        }
        
        // get the squared distances between all atoms
        double dPotEnergyAdded = 0.0;
        
        int lastOffset = 0;
        for(int i = 0; i < atsPerMol.length-1; i++){
            lastOffset += atsPerMol[i];
        }
        
        final int firstLoopAtomNo = (hasRigidEnv) ? lastOffset - 1: iNoOfAtoms-1;
        for (int i = 0; i < firstLoopAtomNo; i++) {

            if(atomNos[i] == 0){continue;} // dummy

            final double dEpsilon1 = eps[i];
            final double dSigma1 = sig[i];
            final double x0 = xyz1D[i];
            final double y0 = xyz1D[i+iNoOfAtoms];
            final double z0 = xyz1D[i+2*iNoOfAtoms];

            for (int j = i + 1; j < iNoOfAtoms; j++) {

                if(atomNos[j] == 0){continue;} // dummy

                // the Lennard-Jones parameters, set two
                final double dEpsilon2 = eps[j];
                final double dSigma2 = sig[j];
                
                final double dEpsilon = Math.sqrt(dEpsilon1 * dEpsilon2);
                final double dSigma = 0.5 * (dSigma1 + dSigma2);

                // the cutoff distance
                final double dSeam = 0.64 * dSigma;
                final double dSeamSquared = dSeam * dSeam;

                final double dDistX = x0 - xyz1D[j];
                final double dDistY = y0 - xyz1D[j+iNoOfAtoms];
                final double dDistZ = z0 - xyz1D[j+2*iNoOfAtoms];
                final double dDistSquared = dDistX * dDistX + dDistY * dDistY + dDistZ * dDistZ;
                if (dDistSquared > dSeamSquared) {
                    final double dInvRPow2 = dSigma * dSigma / dDistSquared;
                    final double dInvRPow6 = dInvRPow2 * dInvRPow2 * dInvRPow2;
                    final double dInvRPow12 = dInvRPow6 * dInvRPow6;

                    /*
                     * The used Lennard-Jones potential is of the form
                     * V(r_{ij})= 4\cdot\epsilon\cdot \left[\left(\dfrac{\sigma}{r_{ij}}\right)^{12}-\left(\dfrac{\sigma}{r_{ij}}^{6}\right)\right].
                     */
                    final double tmp = 4.0 * dEpsilon * (dInvRPow12 - dInvRPow6);
                    energyparts[i] += tmp;
                    energyparts[j] += tmp;
                    dPotEnergyAdded += tmp;
                } else {
                    //System.out.println("Atoms are too close together");
                    
                    // more constants... needed for cutting of the potential
                    final double invSeam = 1.0/dSeam;
                    final double dConst1 = (4.0 * dEpsilon * (t112 - t1Hex) - 10000.0) *invSeam;
                    final double dConst2 = 10000.0;
                    
                    final double dDist = Math.sqrt(dDistSquared);
                    final double tmp = dConst1 * dDist + dConst2;
                    energyparts[i] += tmp;
                    energyparts[j] += tmp;
                    dPotEnergyAdded += tmp;
                }
            }
        }
        
        return dPotEnergyAdded;
    }

    /**
     * same as above in blue...
     * This provides the gradient, derived with exactly the same methods as above.
     */
    @Override
    public void gradientCalculation(final long lID, final int iIteration,
            final double[] xyz1D, final String[] saAtomTypes, final short[] atomNos,
            final int[] atsPerMol, final double[] energyparts,
            final int iNoOfAtoms, final float[] faCharges, final short[] iaSpins, final BondInfo bonds,
            final Gradient gradient, final boolean hasRigidEnv) {

        // zero the partial contributions
        for(int i = 0; i < energyparts.length; i++) energyparts[i] = 0.0;
        
        // the matrix for the gradient
        gradient.zeroGradient();
        final double[][] daGradientMat = gradient.getTotalGradient();

        // some cutoff constants
        final double t1 = 1.0/0.64;
        final double t1Sq = t1*t1;
        final double t1Hex = t1Sq*t1Sq*t1Sq;
        final double t112 = t1Hex*t1Hex;
        
        // get all LJ parameters in O(N)
        final double[] eps = new double[iNoOfAtoms];
        final double[] sig = new double[iNoOfAtoms];
        for (int i = 0; i < iNoOfAtoms; i++) {
            if(atomNos[i] == 0){continue;} // dummy
            eps[i] = AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[i]);
            sig[i] = AtomicProperties.giveLennardJonesSigma(saAtomTypes[i]);
        }
        
        // calculate all pair distances
        double dPotEnergyAdded = 0.0;
        
        int lastOffset = 0;
        for(int i = 0; i < atsPerMol.length-1; i++){
            lastOffset += atsPerMol[i];
        }
        
        final int firstLoopAtomNo = (hasRigidEnv) ? lastOffset - 1: iNoOfAtoms-1;
        for (int i = 0; i < firstLoopAtomNo; i++) {

            if(atomNos[i] == 0){continue;} // dummy
            
            // LJ parameters, set one. getting it here is of course faster
            final double dEpsilon1 = eps[i];
            final double dSigma1 = sig[i];
            final double x0 = xyz1D[i];
            final double y0 = xyz1D[i+iNoOfAtoms];
            final double z0 = xyz1D[i+2*iNoOfAtoms];

            double gradXI = daGradientMat[0][i];
            double gradYI = daGradientMat[1][i];
            double gradZI = daGradientMat[2][i];
            
            for (int j = (i + 1); j < iNoOfAtoms; j++) {

                if(atomNos[j] == 0){continue;} // dummy

                // the Lennard-Jones parameters, set two
                final double dEpsilon2 = eps[j];
                final double dSigma2 = sig[j];

                final double dEpsilon = Math.sqrt(dEpsilon1 * dEpsilon2);
                final double dSigma = 0.5 * (dSigma1 + dSigma2);

                // the cutoff distance
                final double dSeam = 0.64 * dSigma;
                final double dSeamSquared = dSeam * dSeam;
                
                final double dDistX = x0 - xyz1D[j];
                final double dDistY = y0 - xyz1D[j+iNoOfAtoms];
                final double dDistZ = z0 - xyz1D[j+2*iNoOfAtoms];
                final double dDistSquared = dDistX * dDistX + dDistY * dDistY + dDistZ * dDistZ;

                final double dDist = Math.sqrt(dDistSquared);
                final double dDistInv = 1.0 / dDist;
                // divide the distances in all three dimensions by the total distance to cover for the coordinate system
                final double dDivProdX = dDistX * dDistInv;
                final double dDivProdY = dDistY * dDistInv;
                final double dDivProdZ = dDistZ * dDistInv;


                double dTemp;
                if (dDistSquared > dSeamSquared) {
                    final double dInvRPow2 = dSigma * dSigma / dDistSquared;
                    final double dInvRPow6 = dInvRPow2 * dInvRPow2 * dInvRPow2;
                    final double dInvRPow12 = dInvRPow6 * dInvRPow6;

                    /*
                     * The used Lennard-Jones potential is of the form
                     * V(r_{ij})= 4\cdot\epsilon\cdot \left[\left(\dfrac{\sigma}{r_{ij}}\right)^{12}-\left(\dfrac{\sigma}{r_{ij}}^{6}\right)\right].
                     */

                    final double tmp = 4.0 * dEpsilon * (dInvRPow12 - dInvRPow6);
                    energyparts[i] += tmp;
                    energyparts[j] += tmp;
                    dPotEnergyAdded += tmp;
                    dTemp = -48.0 * dEpsilon * dInvRPow12 * dDistInv + 24.0 * dEpsilon * dInvRPow6 * dDistInv;
                } else {
                    if(GlobalConfig.DEBUGLEVEL > 0) System.err.println("WARNING: LJ gradient: Atoms too close together, we take the cutoff potential.");
                    
                    // more constants... needed for cutting of the potential
                    final double dConst1 = (4.0 * dEpsilon * (t112 - t1Hex) - 10000.0) / dSeam;
                    final double dConst2 = 10000.0;
                    
                    final double tmp = dConst1 * dDist + dConst2;
                    energyparts[i] += tmp;
                    energyparts[j] += tmp;
                    dPotEnergyAdded += tmp;
                    dTemp = dConst1;
                }
                gradXI += dTemp * dDivProdX;
                daGradientMat[0][j] -= dTemp * dDivProdX;
                gradYI += dTemp * dDivProdY;
                daGradientMat[1][j] -= dTemp * dDivProdY;
                gradZI += dTemp * dDivProdZ;
                daGradientMat[2][j] -= dTemp * dDivProdZ;
            }
            
            daGradientMat[0][i] = gradXI;
            daGradientMat[1][i] = gradYI;
            daGradientMat[2][i] = gradZI;
        }
        
        gradient.setTotalEnergy(dPotEnergyAdded);
    }
}
