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
 * This provides an implementation of a backend for Lennard-Jones clusters. For
 * the time being it is NOT possible for this force field to calculate
 * heterogeneous clusters. Besides, as described by the Backend interface, it
 * depends on a 1D array of coordinates for calculation.
 * @author Johannes Dieterich
 * @version 2020-03-21
 */
public class LennardJonesFF implements CartesianFullBackend {

    // the ID
    private static final long serialVersionUID = (long) 20200201;
    
    private final boolean cache;
    private double eps;
    private double sig;
    
    public LennardJonesFF(final boolean caching){
        this.cache = caching;
        this.eps = Double.NaN;
        this.sig = Double.NaN;
    }

    @Override
    public LennardJonesFF clone(){
        return new LennardJonesFF(cache);
    }

    @Override
    public String getMethodID(){
        return "LJ Force Field for homogenous cases";
    }

    /**
     * This just provides a VERY simple LJ potential energy calculation.
     * ATTENTION: This works so far just for clusters being build from identical
     * atom types and assumes that the cluster is build from the FIRST atom!
     */
    @Override
    public double energyCalculation(final long lID, final int iIteration,
            final double[] xyz1D, final String[] saAtomTypes, final short[] atomNos,
            final int[] atsPerMol, final double[] energyparts,
            final int iNoOfAtoms, final float[] faCharges, final short[] iaSpins, final BondInfo bonds) {
        
        // zero the partial contributions
        for(int i = 0; i < energyparts.length; i++) energyparts[i] = 0.0;

        // the Lennard-Jones parameters
        final double d4Epsilon;
        final double dSigma;
        if(Double.isNaN(eps)){
            d4Epsilon = 4.0*AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[0]);
            dSigma = AtomicProperties.giveLennardJonesSigma(saAtomTypes[0]);
            if(cache){
                eps = d4Epsilon;
                sig = dSigma;
            }
        } else {
            d4Epsilon = eps;
            dSigma = sig;
        }

        // the cutoff distance
        final double dSeam = 0.64 * dSigma;
        final double dSeamSquared = dSeam * dSeam;

        // more constants... needed for cutting off the potential
        final double t = 1.0/0.64;
        final double tSq = t*t;
        final double t6 = tSq*tSq*tSq;
        final double t12 = t6*t6;
        final double dConst1 = (4.0 * (t12 - t6) - 10000.0) / dSeam;
        final double dConst2 = 10000.0;

        // get the squared distances between all atoms
        double dPotEnergyAdded = 0.0;
        for (int i = 0; i < iNoOfAtoms - 1; i++) {
            final double x0 = xyz1D[i];
            final double y0 = xyz1D[i+iNoOfAtoms];
            final double z0 = xyz1D[i+2*iNoOfAtoms];
            for (int j = i + 1; j < iNoOfAtoms; j++) {
                // calculate pair distance
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

                    final double tmp = d4Epsilon * (dInvRPow12 - dInvRPow6);
                    energyparts[i] += tmp;
                    energyparts[j] += tmp;
                    dPotEnergyAdded += tmp;

                } else {
                    //System.out.println("Atoms are too close together");
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
     * ATTENTION: This works so far just for clusters being build from identical
     * atom types and assumes that the cluster is build from the FIRST atom!
     */
    @Override
    public void gradientCalculation(final long lID, final int iIteration,
            final double[] xyz1D, final String[] saAtomTypes, final short[] atomNos,
            final int[] atsPerMol, final double[] energyparts,
            final int iNoOfAtoms, final float[] faCharges, final short[] iaSpins, final BondInfo bonds,
            final Gradient gradient) {
        
        // zero the partial contributions
        for(int i = 0; i < energyparts.length; i++) energyparts[i] = 0.0;

        // the matrix for the gradient
        gradient.zeroGradient();
        final double[][] daGradientMat = gradient.getTotalGradient();

        // the Lennard-Jones parameters
        final double d4Epsilon;
        final double dSigma;
        if(Double.isNaN(eps)){
            d4Epsilon = 4.0*AtomicProperties.giveLennardJonesEpsilon(saAtomTypes[0]);
            dSigma = AtomicProperties.giveLennardJonesSigma(saAtomTypes[0]);
            if(cache){
                eps = d4Epsilon;
                sig = dSigma;
            }
        } else {
            d4Epsilon = eps;
            dSigma = sig;
        }
        final double d24Epsilon = d4Epsilon*6.0;
        final double d48Epsilon = d24Epsilon*2.0;

        // the cutoff distance
        final double dSeam = 0.64 * dSigma;
        final double dSeamSquared = dSeam*dSeam;

        // more constants... needed for cutting off the potential
        final double t = 1.0/0.64;
        final double tSq = t*t;
        final double t6 = tSq*tSq*tSq;
        final double t12 = t6*t6;
        final double dConst1 = (4.0 * (t12 - t6) - 10000.0) / dSeam;
        final double dConst2 = 10000.0;

        // calculate all pair distances
        double dPotEnergyAdded = 0.0;
        for (int i = 0; i < (iNoOfAtoms - 1); i++) {
            
            final double x0 = xyz1D[i];
            final double y0 = xyz1D[i+iNoOfAtoms];
            final double z0 = xyz1D[i+2*iNoOfAtoms];
            
            double gradXI = daGradientMat[0][i];
            double gradYI = daGradientMat[1][i];
            double gradZI = daGradientMat[2][i];
            
            for (int j = (i + 1); j < iNoOfAtoms; j++) {
                
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

                    final double tmp = d4Epsilon * (dInvRPow12 - dInvRPow6);
                    energyparts[i] += tmp;
                    energyparts[j] += tmp;
                    dPotEnergyAdded += tmp;
                    dTemp = -d48Epsilon * dInvRPow12 * dDistInv + d24Epsilon * dInvRPow6 * dDistInv;
                } else {
                    if(GlobalConfig.DEBUGLEVEL > 0) System.err.println("WARNING: LJ gradient: Atoms too close together, we take the cutoff potential.");
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
