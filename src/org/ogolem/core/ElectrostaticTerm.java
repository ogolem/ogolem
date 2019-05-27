/**
Copyright (c) 2010-2012, J. M. Dieterich
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
 * Coulomb term for usage in force fields. Computes in atomic units. Can be made
 * to correspond to classical MM implementations, where the 1-2 and 1-3
 * contributions are neglected and the 1-4 contribution is scaled. If caching is
 * enabled with scaling, the same object can only compute structures of the same
 * system!
 * @author Johannes Dieterich
 * @version 2016-02-24
 */
public final class ElectrostaticTerm implements InteractionTerm {

    private static final long serialVersionUID = (long) 20110126;
    private final boolean use14Scaling;
    private final double scale14;
    private final boolean useCaching;
    private boolean[] cache13;
    private boolean[] cache14;

    public ElectrostaticTerm(boolean employ14Scaling, double scaleFac14,
            boolean useCache){
        this.use14Scaling = employ14Scaling;
        this.scale14 = scaleFac14;
        this.useCaching = useCache;
        // initialize the caches explicitly
        this.cache13 = null;
        this.cache14 = null;
    }
    
    public ElectrostaticTerm(final ElectrostaticTerm orig){
        this.use14Scaling = orig.use14Scaling;
        this.useCaching = orig.useCaching;
        this.scale14 = orig.scale14;
        
        if(orig.cache13 != null) this.cache13 = orig.cache13.clone();
        else this.cache13 = null;
        
        if(orig.cache14 != null) this.cache14 = orig.cache14.clone();
        else this.cache14 = null;
    }

    @Override
    public double partialInteraction(final Topology topology){

        // get the distances and charges
        final float[] charges = topology.getCharges();
        final double[][] dists = topology.getDistances();
        final BondInfo bonds = topology.getBonds();
        final int iNoOfAtoms = dists.length;

        if(useCaching && cache13 == null){
            // initialize the caches
            final int noOfContr = (iNoOfAtoms*iNoOfAtoms-iNoOfAtoms)/2;
            cache13 = new boolean[noOfContr];
            cache14 = new boolean[noOfContr];
            
            // populate caches on the fly, easiest that way
            int interCounter = 0;
            for (int i = 0; i < iNoOfAtoms - 1; i++) {
                for (int j = i + 1; j < iNoOfAtoms; j++) {
                    cache13[interCounter] = topology.does13exist(i, j);
                    cache14[interCounter] = topology.does14exist(i, j);
                    interCounter++;
                }
            }
        }

        int interCount = -1;
        double dCoulombEnergy = 0.0;
        for(int i = 0; i < iNoOfAtoms-1; i++){
            for(int j = i+1; j < iNoOfAtoms; j++){
                interCount++;
                
                double scale = 1.0;
                if(this.use14Scaling){
                    
                    if(!useCaching){
                        if(bonds.hasBond(i,j)) continue;
                        else if(topology.does13exist(i, j)) continue;
                        else if(topology.does14exist(i, j)) scale = scale14;
                    } else{
                        if(bonds.hasBond(i,j)) continue;
                        else if(cache13[interCount]) continue;
                        else if(cache14[interCount]) scale = scale14;
                    }
                }

                dCoulombEnergy += scale*(charges[i]*charges[j])/dists[i][j];
            }
        }

        return dCoulombEnergy;
    }

    @Override
    public Gradient partialGradient(final Topology topology){
        
        // get the distances, charges and positions
        final float[] charges = topology.getCharges();
        final double[][] dists = topology.getDistances();
        final double[][][] diffs = topology.getPosDiffs();
        final BondInfo bonds = topology.getBonds();
        final int iNoOfAtoms = dists.length;

        if(useCaching && cache13 == null){
            // initialize the caches
            final int noOfContr = (iNoOfAtoms*iNoOfAtoms-iNoOfAtoms)/2;
            cache13 = new boolean[noOfContr];
            cache14 = new boolean[noOfContr];
            
            // populate caches on the fly, easiest that way
            int interCounter = 0;
            for (int i = 0; i < iNoOfAtoms - 1; i++) {
                for (int j = i + 1; j < iNoOfAtoms; j++) {
                    cache13[interCounter] = topology.does13exist(i, j);
                    cache14[interCounter] = topology.does14exist(i, j);
                    interCounter++;
                }
            }
        }

        final double[][] daGradient = new double[3][iNoOfAtoms];

        int interCount = -1;
        double dCoulombEnergy = 0.0;
        for(int i = 0; i < iNoOfAtoms-1; i++){
            for(int j = i+1; j < iNoOfAtoms; j++){
                interCount++;

                double scale = 1.0;
                if(this.use14Scaling){

                    if(!useCaching){
                        if(bonds.hasBond(i,j)) continue;
                        else if(topology.does13exist(i, j)) continue;
                        else if(topology.does14exist(i, j)) scale = scale14;
                    } else{
                        if(bonds.hasBond(i,j)) continue;
                        else if(cache13[interCount]) continue;
                        else if(cache14[interCount]) scale = scale14;
                    }
                }

                final double dDistX = diffs[i][j][0];
                final double dDistY = diffs[i][j][1];
                final double dDistZ = diffs[i][j][2];

                // divide the distances in all three dimensions by the total distance to cover for the coordinate system
                final double dist = dists[i][j];
                final double distInv = 1.0/dist;
                final double dDivProdX = dDistX * distInv;
                final double dDivProdY = dDistY * distInv;
                final double dDivProdZ = dDistZ * distInv;

                final double dTempDeriv = -scale*(charges[i]*charges[j])*distInv*distInv;

                daGradient[0][i] += dTempDeriv * dDivProdX;
                daGradient[0][j] -= dTempDeriv * dDivProdX;
                daGradient[1][i] += dTempDeriv * dDivProdY;
                daGradient[1][j] -= dTempDeriv * dDivProdY;
                daGradient[2][i] += dTempDeriv * dDivProdZ;
                daGradient[2][j] -= dTempDeriv * dDivProdZ;

                // the actual energy
                dCoulombEnergy += scale*(charges[i]*charges[j])* distInv;
            }
        }

        final Gradient partialGrad = new Gradient();
        partialGrad.setGradientTotal(daGradient);
        partialGrad.setTotalEnergy(dCoulombEnergy);

        return partialGrad;
    }
}
