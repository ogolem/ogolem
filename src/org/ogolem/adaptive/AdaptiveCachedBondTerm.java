/**
Copyright (c) 2012, J. M. Dieterich
              2016, J. M. Dieterih and B. Hartke
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
package org.ogolem.adaptive;

import org.ogolem.core.BondInfo;
import org.ogolem.core.Gradient;
import org.ogolem.core.Topology;

/**
 * Simple harmonic potential for bonded interactions. Using specific terms makes
 * little to no sense in case that the reference cartesians are different
 * systems! Also, if one enables caching, the system is not allowed to change
 * for different invocations of the same object!
 * This makes aggressive use of caching in order to optimize performance.
 * @author Johannes Dieterich
 * @version 2016-02-24
 */
public final class AdaptiveCachedBondTerm extends AdaptiveBondedHarmonicTerm {

    private static final long serialVersionUID = (long) 20110617;

    private int[][] allBonds;

    public AdaptiveCachedBondTerm(boolean useTermsForEachPair, double distCut)
            throws Exception{
        super(useTermsForEachPair, distCut, true);
    }

    private AdaptiveCachedBondTerm(AdaptiveCachedBondTerm orig){
        super(orig);
        if(orig.allBonds != null){
            allBonds = new int[orig.allBonds.length][2];
            for(int i = 0; i < orig.allBonds.length; i++){
                allBonds[i][0] = orig.allBonds[i][0];
                allBonds[i][1] = orig.allBonds[i][1];
            }
        }
    }

    @Override
    public AdaptiveCachedBondTerm clone(){
        return new AdaptiveCachedBondTerm(this);
    }

    @Override
    public double partialInteraction(Topology topology, AdaptiveParameters params){

        final BondInfo bonds = topology.getBonds();
        final double[][] dists = topology.getDistances();
        final String[] atoms = topology.getAtomNames();
        final int iNoOfAtoms = dists.length;

        if(paramOffsetCache == null){
            int noOfBonds = 0;
            for(int i = 0; i < iNoOfAtoms-1; i++){
                for(int j = i+1; j < iNoOfAtoms; j++){
                    if(bonds.hasBond(i, j)){
                        noOfBonds++;
                    }
                }
            }
            paramOffsetCache = new int[noOfBonds];
            allBonds = new int[noOfBonds][2];
            
            // fill cache up on the fly
            int counter = -1;
            for(int i = 0; i < iNoOfAtoms-1; i++){
                for(int j = i+1; j < iNoOfAtoms; j++){
                   if(!bonds.hasBond(i, j)) continue;
                   counter++;
                   // fill cache up a little
                   final String s = buildKey(atoms, ids, i, j, params, specialIDs);
                   paramOffsetCache[counter] = params.getStartPointForKey(s);
                   allBonds[counter][0] = i;
                   allBonds[counter][1] = j;
                }
            }
        }

        final double[] p = params.getAllParamters();
        double energy = 0.0;
        for(int counter = 0; counter < allBonds.length; counter++){
                
            // use the cache
            final int offset = paramOffsetCache[counter];
            final double dist = dists[allBonds[counter][0]][allBonds[counter][1]];
            final double distX0 = dist-p[offset+1];
            
            // check wether we are so far apart that we need a cutoff
            if (distX0 >= distCutoff) {
                energy += CUTOFF;
                continue;
            }

            energy += 0.5 * p[offset] * distX0 * distX0;
        }

        // check for NaN/Infinity cases and catch them
        if(Double.isInfinite(energy) || Double.isNaN(energy)) return FixedValues.NONCONVERGEDENERGY;

        return energy;
    }

    @Override
    public Gradient partialCartesianGradient(Topology topology,
            AdaptiveParameters params){

        final int noOfAtoms = topology.getNumberOfAtoms();
        final String[] atoms = topology.getAtomNames();
        final BondInfo bonds = topology.getBonds();
        final double[][][] diffs = topology.getPosDiffs();
        final double[][] dists = topology.getDistances();

        final Gradient analyticalGradient = new Gradient();
        final double[][] grad = new double[3][noOfAtoms];

        if(paramOffsetCache == null){
            int noOfBonds = 0;
            for(int i = 0; i < noOfAtoms-1; i++){
                for(int j = i+1; j < noOfAtoms; j++){
                    if(bonds.hasBond(i, j)){
                        noOfBonds++;
                    }
                }
            }
            paramOffsetCache = new int[noOfBonds];
            allBonds = new int[noOfBonds][2];
            
            // fill cache up on the fly
            int counter = -1;
            for(int i = 0; i < noOfAtoms-1; i++){
                for(int j = i+1; j < noOfAtoms; j++){
                   if(!bonds.hasBond(i, j)) continue;
                   counter++;
                   // fill cache up a little
                   final String s = buildKey(atoms, ids, i, j, params, specialIDs);
                   paramOffsetCache[counter] = params.getStartPointForKey(s);
                   allBonds[counter][0] = i;
                   allBonds[counter][1] = j;
                }
            }
        }

        final double[] p = params.getAllParamters();
        double energy = 0.0;
        for(int counter = 0; counter < allBonds.length; counter++){
                
            // use the cache
            final int offset = paramOffsetCache[counter];
            final int i = allBonds[counter][0];
            final int j = allBonds[counter][1];
            final double dist = dists[i][j];
            final double distX0 = dist-p[offset+1];
            
            // check wether we are so far apart that we need a cutoff
            if (distX0 >= distCutoff) {
                energy += CUTOFF;
                continue;
            }

            double tempDeriv = p[offset] * (p[offset + 1] - dist) / dist;

            // normalize NaNs and infinity
            if (Double.isNaN(tempDeriv) || Double.isInfinite(tempDeriv)) {
                tempDeriv = FixedValues.NONCONVERGEDGRADIENT;
            }

            grad[0][i] -= tempDeriv * (diffs[i][j][0]);
            grad[0][j] += tempDeriv * (diffs[i][j][0]);
            grad[1][i] -= tempDeriv * (diffs[i][j][1]);
            grad[1][j] += tempDeriv * (diffs[i][j][1]);
            grad[2][i] -= tempDeriv * (diffs[i][j][2]);
            grad[2][j] += tempDeriv * (diffs[i][j][2]);
            
            energy += 0.5 * p[offset] * distX0 * distX0;
        }

        // normalize NaNs and infinity
        if(Double.isNaN(energy) || Double.isInfinite(energy)) energy = FixedValues.NONCONVERGEDENERGY;

        analyticalGradient.setGradientTotal(grad);
        analyticalGradient.setTotalEnergy(energy);

        return analyticalGradient;
    }
}
