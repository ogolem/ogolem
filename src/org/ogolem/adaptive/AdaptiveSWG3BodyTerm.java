/**
Copyright (c) 2016, J. M. Dieterich and B. Hartke
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

import java.util.ArrayList;
import java.util.List;
import org.apache.commons.math3.util.FastMath;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.Topology;
import org.ogolem.helpers.Tuple3D;

/**
 * A Stillinger-Weber-Gong 3 body term.
 * @author Johannes Dieterich
 * @version 2016-02-23
 */
class AdaptiveSWG3BodyTerm implements AdaptiveInteractionTerm {

    private static final long serialVersionUID = (long) 20160219;
    
    private final boolean useCaching;
    private int[] paramOffsetCache2b1;
    private int[][] paramOffsetCache3b;
    
    AdaptiveSWG3BodyTerm(final boolean useCaching){
        this.useCaching = useCaching;
    }
    
    AdaptiveSWG3BodyTerm(final AdaptiveSWG3BodyTerm orig){
        this.useCaching = orig.useCaching;
        
        if(orig.paramOffsetCache2b1 != null){this.paramOffsetCache2b1 = orig.paramOffsetCache2b1.clone();}
        else {this.paramOffsetCache2b1 = null;}
        
        
        if(orig.paramOffsetCache3b != null){
            this.paramOffsetCache3b[0] = orig.paramOffsetCache3b[0].clone();
            this.paramOffsetCache3b[1] = orig.paramOffsetCache3b[1].clone();
            this.paramOffsetCache3b[2] = orig.paramOffsetCache3b[2].clone();
        }
        else {this.paramOffsetCache3b = null;}
    }
    
    @Override
    public AdaptiveSWG3BodyTerm clone() {
        return new AdaptiveSWG3BodyTerm(this);
    }

    @Override
    public double partialInteraction(final Topology topology, final AdaptiveParameters params) {
        
        final int noOfAtoms = topology.getNumberOfAtoms();
        
        if(noOfAtoms < 3){return 0.0;} // 3 body!
                
        final double[][] dists = topology.getDistances();
        final String[] atomNames = topology.getAtomNames();
        final double[] p = params.getAllParamters();
        final short[] atomNos = topology.getAtomicNumbers();
        
        if(useCaching && paramOffsetCache2b1 == null){
            // initialize caches
            initializeCaches(params, atomNames);
        }
        
        final double oneThird = 1.0/3.0;
        
        int counter;
        double energy = 0.0;
        
        int off2b1;
        int off2b2;
        int off2b3;
        int off3b;
        counter = -1;
        int c = -1;

        for (int k = 0; k < noOfAtoms; k++) {
            if(atomNos[k] == 0){continue;} // dummy
            for (int j = 0; j < k; j++) {

                if(atomNos[j] == 0){continue;} // dummy
                
                final double distJK = dists[j][k];
                
                c++;
                if (!useCaching) {
                    off2b1 = posOfKey2sin3(atomNames, k, j, params);
                } else {
                    off2b1 = paramOffsetCache2b1[c];
                }

                // ok, check whether the distance is above the equilibrium distance for this pair
                final double r02b1 = p[off2b1];
                if (distJK > r02b1) {
                    continue;
                }

                final double dRJKRCC = 1.0 / (distJK - r02b1);
                final double gam2b1 = p[off2b1 + 1];
                final double dGRJKR = gam2b1 * dRJKRCC;
                final double dFJK = FastMath.exp(dGRJKR);

                for (int i = 0; i < j; i++) {

                    if(atomNos[i] == 0){continue;} // dummy
                    
                    counter++;
                    
                    final double distIK = dists[i][k];
                    final double distIJ = dists[i][j];

                    if (!useCaching) {
                        off3b = posOfKey3in3(atomNames, i, j, k, params);
                        off2b2 = posOfKey2sin3(atomNames, i, j, params);
                        off2b3 = posOfKey2sin3(atomNames, i, k, params);
                    } else {
                        off3b = paramOffsetCache3b[2][counter];
                        off2b2 = paramOffsetCache3b[0][counter];
                        off2b3 = paramOffsetCache3b[1][counter];
                    }

                    // ok, check whether the distance is above the equilibrium distance for this pair
                    final double r02b2 = p[off2b2];
                    final double r02b3 = p[off2b3];
                    if (distIJ > r02b2) {
                        continue;
                    }
                    if (distIK > r02b3) {
                        continue;
                    }

                    final double gam2b2 = p[off2b2 + 1];
                    final double gam2b3 = p[off2b3 + 1];
                    
                    // now we can finally start to evaluate the contributions
                    final double dFIJ = FastMath.exp(gam2b2 / (distIJ - r02b2));
                    final double dFIK = FastMath.exp(gam2b3 / (distIK - r02b3));

                    final double dCTIJK = (distIJ * distIJ + distJK * distJK
                            - distIK * distIK) / (2.0 * distIJ * distJK);
                    final double dCTIKJ = (distIK * distIK + distJK * distJK
                            - distIJ * distIJ) / (2.0 * distIK * distJK);
                    final double dCTJIK = (distIJ * distIJ + distIK * distIK
                            - distJK * distJK) / (2.0 * distIJ * distIK);

                    final double c0 = p[off3b];
                    final double c1 = p[off3b+1];
                    final double lam = p[off3b+2];
                    
                    final double t1 = dCTJIK + oneThird;
                    final double dXJIK = t1 * t1;
                    final double t2 = dCTJIK + c0;
                    final double dYJIK = t2 * t2 + c1;

                    final double t3 = dCTIJK + oneThird;
                    final double dXIJK = t3 * t3;
                    final double t4 = dCTIJK + c0;
                    final double dYIJK = t4 * t4 + c1;

                    final double t5 = dCTIKJ + oneThird;
                    final double dXIKJ = t5 * t5;
                    final double t6 = dCTIKJ + c0;
                    final double dYIKJ = t6 * t6 + c1;

                    energy += lam * (dFIJ * dFIK * dXJIK * dYJIK + dFIJ * dFJK * dXIJK * dYIJK + dFIK * dFJK * dXIKJ * dYIKJ);
                }
            }
        }

        return (energy >= FixedValues.NONCONVERGEDENERGY || Double.isNaN(energy)) ? FixedValues.NONCONVERGEDENERGY : energy;
    }

    @Override
    public double partialParamGradient(final Topology topology, final AdaptiveParameters params, final double[] grad) {
        
        // XXX numerical
        final int[] startEnd = params.getStartAndEndForPrefix("adaptiveswg3b:");
        final int start = startEnd[0];
        final int end = startEnd[1];
        
        final ParameterGradient gr = NumericalGradients.calculateParamGrad(topology, params, this, start, end);
        final double[] grData = gr.getTotalGradient();
        for(int i = start; i < end; i++){
            grad[i] += grData[i];
        }
        
        return gr.getTotalEnergy();
    }

    @Override
    public Gradient partialCartesianGradient(final Topology topology, final AdaptiveParameters params) {
        
        final int noOfAtoms = topology.getNumberOfAtoms();
        
        final Gradient gradient = new Gradient(3,noOfAtoms);
        
        if(noOfAtoms < 3){
            // 3 body!
            return gradient;
        }

        final double[][] xyz = topology.getPositions();
        final double[][] dists = topology.getDistances();
        final String[] atomNames = topology.getAtomNames();
        final double[] p = params.getAllParamters();
        final short[] atomNos = topology.getAtomicNumbers();
        
        if(useCaching && paramOffsetCache2b1 == null){
            // initialize caches
            initializeCaches(params, atomNames);
        }
        
        final double[][] gradMat = gradient.getTotalGradient();
        
        final double oneThird = 1.0/3.0;
        
        int counter;
        double energy = 0.0;
        
        int off2b1;
        int off2b2;
        int off2b3;
        int off3b;
        counter = -1;
        int c = -1;

        for (int k = 0; k < noOfAtoms; k++) {
            
            if(atomNos[k] == 0){continue;} // dummy
            
            final double posKX = xyz[0][k];
            final double posKY = xyz[1][k];
            final double posKZ = xyz[2][k];
            
            for (int j = 0; j < k; j++) {
                
                if(atomNos[j] == 0){continue;} // dummy

                final double distJK = dists[j][k];
                final double posJX = xyz[0][j];
                final double posJY = xyz[1][j];
                final double posJZ = xyz[2][j];
                final double diffJKX = posJX-posKX;
                final double diffJKY = posJY-posKY;
                final double diffJKZ = posJZ-posKZ;
                
                c++;
                if (!useCaching) {
                    off2b1 = posOfKey2sin3(atomNames, k, j, params);
                } else {
                    off2b1 = paramOffsetCache2b1[c];
                }

                // ok, check whether the distance is above the equilibrium distance for this pair
                final double r02b1 = p[off2b1];
                if (distJK > r02b1) {
                    continue;
                }

                final double dRJKRCC = 1.0 / (distJK - r02b1);
                final double gam2b1 = p[off2b1 + 1];
                final double dGRJKR = gam2b1 * dRJKRCC;
                final double dFJK = FastMath.exp(dGRJKR);
                
                for (int i = 0; i < j; i++) {

                    if(atomNos[i] == 0){continue;} // dummy
                    
                    counter++;
                    
                    final double distIK = dists[i][k];
                    final double distIJ = dists[i][j];
                    
                    final double posIX = xyz[0][i];
                    final double posIY = xyz[1][i];
                    final double posIZ = xyz[2][i];
                    final double diffIKX = posIX-posKX;
                    final double diffIKY = posIY-posKY;
                    final double diffIKZ = posIZ-posKZ;
                    final double diffIJX = posIX-posJX;
                    final double diffIJY = posIY-posJY;
                    final double diffIJZ = posIZ-posJZ;

                    if (!useCaching) {
                        off3b = posOfKey3in3(atomNames, i, j, k, params);
                        off2b2 = posOfKey2sin3(atomNames, i, j, params);
                        off2b3 = posOfKey2sin3(atomNames, i, k, params);
                    } else {
                        off3b = paramOffsetCache3b[2][counter];
                        off2b2 = paramOffsetCache3b[0][counter];
                        off2b3 = paramOffsetCache3b[1][counter];
                    }

                    // ok, check whether the distance is above the equilibrium distance for this pair
                    final double r02b2 = p[off2b2];
                    final double r02b3 = p[off2b3];
                    if (distIJ > r02b2) {
                        continue;
                    }
                    if (distIK > r02b3) {
                        continue;
                    }

                    final double gam2b2 = p[off2b2 + 1];
                    final double gam2b3 = p[off2b3 + 1];
                    
                    // now we can finally start to evaluate the contributions
                    final double dRIJRCC = 1.0 / (distIJ-r02b2);
                    final double dGRIJR = gam2b2 * dRIJRCC;
                    final double dFIJ = FastMath.exp(dGRIJR);

                    final double dRIKRCC = 1.0 / (distIK-r02b3);
                    final double dGRIKR = gam2b3 * dRIKRCC;
                    final double dFIK = FastMath.exp(dGRIKR);

                    final double dCTIJK = (distIJ * distIJ + distJK * distJK
                            - distIK * distIK) / (2.0 * distIJ * distJK);
                    final double dCTIKJ = (distIK * distIK + distJK * distJK
                            - distIJ * distIJ) / (2.0 * distIK * distJK);
                    final double dCTJIK = (distIJ * distIJ + distIK * distIK
                            - distJK * distJK) / (2.0 * distIJ * distIK);

                    final double c0 = p[off3b];
                    final double c1 = p[off3b+1];
                    final double lam = p[off3b+2];
                    
                    final double dXJIK = dCTJIK + oneThird;
                    final double dXJIK2 = dXJIK * dXJIK;

                    final double dYJIK = dCTJIK + c0;
                    final double dYJIK2 = dYJIK * dYJIK + c1;

                    final double dXIJK = dCTIJK + oneThird;
                    final double dXIJK2 = dXIJK * dXIJK;

                    final double dYIJK = dCTIJK + c0;
                    final double dYIJK2 = dYIJK * dYIJK + c1;

                    final double dXIKJ = dCTIKJ + oneThird;
                    final double dXIKJ2 = dXIKJ * dXIKJ;

                    final double dYIKJ = dCTIKJ + c0;
                    final double dYIKJ2 = dYIKJ * dYIKJ + c1;

                    final double dV31 = dFIJ * dFIK * dXJIK2 * dYJIK2;
                    final double dV32 = dFIJ * dFJK * dXIJK2 * dYIJK2;
                    final double dV33 = dFIK * dFJK * dXIKJ2 * dYIKJ2;

                    final double dPJIK = lam * dFIJ * dFIK * 2.0 * dXJIK * (dYJIK2 + dXJIK * dYJIK) * distJK / distIJ / distIK;
                    final double dPIJK = lam * dFIJ * dFJK * 2.0 * dXIJK * (dYIJK2 + dXIJK * dYIJK) * distIK / distIJ / distJK;
                    final double dPIKJ = lam * dFIK * dFJK * 2.0 * dXIKJ * (dYIKJ2 + dXIKJ * dYIKJ) * distIJ / distIK / distJK;

                    final double dDJIKIJ = -lam * dV31 * dGRIJR * dRIJRCC + dPJIK * dCTIJK;
                    final double dDJIKIK = -lam * dV31 * dGRIKR * dRIKRCC + dPJIK * dCTIKJ;
                    final double dDJIKJK = -dPJIK;

                    final double dDIJKIJ = -lam * dV32 * dGRIJR * dRIJRCC + dPIJK * dCTJIK;
                    final double dDIJKIK = -dPIJK;
                    final double dDIJKJK = -lam * dV32 * dGRJKR * dRJKRCC + dPIJK * dCTIKJ;

                    final double dDIKJIJ = -dPIKJ;
                    final double dDIKJIK = -lam * dV33 * dGRIKR * dRIKRCC + dPIKJ * dCTJIK;
                    final double dDIKJJK = -lam * dV33 * dGRJKR * dRJKRCC + dPIKJ * dCTIJK;

                    final double dDSumIJ = dDJIKIJ + dDIJKIJ + dDIKJIJ;
                    final double dDSumIK = dDJIKIK + dDIJKIK + dDIKJIK;
                    final double dDSumJK = dDJIKJK + dDIJKJK + dDIKJJK;

                    final double xxrIJ = diffIJX/distIJ;
                    final double xxrIK = diffIKX/distIK;
                    final double xxrJK = diffJKX/distJK;
                    
                    gradMat[0][i] +=  dDSumIJ * xxrIJ + dDSumIK * xxrIK;
                    gradMat[0][j] += -dDSumIJ * xxrIJ + dDSumJK * xxrJK;
                    gradMat[0][k] += -dDSumIK * xxrIK - dDSumJK * xxrJK;

                    final double yyrIJ = diffIJY/distIJ;
                    final double yyrIK = diffIKY/distIK;
                    final double yyrJK = diffJKY/distJK;
                    
                    gradMat[1][i] +=  dDSumIJ * yyrIJ + dDSumIK * yyrIK;
                    gradMat[1][j] += -dDSumIJ * yyrIJ + dDSumJK * yyrJK;
                    gradMat[1][k] += -dDSumIK * yyrIK - dDSumJK * yyrJK;
                    
                    final double zzrIJ = diffIJZ/distIJ;
                    final double zzrIK = diffIKZ/distIK;
                    final double zzrJK = diffJKZ/distJK;
                    
                    gradMat[2][i] +=  dDSumIJ * zzrIJ + dDSumIK * zzrIK;
                    gradMat[2][j] += -dDSumIJ * zzrIJ + dDSumJK * zzrJK;
                    gradMat[2][k] += -dDSumIK * zzrIK - dDSumJK * zzrJK;

                    energy += lam * (dV31 + dV32 + dV33);                    
                }
            }
        }

        if(energy < FixedValues.NONCONVERGEDENERGY && !Double.isNaN(energy)){
            gradient.setTotalEnergy(energy);
        } else {
            gradient.setTotalEnergy(FixedValues.NONCONVERGEDENERGY);
        }
        
        return gradient;
    }

    @Override
    public Tuple3D<String[], int[], Integer> requiredParams(final ArrayList<CartesianCoordinates> cartesians, final ArrayList<Topology> topologies, final String sMethod) {
        
        // make a list of all the atom types
        final List<String> allAtoms = new ArrayList<>();
        for(final Topology topo : topologies){
            final String[] atoms = topo.getAtomNames();
            for(int i = 0; i < atoms.length; i++){
                if(!allAtoms.contains(atoms[i]) && !atoms[i].equalsIgnoreCase("XX")){
                    allAtoms.add(atoms[i]);
                }
            }
        }
        
        int noPairs = 0;
        int noTriples = 0;
        for(int k = 0; k < allAtoms.size(); k++){
            for(int j = k; j < allAtoms.size(); j++){
                noPairs++;
                for(int i = j; i < allAtoms.size(); i++){
                    noTriples++;
                }
            }
        }
        
        // for each pair, we need 2 parameters, for each triple we need 3 parameters
        final int noParams = 2*noPairs+3*noTriples;
        
        final int[] noPerKey = new int[noPairs+noTriples];
        final String[] keys = new String[noPairs+noTriples];
        
        int c = 0;
        for(int k = 0; k < allAtoms.size(); k++){
            for(int j = k; j < allAtoms.size(); j++){
                noPerKey[c] = 2;
                keys[c] = "adaptiveswg3b2b:" + allAtoms.get(k) + allAtoms.get(j);
                c++;
                for(int i = j; i < allAtoms.size(); i++){
                    noPerKey[c] = 3;
                    keys[c] = "adaptiveswg3b:" + allAtoms.get(k) + allAtoms.get(j) + allAtoms.get(i);
                    c++;
                }
            }
        }
        
        return new Tuple3D<>(keys,noPerKey,noParams);
    }

    @Override
    public double[][] bordersForMyParams(final AdaptiveParameters params) {
        
        // not all of the parameters belong to us
        final List<String> corrKeys2B = params.getKeysStartingWith("adaptiveswg3b2b:");

        final int noKeys2B = corrKeys2B.size();
        final int noParams2B = noKeys2B*2;
        
        final List<String> corrKeys3B = params.getKeysStartingWith("adaptiveswg3b:");

        final int noKeys3B = corrKeys3B.size();
        final int noParams3B = noKeys3B*3;

        final double[][] borders = new double[2][noParams2B+noParams3B];

        int counter = 0;
        for(int i = 0; i < noKeys2B; i++){

            // first one: r0 XXX NOT OPTIMAL
            borders[0][counter] = 0.0;
            borders[1][counter] = 20.0;
            counter++;

            // second one: gamma XXX NOT OPTIMAL
            borders[0][counter] = 0.0;
            borders[1][counter] = 20.0;
            counter++;
        }
        
        for(int i = 0; i < noKeys2B; i++){

            // first one: c0 XXX NOT OPTIMAL
            borders[0][counter] = -10.0;
            borders[1][counter] = 10.0;
            counter++;

            // second one: c1 XXX NOT OPTIMAL
            borders[0][counter] = -10.0;
            borders[1][counter] = 10.0;
            counter++;
            
            // second one: lambda XXX NOT OPTIMAL
            borders[0][counter] = 0.0;
            borders[1][counter] = 1000.0;
            counter++;
        }

        return borders;
    }
    
    private int posOfKey2sin3(final String[] atoms, final int i, final int j, final AdaptiveParameters params){

        String key = "adaptiveswg3b2b:" + atoms[i] + atoms[j];
        int pos = params.getStartPointForKey(key);
        if(pos >= 0){return pos;}

        key = "adaptiveswg3b2b:" + atoms[j] + atoms[i];
        pos = params.getStartPointForKey(key);
        if(pos >= 0){return pos;}

        throw new RuntimeException("WARNING: No parameters for key " + key + ".");
    }

    private int posOfKey3in3(final String[] atoms, final int i, final int j, final int k, final AdaptiveParameters params){

        String key = "adaptiveswg3b:" + atoms[i] + atoms[j] + atoms[k];
        int pos = params.getStartPointForKey(key);
        if(pos >= 0){return pos;}

        key = "adaptiveswg3b:" + atoms[i] + atoms[k] + atoms[j];
        pos = params.getStartPointForKey(key);
        if(pos >= 0){return pos;}

        key = "adaptiveswg3b:" + atoms[j] + atoms[i] + atoms[k];
        pos = params.getStartPointForKey(key);
        if(pos >= 0){return pos;}

        key = "adaptiveswg3b:" + atoms[j] + atoms[k] + atoms[i];
        pos = params.getStartPointForKey(key);
        if(pos >= 0){return pos;}

        key = "adaptiveswg3b:" + atoms[k] + atoms[i] + atoms[j];
        pos = params.getStartPointForKey(key);
        if(pos >= 0){return pos;}

        key = "adaptiveswg3b:" + atoms[k] + atoms[j] + atoms[i];
        pos = params.getStartPointForKey(key);
        if(pos >= 0){return pos;}

        // we are through with all possibilities :-(
        throw new RuntimeException("WARNING: No parameters for key " + key + ".");
    }
    
    private void initializeCaches(final AdaptiveParameters params, final String[] atoms){

        // how many 3b? I am too bored to calculate this now...
        int counter = 0;
        int c = 0;
        for (int i = 0; i < atoms.length; i++) {
            if(atoms[i].equalsIgnoreCase("XX")){continue;} // dummy
            for (int j = 0; j < i; j++) {
                if(atoms[j].equalsIgnoreCase("XX")){continue;} // dummy
                c++;
                for (int k = 0; k < j; k++) {
                    if(atoms[k].equalsIgnoreCase("XX")){continue;} // dummy
                    counter++;
                }
            }
        }

        this.paramOffsetCache3b = new int[3][counter];
        this.paramOffsetCache2b1 = new int[c];

        counter = 0;
        c = 0;
        for (int i = 0; i < atoms.length; i++) {
            if(atoms[i].equalsIgnoreCase("XX")){continue;} // dummy
            for (int j = 0; j < i; j++) {
                if(atoms[j].equalsIgnoreCase("XX")){continue;} // dummy
                paramOffsetCache2b1[c] = posOfKey2sin3(atoms, i, j, params);
                c++;
                for (int k = 0; k < j; k++) {
                    if(atoms[k].equalsIgnoreCase("XX")){continue;} // dummy
                    // 2b part
                    paramOffsetCache3b[0][counter] = posOfKey2sin3(atoms, k, j, params);
                    paramOffsetCache3b[1][counter] = posOfKey2sin3(atoms, k, i, params);
                    // 3b part
                    paramOffsetCache3b[2][counter] = posOfKey3in3(atoms, i, j, k, params);
                    counter++;
                }
            }
        }
    }
}
