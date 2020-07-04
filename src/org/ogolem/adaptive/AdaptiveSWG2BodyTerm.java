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
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.Topology;
import org.ogolem.helpers.Tuple3D;

/**
 * The Stillinger-Weber-Gong 2 body term.
 * @author Johannes Dieterich
 * @version 2016-02-25
 */
class AdaptiveSWG2BodyTerm  implements AdaptiveInteractionTerm {

    private static final long serialVersionUID = (long) 20160219;

    private final boolean useCaching;
    private final double blowFacClose;
    
    private int[] paramOffsetCache2b;
    
    AdaptiveSWG2BodyTerm(final boolean useCaching, final double blowFacClose){
        this.useCaching = useCaching;
        this.blowFacClose = blowFacClose;
    }
    
    AdaptiveSWG2BodyTerm(final AdaptiveSWG2BodyTerm orig){
        this.useCaching = orig.useCaching;
        if(orig.paramOffsetCache2b != null){this.paramOffsetCache2b = orig.paramOffsetCache2b.clone();}
        else {this.paramOffsetCache2b = null;}
        this.blowFacClose = orig.blowFacClose;
    }
    
    @Override
    public AdaptiveSWG2BodyTerm clone() {
        return new AdaptiveSWG2BodyTerm(this);
    }

    @Override
    public double partialInteraction(final Topology topology, final AdaptiveParameters params) {
        
        final int noOfAtoms = topology.getNumberOfAtoms();
        final double[][] dists = topology.getDistances();
        final String[] atomNames = topology.getAtomNames();
        final double[] p = params.getAllParamters();
        final short[] atomNos = topology.getAtomicNumbers();
        
        if(useCaching && paramOffsetCache2b == null){
            // initialize caches
            initializeCaches(params, atomNames);
        }
        
        int c = -1;
        double energy = 0.0;
        for(int i = 0; i < noOfAtoms-1; i++){
            if(atomNos[i] == 0){continue;} // dummy
            final double rad1 = AtomicProperties.giveRadius(atomNos[i]);
            for(int j = i+1; j < noOfAtoms; j++){
                
                if(atomNos[j] == 0){continue;} // dummy

                final double rad2 = AtomicProperties.giveRadius(atomNos[j]);
                final double addedFudged = blowFacClose*(rad1+rad2);
                
                c++;
                final double distIJ = dists[i][j];
                if(distIJ <= addedFudged){
                    energy += FixedValues.MAXTOEMERGENCY;
                    continue;
                }
                
                final int offP = (useCaching) ? paramOffsetCache2b[c] : posOfKey(atomNames, i, j, params);
                final double r0 = p[offP];
                
                if(distIJ >= r0){continue;}
                
                final double a = p[offP+1];
                final double b = p[offP+2];
                final double alpha = p[offP+3];
                
                final double rSq = distIJ*distIJ;
                final double r4 = rSq*rSq;
                
                final double e = a*(b/r4-1)*FastMath.exp(alpha/(distIJ-r0));
                
                energy += e;
            }
        }
        
        return (energy >= FixedValues.NONCONVERGEDENERGY || Double.isNaN(energy)) ? FixedValues.NONCONVERGEDENERGY : energy;
    }

    @Override
    public double partialParamGradient(final Topology topology, final AdaptiveParameters params, final double[] grad) {
        
        // XXX numerical
        final int[] startEnd = params.getStartAndEndForPrefix("adaptiveswg2b:");
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
        final double[][] xyz = topology.getPositions();
        final double[][] dists = topology.getDistances();
        final String[] atomNames = topology.getAtomNames();
        final double[] p = params.getAllParamters();
        final short[] atomNos = topology.getAtomicNumbers();
        
        if(useCaching && paramOffsetCache2b == null){
            // initialize caches
            initializeCaches(params, atomNames);
        }
        
        final Gradient gradient = new Gradient(3,noOfAtoms);
        final double[][] gradMat = gradient.getTotalGradient();
        
        int c = -1;
        double energy = 0.0;
        for(int i = 0; i < noOfAtoms-1; i++){
            
            if(atomNos[i] == 0){continue;} // dummy
            
            final double rad1 = AtomicProperties.giveRadius(atomNos[i]);
            final double posIX = xyz[0][i];
            final double posIY = xyz[1][i];
            final double posIZ = xyz[2][i];
            
            for(int j = i+1; j < noOfAtoms; j++){

                if(atomNos[j] == 0){continue;} // dummy
                
                final double rad2 = AtomicProperties.giveRadius(atomNos[j]);
                final double addedFudged = blowFacClose*(rad1+rad2);
                
                c++;
                final double distIJ = dists[i][j];
                
                if(distIJ <= addedFudged){
                    energy += FixedValues.MAXTOEMERGENCY;
                    
                    gradMat[0][i] += FixedValues.NONCONVERGEDGRADIENT;
                    gradMat[0][j] -= FixedValues.NONCONVERGEDGRADIENT;
                
                    gradMat[1][i] += FixedValues.NONCONVERGEDGRADIENT;
                    gradMat[1][j] -= FixedValues.NONCONVERGEDGRADIENT;
                
                    gradMat[2][i] += FixedValues.NONCONVERGEDGRADIENT;
                    gradMat[2][j] -= FixedValues.NONCONVERGEDGRADIENT;
                    
                    continue;
                }
                
                final int offP = (useCaching) ? paramOffsetCache2b[c] : posOfKey(atomNames, i, j, params);
                final double r0 = p[offP];
                
                if(distIJ >= r0){continue;}
                
                final double a = p[offP+1];
                final double b = p[offP+2];
                final double alpha = p[offP+3];
                
                final double rSq = distIJ*distIJ;
                final double r4 = rSq*rSq;
                
                final double rjkrrc = 1.0/(distIJ-r0);
                final double arjkr = alpha*rjkrrc;
                final double fjk = FastMath.exp(arjkr);
                
                final double bbr = b/r4-1.0;
                
                //dvdr=-aa*fjk*(4.d0*bb*r4/rjk+arjkr*rjkrrc*bbr)
                final double dvdr = -a*fjk*(4*b/r4/distIJ + arjkr*rjkrrc*bbr);
                
                final double diffX = posIX - xyz[0][j];
                final double diffY = posIY - xyz[1][j];
                final double diffZ = posIZ - xyz[2][j];
                
                final double invDist = 1/distIJ;
                
                gradMat[0][i] += dvdr*diffX*invDist;
                gradMat[0][j] -= dvdr*diffX*invDist;
                
                gradMat[1][i] += dvdr*diffY*invDist;
                gradMat[1][j] -= dvdr*diffY*invDist;
                
                gradMat[2][i] += dvdr*diffZ*invDist;
                gradMat[2][j] -= dvdr*diffZ*invDist;
                
                final double e = a*(b/r4-1)*fjk;
                
                energy += e;
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
        
        // per pair of atom types, we need 4 parameters: r0, a, b, alpha
        final int noPairs = (allAtoms.size()+1)*allAtoms.size()/2;
        final int noParams = 4*noPairs;
        
        final int[] noPerKey = new int[noPairs];
        final String[] keys = new String[noPairs];
        
        int c = 0;
        for(int i = 0; i < allAtoms.size(); i++){
            for(int j = 0; j < allAtoms.size(); j++){
                keys[c] = "adaptiveswg2b:" + allAtoms.get(i) + allAtoms.get(j);
                noPerKey[c] = 4;
                c++;
            }
        }
        
        return new Tuple3D<>(keys,noPerKey,noParams);
    }

    @Override
    public double[][] bordersForMyParams(final AdaptiveParameters params) {
        
        // not all of the parameters belong to us
        final List<String> corrKeys = params.getKeysStartingWith("adaptiveswg2b:");

        final int noKeys = corrKeys.size();
        final int noParams = noKeys*4;

        final double[][] borders = new double[2][noParams];

        int counter = 0;
        for(int i = 0; i < noKeys; i++){

            // first one: r0 XXX NOT OPTIMAL
            borders[0][counter] = 0.0;
            borders[1][counter] = 20.0;
            counter++;

            // second one: a XXX NOT OPTIMAL
            borders[0][counter] = 0.0;
            borders[1][counter] = 1000.0;
            counter++;

            // third one: b XXX NOT OPTIMAL
            borders[0][counter] = 0.0;
            borders[1][counter] = 1000.0;
            counter++;
            
            // fourth one: alpha XXX NOT OPTIMAL
            borders[0][counter] = 0.0;
            borders[1][counter] = 20.0;
            counter++;
        }

        return borders;
    }
    
    private int posOfKey(final String[] atoms, final int i, final int j, final AdaptiveParameters params){

        String key = "adaptiveswg2b:" + atoms[i] + atoms[j];
        int pos = params.getStartPointForKey(key);
        if(pos >= 0){return pos;}

        key = "adaptiveswg2b:" + atoms[j] + atoms[i];
        pos = params.getStartPointForKey(key);
        if(pos >= 0){return pos;}

        throw new RuntimeException("WARNING: No parameters for key " + key + ".");
    }
    
    private void initializeCaches(final AdaptiveParameters params, final String[] atoms){

        // how many 2b? I am too bored to calculate this now...
        int c = 0;
        for (int i = 0; i < atoms.length; i++) {
            if(atoms[i].equalsIgnoreCase("XX")){continue;} // dummy
            for (int j = 0; j < i; j++) {
                if(atoms[j].equalsIgnoreCase("XX")){continue;} // dummy
                c++;
            }
        }

        this.paramOffsetCache2b = new int[c];

        c = 0;
        for (int i = 0; i < atoms.length; i++) {
            if(atoms[i].equalsIgnoreCase("XX")){continue;} // dummy
            for (int j = 0; j < i; j++) {
                if(atoms[i].equalsIgnoreCase("XX")){continue;} // dummy
                paramOffsetCache2b[c] = posOfKey(atoms, i, j, params);
                c++;
            }
        }
    }
}
