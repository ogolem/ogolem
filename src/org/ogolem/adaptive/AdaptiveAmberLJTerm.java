/**
Copyright (c) 2011-2014, J. M. Dieterich
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
package org.ogolem.adaptive;

import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.LinkedList;
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.Topology;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.io.InputPrimitives;

/**
 * A term for the slightly modified AMBER version of the classic LJ(6,12,6)
 * potential. In case that caching should be enabled, the term may only be invoked
 * on the exact same systems!
 * Also calculates the electrostatic contribution as this safes some CPU cycles.
 * @author Johannes Dieterich
 * @version 2016-02-24
 */
final class AdaptiveAmberLJTerm implements AdaptiveInteractionTerm {

    private static final long serialVersionUID = (long) 20110126;
    private final static double CLOSECUT = 10;
    private final boolean useSpecialIDs;
    private final double blowDist;
    private final double blowClose;
    private final double scaleFac14;
    private final double scaleFac14E;
    private final String[] ids;
    private final boolean useCaching;
    private int[] paramPosCache;
    private boolean[] is13Cache;
    private boolean[] is14Cache;

    AdaptiveAmberLJTerm(boolean specialIDs, double distBlow, double closeBlow,
            double scale14, double scale14E, boolean useCache) throws Exception{
        this.useSpecialIDs = specialIDs;
        this.blowDist = distBlow;
        this.blowClose = closeBlow;
        this.scaleFac14 = scale14;
        this.scaleFac14E = scale14;
        this.useCaching = useCache;

        if(specialIDs){
            // read special IDs in
            final String aux = "ids-amberlj.aux";
            try{
                final String[] sa = InputPrimitives.readFileIn(aux);
                ids = new String[sa.length];
                for(int i = 0; i < ids.length; i++){
                    final String[] sa2 = sa[i].split("\\s+");
                    ids[i] = sa2[1];
                }
            } catch(Exception e){
                throw new Exception("Couldn't use ids-amberlj.aux.",e);
            }
        } else{
            this.ids = null;
        }

        // initialize the caches anyways
        this.paramPosCache = null;
        this.is13Cache = null;
        this.is14Cache = null;
    }

    private AdaptiveAmberLJTerm(AdaptiveAmberLJTerm orig){
        this.useSpecialIDs = orig.useSpecialIDs;
        this.blowDist = orig.blowDist;
        this.blowClose = orig.blowClose;
        this.scaleFac14 = orig.scaleFac14;
        this.scaleFac14E = orig.scaleFac14E;
        this.useCaching = orig.useCaching;

        if(orig.ids != null) this.ids = orig.ids.clone();
        else this.ids = null;

        if(orig.paramPosCache != null) this.paramPosCache = orig.paramPosCache.clone();
        else this.paramPosCache = null;

        if(orig.is13Cache != null) this.is13Cache = orig.is13Cache.clone();
        else this.is13Cache = null;

        if(orig.is14Cache != null) this.is14Cache = orig.is14Cache.clone();
        else this.is14Cache = null;
    }

    @Override
    public AdaptiveAmberLJTerm clone(){
        return new AdaptiveAmberLJTerm(this);
    }

    @Override
    public double partialInteraction(Topology topology, AdaptiveParameters params){

        final String[] atoms = topology.getAtomNames();
        final short[] nos = topology.getAtomicNumbers();
        final double[][] dists = topology.getDistances();
        final BondInfo bonds = topology.getBonds();
        final float[] charges = topology.getCharges();
        final int noOfAtoms = atoms.length;
        final double[] allParams = params.getAllParamters();
        double[] dscr = new double[2];

        boolean freshCache = false;
        if(useCaching && paramPosCache == null){
            final int noOfContr = (noOfAtoms*noOfAtoms-noOfAtoms)/2;

            // make new caches
            paramPosCache = new int[noOfAtoms];
            is13Cache = new boolean[noOfContr];
            is14Cache = new boolean[noOfContr];
            freshCache = true;
        }

        // populate caches, easiest that way
        if(useCaching && freshCache){
            for (int i = 0; i < noOfAtoms; i++) {
                final String s1;
                if (useSpecialIDs) s1 = ids[i];
                else s1 = atoms[i];

                paramPosCache[i] = params.getStartPointForKey("amberlj:" + s1);
            }

            int interCounter = 0;
            for (int i = 0; i < noOfAtoms - 1; i++) {
                for (int j = i + 1; j < noOfAtoms; j++) {
                    is13Cache[interCounter] = topology.does13exist(i, j);
                    is14Cache[interCounter] = topology.does14exist(i, j);
                    interCounter++;
                }
            }
        }

        int interCounter = -1;
        double energy = 0.0;
        for(int i = 0; i < noOfAtoms-1; i++){

            final double rad1 = AtomicProperties.giveRadius(nos[i]);

            for(int j = i+1; j < noOfAtoms; j++){

                interCounter++;
                final double dist = dists[i][j];
                
                final double rad2 = AtomicProperties.giveRadius(nos[j]);
                
                // check if we even need to compute (or scale)
                double scale = 1.0;
                double scaleE = 1.0;
                if(useCaching){
                    // 1-2 contribution: continue
                    if(bonds.hasBond(i, j)) continue;
                    // 1-3 contribution: continue
                    else if(is13Cache[interCounter]) continue;
                    // 1-4 contribution: scale
                    else if(is14Cache[interCounter]){
                        scale = scaleFac14;
                        scaleE = scaleFac14E;
                    }
                } else{
                    if(bonds.hasBond(i, j)) continue;
                    else if(topology.does13exist(i, j)) continue;
                    else if(topology.does14exist(i, j)){
                        scale = scaleFac14;
                        scaleE = scaleFac14E;
                    }
                }
                
                // always calculate the electrostatics
                final double distInv = 1.0/dist;
                energy += scaleE*(charges[i]*charges[j])*distInv;

                // cutoffs for the LJ
                if(dist > blowDist*(rad1+rad2)) continue;
                else if(dist < blowClose*(rad1+rad2)){
                    energy += CLOSECUT;
                    continue;
                }

                // get params
                if(useCaching) getParams(allParams, paramPosCache[i], paramPosCache[j], dscr);
                else getParams(params, atoms, ids, i, j, useSpecialIDs, dscr);
                
                if(dscr == null){
                    //System.err.println("WARNING: No (reasonable) AmberLJ parameters "
                    //        + "for " + atoms[i] + i + " and " + atoms[j] + j + ". Adding nonconverged energy.");
                    energy += FixedValues.NONCONVERGEDENERGY;
                    dscr = new double[2];
                    continue;
                }

                // compute
                final double radTerm = dscr[1]*distInv;
                final double radTerm2 = radTerm*radTerm;
                final double radTerm6 = radTerm2*radTerm2*radTerm2;
                final double radTerm12 = radTerm6*radTerm6;
                energy += scale*dscr[0]*(radTerm12-2*radTerm6);
            }
        }

        return energy;
    }

    @Override
    public double partialParamGradient(Topology topology, AdaptiveParameters params, final double[] grad){

        final int[] startEnd = params.getStartAndEndForPrefix("amberlj:");
        if(startEnd == null) return FixedValues.NONCONVERGEDENERGY;

        final BondInfo bonds = topology.getBonds();
        final double[][] dists = topology.getDistances();
        final String[] atoms = topology.getAtomNames();
        final short[] nos = topology.getAtomicNumbers();
        final int noOfAtoms = topology.getNumberOfAtoms();

        boolean freshCache = false;
        if(useCaching && paramPosCache == null){
            final int noOfContr = (noOfAtoms*noOfAtoms-noOfAtoms)/2;

            // make new caches
            paramPosCache = new int[noOfAtoms];
            is13Cache = new boolean[noOfContr];
            is14Cache = new boolean[noOfContr];
            freshCache = true;
        }

        // populate caches, easiest that way
        if(useCaching && freshCache){
            for (int i = 0; i < noOfAtoms; i++) {
                final String s1;
                if (useSpecialIDs) s1 = ids[i];
                else s1 = atoms[i];

                paramPosCache[i] = params.getStartPointForKey("amberlj:" + s1);
            }

            int interCounter = 0;
            for (int i = 0; i < noOfAtoms - 1; i++) {
                for (int j = i + 1; j < noOfAtoms; j++) {
                    is13Cache[interCounter] = topology.does13exist(i, j);
                    is14Cache[interCounter] = topology.does14exist(i, j);
                    interCounter++;
                }
            }
        }

        // loop over all interactions
        double e = 0.0;
        double[] all = params.getAllParamters();
        int off1;
        int off2;
        int counter = -1;
        for (int i = 0; i < noOfAtoms - 1; i++) {

            // get params and offset for this atom
            if(useCaching){
                off1 = paramPosCache[i];
            } else{
                String key;
                if(useSpecialIDs) key = "amberlj:" + ids[i];
                else key = "amberlj:" + atoms[i];
                off1 = params.getStartPointForKey(key);
            }

            final double rad1 = AtomicProperties.giveRadius(nos[i]);

            for (int j = i + 1; j < noOfAtoms; j++) {
                counter++;
                
                final double dist = dists[i][j];

                // get params and offset for this atom
                if (useCaching) {
                    off2 = paramPosCache[j];
                } else {
                    String key;
                    if (useSpecialIDs) key = "amberlj:" + ids[j];
                    else key = "amberlj:" + atoms[j];
                    off2 = params.getStartPointForKey(key);
                }

                final double rad2 = AtomicProperties.giveRadius(nos[j]);

                // cutoffs: first part
                if (dist > (blowDist * (rad1 + rad2))) continue;
                // cutoffs: second part
                else if(dist < blowClose*(rad1+rad2)){
                    grad[off1    ] += FixedValues.NONCONVERGEDGRADIENT;
                    grad[off1 + 1] += FixedValues.NONCONVERGEDGRADIENT;
                    grad[off2    ] += FixedValues.NONCONVERGEDGRADIENT;
                    grad[off2 + 1] += FixedValues.NONCONVERGEDGRADIENT;
                    e += FixedValues.NONCONVERGEDENERGY;
                    continue;
                }

                // check if we even need to compute (or scale)
                double scale = 1.0;
                if (useCaching) {
                    // 1-2 contribution: continue
                    if (bonds.hasBond(i, j)) continue;
                    // 1-3 constribution: continue
                    else if(is13Cache[counter]) continue;
                    // 1-4 constribution: scale
                    else if(is14Cache[counter]) scale = scaleFac14;
                } else {
                    if (bonds.hasBond(i, j)) continue;
                    else if(topology.does13exist(i, j)) continue;
                    else if(topology.does14exist(i, j))  scale = scaleFac14;
                }

                // check whether the eps are negative, might happen, unfortunately
                boolean cont = false;
                if(all[off1] < 0.0){
                    //System.err.println("WARNING: eps1 < 0 in Amber-LJ. Forcing gradient.");
                    grad[off1] += FixedValues.NONCONVERGEDGRADIENT;
                    e += FixedValues.NONCONVERGEDENERGY;
                    cont = true;
                }

                if(all[off2] < 0.0){
                    //System.err.println("WARNING: eps2 < 0 in Amber-LJ. Forcing gradient.");
                    grad[off2] += FixedValues.NONCONVERGEDGRADIENT;
                    e += FixedValues.NONCONVERGEDENERGY;
                    cont = true;
                }

                if(cont) continue;

                final double epsMix = sqrt(all[off1]*all[off2]);
                final double sigMix = 0.5*(all[off1+1]+all[off2+1]);

                final double radTerm  = sigMix/dist;
                final double radTerm2 = radTerm*radTerm;
                final double radTerm6 = radTerm2*radTerm2*radTerm2;
                final double radTerm12= radTerm6*radTerm6;
                
                final double term1 = 2*radTerm6-radTerm12;
                final double dist2 = 1.0/(dist*dist);
                final double dist6 = dist2*dist2*dist2;
                final double dist12 = dist6*dist6;
                final double sig2 = sigMix*sigMix;
                final double sig4 = sig2*sig2;
                final double sig5 = sig4*sigMix;
                final double sig6 = sig5*sigMix;
                final double sig11 = sig5*sig6;
                final double term2 = -scale*epsMix*(6*sig5*dist6-6*sig11*dist12);

                grad[off1    ] += -scale*all[off2]*term1/(2*epsMix);
                grad[off1 + 1] += term2;
                grad[off2    ] += -scale*all[off1]*term1/(2*epsMix);
                grad[off2 + 1] += term2;
                e += scale*epsMix*(radTerm12-2*radTerm6);
            }
        }

        return e;
    }

    @Override
    public Gradient partialCartesianGradient(Topology topology,
            AdaptiveParameters params){

        final String[] atoms = topology.getAtomNames();
        final short[] nos = topology.getAtomicNumbers();
        final double[][] dists = topology.getDistances();
        final double[][][] diffs = topology.getPosDiffs();
        final BondInfo bonds = topology.getBonds();
        final float[] charges = topology.getCharges();
        final int noOfAtoms = atoms.length;
        final double[] allParams = params.getAllParamters();

        boolean freshCache = false;
        if(useCaching && paramPosCache == null){
            int noOfContr = 0;
            for(int i = 0; i < noOfAtoms-1; i++){
                for(int j = i+1; j < noOfAtoms; j++){
                    noOfContr++;
                }
            }
            // make new caches
            paramPosCache = new int[noOfAtoms];
            is13Cache = new boolean[noOfContr];
            is14Cache = new boolean[noOfContr];
            freshCache = true;
        }

        // populate caches, easiest that way
        if(useCaching && freshCache){
            for (int i = 0; i < noOfAtoms; i++) {
                final String s1;
                if (useSpecialIDs) s1 = ids[i];
                else s1 = atoms[i];
                
                paramPosCache[i] = params.getStartPointForKey("amberlj:" + s1);
            }

            int interCounter = 0;
            for (int i = 0; i < noOfAtoms - 1; i++) {
                for (int j = i + 1; j < noOfAtoms; j++) {
                    is13Cache[interCounter] = topology.does13exist(i, j);
                    is14Cache[interCounter] = topology.does14exist(i, j);
                    interCounter++;
                }
            }
        }

        final Gradient gradient = new Gradient();
        final double[][] grad = new double[3][noOfAtoms];
        double[] dscr = new double[2];
        int interCounter = -1;
        double energy = 0.0;
        for(int i = 0; i < noOfAtoms-1; i++){

            final double rad1 = AtomicProperties.giveRadius(nos[i]);

            for(int j = i+1; j < noOfAtoms; j++){        
                interCounter++;
                
                final double dist = dists[i][j];
                
                final double rad2 = AtomicProperties.giveRadius(nos[j]);

                                // check if we even need to compute (or scale)
                double scale = 1.0;
                double scaleE = 1.0;
                if(useCaching){
                    if(bonds.hasBond(i, j)) continue;
                    if(is13Cache[interCounter]) continue;
                    if(is14Cache[interCounter]){
                        scale = scaleFac14;
                        scaleE = scaleFac14E;
                    }
                } else{
                    if(bonds.hasBond(i, j)) continue;
                    if(topology.does13exist(i, j)) continue;
                    if(topology.does14exist(i, j)){
                        scale = scaleFac14;
                        scaleE = scaleFac14E;
                    }
                }
                
                final double distInv = 1.0/dist;
                final double dDistX = diffs[i][j][0];
                final double dDistY = diffs[i][j][1];
                final double dDistZ = diffs[i][j][2];

                // divide the distances in all three dimensions by the total distance to cover for the coordinate system
                final double dDivProdX = dDistX * distInv;
                final double dDivProdY = dDistY * distInv;
                final double dDivProdZ = dDistZ * distInv;

                double tmpDeriv = -scaleE*(charges[i]*charges[j])*distInv*distInv;

                // the actual Coulomb energy
                energy += scaleE*(charges[i]*charges[j])* distInv;

                // cutoffs for the LJ term
                if(dist > blowDist*(rad1+rad2)){}
                else if(dist < blowClose*(rad1+rad2)){
                    energy += CLOSECUT;
                    tmpDeriv += CLOSECUT;
                } else {
                    

                    // get params, be noisy if there are none
                    if(useCaching) getParams(allParams, paramPosCache[i], paramPosCache[j], dscr);
                    else getParams(params, atoms, ids, i, j, useSpecialIDs, dscr);
                
                    if(dscr == null){
                        System.err.println("WARNING: No (reasonable) AmberLJ parameters for " + atoms[i] + i + " and " + atoms[j] + j);
                        dscr = new double[2];
                        continue;
                    }
                
                    // compute
                    final double radTerm = dscr[1]*distInv;
                    final double radTerm2 = radTerm*radTerm;
                    final double radTerm6 = radTerm2*radTerm2*radTerm2;
                    final double radTerm12 = radTerm6*radTerm6;
                    final double t1 = scale*dscr[0];
                
                    tmpDeriv +=  t1 * (-12.0 * radTerm12 * distInv + 12.0 * radTerm6 * distInv);
                
                    energy += t1*(radTerm12-2*radTerm6);
                }
                
                final double derivX = tmpDeriv * dDivProdX;
                final double derivY = tmpDeriv * dDivProdY;
                final double derivZ = tmpDeriv * dDivProdZ;

                grad[0][i] += derivX;
                grad[0][j] -= derivX;
                grad[1][i] += derivY;
                grad[1][j] -= derivY;
                grad[2][i] += derivZ;
                grad[2][j] -= derivZ;    
            }
        }

        gradient.setGradientTotal(grad);
        gradient.setTotalEnergy(energy);

        return gradient;
    }

    @Override
    public Tuple3D<String[],int[],Integer> requiredParams(ArrayList<CartesianCoordinates>
            cartesians, ArrayList<Topology> topologies, String sMethod){

        System.out.println("INFO: We are using Lorentz-Berthelot mixing rules "
                + "and assuming all geometries to be of the same system!");

        String[] sa;
        if(useSpecialIDs){
            sa = ids;
        } else{
            sa = cartesians.get(0).getAllAtomTypes();
        }

        // figure all non-redundant keys out
        final LinkedList<String> nonRedKeys = new LinkedList<>();
        for(String s : sa){
            if(!nonRedKeys.contains(s)){
                nonRedKeys.add(s);
            }
        }

        final String[] keys = new String[nonRedKeys.size()];
        final int[] paramsPerKey = new int[keys.length];
        int sum = 0;
        for(int i = 0; i < keys.length; i++){
            keys[i] = "amberlj:" + nonRedKeys.get(i);
            sum += 2;
            paramsPerKey[i] = 2;
        }

        final Tuple3D<String[],int[],Integer> result = new Tuple3D<>(keys,paramsPerKey,sum);

        return result;
    }

    @Override
    public double[][] bordersForMyParams(AdaptiveParameters params){

        final int[] startEnd = params.getStartAndEndForPrefix("amberlj:");
        if(startEnd == null) return null;
        final int size = startEnd[1]-startEnd[0];
        final double[][] borders = new double[2][size+1];

        int i = 0;
        while(i < size){
            borders[0][i] = -0.01;
            borders[1][i] = 0.01;
            i++;

            borders[0][i] = 0.0;
            borders[1][i] = 15.0;
            i++;
        }
        
        return borders;
    }

    private static void getParams(final AdaptiveParameters params, final String[] atoms,
            final String[] ids, final int i, final int j, final boolean useIDs,
            double[] d){

        String s1;
        String s2;
        if(useIDs){
            s1 = ids[i];
            s2 = ids[j];
        } else{
            s1 = atoms[i];
            s2 = atoms[j];
        }

        if(s1.equalsIgnoreCase(s2)){
            d = params.getParametersForKey("amberlj:" + s1);
            return;
        }

        final double[] d1 = params.getParametersForKey("amberlj:" + s1);
        final double[] d2 = params.getParametersForKey("amberlj:" + s2);
        if(d1 == null || d2 == null){
            d = null;
            return;
        }

        // Lorentz-Berthelot mixing
        d[0] = Math.sqrt(d1[0] * d2[0]);
        d[1] = 0.5 * (d1[1] + d2[1]);
    }

    private static void getParams(final double[] allParams, final int first1,
            final int first2, double[] d){

        // no params for this pair
        if(first1 < 0 || first2 < 0){
            d[0] = 0.0; d[1] = 0.0;
            return;
        }
        
        // one below 0, will cause problems!
        if(allParams[first1] < 0.0 || allParams[first2] < 0.0){
            d[0] = 0.0; d[1] = 0.0;
            return;
        }

        // Lorentz-Berthelot mixing
        d[0] = sqrt(allParams[first1]*allParams[first2]);
        d[1] = 0.5 * (allParams[first1+1]+allParams[first2+1]);
    }
}
