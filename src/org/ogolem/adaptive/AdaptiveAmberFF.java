/**
Copyright (c) 2011-2014, J. M. Dieterich
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
package org.ogolem.adaptive;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.Topology;
import org.ogolem.helpers.Tuple;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.io.InputPrimitives;
import org.ogolem.math.AcosLookup;
import org.ogolem.math.CosLookup;
import org.ogolem.math.SinLookup;

/**
 * An adaptivable AMBER-style force field. Please note that this FF ONLY
 * guarantees to work with one system at a time!
 * @author Johannes Dieterich
 * @version 2016-02-24
 */
public class AdaptiveAmberFF extends AbstractAdaptiveBackend {

    private static final long serialVersionUID = (long) 20150727;
    private static final boolean DEBUG = false;

    private final AdaptiveInteractionTerm[] terms;
    private final AdaptiveParameters params;
    private final List<int[]> contr13;
    private final List<int[]> contr14;
    private final boolean useCaches;
    private final boolean useShifter;
    private ArrayList<Topology> topoCache;
    private final AmberMath math;
    private double[][] xyz; // used as an object cache when in backend mode.

    public AdaptiveAmberFF(boolean bIsInAdaptive, AdaptiveParameters parameters,
            boolean useCaching, int whichAcos, int whichCos, int whichSin, boolean useTotalEnergyShifter)
            throws Exception{
        this.useCaches = useCaching;
        this.useShifter = useTotalEnergyShifter;
        this.terms = (useShifter) ? new AdaptiveInteractionTerm[5] : new AdaptiveInteractionTerm[4];
        // add individual terms
        if(useShifter){
            final AdaptiveInteractionTerm shift = new AdaptiveTotalEnergyShifter(true);
            terms[4] = shift;
        }
        
        this.math = new AmberMath(whichAcos, whichCos, whichSin);
        if(!useCaching) terms[0] = new AdaptiveBondedHarmonicTerm(true,20.0,useCaching);
        else terms[0] = new AdaptiveCachedBondTerm(true,20.0);
        terms[1] = new AdaptiveAmberAngleTerm(true, useCaching, math);
        terms[2] = new AdaptiveAmberDihedralTerm(true, useCaching, math);
        terms[3] = new AdaptiveAmberLJTerm(true,20.0,0.8,0.5, 1/1.2, useCaching);
        // parameters
        if(!bIsInAdaptive){
            if(parameters != null){
                params = parameters;
            } else {
                System.out.println("INFO: Performance penalty encountered, by " +
                        "needing to read the parameters in again and again. Consider changing. :-)");
                final String sFile = "adaptive-amberff.param";
                String[] saData = null;
                try {
                    saData = Input.ReadFile(sFile);
                } catch (Exception e) {
                    System.err.println("ERROR: Couldn't read parameter file in: " + e.toString());
                    throw e;
                }
                params = new AdaptiveParameters(saData, 0);
            }
        } else {
            // we simply do not need them
            params = null;
            topoCache = new ArrayList<>();
        }

        // create topology: read and parse aux files, might throw various exceptions
        final String[] contrib13 = InputPrimitives.readFileIn("amberff-list13.aux");
        final String[] contrib14 = InputPrimitives.readFileIn("amberff-list14.aux");

        this.contr13 = new ArrayList<>(contrib13.length);
        this.contr14 = new ArrayList<>(contrib14.length);

        // 1-3 contributions
        for(final String s : contrib13){
            final String[] sa = s.trim().split("\\s+");
            final int[] ia = new int[3];
            for(int i = 0; i < 3; i++){
                ia[i] = Integer.parseInt(sa[i]);
            }

            // add
            contr13.add(ia);
        }

        // 1-4 contributions
        for(final String s : contrib14){

            if(s.contains("null")) break;

            final String[] sa = s.trim().split("\\s+");
            final int[] ia = new int[4];
            for(int i = 0; i < 4; i++){
                ia[i] = Integer.parseInt(sa[i]);
            }

            // add
            contr14.add(ia);
        }        
    }

    private AdaptiveAmberFF(final AdaptiveAmberFF orig){
        this.math = orig.math.clone();
        this.useCaches = orig.useCaches;
        this.useShifter = orig.useShifter;
        this.terms = new AdaptiveInteractionTerm[orig.terms.length];
        int c = 0;
        for(AdaptiveInteractionTerm term : orig.terms){
            terms[c] = term.clone();
            c++;
        }
        // initialize the electrostatic
        // parameters
        if(orig.params != null )this.params = new AdaptiveParameters(orig.params);
        else{
            this.params = null;
            this.topoCache = new ArrayList<>();
            orig.topoCache.forEach((topo) -> {
                this.topoCache.add(new Topology(topo));
            });
        }
        // 1-3 and 1-4 contributions, not allowed to be null
        this.contr13 = new ArrayList<>(orig.contr13);
        this.contr14 = new ArrayList<>(orig.contr14);
    }

    @Override
    public AdaptiveAmberFF clone(){
        return new AdaptiveAmberFF(this);
    }

    @Override
    public String getMethodID(){
        return "Adaptive Amber FF";
    }

    @Override
    public void gradientCalculation(long lID, int iIteration, double[] xyz1D,
            String[] saAtomTypes, short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges, short[] iaSpins,
            final BondInfo bonds, final Gradient gradient){
        
        if(xyz == null || !useCaches) xyz = new double[3][iNoOfAtoms];
        System.arraycopy(xyz1D, 0, xyz[0], 0, iNoOfAtoms);
        System.arraycopy(xyz1D, iNoOfAtoms, xyz[1], 0, iNoOfAtoms);
        System.arraycopy(xyz1D, iNoOfAtoms * 2, xyz[2], 0, iNoOfAtoms);

        // then the topology
        final Topology topology = new Topology(saAtomTypes, xyz, bonds, faCharges, iaSpins, atomNos,
                contr13, contr14, true);

        // now compute all contributions
        final ArrayList<Gradient> gradContrib = new ArrayList<>();
        double energy = 0.0;
        for (final AdaptiveInteractionTerm term : terms) {
            final Gradient grad = term.partialCartesianGradient(topology, params);
            energy += grad.getTotalEnergy();
            gradContrib.add(grad);
        }

        // add the gradient up
        final int[] dims = {3,iNoOfAtoms};
        final Gradient gradientTmp = new Gradient(gradContrib, dims);
        gradient.setTotalEnergy(energy);

        if(DEBUG){
            // numerical gradient to check
            Gradient numGrad = NumericalGradients.numericalGradient(lID, iIteration, xyz1D,
                saAtomTypes, atomNos, atsPerMol, energyparts, iNoOfAtoms, faCharges, iaSpins,
                bonds, this);
        }
        
        gradient.copyDataIn(gradientTmp);
    }

    @Override
    public double energyCalculation(long lID, int iIteration, double[] xyz1D,
            String[] saAtomTypes, short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms,
            float[] faCharges, short[] iaSpins, final BondInfo bonds){
        
        if(xyz == null || !useCaches) xyz = new double[3][iNoOfAtoms];
        System.arraycopy(xyz1D, 0, xyz[0], 0, iNoOfAtoms);
        System.arraycopy(xyz1D, iNoOfAtoms, xyz[1], 0, iNoOfAtoms);
        System.arraycopy(xyz1D, iNoOfAtoms * 2, xyz[2], 0, iNoOfAtoms);
        
        // the topology
        final Topology topology = new Topology(saAtomTypes, xyz,
                bonds, faCharges, iaSpins, atomNos, contr13,
                contr14, false);

        // now all interaction terms
        double energy = 0.0;
        for (final AdaptiveInteractionTerm term : terms) {
            energy += term.partialInteraction(topology, params);
        }

        return energy;
    }

    @Override
    public double energyOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds){

        Topology topology;
        if(topoCache.size() < geomID+1){
            // first the topology 
            topology = new Topology(cartes.getAllAtomTypes(), cartes.getAllXYZCoord(),
                    bonds, cartes.getAllCharges(), cartes.getAllSpins(), cartes.getAllAtomNumbers(), contr13,
                    contr14, false);
            topoCache.add(topology);
        } else{
            topology = topoCache.get(geomID);
        }

        // now all interaction terms
        double energy = 0.0;
        for (final AdaptiveInteractionTerm term : terms) {
            final double e = term.partialInteraction(topology, params);
            energy += e;
        }

        return energy;
    }

    @Override
    public double gradientOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds,
            final double[] grad){
        
        Topology topology;
        if(topoCache.size() < geomID+1){
            // first the topology 
            topology = new Topology(cartes.getAllAtomTypes(), cartes.getAllXYZCoord(),
                    bonds, cartes.getAllCharges(), cartes.getAllSpins(), cartes.getAllAtomNumbers(), contr13,
                    contr14, false);
            topoCache.add(topology);
        } else{
            topology = topoCache.get(geomID);
        }

        // now all interaction terms
        double e = 0.0;
        for (final AdaptiveInteractionTerm term : terms) {
            e += term.partialParamGradient(topology, params, grad);
        }

        if(DEBUG){
            final double[] num = new double[grad.length];
            final double numE = NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, num);
            assert(Math.abs(numE-e)< 1e-8);
            final double[] analy = grad;
            for(int i = 0; i < num.length; i++){
                if(Math.abs(num[i]-analy[i]) >= 1e-6){
                    System.out.println("DEBUG: Difference in param grad at pos " + i + " num: " + num[i] + " analytical: " + analy[i]);
                }
            }
        }

        return e;
    }

    @Override
    public double[][] minMaxBordersForParams(final AdaptiveParameters params){

        ArrayList<double[][]> borders = new ArrayList<>();
        for(AdaptiveInteractionTerm term : terms){
            double[][] daBorders = term.bordersForMyParams(params);
            if(daBorders != null) borders.add(daBorders);
        }

        // copy around
        double[][] daBorders = new double[2][params.getNumberOfParamters()];
        int iOffset = 0;
        for(double[][] daPartBorders : borders){
            int iLength = daPartBorders[0].length;
            System.arraycopy(daPartBorders[0], 0, daBorders[0], iOffset, iLength);
            System.arraycopy(daPartBorders[1], 0, daBorders[1], iOffset, iLength);
            iOffset += iLength;
        }

        return daBorders;
    }

    /**
     * Please note a couple of special issues with this particular implementation
     * 1. The used FF terms assume all structures to be of the same system.
     * 2. The used reference is the first structure.
     * 3. Bonds are once assigned here based on a blow factor of 1.3 and the first reference structure.
     * 4. Use with care. A custom parameter stub might be the better idea in most cases!
     * @param refCartes
     * @param sMethod
     * @return An intial parameter stub.
     */
    @Override
    public AdaptiveParameters createInitialParameterStub(final ArrayList<CartesianCoordinates> refCartes,
            final String sMethod){

        int iParamSum = 0;
        ArrayList<Tuple<String,Integer>> paramsPerKey = new ArrayList<>();

        final CartesianCoordinates refCart = refCartes.get(0);
        final BondInfo refBonds = org.ogolem.core.CoordTranslation.checkForBonds(refCart, 1.3);
        final Topology topology = new Topology(refCart.getAllAtomTypes(),refCart.getAllXYZCoord(),
                refBonds, refCart.getAllCharges(), refCart.getAllSpins(), refCart.getAllAtomNumbers(),
                contr13, contr14);
        final ArrayList<Topology> tops = new ArrayList<>();
        tops.add(topology);

        for(AdaptiveInteractionTerm term : terms){

            Tuple3D<String[],int[],Integer> tupel = term.requiredParams(refCartes, tops, sMethod);

            String[] saKeys = tupel.getObject1();
            int[] iaParamsKey = tupel.getObject2();
            for(int i = 0; i < saKeys.length; i++){
                Tuple<String,Integer> tup = new Tuple<>(saKeys[i],iaParamsKey[i]);
                paramsPerKey.add(tup);
            }

            iParamSum += tupel.getObject3();
        }

        String[] saKeys = new String[paramsPerKey.size()];
        int[] iaParamsPerKey = new int[paramsPerKey.size()];
        for(int i = 0; i < saKeys.length; i++){
            Tuple<String, Integer> tupel = paramsPerKey.get(i);
            saKeys[i] = tupel.getObject1();
            iaParamsPerKey[i] = tupel.getObject2();
        }

        final AdaptiveParameters paramStub = new AdaptiveParameters(iParamSum, -1, saKeys, iaParamsPerKey, sMethod);

        return paramStub;
    }
    
    public static class AmberMath implements Serializable, Cloneable {
        
        private static final long serialVersionUID = (long) 20111124;
        
        private final int wAcos;
        private final int wCos;
        private final int wSin;
        private final AcosLookup acosLook;
        private final CosLookup cosLook;
        private final SinLookup sinLook;
        
        public AmberMath(final int whichAcos, final int whichCos, final int whichSin){
            
            switch(whichAcos){
                case 0: acosLook = null; wAcos = 0; break;
                case 1: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[0]); wAcos = 1; break;
                case 2: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[1]); wAcos = 1; break;
                case 3: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[2]); wAcos = 1; break;
                case 4: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[3]); wAcos = 1; break;
                case 5: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[4]); wAcos = 1; break;
                case 11: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[0]); wAcos = 2; break;
                case 12: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[1]); wAcos = 2; break;
                case 13: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[2]); wAcos = 2; break;
                case 14: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[3]); wAcos = 2; break;
                case 15: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[4]); wAcos = 2; break;
                case 20: acosLook = null; wAcos = 3; break;
                default: acosLook = null; wAcos = 0; break;
            }
            
            switch(whichCos){
                case 0: cosLook = null; wCos = 0; break;
                case 1: cosLook = new CosLookup(CosLookup.POSSIBLE_SIZES[0]); wCos = 1; break;
                case 2: cosLook = new CosLookup(CosLookup.POSSIBLE_SIZES[1]); wCos = 1; break;
                case 3: cosLook = new CosLookup(CosLookup.POSSIBLE_SIZES[2]); wCos = 1; break;
                case 4: cosLook = new CosLookup(CosLookup.POSSIBLE_SIZES[3]); wCos = 1; break;
                case 5: cosLook = new CosLookup(CosLookup.POSSIBLE_SIZES[4]); wCos = 1; break;
                case 11: cosLook = new CosLookup(CosLookup.POSSIBLE_SIZES[0]); wCos = 2; break;
                case 12: cosLook = new CosLookup(CosLookup.POSSIBLE_SIZES[1]); wCos = 2; break;
                case 13: cosLook = new CosLookup(CosLookup.POSSIBLE_SIZES[2]); wCos = 2; break;
                case 14: cosLook = new CosLookup(CosLookup.POSSIBLE_SIZES[3]); wCos = 2; break;
                case 15: cosLook = new CosLookup(CosLookup.POSSIBLE_SIZES[4]); wCos = 2; break;
                case 20: cosLook = null; wCos = 3; break;
                default: cosLook = null; wCos = 0; break;
            }
            
            switch(whichSin){
                case 0: sinLook = null; wSin = 0; break;
                case 1: sinLook = new SinLookup(SinLookup.POSSIBLE_SIZES[0]); wSin = 1; break;
                case 2: sinLook = new SinLookup(SinLookup.POSSIBLE_SIZES[1]); wSin = 1; break;
                case 3: sinLook = new SinLookup(SinLookup.POSSIBLE_SIZES[2]); wSin = 1; break;
                case 4: sinLook = new SinLookup(SinLookup.POSSIBLE_SIZES[3]); wSin = 1; break;
                case 5: sinLook = new SinLookup(SinLookup.POSSIBLE_SIZES[4]); wSin = 1; break;
                case 11: sinLook = new SinLookup(SinLookup.POSSIBLE_SIZES[0]); wSin = 2; break;
                case 12: sinLook = new SinLookup(SinLookup.POSSIBLE_SIZES[1]); wSin = 2; break;
                case 13: sinLook = new SinLookup(SinLookup.POSSIBLE_SIZES[2]); wSin = 2; break;
                case 14: sinLook = new SinLookup(SinLookup.POSSIBLE_SIZES[3]); wSin = 2; break;
                case 15: sinLook = new SinLookup(SinLookup.POSSIBLE_SIZES[4]); wSin = 2; break;
                case 20: sinLook = null; wSin = 3; break;
                default: sinLook = null; wSin = 0; break;
            }
            
            System.out.println(" " + wSin + "  " + wCos + " " + wAcos);
        }
        
        AmberMath(final AmberMath orig){
            this.wAcos = orig.wAcos;
            this.wCos = orig.wCos;
            this.wSin = orig.wSin;
            
            if(orig.acosLook != null) this.acosLook = orig.acosLook.clone();
            else this.acosLook = null;
            
            if(orig.cosLook != null) this.cosLook = orig.cosLook.clone();
            else this.cosLook = null;
            
            if(orig.sinLook != null) this.sinLook = orig.sinLook.clone();
            else this.sinLook = null;
        }
        
        @Override
        public AmberMath clone(){
            return new AmberMath(this);
        }
        
        final double acos(final double x){
            switch(wAcos){
                case 0: return Math.acos(x);
                case 1: return acosLook.acosInter(x);
                case 2: return acosLook.acosNonInter(x);
                case 3: return org.apache.commons.math3.util.FastMath.acos(x);
                default: return Math.acos(x);
            }
        }
        
        final double cos(final double x){
            switch(wCos){
                case 0: return Math.cos(x);
                case 1: return cosLook.cosInter(x);
                case 2: return cosLook.cosNonInter(x);
                case 3: return org.apache.commons.math3.util.FastMath.cos(x);
                default: return Math.cos(x);
            }
        }
        
        final double sin(final double x){
            switch(wSin){
                case 0: return Math.sin(x);
                case 1: return sinLook.sinInter(x);
                case 2: return sinLook.sinNonInter(x);
                case 3: return org.apache.commons.math3.util.FastMath.sin(x);
                default: return Math.sin(x);
            }
        }
    }
}
