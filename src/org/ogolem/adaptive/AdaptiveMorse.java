/**
Copyright (c) 2010-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
              2017, J. M. Dieterich and B. Hartke
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
import java.util.Iterator;
import java.util.LinkedList;
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;

/**
 * A morse potential.
 * @author Johannes Dieterich
 * @version 2017-03-03
 */
public final class AdaptiveMorse extends AbstractAdaptiveBackend{

    private static final long serialVersionUID = (long) 20150727;
    private final AdaptiveParameters params;
    private final double closeBlow;
    private final double distBlow;
    private final boolean useCaching;
    private int[] paramOffsetCache;
    private double[] radiiCache;

    public AdaptiveMorse(final boolean bIsInAdaptive, final String sMethod,
            final AdaptiveParameters parameters){

        // extract the cutoffs from the string
        final String s = sMethod.substring(14).trim();
        final String[] sa = s.split("\\,");
        double d1 = 0.2;
        double d2 = 20.0;
        boolean b = false;
        try{
            d1 = Double.parseDouble(sa[0].trim());
            d2 = Double.parseDouble(sa[1].trim());
            b = Boolean.parseBoolean(sa[2].trim());
        } catch(Exception e){
            System.err.println("ERROR: Couldn't extract Morse cutoffs from method string. Using defaults. " + e.toString());
        }

        this.closeBlow = d1;
        this.distBlow = d2;
        this.useCaching = b;

        if(!bIsInAdaptive){
            if(parameters != null){
                params = parameters;
            } else {
                System.out.println("INFO: Performance penalty encountered, by " +
                        "needing to read the parameters in again and again. Consider changing. :-)");
                final String sFile = "adaptive-morse.param";
                String[] saData;
                try {
                    saData = Input.ReadFile(sFile);
                } catch (Exception e) {
                    System.err.println("ERROR: Couldn't read parameter file in: " + e.toString());
                    saData = new String[0];
                }
                params = new AdaptiveParameters(saData,0);
            }
        } else {
            // we simply do not need them
            params = null;
        }
        // initialize a null'd cache anyhow
        paramOffsetCache = null;
    }

    private AdaptiveMorse(AdaptiveMorse orig){
        if(orig.params != null){
            this.params = new AdaptiveParameters(orig.params);
        } else{
            this.params = null;
        }
        this.closeBlow = orig.closeBlow;
        this.distBlow = orig.distBlow;
        this.useCaching = orig.useCaching;
        if(orig.paramOffsetCache != null) {this.paramOffsetCache = orig.paramOffsetCache.clone();}
        else {this.paramOffsetCache = null;}
    }

    @Override
    public AdaptiveMorse clone(){
        return new AdaptiveMorse(this);
    }

    @Override
    public String getMethodID(){
        return "Adaptive Morse";
    }

    @Override
    public void gradientCalculation(long lID, int iIteration, double[] daXYZ1D,
            String[] saAtomTypes, short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges, short[] iaSpins,
            final BondInfo bonds, final Gradient gradient, final boolean hasRigidEnv){

        final double[][] xyz = new double[3][iNoOfAtoms];
        System.arraycopy(daXYZ1D, 0, xyz[0], 0, iNoOfAtoms);
        System.arraycopy(daXYZ1D, iNoOfAtoms, xyz[1], 0, iNoOfAtoms);
        System.arraycopy(daXYZ1D, 2*iNoOfAtoms, xyz[2], 0, iNoOfAtoms);

        gradient.zeroGradient();
        final double[][] grad = gradient.getTotalGradient();

        if(useCaching && paramOffsetCache == null){
            initializeCaches(params, saAtomTypes);
        }

        double[] p = params.getAllParamters();
        double energy = 0.0;
        int offset = 0;
        int counter = -1;
        
        int lastOffset = 0;
        for(int i = 0; i < atsPerMol.length-1; i++){
            lastOffset += atsPerMol[i];
        }
        
        final int firstLoopAtomNo = (hasRigidEnv) ? lastOffset - 1: iNoOfAtoms-1;
        for (int i = 0; i < firstLoopAtomNo; i++) {

            double rad1;
            if(useCaching) {rad1 = radiiCache[i];}
            else {rad1 = AtomicProperties.giveRadius(saAtomTypes[i]);}

            for(int j = i+1; j < iNoOfAtoms; j++){

                counter++;

                double rad2;
                if(useCaching) {rad2 = radiiCache[j];}
                else {rad2 = AtomicProperties.giveRadius(saAtomTypes[j]);}

                final double dX = xyz[0][i]-xyz[0][j];
                final double dY = xyz[1][i]-xyz[1][j];
                final double dZ = xyz[2][i]-xyz[2][j];
                final double dist = Math.sqrt(dX*dX + dY*dY + dZ*dZ);

                // distance cutoff
                if(dist > distBlow*(rad1+rad2)){
                    // we assume no contribution
                    continue;
                }

                if(!useCaching){
                    p = getParamsForDouble(saAtomTypes[i], saAtomTypes[j], params);
                    if(p == null) {continue;}
                } else{
                    offset = paramOffsetCache[counter];
                }

                final double expTerm = Math.exp(-p[offset+1]*(dist-p[offset+2]));
                final double expPart = 2.0*p[offset]*p[offset+1]*expTerm*(1-expTerm);

                final double distInv = 1/dist;
                final double divProdX = (xyz[0][i] - xyz[0][j]) * distInv;
                final double divProdY = (xyz[1][i] - xyz[1][j]) * distInv;
                final double divProdZ = (xyz[2][i] - xyz[2][j]) * distInv;

                grad[0][i] += expPart*divProdX;
                grad[0][j] -= expPart*divProdX;
                grad[1][i] += expPart*divProdY;
                grad[1][j] -= expPart*divProdY;
                grad[2][i] += expPart*divProdZ;
                grad[2][j] -= expPart*divProdZ;

                final double t = 1-expTerm;
                energy += p[offset] * t*t;
            }
        }

        // put everything in
        gradient.setTotalEnergy(energy);
    }

    @Override
    public double energyCalculation(long lID, int iIteration, double[] daXYZ1D,
            String[] saAtomTypes, short[] atomNos, int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges, short[] iaSpins,
            final BondInfo bonds, final boolean hasRigidEnv){

        // this is a fake and just allowed since it is intermediate
        final int[] iaAtsPerMol = {iNoOfAtoms};
        final CartesianCoordinates cartes = new CartesianCoordinates(iNoOfAtoms, 1, iaAtsPerMol);
        cartes.setAllCharges(faCharges);
        cartes.setAllSpins(iaSpins);
        cartes.setAllAtomTypes(saAtomTypes);
        cartes.setAll1DCartes(daXYZ1D, iNoOfAtoms);

        final double dEnergy = energyOfStructWithParams(cartes, params, (int) lID, bonds, hasRigidEnv);

        return dEnergy;
    }

    @Override
    public double energyOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds){
        return energyOfStructWithParams(cartes, params, geomID, bonds, cartes.containedEnvType() == CartesianCoordinates.ENVTYPE.RIGID);
    }
        
    public double energyOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds,
            final boolean hasRigidEnv){

        final int noOfAtoms = cartes.getNoOfAtoms();
        final String[] atoms = cartes.getAllAtomTypes();
        final double[][] xyz = cartes.getAllXYZCoord();

        if(useCaching && paramOffsetCache == null){
            initializeCaches(params, atoms);
        }

        double[] p = params.getAllParamters();
        double energy = 0.0;
        int counter = -1;
        int offset = 0;
        
        final int[] atsPerMol = cartes.getAllAtomsPerMol();
        int lastOffset = 0;
        for(int i = 0; i < atsPerMol.length-1; i++){
            lastOffset += atsPerMol[i];
        }
        
        final int firstLoopAtomNo = (hasRigidEnv) ? lastOffset - 1: noOfAtoms-1;
        for (int i = 0; i < firstLoopAtomNo; i++) {

            double rad1;
            if(useCaching) {rad1 = radiiCache[i];}
            else {rad1 = AtomicProperties.giveRadius(atoms[i]);}

            for(int j = i+1; j < noOfAtoms; j++){

                counter++;

                double rad2;
                if(useCaching) {rad2 = radiiCache[j];}
                else {rad2 = AtomicProperties.giveRadius(atoms[j]);}

                if (!useCaching) {
                    // get params
                    p = getParamsForDouble(atoms[i], atoms[j], params);
                    if (p == null) {continue;}
                } else {
                    // use the cache
                    offset = paramOffsetCache[counter];
                }

                final double dX = xyz[0][i]-xyz[0][j];
                final double dY = xyz[1][i]-xyz[1][j];
                final double dZ = xyz[2][i]-xyz[2][j];
                final double dist = Math.sqrt(dX*dX + dY*dY + dZ*dZ);

                // close distance cutoff
                if(dist < closeBlow*(rad1+rad2)){
                    energy += 1000.0;
                    continue;
                }

                // distance cutoff
                if(dist > distBlow*(rad1+rad2)){
                    // we assume no contribution
                    continue;
                }

                final double t = 1-Math.exp(-p[offset+1]*(dist-p[offset+2]));
                energy += p[offset] * t*t;
            }
        }

        return energy;
    }

    @Override
    public double gradientOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, final int geomID, final BondInfo bonds,
            final double[] daGrad){

        final double[] daParams = params.getAllParamters();
        final String[] saAtoms = cartes.getAllAtomTypes();
        final double[][] xyz = cartes.getAllXYZCoord();
        final int iNoOfAtoms = saAtoms.length;

        //XXX enable caching and rewrite, this implementation is very inefficient!
        for(int param = 0; param < daGrad.length; param += 3){

            // loop over all interactions
            for(int i = 0; i < iNoOfAtoms-1; i++){
                for(int j = i+1; j < iNoOfAtoms; j++){
                    
                    final String sPair1 = "adaptivemorse:" + saAtoms[i] + "_" + saAtoms[j];
                    final String sPair2 = "adaptivemorse:" + saAtoms[j] + "_" + saAtoms[i];

                    if(   !sPair1.equalsIgnoreCase(params.getForWhichKey(param))
                       && !sPair2.equalsIgnoreCase(params.getForWhichKey(param))
                      ) {continue;}
                    
                    final double dX = xyz[0][i]-xyz[0][j];
                    final double dY = xyz[1][i]-xyz[1][j];
                    final double dZ = xyz[2][i]-xyz[2][j];

                    final double dist = Math.sqrt(dX*dX+dY*dY+dZ*dZ);

                    // please note that the missing - is on purpose in the first term!
                    final double inner = daParams[param+1]*(dist-daParams[param+2]);
                    final double expTermInv = 1/Math.exp(inner);
                    final double tmp = 1-Math.exp(-inner);
                    final double expT = tmp*tmp;

                    daGrad[param    ] += expT*expT;
                    daGrad[param + 1] += 2*daParams[param]*expT*(dist-daParams[param+2]) * expTermInv;
                    daGrad[param + 2] += -2*daParams[param]*daParams[param+1]*expT * expTermInv;
                }
            }

        }

        /*
         *  since the gradient class wants an array of arrays, we adjust a little
         */
        // XXX energy evalution
        final double energy = energyOfStructWithParams(cartes, params, geomID, bonds);

        return energy;
     }

    @Override
    public double[][] minMaxBordersForParams(final AdaptiveParameters params){

        final int noOfParams = params.getNumberOfParamters();
        final double[][] borders = new double[2][noOfParams];

        for(int i = 0; i < noOfParams; i += 3){
            // first one: D
            borders[0][i] = 0.0;
            borders[1][i] = 1.0;

            // second one: a
            borders[0][i+1] = 0.0;
            borders[1][i+1] = 10.0;

            // third one: r_e
            borders[0][i+2] = 0.0;
            borders[1][i+2] = 10.0;
        }

        return borders;
    }

    @Override
    public AdaptiveParameters createInitialParameterStub(final ArrayList<CartesianCoordinates> refCartes,
            final String sMethod){

        // loop through all reference geometries and figure a list of non-redundant atom tyes out
        final LinkedList<String> llAtoms = new LinkedList<>();

        final Iterator<CartesianCoordinates> itRefGeoms = refCartes.iterator();

        while(itRefGeoms.hasNext()){
            final CartesianCoordinates cartesTemp = itRefGeoms.next();
            final String[] saAtoms = cartesTemp.getAllAtomTypes();

            for (final String saAtom : saAtoms) {
                if (!llAtoms.contains(saAtom)) {
                    // add it to the list
                    llAtoms.add(saAtom);
                }
            }
        }

        int iNoOfPairs = 0;
        for (int i = 0; i < llAtoms.size(); i++) {
            for (int j = 0; j <= i; j++) {
                iNoOfPairs++;
            }
        }

        // we need a set of parameters for each pair
        final int iNoOfAtoms = llAtoms.size();
        final String[] saAtoms = new String[iNoOfPairs];
        final int[] iaParamsPerAt = new int[iNoOfPairs];

        int iCounter = 0;
        int iParamSum = 0;
        for (int i = 0; i < iNoOfAtoms; i++) {
            for (int j = i; j < iNoOfAtoms; j++) {
                final String sPair = "adaptivemorse:" + llAtoms.get(i) + "_" + llAtoms.get(j);
                saAtoms[iCounter] = sPair;
                // independent of which pair, always just 3 parameters
                iaParamsPerAt[iCounter] = 3;
                iParamSum += iaParamsPerAt[iCounter];
                iCounter++;
            }
        }

        final AdaptiveParameters paramStub = new AdaptiveParameters(iParamSum, -1, saAtoms, iaParamsPerAt, sMethod);

        return paramStub;
    }

    private static double[] getParamsForDouble(String s1, String s2, AdaptiveParameters params){

        final String pair = "adaptivemorse:" + s1 + "_" + s2;
        double [] p = params.getParametersForKey(pair);
        if(p != null) {return p;}

        final String pair2 = "adaptivemorse:" + s2 + "_" + s1;
        p = params.getParametersForKey(pair2);

        return p;
    }
    
    private void initializeCaches(final AdaptiveParameters params, final String[] atoms){

        final int noOfAtoms = atoms.length;
        this.radiiCache = new double[noOfAtoms];
        this.paramOffsetCache = new int[(noOfAtoms*noOfAtoms-noOfAtoms)/2];

        for(int i = 0; i < noOfAtoms; i++){
            radiiCache[i] = AtomicProperties.giveRadius(atoms[i]);
        }

        int counter = 0;
        for(int i = 0; i < noOfAtoms -1; i++){
            for(int j = i+1; j < noOfAtoms; j++){
                paramOffsetCache[counter] = posOfKey(atoms, i, j, params);
                counter++;
            }
        }
    }

    private int posOfKey(final String[] atoms, final int i, final int j, final AdaptiveParameters params){

        // get parameter position for this pair
        String pair = "adaptivemorse:" + atoms[i] + "_" + atoms[j];
        int pos = params.getStartPointForKey(pair);
        if(pos >= 0) {return pos;}

        pair = "adaptivemorse:" + atoms[j] + "_" + atoms[i];
        pos = params.getStartPointForKey(pair);
        if(pos >= 0) {return pos;}

        System.err.println("WARNING: No parameters for key " + pair + ". Returning wrong ones!");
        return 0;
    }
}
