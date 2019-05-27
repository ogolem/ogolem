/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010     , J. M. Dieterich
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
package org.ogolem.switches;

import java.util.ArrayList;
import org.ogolem.core.CartesianCoordinates;

/**
 * The fitness function mashing together all possible criteria.
 * @author Johannes Dieterich
 * @version 2010-10-14
 */
public final class FitnessFunction {

    private final double dCisAngle;

    private final double dTransAngle;

    private final double dBlowBondDetect;

    private final double dMinEnergyDiff;

    private final double dExactEnergyFac;

    private final double dExactEnergyCis;

    private final double dExactEnergyTrans;

    private final double dSigFacA;

    private final double dSigFacB;

    private final double dThreshDihedralCis;

    private final double dThreshDihedralTrans;

    private final int[] iaWhichAtsDihedralCis;

    private final int[] iaWhichAtsDihedralTrans;

    private final boolean bDebug;

    // the objects that we need to accomplish the task
    private final LocalOpt locopt;
    private final ExcitorGateway excitor;


    public FitnessFunction(final SwitchesConfig swConfig){
        this.dMinEnergyDiff = swConfig.dMinEnergyDiff;
        this.dExactEnergyFac = swConfig.dExactEnergyFactor;
        this.dExactEnergyCis = swConfig.dExactEnergyCis;
        this.dExactEnergyTrans = swConfig.dExactEnergyTrans;
        this.dSigFacA = swConfig.dSigFactorA;
        this.dSigFacB = swConfig.dSigFactorB;
        this.dThreshDihedralCis = swConfig.dThreshDihedralCis;
        this.dThreshDihedralTrans = swConfig.dThreshDihedralTrans;
        this.dCisAngle = swConfig.dOptimalCisVal;
        this.dTransAngle = swConfig.dOptimalTransVal;
        this.dBlowBondDetect = swConfig.dBlowBondsFac;
        this.iaWhichAtsDihedralCis = swConfig.iaWhichAtsForDihedralCis;
        this.iaWhichAtsDihedralTrans = swConfig.iaWhichAtsForDihedralTrans;
        this.excitor = new ExcitorGateway(swConfig, swConfig.iWhichExcitor);
        this.locopt = new LocalOpt(swConfig, swConfig.iWhichLocAlgo);
        this.bDebug = swConfig.bDebug;
    }

    /**
     * Evaluates the fitness of a given switch. Will exit if something goes wrong.
     * @param sw
     * @return a tupel containing the fitness and the S0 S1 excitations energies for cis and trans
     */
    public Tupel<Double,Double,Double> fitnessOfSwitch(final Switch sw){
        
        double dFitness = 0;

        if(bDebug){
            System.out.println("DEBUG: Backbone (cis) before glueing.");
            final String[] saCis = sw.getBackbone().getCartesCopy(true).createPrintableCartesians();
            final String[] saTrans = sw.getBackbone().getCartesCopy(false).createPrintableCartesians();

            for(final String s : saCis){
                System.out.println(s);
            }

            System.out.println("DEBUG: Backbone (trans) before glueing.");
            for(final String s : saTrans){
                System.out.println(s);
            }
        }

        // first we try to magically glue the sidechains and backbones together
        final ArrayList<CartesianCoordinates> alCartes;

        /*
         * For some reason this is very error (aka exception) prone. Since speed in this function
         * isn't absolutely important, we leave it in a try/catch block for the time being
         */
        if(true){
            try{
                alCartes = MagicGlue.glueTogether(sw.getBackbone(),
                    sw.getSidechains(), dBlowBondDetect);
            } catch(Exception e){
                System.err.println("DEBUG: Exception occured in the glueing. Stack coming.");
                e.printStackTrace(System.err);
                final Tupel<Double,Double,Double> retVals = new Tupel<>(FixedValues.UNGLUEABLEENERGY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
                return retVals;
            }
        } else{
            alCartes = MagicGlue.glueTogether(sw.getBackbone(),
                sw.getSidechains(), dBlowBondDetect);
        }

        // put it into the switch
        sw.setCartesians(alCartes);

        // check whether any of the two are null (so glueing didn't work)
        boolean bDoReturn = false;
        if(alCartes.get(0) == null){
            dFitness += FixedValues.UNGLUEABLEENERGY;
            bDoReturn = true;
        } else if(alCartes.get(1) == null){
            dFitness += FixedValues.UNGLUEABLEENERGY;
            bDoReturn = true;
        }

        if(bDoReturn){
            final Tupel<Double,Double,Double> retVals = new Tupel<>(dFitness, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
            return retVals;
        }

        // then we try to locally optimize the things
        final CartesianCoordinates cartesCis = locopt.locOptThis(alCartes.get(0), (int) sw.getID());
        final CartesianCoordinates cartesTrans = locopt.locOptThis(alCartes.get(1), (int) sw.getID());

        if(cartesCis == null){
            dFitness += FixedValues.UNSTABLEENERGY;
            bDoReturn = true;
        } else if(cartesTrans == null){
            dFitness += FixedValues.UNSTABLEENERGY;
            bDoReturn = true;
        }

        if(bDoReturn){
            final Tupel<Double,Double,Double> retVals = new Tupel<>(dFitness, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
            return retVals;
        }

        if(bDebug){

            System.out.println("DEBUG: Attempting dihedral check.");

            final String[] saCartesCis = cartesCis.createPrintableCartesians();
            final String[] saCartesTrans = cartesTrans.createPrintableCartesians();

            System.out.println("DEBUG: Cis cartes coming:");
            for(final String sTempCis : saCartesCis){
                System.out.println(sTempCis);
            }

            System.out.println("DEBUG: Trans cartes coming:");
            for(final String sTempTrans : saCartesTrans){
                System.out.println(sTempTrans);
            }
        }

        /*
         * get the coordinates together for the cis/trans angle value check
         */
        final double[][] daRespCoordsCis = new double[4][3];
        final double[][] daRespCoordsTrans = new double[4][3];

        for(int i = 0; i < 4; i++){
            daRespCoordsCis[i] = cartesCis.getXYZCoordinatesOfAtom(iaWhichAtsDihedralCis[i]);
            daRespCoordsTrans[i] = cartesTrans.getXYZCoordinatesOfAtom(iaWhichAtsDihedralTrans[i]);
        }

        final boolean bCisAngle = LittleHelpers.checkDihedral(dCisAngle,
                dThreshDihedralCis, daRespCoordsCis);
        final boolean bTransAngle = LittleHelpers.checkDihedral(dTransAngle,
                dThreshDihedralTrans, daRespCoordsTrans);

        if(!bCisAngle || !bTransAngle){
            dFitness += FixedValues.UNSTABLEENERGY;
            bDoReturn = true;
            if(bDebug){
                System.out.println("DEBUG: The cis angle was " + bCisAngle + " and the trans was " + bTransAngle);
            }
        }

        if(bDoReturn){
            final Tupel<Double,Double,Double> retVals = new Tupel<>(dFitness, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
            return retVals;
        } else{
            // both structs are fine, set them in
            final ArrayList<CartesianCoordinates> alCartesTemp = new ArrayList<>(2);
            alCartesTemp.add(cartesCis);
            alCartesTemp.add(cartesTrans);
            sw.setCartesians(alCartesTemp);
        }

        // now we calculate the s0->s1 excitation energy
        final double dS0S1Cis = excitor.s0s1TransitionEnergy(cartesCis, (int) sw.getID());
        final double dS0S1Trans = excitor.s0s1TransitionEnergy(cartesTrans, (int) sw.getID());

        if(dS0S1Cis >= FixedValues.UNEXCITABLEENERGY || dS0S1Trans >= FixedValues.UNEXCITABLEENERGY){
            // at least one of the two failed for some reason
            dFitness += FixedValues.UNEXCITABLEENERGY;
            final Tupel<Double,Double,Double> retVals = new Tupel<>(dFitness, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
            return retVals;
        }

        // calculate the first term: how different are the excitation energies
        final double dDelta = Math.abs(dS0S1Cis - dS0S1Trans);
        final double dDiff = Math.abs(dDelta - dMinEnergyDiff);

        /*
         * the sigmoid function has the principle form
         * P(x)= 1 / (1+Math.exp(-x))
         * we need to adjust that in three spots of course
         * 1) min/max values
         * 2) shifting in y
         * 3) steepness
         * those are tunables
         */

        final double dSigmoid = dSigFacA * (1.0/(1.0+Math.exp(dDiff*dSigFacB)));

        if(dSigFacA != 0.0){
            dFitness += dSigmoid;
        }

        // now the term describing how far off the excitation energies are from our optimum
        final double dDiffCis = Math.abs(dS0S1Cis - dExactEnergyCis);
        final double dDiffTrans = Math.abs(dS0S1Trans - dExactEnergyTrans);

        // add that term to the fitness
        dFitness += dExactEnergyFac * (dDiffCis + dDiffTrans);

        /*
         * check for NaN and infinity cases
         */
        if(Double.isInfinite(dFitness) || Double.isNaN(dFitness)){
            // this is NOT good
            System.err.println("WARNING: Fitness is either NaN or infinite. Making it sensible. " +
                    "Probably your fitness function parameters are wrong?");
            dFitness = FixedValues.UNGLUEABLEENERGY;
        }

        /*
         * check whether the fitness is EXACTLY 0.0. This is EXTREMELY unlikely!
         */
        if(dFitness == 0.0){
            System.err.println("WARNING: We detect a fitness of exactly 0.0 . " +
                    "This is an extremely unlikely case and probably indicates " +
                    "a problem in the code. Making the fitness sensible now...");
            dFitness = FixedValues.UNGLUEABLEENERGY;
        }

        final Tupel<Double,Double,Double> retVals = new Tupel<>(dFitness, dS0S1Cis, dS0S1Trans);

        return retVals;
    }
}
