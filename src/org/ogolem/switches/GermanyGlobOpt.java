/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2012, J. M. Dieterich
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

/**
 * Standard one-point genotype crossover.
 * @author Johannes Dieterich
 * @version 2010-04-17
 */
final class GermanyGlobOpt implements SwitchesDarwin{

    private final boolean bMoreMut;

    private final Taboos taboos;

    private final FitnessFunction fitness;

    GermanyGlobOpt(final SwitchesConfig swConfig){
        this.bMoreMut = swConfig.bMoreMutation;
        this.taboos = Taboos.getReference();
        this.fitness = new FitnessFunction(swConfig);
    }

    @Override
    public Switch doTheGlobOpt(final int iID, final Switch sMother, final Switch sFather){

        // let the switches cross
        final ArrayList<Switch> vChilds = Cross(sMother, sFather);

        // mutate
        final Switch swOne = Mutate(vChilds.get(0));
        final Switch swTwo = Mutate(vChilds.get(1));

        // put the correct IDs in
        swOne.setID(iID);
        swTwo.setID(iID);

        // check for taboos
        final boolean bKnownOne = taboos.isThisKnown(swOne);
        final boolean bKnownTwo = taboos.isThisKnown(swTwo);

        // just optimize if there is no taboo
        double dFitnessOne = FixedValues.UNGLUEABLEENERGY;
        double dFitnessTwo = FixedValues.UNGLUEABLEENERGY;

        if(!bKnownOne){
            final Tupel<Double, Double, Double> energies = fitness.fitnessOfSwitch(swOne);
            dFitnessOne = energies.getObject1();
            swOne.setFitness(dFitnessOne);
            swOne.setS0S1EnergyCis(energies.getObject2());
            swOne.setS0S1EnergyTrans(energies.getObject3());
        }

        if(!bKnownTwo){
            final Tupel<Double, Double, Double> energies = fitness.fitnessOfSwitch(swTwo);
            dFitnessTwo = energies.getObject1();
            swTwo.setFitness(dFitnessOne);
            swTwo.setS0S1EnergyCis(energies.getObject2());
            swTwo.setS0S1EnergyTrans(energies.getObject3());
        }

        if(bKnownOne && bKnownTwo){
            return null;
        } else if(bKnownOne && !bKnownTwo){
            return swTwo;
        } else if(bKnownTwo && !bKnownOne){
            return swOne;
        } else if(!bKnownTwo && !bKnownOne && dFitnessOne <= dFitnessTwo){
            return swOne;
        } else if(!bKnownTwo && !bKnownOne && dFitnessTwo < dFitnessOne){
            return swTwo;
        } else{
            System.err.println("ERROR: Unknown case in globopt. Contact the author please.");
            return null;
        }
    }

    @Override
    public ArrayList<Switch> Cross(final Switch sMother, final Switch sFather){

        final ArrayList<Switch> alChildsFinal = new ArrayList<>(2);

        final ArrayList<Color> alMother = sMother.getCopyOfColors();
        final ArrayList<Color> alFather = sFather.getCopyOfColors();

        final ArrayList<ArrayList<Color>> alChilds = GlobOptAtomics.genotypeCross(alMother, alFather);

        final Switch sChildOne = new Switch(sMother);
        final Switch sChildTwo = new Switch(sFather);

        // set the new colors in
        sChildOne.setAllColors(alChilds.get(0));
        sChildTwo.setAllColors(alChilds.get(1));

        // return the children
        alChildsFinal.add(sChildOne);
        alChildsFinal.add(sChildTwo);
        
        return alChildsFinal;
    }

    @Override
    public Switch Mutate(final Switch switchStart){

        final Switch switchEnd = new Switch(switchStart);

        ArrayList<Color> alColors = switchStart.getCopyOfColors();
        
        alColors = GlobOptAtomics.genotypeMutation(alColors, bMoreMut);

        switchEnd.setAllColors(alColors);

        return switchEnd;
    }
}
