/**
Copyright (c) 2010-2012, J. M. Dieterich
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
package org.ogolem.rmi;

import org.ogolem.switches.ColorPalette;
import org.ogolem.switches.FitnessFunction;
import org.ogolem.switches.FixedValues;
import org.ogolem.switches.Switch;
import org.ogolem.switches.SwitchesConfig;
import org.ogolem.switches.Tupel;

/**
 * Initialization of a random switch.
 * @author Johannes Dieterich
 * @version 2012-06-24
 */
final class SwitchInitTask implements Task<Switch>{

    private static final long serialVersionUID = (long) 20101114;
    private final ColorPalette palette;
    private final SwitchesConfig config;
    private final Switch sw;
    private final int id;

    SwitchInitTask(final Switch refSwitch, final int futureID,
            final SwitchesConfig swConfig, final ColorPalette cp){
        this.sw = new Switch(refSwitch);
        this.palette = cp;
        this.id = futureID;
        this.config = swConfig;
    }

    @Override
    public Result<Switch> executeTask(int onClient) {

        // set random colors
        sw.changeToRandomColors();

        sw.setID(id);
        
        // we need to figure the fitness out
        final FitnessFunction fitness = new FitnessFunction(config);
        final Tupel<Double, Double, Double> energies = fitness.fitnessOfSwitch(sw);
        sw.setFitness(energies.getObject1());
        sw.setS0S1EnergyCis(energies.getObject2());
        sw.setS0S1EnergyTrans(energies.getObject3());
        
        return new Result<>(sw,true,onClient);
    }

    @Override
    public Result<Switch> getDummyAnswer(int onClient){
        sw.changeToRandomColors();
        sw.setID(id);
        sw.setFitness(FixedValues.UNGLUEABLEENERGY);

        return new Result<>(sw,true,onClient);
    }
}
