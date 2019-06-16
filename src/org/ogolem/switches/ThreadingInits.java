/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
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

import java.util.concurrent.*;
import org.ogolem.generic.genericpool.GenericPool;

/**
 * A threading initial fill of the switch population.
 * @author Johannes Dieterich
 * @version 2010-08-03
 */
final class ThreadingInits {

    private final SwitchesConfig swConfig;
    private final GenericPool<Color,Switch> pool;
    private final int iNoOfThreads;

    ThreadingInits(final SwitchesConfig switchesConfig, final int iNumberOfThreads,
            final GenericPool<Color,Switch> pool){
        this.iNoOfThreads = iNumberOfThreads;
        this.swConfig = switchesConfig;
        this.pool = pool;
    }

    void doTheInitialPoolFill(final int iPoolSize, final Switch refSwitch){
        final ExecutorService threadpool = Executors.newFixedThreadPool(iNoOfThreads);

        // do it for all
        for(int i = 0; i < iPoolSize; i++){
            threadpool.submit(createInitialFillTask(refSwitch, i, swConfig,
                    pool, Taboos.getReference()));
        }

        // shut the pool down
        threadpool.shutdown();

        // wait for all threads to finish
        try{
            threadpool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        } catch(InterruptedException e){
            System.err.println("Threadpool reached wallclock limit. This should really NEVER happen! " + e.toString());
        }
    }

    private static Runnable createInitialFillTask(final Switch refSwitch,
            final int i, final SwitchesConfig switchesConfig,
            final GenericPool<Color,Switch> pool, final Taboos taboos){

        return () -> {
            final Switch sSwitch = new Switch(refSwitch);
            
            // we need to set random colors
            sSwitch.changeToRandomColors();
            
            sSwitch.setID(i);
            
            // we need to figure the fitness out
            final FitnessFunction fitness = new FitnessFunction(switchesConfig);
            final Tupel<Double, Double, Double> energies = fitness.fitnessOfSwitch(sSwitch);
            sSwitch.setFitness(energies.getObject1());
            sSwitch.setS0S1EnergyCis(energies.getObject2());
            sSwitch.setS0S1EnergyTrans(energies.getObject3());
            
            // now add the done switch to the pool and the taboo list
            pool.addIndividualForced(sSwitch, sSwitch.getFitness());
            taboos.addTaboo(sSwitch);
        };
    }
}
