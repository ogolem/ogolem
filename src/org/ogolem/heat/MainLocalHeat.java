/**
Copyright (c) 2017, J. M. Dieterich and B. Hartke
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
package org.ogolem.heat;

import org.ogolem.core.Geometry;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Input;
import org.ogolem.core.Molecule;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericInitializer;
import java.util.Random;

/**
 * Entry point to use local heat pulses as a global optimization (as they were
 * originally intended by Arnulf Moebius.
 * @author Johannes Dieterich
 * @version 2017-06-18
 */
public class MainLocalHeat {
    
    public static void run(final String[] args){
        
        //final Random r = new Random();
        //Lottery.setGenerator(new org.ogolem.random.StandardRNG(r.nextLong()));
        final String inpFile = args[0];
        
        GlobalConfig gc = null;
        try{
            gc = Input.ConfigureMe(inpFile);
        } catch(Exception e){
            System.err.println("Failure to configure.");
            e.printStackTrace(System.err);
            System.exit(1);
        }
        
        GenericFitnessFunction<Molecule,Geometry> fitFunc = null;
        try{
            fitFunc = gc.getFitnessFunction();
        } catch(Exception e){
            System.err.println("Can't get fitness function.");
            e.printStackTrace(System.err);
            System.exit(2);
        }
        
        System.out.println("Using fitness function: " + fitFunc.getMyID());
        System.out.println("Please make sure this is local heat!");
        
        GenericInitializer<Molecule,Geometry> initer = null;
        try{
            initer = gc.getInitializer();
        } catch(Exception e){
            System.err.println("Can't get initer!");
            e.printStackTrace(System.err);
            System.exit(3);
        }
        
        Geometry example = null;
        try{
            example = gc.getExample();
        } catch(Exception e){
            System.err.println("Can't get example.");
            e.printStackTrace(System.err);
            System.exit(4);
        }
        
        // now get the random start
        final Geometry g = initer.initializeOnly(example, 1);
        
        // do the local heat pulsing
        final Geometry fin = fitFunc.fitness(g, false);
        
        System.out.println("Best geometry found:");
        final String[] forma = fin.makePrintableAbsoluteCoord(true);
        for(final String s: forma){
            System.out.println(s);
        }
    }
}
