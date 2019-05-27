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
package org.ogolem.rmi;

import org.ogolem.adaptive.AdaptiveConf;
import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.core.Geometry;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Input;
import org.ogolem.core.InputSanityCheck;
import org.ogolem.core.Molecule;
import org.ogolem.generic.Optimizable;
import org.ogolem.helpers.Tuple;

/**
 * A factory for threading RMI backends.
 * @author Johannes Dieterich
 * @version 2016-01-09
 */
class ThreadingClientFactory {

    @SuppressWarnings("rawtypes")
    static GenericThreadingClientBackend<?,? extends Optimizable> constructBackend(final String jobType,
            final String inputFile, final int threads, final int myID, final RMICommunication<?> comm,
            final int maxTasks, final boolean doMaxStruct, final int indsToMerge, final long sleepTime) throws Exception {
        
        if(jobType.equalsIgnoreCase("geom")){
            
            final GlobalConfig config = Input.ConfigureMe(inputFile);
            
            final Tuple<Boolean, String[]> sanity = InputSanityCheck.isInputSane(config);
            if (!sanity.getObject1()) {
                System.err.println("ERROR: Your configuration is insane, I am afraid. Exiting.");
                for (String s : sanity.getObject2()) {
                    System.err.println(s);
                }
                throw new Exception("Insane configuration.");
            }
            
            @SuppressWarnings("unchecked")
            final GenericThreadingClientBackend<Molecule,Geometry> threader = new GenericThreadingClientBackend<>(config,threads,
                    myID,(RMICommunication<Geometry>)comm,maxTasks,doMaxStruct,indsToMerge,sleepTime);
            
            return threader;
        } else if(jobType.equalsIgnoreCase("param")){
            
            GlobalConfig globConf = null;
            try {
                globConf = Input.ConfigureMe(inputFile);
            } catch (Exception e) {
                // catch also the whole rest of exceptions, out of safety
                System.err.println("Failure in creating the global configuration!");
                throw e;
            }

            // check the input for sanity
            final Tuple<Boolean, String[]> sanity = InputSanityCheck.isInputSane(globConf);
            if (!sanity.getObject1()) {
                System.err.println("ERROR: Your configuration is insane, I am afraid. Exiting.");
                for (String s : sanity.getObject2()) {
                    System.err.println(s);
                }
                throw new Exception("Insane configuration.");
            }

            final AdaptiveConf config = globConf.getAdaptiveConf();
            
            @SuppressWarnings("unchecked")
            final GenericThreadingClientBackend<Double,AdaptiveParameters> threader = new GenericThreadingClientBackend<>(config,threads,
                    myID,(RMICommunication<AdaptiveParameters>)comm,maxTasks,doMaxStruct,indsToMerge,sleepTime);
            
            return threader;
        } /*else if(jobType.equalsIgnoreCase("switch")){
            
            // read the configuration file in
            SwitchesConfig config = SwitchesInput.readConfigIn(inputFile);
            // TODO since SwitchesConfig does not implement Configuration, this does not work yet.
            try {
                config = SwitchesInput.readConfigIn(inputFile);
            } catch (Exception e) {
                System.err.println("ERROR: Can't configure myself. Aborting.");
                throw e;
            }
            
            @SuppressWarnings("unchecked")
            final GenericThreadingClientBackend<Color,Switch> threader = new GenericThreadingClientBackend<Color,Switch>(config,threads,
                    myID,(RMICommunication<Switch>)comm,maxTasks,doMaxStruct,indsToMerge,sleepTime);
            
            return threader;
        }*/
        
        throw new RuntimeException("No optimization job " + jobType + " supported by threading client backend for RMI.");
    }
}
