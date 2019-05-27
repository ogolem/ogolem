/**
Copyright (c) 2010-2014, J. M. Dieterich
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
package org.ogolem.rmi;

import java.io.IOException;
import java.util.ArrayList;
import org.ogolem.adaptive.AdaptiveConf;
import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.core.Geometry;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Input;
import org.ogolem.core.InputSanityCheck;
import org.ogolem.core.Molecule;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.GenericPoolConfig;
import org.ogolem.helpers.Tuple;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.switches.Backbone;
import org.ogolem.switches.Color;
import org.ogolem.switches.ColorPalette;
import org.ogolem.switches.Switch;
import org.ogolem.switches.SwitchesConfig;
import org.ogolem.switches.SwitchesInput;

/**
 * Constructs all jobs from input.
 * @author Johannes Dieterich
 * @version 2016-04-10
 */
final class JobFactory {
    
    static RMICommunication<?> constructJob(final RMICommunication<?> serverComm,
            final String optType, final String inputFile, final long timeOut, final long timeOutJob,
            final int maxTasks, final int maxIndiesExchange, final AliveCheck alive, final int myID){
        return constructJob(serverComm, true, optType, inputFile, timeOut, timeOutJob, -1, maxTasks, maxIndiesExchange, alive, myID);
    }
    
    static RMICommunication<?> constructJob(final String optType,
            final String inputFile, final long timeOut, final long timeOutJob,
            final int noProxies, final int maxIndiesExchange){
        return constructJob(null, false, optType, inputFile, timeOut, timeOutJob, noProxies, -1, maxIndiesExchange, null, -1);
    }
    
    @SuppressWarnings("unchecked") // warnigns come from the cast <?> to, e.g., <Geometry> which we know will be fine
    private static RMICommunication<?> constructJob(final RMICommunication<?> serverComm,
            final boolean isProxy, final String optType,
            final String inputFile, final long timeOut, final long timeOutJob,
            final int noProxies, final int maxTasks, final int maxIndiesExchange,
            final AliveCheck alive, final int myID){
        
        // unconditionally, we will always need an output directory. create it.
        final Tuple<String,String> dirs = ManipulationPrimitives.outDirAndBaseName(inputFile);
        final String outFolder = dirs.getObject1();
        
        try{
            OutputPrimitives.createAFolder(outFolder);
        } catch(IOException e){
            System.err.println("ERROR: Couldn't create output folder on RMI SERVER! Folder name is " + outFolder);
            e.printStackTrace(System.err);
            System.exit(2);
        }

        if(optType.equalsIgnoreCase("geom")){
            try{
                return setupGeometryOptimization((RMICommunication<Geometry>) serverComm,
                        isProxy, inputFile, timeOut, timeOutJob, noProxies, maxTasks, maxIndiesExchange,
                        alive, myID);
            } catch(Exception e){
                e.printStackTrace(System.err);
                return null;
            }
        } else if(optType.equalsIgnoreCase("param")){
            try{
                return setupParameterOptimization((RMICommunication<AdaptiveParameters>) serverComm,
                        isProxy, inputFile, timeOut, timeOutJob, noProxies, maxTasks, maxIndiesExchange,
                        alive, myID);
            } catch(Exception e){
                e.printStackTrace(System.err);
                return null;
            }
        } else if(optType.equalsIgnoreCase("switch")){
            
            if(isProxy){
                System.err.println("ERROR: switch optimization and proxy-ing does not work!");
                return null;
            }
            
            try{
                return setupSwitchOptimization(inputFile, timeOut, timeOutJob, noProxies);
            } catch(Exception e){
                e.printStackTrace(System.err);
                return null;
            }
        } else{
            System.err.println("ERROR: Unknown request " + optType + ". Using empty job now.");
            try{
                return new RMICommImpl<>(new TaskQueue<>(new EmptyJob()),true, 0, 0, noProxies);
            } catch(Exception e){
                e.printStackTrace(System.err);
                return null;
            }
        }
    }

    private static RMICommImpl<Geometry> setupGeometryOptimization(final RMICommunication<Geometry> serverComm,
            final boolean isProxy, final String inputFile,
            final long timeOut, final long timeOutJob, final int noProxies, final int maxTasks,
            final int maxIndiesExchange, final AliveCheck alive, final int myID)
        throws Exception{

        // check the input file
        if (inputFile.equalsIgnoreCase("")) {
            System.err.println("Please specify an input file!");
            throw new Exception("No input file specified.");
        } else if (inputFile.endsWith(".ogo") != true) {
            throw new Exception("No correct input file specified.");
        }

        // get the configuration
        GlobalConfig globConf = null;
        try {
            globConf = Input.ConfigureMe(inputFile);
        } catch (Exception e){
            // catch also the whole rest of exceptions, out of safety
            System.err.println("Failure in creating the global configuration!");
            throw e;
        }

        // check the input for sanity
        final Tuple<Boolean,String[]> sanity = InputSanityCheck.isInputSane(globConf);
        if(!sanity.getObject1()){
            System.err.println("ERROR: Your configuration is insane, I am afraid. Exiting.");
            for(String s : sanity.getObject2()){
                System.err.println(s);
            }
            throw new Exception("Insane configuration.");
        }

        final GenericPoolConfig<Molecule,Geometry> poolC = globConf.getGenericPoolConfig();
        if(isProxy){poolC.beSilent();}
        final GenericPool<Molecule,Geometry> pool = GenericPool.getInstance(poolC,
                globConf.getExample());
        assert(pool != null);
        final GenericHistory<Molecule,Geometry> history = GenericHistory.getReference(globConf.getGenericHistoryConfig());
        assert(history != null);
        
        final Job<Geometry> job;
        if(isProxy){
            job = new GenericProxyJob<>(globConf, pool, history, serverComm, maxTasks, maxIndiesExchange, alive, myID);
        } else {
            job = new GenericGlobOptJob<>(globConf,
                pool,history, globConf.getReader(),globConf.getWriter(globConf.getOutputFolder()),
                globConf.getOutputFolder(),globConf.getOutputFile());
        }
        final TaskQueue<Geometry> queue = new TaskQueue<>(job);

        return new RMICommImpl<>(queue,false, timeOut, timeOutJob, noProxies);
    }
    
    private static RMICommImpl<AdaptiveParameters> setupParameterOptimization(
            final RMICommunication<AdaptiveParameters> serverComm, final boolean isProxy,
            final String inputFile, final long timeOut, final long timeOutJob, final int noProxies,
            final int maxTasks, final int maxIndiesExchange, final AliveCheck alive, final int myID)
        throws Exception{

        // check the input file
        if (inputFile.equalsIgnoreCase("")) {
            System.err.println("Please specify an input file!");
            throw new Exception("No input file specified.");
        } else if (inputFile.endsWith(".ogo") != true) {
            throw new Exception("No correct input file specified.");
        }

        // get the configuration
        GlobalConfig globConf = null;
        try {
            globConf = Input.ConfigureMe(inputFile);
        } catch (Exception e){
            // catch also the whole rest of exceptions, out of safety
            System.err.println("Failure in creating the global configuration!");
            throw e;
        }

        // check the input for sanity
        final Tuple<Boolean,String[]> sanity = InputSanityCheck.isInputSane(globConf);
        if(!sanity.getObject1()){
            System.err.println("ERROR: Your configuration is insane, I am afraid. Exiting.");
            for(String s : sanity.getObject2()){
                System.err.println(s);
            }
            throw new Exception("Insane configuration.");
        }

        final AdaptiveConf adapConf = globConf.getAdaptiveConf();
        final GenericPoolConfig<Double,AdaptiveParameters> poolC = adapConf.getGenericPoolConfig();
        if(isProxy){poolC.beSilent();}
        final GenericPool<Double,AdaptiveParameters> pool = GenericPool.getInstance(poolC,
                adapConf.getExample());
        assert(pool != null);
        final GenericHistory<Double,AdaptiveParameters> history = GenericHistory.getReference(adapConf.getGenericHistoryConfig());
        assert(history != null);
        
        final Job<AdaptiveParameters> job;
        if(isProxy){
            job = new GenericProxyJob<>(adapConf, pool, history, serverComm, maxTasks, maxIndiesExchange, alive, myID);
        } else {
            job = new GenericGlobOptJob<>(adapConf,
                pool,history, adapConf.getReader(),adapConf.getWriter(globConf.getOutputFolder()),
                globConf.getOutputFolder(),globConf.getOutputFile());
        }
        final TaskQueue<AdaptiveParameters> queue = new TaskQueue<>(job);

        return new RMICommImpl<>(queue,false, timeOut, timeOutJob, noProxies);
    }

    private static RMICommImpl<Switch> setupSwitchOptimization(final String inputFile,
            final long timeOut, final long timeOutJob, final int noProxies)
        throws Exception{

        // read the configuration file in
        SwitchesConfig swConfig = null;
        try{
            swConfig = SwitchesInput.readConfigIn(inputFile);
        } catch(Exception e){
            System.err.println("ERROR: Can't configure myself. Aborting.");
            throw e;
        }

        // now do some "afterconfiguring"
        final ArrayList<Color> alColors = new ArrayList<>();

        final ColorPalette palette = ColorPalette.getReference();

        // yeah, this is pretty nasty. but this is no performance penalty
        for(int i = 0; i < SwitchesConfig.backbone.getConnectsCopy(true).length; i++){
            alColors.add(palette.getRandomColor());
        }

        final Switch refSwitch = new Switch(new Backbone(SwitchesConfig.backbone), alColors);

        final SwitchGlobOptJob job = new SwitchGlobOptJob(swConfig, refSwitch);
        final TaskQueue<Switch> queue = new TaskQueue<>(job);

        return new RMICommImpl<>(queue,false, timeOut, timeOutJob, noProxies);
    }

}
