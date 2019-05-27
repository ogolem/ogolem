/**
Copyright (c) 2014, J. M. Dieterich
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
package org.ogolem.clusters;

import java.io.File;
import java.io.IOException;
import org.ogolem.core.Geometry;
import org.ogolem.core.GeometryReader;
import org.ogolem.core.GeometryWriter;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Input;
import org.ogolem.core.InputSanityCheck;
import org.ogolem.core.Molecule;
import org.ogolem.generic.threading.GenericOGOLEMOptimization;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.generichistory.GenericHistoryConfig;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.GenericPoolConfig;
import org.ogolem.helpers.Tuple;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * The entry point for cluster global optimization, threading style.
 * @author Johannes Dieterich
 * @version 2016-04-10
 */
public class MainThreadingClusterGlobOpt {
    
    public static void execute(final String[] args){
        
        if(args[0].equalsIgnoreCase("help")){
            System.out.println("This is the threading cluster global optimization functionality.");
            System.out.println("Required arguments:");
            System.out.println(" * the input file (.ogo format)");
            System.out.println(" * the number of threads to be used");
            System.exit(0);
        }
        
        final String configFile = args[0];
        if(!configFile.endsWith(".ogo")){
            System.err.println("ERROR: Specified input file not an .ogo one. " + configFile);
            System.exit(1);
        }
        final int noThreads = Integer.parseInt(args[1]);
        
        final Tuple<String,String> dirs = ManipulationPrimitives.outDirAndBaseName(configFile);
        final String outFolder = dirs.getObject1();
        final String baseName = dirs.getObject2();
        
        final String outFile = outFolder + File.separator + baseName + ".out";
        
        // parse config
        GlobalConfig conf = null;
        try{
            conf = Input.ConfigureMe(configFile);
        } catch(Exception e){
            System.err.append("ERROR: Failure configure ogolem from input.");
            e.printStackTrace(System.err);
            System.exit(1);
        }
        
        conf.setOutputFolder(outFolder);
        conf.setOutputFile(outFile);
        
        // check the input for sanity
        final Tuple<Boolean,String[]> sanity = InputSanityCheck.isInputSane(conf);
        if(!sanity.getObject1()){
            System.err.println("ERROR: Your configuration is insane, I am afraid. Exiting.");
            for(String s : sanity.getObject2()){
                System.err.println(s);
            }
            System.exit(1);
        }
        
        // serialize the globconf for restart purposes
	final String configBin = outFolder + "-config.bin";
        try{
            OutputPrimitives.writeObjToBinFile(configBin, conf);
        } catch(Exception e){
            System.err.println("ERROR: Failure in binary writing the global configuration. " + e.toString());
            e.printStackTrace(System.err);
            System.exit(37);
        }


        final boolean checkMove = !conf.doesRestart();
        if(checkMove){
            try{
                ManipulationPrimitives.moveOldFolders(outFolder, 9, System.getProperty("user.dir"));
            } catch (Exception e) {
                System.err.println("ERROR: Failure in initial file/folder operations. " + e.toString());
                e.printStackTrace(System.err);
                System.exit(1);
            }
        }
        
        // create folder
        try{
            OutputPrimitives.createAFolder(outFolder);
        } catch(IOException e){
            System.err.println("ERROR: Couldn't create output folder!");
            e.printStackTrace(System.err);
            System.exit(2);
        }
                
        // assemble some helper objects
        final GeometryWriter writer = new GeometryWriter(outFolder);
        final GeometryReader reader = new GeometryReader();
        GenericHistory<Molecule,Geometry> history = null;
        GenericPool<Molecule,Geometry> pool = null;
        try{
            assert(conf != null);
            final GenericHistoryConfig hisConf = conf.getGenericHistoryConfig();
            history = GenericHistory.getReference(hisConf);
            final GenericPoolConfig<Molecule,Geometry> poolConf = conf.getGenericPoolConfig();
            final Geometry example = new Geometry(conf.geoConfCopy());
            pool = GenericPool.getInstance(poolConf,example);
        } catch(Exception e){
            System.err.println("ERROR: Couldn't instantiate helper objects in the cluster optimizer.");
            e.printStackTrace(System.err);
            System.exit(4);
        }
        
        assert(pool != null);
        assert(history != null);
        
        // start our generic global optimization
        final GenericOGOLEMOptimization<Molecule,Geometry> opter
                = new GenericOGOLEMOptimization<>(conf, pool, history, outFile, 
                        outFolder, writer, reader, conf.getTasksToSubmit());
        
        // run it
        opter.globOpt(noThreads);
        
        // be done ;-)
    }
}
