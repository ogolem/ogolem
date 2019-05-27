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

import java.io.IOException;
import java.rmi.registry.LocateRegistry;
import java.rmi.registry.Registry;
import org.ogolem.generic.Optimizable;
import org.ogolem.helpers.Tuple;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * A threading client for the RMI setup.
 * @author Johannes Dieterich
 * @version 2016-04-14
 */
public class MainRMIThreader {
   
    private static final boolean DEBUG = false;
    
    public static void main(final String args[]) {
        
        //XXX fix...
        if (args == null || args[0].equalsIgnoreCase("--help")) {
            System.out.println("Sorry, no help yet...");
            System.exit(0);
        }

        String server = null;
        String name = "RMICommunication";
        long sleepTime = 5000;
        int port = 1099;
        String jobType = null;
        String inputFile = null;
        int noThreads = Runtime.getRuntime().availableProcessors();
        boolean doMaxStruct = false;
        int indsToMerge = -1;
        int maxTasks = 1000;
        
        for(int argID = 0; argID < args.length; argID++){
            final String arg = args[argID];
            if(arg.equalsIgnoreCase("--server")){
                argID++;
                server = args[argID];
            } else if(arg.equalsIgnoreCase("--jobtype")){
                argID++;
                jobType = args[argID];
            } else if(arg.equalsIgnoreCase("--inputfile")){
                argID++;
                inputFile = args[argID];
            } else if(arg.equalsIgnoreCase("--commname")){
                argID++;
                name = args[argID];
            } else if(arg.equalsIgnoreCase("--sleeptime")){
                argID++;
                sleepTime = Long.parseLong(args[argID])*1000;
            } else if(arg.equalsIgnoreCase("--serverregistryport")){
                argID++;
                port = Integer.parseInt(args[argID]);
            } else if(arg.equalsIgnoreCase("--threads")){
                argID++;
                noThreads = Integer.parseInt(args[argID]);
            } else if(arg.equalsIgnoreCase("--mergeinds")){
                argID++;
                indsToMerge = Integer.parseInt(args[argID]);
                if(indsToMerge > 0){
                    doMaxStruct = true;
                }
            } else if(arg.equalsIgnoreCase("--maxchunk")){
                argID++;
                maxTasks = Integer.parseInt(args[argID]);
            } else {
                System.err.println("Unknown argument " + arg + ". Shutting down.");
                System.exit(1);
            }
        }
        
        String keySuffix = System.getenv("OGO_RMIKEY");
        if(keySuffix == null){
            keySuffix = "Super secret password: Nakatomi Socrates";
        }
        
        // check for mandatory settings
        if(server == null){System.err.println("No server specified."); System.exit(2);}
        if(jobType == null){System.err.println("No jobType specified."); System.exit(3);}
        if(inputFile == null){System.err.println("No inputFile specified."); System.exit(4);}
                
        final String key = "Client speaking, I am here. " + keySuffix;
        
        // first we try to connect to our master and exchange keys
        final long waittimeClientMS = 2000; // two seconds.
        RMICommunication<?> serverComm = null;
        int myID = -1;
        for(int attempt = 0; attempt < 3; attempt++){
            try{
                final Registry registry = LocateRegistry.getRegistry(server,port);
            
                if(DEBUG){
                    final String[] sa = registry.list();
                    System.out.println("DEBUG: List of remote objects in registry:");
                    for(final String s : sa){
                        System.out.println("DEBUG: " + s);
                    }
                }
                serverComm = (RMICommunication<?>) registry.lookup(name);
            
                if(serverComm == null){
                    System.err.println("WARNING: Communication object is null. Trying again....");
                }
            
                // exchange keys
                final Tuple<String,Integer> answer = serverComm.registerWithMaster(key);
                if(!answer.getObject1().equalsIgnoreCase("Master speaking, everything fine. " + keySuffix)){
                    // not good, exit (we want to avoid double registration...)
                    System.err.println("ERROR: Wrong answer from master.s Aborting client. Server was " + server + ". Answer was " + answer.getObject1());
                    System.exit(27);
                } else{
                    myID = answer.getObject2();
                    // yup, ID assigned!
                    break;
                }
            } catch(Exception e){
                System.err.println("RMI exception: " + e.getMessage() + ". Trying again... Server is " + server);
                e.printStackTrace(System.err);
            }
            
            // sleep for a bit...
            try{
                Thread.sleep(waittimeClientMS);
            } catch (Exception e){
                // whatever...
                e.printStackTrace(System.err);
            }
            
            if(attempt == 3){
                System.err.println("ERROR: Failure to connect to initial server happened three times. Giving up...");
                System.exit(42);
            }
        }
        
        // create, if necessary, the output directory
        final String outFolder = inputFile.substring(0,inputFile.indexOf(".ogo"));
        try {
            ManipulationPrimitives.moveOldFolders(outFolder, 9, System.getProperty("user.dir"));
        } catch (Exception e) {
            System.err.println("ERROR: Failure in initial file/folder operations. " + e.toString());
            e.printStackTrace(System.err);
            System.exit(1);
        }
        
        // create folder
        try{
            OutputPrimitives.createAFolder(outFolder);
        } catch(IOException e){
            System.err.println("ERROR: Couldn't create output folder!");
            e.printStackTrace(System.err);
            System.exit(2);
        }
        
        final long offsetTime = System.currentTimeMillis();
        
        // the check for the aliveness of the server
        final AliveCheck aliveCheck = new AliveCheck(serverComm, sleepTime, myID);        
        
        final long setupTime = System.currentTimeMillis();
        // now do the threading optimization one chunk at a time
        int noInitChunks = -42;
        int noGlobOptChunks = -42;
        long initFetchS = -42;
        long initWorkS = -42;
        long optWorkS = -42;
        long optFetchS = -42;
        long waitedS = -42;
        try{
            @SuppressWarnings("rawtypes")            
            final GenericThreadingClientBackend<?,? extends Optimizable> threader = ThreadingClientFactory.constructBackend(jobType,
                    inputFile, noThreads, myID, serverComm, maxTasks,doMaxStruct,indsToMerge,sleepTime);
            threader.doAllTasks();
            noInitChunks = threader.getNoOfInitChunks();
            noGlobOptChunks = threader.getNoOfGlobOptChunks();
            initFetchS = threader.getTimeInSInitFetch();
            initWorkS = threader.getTimeInSInitWork();
            optWorkS = threader.getTimeInSOptWork();
            optFetchS = threader.getTimeInSOptFetch();
            waitedS = threader.getTimeInSWaited();
        } catch(Exception e){
            System.err.println("Exception in threading RMI client.");
            e.printStackTrace(System.err);
        }
        
        aliveCheck.done();
        
        // be done
        final long totalTime = System.currentTimeMillis();
        System.out.println("THREADING CLIENT ACCOUNTING INFORMATION:");
        System.out.println("# init chunks worked on: " + noInitChunks);
        System.out.println("# globopt chunks worked on: " + noGlobOptChunks);
        System.out.println("Time for setup: " + (setupTime-offsetTime)/1000 + "s.");
        System.out.println("Time for job: " + (totalTime-setupTime)/1000 + "s.");
        System.out.println("   * time for initial fetching: " + initFetchS + "s");
        System.out.println("   * time for initial work: " + initWorkS + "s");
        System.out.println("   * time for optimization fetching: " + optFetchS + "s");
        System.out.println("              - waited: " + waitedS + "s");
        System.out.println("   * time for optimization work: " + optWorkS + "s");
        System.out.println("Total time: " + (totalTime-offsetTime)/1000 + "s.");
        System.out.println("Thank you for running OGOLEM. We hope it was successfull!");
        System.out.println("");
        System.out.println(org.ogolem.helpers.Fortune.randomFortune());
        // this needs to be here since it doesn't return on itself (very robust server ;-) )
        System.exit(0);
    }
}
