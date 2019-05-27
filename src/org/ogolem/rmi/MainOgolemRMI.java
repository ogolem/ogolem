/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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
package org.ogolem.rmi;

import java.rmi.Remote;
import java.rmi.registry.Registry;
import java.rmi.registry.LocateRegistry;
import java.rmi.server.UnicastRemoteObject;

/**
 * Entry point for a purely RMI parallelized version of the OGOLEM framework.
 * @author Johannes Dieterich
 * @version 2016-02-06
 */
public class MainOgolemRMI {

    // to disallow the GC to collect this object, it is made a class variable
    private static RMICommunication<?> comm;
    private static Registry registry;

    /**
     * The server side of the RMI part. Requires as input what should be
     * optimized, followed by the input file name, the timeout and the job timeout (both in seconds),
     * and the port on which the server is supposed to export the remote object and the port on which the
     * rmiregistry is running.
     * @param args
     */
    public static void main(String[] args){
        
        //XXX fix...
        if (args == null || args[0].equalsIgnoreCase("--help")) {
            System.out.println("Sorry, no help yet...");
            System.exit(0);
        }

        // register the security manager
        if (System.getSecurityManager() == null) {
            System.setSecurityManager(new SecurityManager());
        }

        int port = 2500;
        int reg = 1099;
        int noProxies = -1;
        String jobType = null;
        String inputFile = null;
        long timeOut = 10000l;
        long timeOutJob = 10000l;
        String name = "RMICommunication";
        
        for(int argID = 0; argID < args.length; argID++){
            final String arg = args[argID];
            if (arg.equalsIgnoreCase("--jobtype")){
                argID++;
                jobType = args[argID];
            } else if(arg.equalsIgnoreCase("--inputfile")){
                argID++;
                inputFile = args[argID];
            } else if(arg.equalsIgnoreCase("--commname")){
                argID++;
                name = args[argID];
            } else if(arg.equalsIgnoreCase("--myregistryport")){
                argID++;
                reg = Integer.parseInt(args[argID]);
            } else if(arg.equalsIgnoreCase("--myserverport")){
                argID++;
                port = Integer.parseInt(args[argID]);
            } else if(arg.equalsIgnoreCase("--noproxies")){
                argID++;
                noProxies = Integer.parseInt(args[argID]);
            } else if(arg.equalsIgnoreCase("--timeout")){
                argID++;
                timeOut = Long.parseLong(args[argID])*1000;
            } else if(arg.equalsIgnoreCase("--timeoutjob")){
                argID++;
                timeOutJob = Long.parseLong(args[argID])*1000;
            } else {
                System.err.println("ERROR: Can't parse argument " + arg + ". Shutting down.");
                System.exit(222);
            }
        }
        
        // check mandatory settings
        if(jobType == null){System.err.println("No jobType specified."); System.exit(3);}
        if(inputFile == null){System.err.println("No inputFile specified."); System.exit(4);}
        
        // set everything up
        final long offsetTime = System.currentTimeMillis();
        
        comm = JobFactory.constructJob(jobType, inputFile, timeOut, timeOutJob, noProxies, -1);
        try {
            registry = LocateRegistry.createRegistry(reg);
	    final Remote stub =  UnicastRemoteObject.exportObject(comm, port);
            registry.bind(name, stub);
        } catch (Exception e) {
            System.out.println("OGOLEM RMI server exception: " + e.getMessage());
            e.printStackTrace(System.err);
            System.exit(333);
        }

        final long setupTime = System.currentTimeMillis();

        boolean everythingDone = false;
        while(!everythingDone){
            try{
                //XXX ugly!
                Thread.sleep(5000);
            } catch(Exception e){
                System.err.println("WARNING: Can't wait before checking again the server." + e.toString());
                e.printStackTrace(System.err);
            }
            try{
                everythingDone = comm.isEverythingDone();
                if(everythingDone){
                    comm.stopClientList();
                }
            } catch(Exception e){
                System.err.println("ERROR: Can't ask the server how it is doing. "
                        + "Breaking loop and exiting. " + e.toString());
                e.printStackTrace(System.err);
                break;
            }
        }


        try{
            registry.unbind("RMICommunication");
        } catch(Exception e){
            System.err.println("ERROR: Couldn't unexport the remote object. This "
                    + "needs to be fixed by hand.");
            e.printStackTrace(System.err);
        }

        final long totalTime = System.currentTimeMillis();
        System.out.println("SERVER ACCOUNTING INFORMATION:");
        System.out.println("Time for setup: " + (setupTime-offsetTime)/1000 + "s.");
        System.out.println("Time for job: " + (totalTime-setupTime)/1000 + "s.");
        System.out.println("Total time: " + (totalTime-offsetTime)/1000 + "s.");
        System.out.println("Thank you for running OGOLEM. We hope it was successfull!");
        System.out.println("");
        System.out.println(org.ogolem.helpers.Fortune.randomFortune());
        // this needs to be here since it doesn't return on itself (very robust server ;-) )
        System.exit(0);
    }
}
