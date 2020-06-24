/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2011, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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

import java.rmi.registry.LocateRegistry;
import java.rmi.registry.Registry;
import org.ogolem.helpers.Tuple;
import static org.ogolem.rmi.RMICodes.JOBSTATE.*;

/**
 * The entry to the client part of OGOLEM's RMI.
 * @author Johannes Dieterich
 * @version 2020-06-23
 */
public class MainRMIClient {

    private static final boolean DEBUG = false;
    private static final int DEFAULTSLEEPINGTIMEFORNEWTASK = 5000;

    @SuppressWarnings({"rawtypes","unchecked"}) // Generics suck sometimes...
    public static void main(String args[]) {
        
        //XXX fix...
        if (args == null || args[0].equalsIgnoreCase("--help")) {
            System.out.println("Sorry, no help yet...");
            System.exit(0);
        }

        String server = null;
        String name = "RMICommunication";
        long sleepTime = 5000;
        int port = 1099;
        int taskWaitTime = DEFAULTSLEEPINGTIMEFORNEWTASK;
        
        for(int argID = 0; argID < args.length; argID++){
            final String arg = args[argID];
            if(arg.equalsIgnoreCase("--server")){
                argID++;
                server = args[argID];
            } else if(arg.equalsIgnoreCase("--commname")){
                argID++;
                name = args[argID];
            } else if(arg.equalsIgnoreCase("--sleeptime")){
                argID++;
                sleepTime = Long.parseLong(args[argID])*1000;
            } else if(arg.equalsIgnoreCase("--taskwaittime")){
                argID++;
                taskWaitTime = Integer.parseInt(args[argID])*1000;
            } else if(arg.equalsIgnoreCase("--serverregistryport")){
                argID++;
                port = Integer.parseInt(args[argID]);
            } else {
                System.err.println("Unknown argument " + arg + ". Shutting down.");
                System.exit(1);
            }
        }

        // check for mandatory settings
        if(server == null){System.err.println("No server specified."); System.exit(2);}
        
        String keySuffix = System.getenv("OGO_RMIKEY");
        if(keySuffix == null){
            keySuffix = "Super secret password: Nakatomi Socrates";
        }
        
        final String key = "Client speaking, I am here. " + keySuffix;

        // first we try to connect to our server and exchange keys
        final long waittimeClientMS = 2000; // two seconds.
        RMICommunication<?> comm = null;
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
                comm = (RMICommunication<?>) registry.lookup(name);
            
                if(comm == null){
                    System.err.println("WARNING: Communication object is null. Trying again....");
                }
            
                // exchange keys
                final Tuple<String,Integer> answer = comm.registerWithServer(key);
                if(!answer.getObject1().equalsIgnoreCase("Server speaking, everything fine. " + keySuffix)){
                    // not good, exit (we want to avoid double registration...)
                    System.err.println("ERROR: Wrong answer from server. Aborting client. Server was " + server + ". Answer was " + answer.getObject1());
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

        // explicit check whether either the communication or the global configuration are null not needed
        assert(comm != null);

        // the check for the aliveness of the server
        final AliveCheck aliveCheck = new AliveCheck(comm, sleepTime, myID);

        WorkLoop:
        for(;;){

            Task<?> task = null;
            try{
                task = comm.getWorkTask(myID);
            } catch(Exception e){
                System.err.println("ERROR: Can't get the task for working. "
                        + "Aborting." + e.toString());
                e.printStackTrace(System.err);
                System.exit(11);
            }
            
            if(DEBUG && task == null){System.out.println("DEBUG: Client " + myID + " task was null.");}

            RMICodes.JOBSTATE whatNext = WAITING;
            if(task != null){
                // raw type used on purpose, <?> can't be used and we do not know what we get
                Result res = task.executeTask(myID);

                if(res != null){
                    // try to return the result
                    try{
                        whatNext = comm.returnResult(res);
                    } catch(Exception e){
                        System.err.println("ERROR: Can't hand solution over. "  + e.toString() + " Sleeping and trying again. My ID is " + myID + ".");
                        try{
                            //XXX ugly!
                            Thread.sleep(taskWaitTime);
                            whatNext = comm.returnResult(res);
                        } catch(Exception e2){
                            System.err.println("ERROR: Can't hand solution over for the second time. Aborting. "  + e2.toString() + " My ID is " + myID + ".");
                            e2.printStackTrace(System.err);
                            System.exit(9);
                        }
                    }
                } else{
                    System.err.println("ERROR: The result is null. This is a problem. Aborting. My ID is " + myID + ".");
                    System.exit(10);
                }
            }

            if(DEBUG){System.out.println("DEBUG: Client " + myID + " got as next task: " + whatNext);}
            
            if(whatNext == CONTINUE){
                // we loop on
                continue WorkLoop;
            } else if(whatNext == WAITING){
                // waiting for something on the server side
                while(true){

                    try{
                        whatNext = comm.whatNext(myID);
                    } catch(Exception e){
                        System.err.println("ERROR: Can't hand solution over. Aborting. " + e.toString() + " My ID is " + myID + ".");
                        e.printStackTrace(System.err);
                        System.exit(8);
                    }

                    if(DEBUG){System.out.println("DEBUG: Client " + myID + " got as next task (II): " + whatNext);}
                    
                    if(whatNext == CONTINUE){
                        continue WorkLoop;
                    } else if(whatNext == WAITING){
                        // no new information. sleep and then loop on.
                    } else if(whatNext == FINISH){
                        System.out.println("INFO: Shutting down gracefully. Thanks and bye. :-) My ID is " + myID + ".");
                        try {
                            aliveCheck.done();
                            break WorkLoop;
                        } catch (Exception e) {
                            System.err.println("WARNING: Couldn't properly shut the "
                                    + "AliveCheck down. Exiting now. " + e.toString() + " My ID is " + myID + ".");
                            e.printStackTrace(System.err);
                            System.exit(0);
                        }
                    } else if (whatNext == SNAFU) {
                        System.err.println("ERROR: Shutting down ungracefully, server required so. My ID is " + myID + ".");
                        System.exit(33);
                    } else {
                        // unknown state, once more
                        System.err.println("ERROR: Unknown return code after handing "
                                + "result over. Aborting. My ID is " + myID + ".");
                        System.exit(55);
                    }

                    try{
                        //XXX ugly!
                        Thread.sleep(taskWaitTime);
                    } catch(Exception e){
                        System.err.println("ERROR: Client doesn't go to sleep. Aborting. " + e.toString() + " My ID is " + myID + ".");
                        e.printStackTrace(System.err);
                        System.exit(390);
                    }
                }
            } else if(whatNext == FINISH){
                System.out.println("INFO: Shutting down gracefully. Thanks and bye. :-) My ID is " + myID + ".");
                try{
                    aliveCheck.done();
                    break WorkLoop;
                } catch(Exception e){
                    System.err.println("WARNING: Couldn't properly shut the " +
                            "AliveCheck down. Exiting now. " + e.toString() + " My ID is " + myID + ".");
                    e.printStackTrace(System.err);
                    System.exit(0);
                }
            } else if(whatNext == SNAFU){
                System.err.println("ERROR: Shutting down ungracefully, server required so. My ID is " + myID + ".");
                System.exit(33);
            } else{
                // unknown state, once more
                System.err.println("ERROR: Unknown return code after handing " +
                        "result over. Aborting. My ID is " + myID + ".");
                System.exit(55);
            }
        }
     }
}
