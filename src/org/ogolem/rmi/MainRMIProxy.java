/*
Copyright (c) 2014, J. M. Dieterich
              2015-2022, J. M. Dieterich and B. Hartke
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
import java.rmi.Remote;
import java.rmi.registry.LocateRegistry;
import java.rmi.registry.Registry;
import java.rmi.server.UnicastRemoteObject;
import org.ogolem.helpers.Tuple;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Proxy function for the RMI parallelization.
 *
 * @author Johannes Dieterich
 * @version 2022-12-16
 */
public class MainRMIProxy {

  private static final boolean DEBUG = false;

  // to disallow the GC to collect this object, it is made a class variable
  private static RMICommunication<?> comm;
  private static Registry registry;

  /**
   * The proxy side of the RMI part. Requires as input what should be optimized, followed by the
   * input file name, the timeout and the job timeout (both in seconds), and the port on which the
   * proxy is supposed to export the remote object and the port on which the rmiregistry is running.
   * Additionally, it needs the server name, a sleeping time and the port of the server.
   *
   * @param args
   */
  public static void main(String[] args) {

    // XXX fix...
    if (args == null || args[0].equalsIgnoreCase("--help")) {
      System.out.println("Sorry, no help yet...");
      System.exit(0);
    }

    String exportName = "RMICommunication";
    int port = 2500;
    int reg = 1099;
    String serverName = null;
    long sleepTime = 5000;
    int serverPort = 1099;
    int maxTasks = 1000;
    int maxIndiesExchange = -1;
    String jobType = null;
    String inputFile = null;
    long timeOut = 10000l;
    long timeOutJob = 10000l;

    for (int argID = 0; argID < args.length; argID++) {
      final String arg = args[argID];
      if (arg.equalsIgnoreCase("--server")) {
        argID++;
        serverName = args[argID];
      } else if (arg.equalsIgnoreCase("--jobtype")) {
        argID++;
        jobType = args[argID];
      } else if (arg.equalsIgnoreCase("--inputfile")) {
        argID++;
        inputFile = args[argID];
      } else if (arg.equalsIgnoreCase("--commname")) {
        argID++;
        exportName = args[argID];
      } else if (arg.equalsIgnoreCase("--sleeptime")) {
        argID++;
        sleepTime = Long.parseLong(args[argID]) * 1000;
      } else if (arg.equalsIgnoreCase("--serverregistryport")) {
        argID++;
        serverPort = Integer.parseInt(args[argID]);
      } else if (arg.equalsIgnoreCase("--myserverport")) {
        argID++;
        port = Integer.parseInt(args[argID]);
      } else if (arg.equalsIgnoreCase("--registryport")) {
        argID++;
        reg = Integer.parseInt(args[argID]);
      } else if (arg.equalsIgnoreCase("--mergeinds")) {
        argID++;
        maxIndiesExchange = Integer.parseInt(args[argID]);
      } else if (arg.equalsIgnoreCase("--maxchunk")) {
        argID++;
        maxTasks = Integer.parseInt(args[argID]);
      } else if (arg.equalsIgnoreCase("--timeout")) {
        argID++;
        timeOut = Long.parseLong(args[argID]) * 1000;
      } else if (arg.equalsIgnoreCase("--timeoutjob")) {
        argID++;
        timeOutJob = Long.parseLong(args[argID]) * 1000;
      } else {
        System.err.println("Unknown argument " + arg + ". Shutting down.");
        System.exit(1);
      }
    }

    // check for mandatory settings
    if (serverName == null) {
      System.err.println("No server specified.");
      System.exit(2);
    }
    if (jobType == null) {
      System.err.println("No jobType specified.");
      System.exit(3);
    }
    if (inputFile == null) {
      System.err.println("No inputFile specified.");
      System.exit(4);
    }

    String keySuffix = System.getenv("OGO_RMIKEY");
    if (keySuffix == null) {
      keySuffix = "Super secret password: Nakatomi Socrates";
    }

    final String key = "Client speaking, I am here. " + keySuffix;

    // first we try to connect to our server and exchange keys
    final long waittimeClientMS = 2000; // two seconds.
    RMICommunication<?> serverComm = null;
    int myID = -1;
    for (int attempt = 0; attempt < 3; attempt++) {
      try {
        final Registry serverReg = LocateRegistry.getRegistry(serverName, serverPort);

        if (DEBUG) {
          final String[] sa = registry.list();
          System.out.println("DEBUG: List of remote objects in registry:");
          for (final String s : sa) {
            System.out.println("DEBUG: " + s);
          }
        }
        serverComm = (RMICommunication<?>) serverReg.lookup(exportName);

        if (serverComm == null) {
          System.err.println("WARNING: Communication object is null. Trying again....");
        }

        // exchange keys
        final Tuple<String, Integer> answer = serverComm.registerWithServer(key);
        if (!answer
            .getObject1()
            .equalsIgnoreCase("Server speaking, everything fine. " + keySuffix)) {
          // not good, exit (we want to avoid double registration...)
          System.err.println(
              "ERROR: Wrong answer from server. Aborting client. Server was "
                  + serverName
                  + ". Answer was "
                  + answer.getObject1());
          System.exit(27);
        } else {
          myID = answer.getObject2();
          // yup, ID assigned!
          break;
        }
      } catch (Exception e) {
        System.err.println(
            "RMI exception: " + e.getMessage() + ". Trying again... Server is " + serverName);
        e.printStackTrace(System.err);
      }

      // sleep for a bit...
      try {
        Thread.sleep(waittimeClientMS);
      } catch (Exception e) {
        // whatever...
        e.printStackTrace(System.err);
      }

      if (attempt == 3) {
        System.err.println(
            "ERROR: Failure to connect to initial server happened three times. Giving up...");
        System.exit(42);
      }
    }

    // create, if necessary, the output directory
    final Tuple3D<String, String, String> dirs =
        ManipulationPrimitives.outDirAndBaseName(inputFile);
    final String outFolder = dirs.getObject2();

    try {
      ManipulationPrimitives.moveOldFolders(outFolder, 9, System.getProperty("user.dir"));
    } catch (Exception e) {
      System.err.println("ERROR: Failure in initial file/folder operations. " + e.toString());
      e.printStackTrace(System.err);
      System.exit(1);
    }

    // create folder
    try {
      OutputPrimitives.createAFolder(outFolder);
    } catch (IOException e) {
      System.err.println("ERROR: Couldn't create output folder!");
      e.printStackTrace(System.err);
      System.exit(2);
    }

    // set everything up
    final long offsetTime = System.currentTimeMillis();

    // the check for the aliveness of the server
    final AliveCheck aliveCheck = new AliveCheck(serverComm, sleepTime, myID);

    comm =
        JobFactory.constructJob(
            serverComm,
            jobType,
            inputFile,
            timeOut,
            timeOutJob,
            maxTasks,
            maxIndiesExchange,
            aliveCheck,
            myID);

    try {
      registry = LocateRegistry.createRegistry(reg);
      final Remote stub = UnicastRemoteObject.exportObject(comm, port);
      registry.bind("RMICommunication", stub);
    } catch (Exception e) {
      System.out.println("OGOLEM RMI proxy exception: " + e.getMessage());
      e.printStackTrace(System.err);
      System.exit(333);
    }

    final long setupTime = System.currentTimeMillis();

    boolean everythingDone = false;
    while (!everythingDone) {
      try {
        // XXX ugly!
        Thread.sleep(5000);
      } catch (Exception e) {
        System.err.println("WARNING: Can't wait before checking again the server." + e.toString());
        e.printStackTrace(System.err);
      }
      try {
        everythingDone = comm.isEverythingDone();
        if (everythingDone) {
          comm.stopClientList();
        }
      } catch (Exception e) {
        System.err.println(
            "ERROR: Can't ask the server how it is doing. "
                + "Breaking loop and exiting. "
                + e.toString());
        e.printStackTrace(System.err);
        break;
      }
    }

    try {
      registry.unbind("RMICommunication");
    } catch (Exception e) {
      System.err.println(
          "ERROR: Couldn't unexport the remote object. This " + "needs to be fixed by hand.");
      e.printStackTrace(System.err);
    }

    final long totalTime = System.currentTimeMillis();
    System.out.println("PROXY ACCOUNTING INFORMATION:");
    System.out.println("Time for setup: " + (setupTime - offsetTime) / 1000 + "s.");
    System.out.println("Time for job: " + (totalTime - setupTime) / 1000 + "s.");
    System.out.println("Total time: " + (totalTime - offsetTime) / 1000 + "s.");
    System.out.println("Thank you for running OGOLEM. We hope it was successfull!");
    System.out.println("");
    System.out.println(org.ogolem.helpers.Fortune.randomFortune());
    // this needs to be here since it doesn't return on itself (very robust server ;-) )
    System.exit(0);
  }
}
