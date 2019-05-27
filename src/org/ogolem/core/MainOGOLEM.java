/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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
package org.ogolem.core;

import java.io.File;
import java.text.DecimalFormat;
import mpi.MPI;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.mpi.GenericMPIOptimization;
import org.ogolem.helpers.Tuple;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.rmi.GenericGlobOptJob;

/**
 * The entry point for any batch-mode, massively parallel usage of the OGOLEM
 * program suite in cluster global optimization mode.
 * @author Johannes Dieterich
 * @version 2014-04-06
 */
public class MainOGOLEM {

    /**
     * The alpha and omega of the OGOLEM program suite. ;-)
     * @param args Arguments, where the first parameter is the input file with {@code .ogo} as an extension.
     */
    public static void main(String[] args) {
        /*
         * The MPI standard does NOT define, how a MPI-program handles any statements BEFORE
         * the MPI.Init() call.
         * Therefore, out of safety, nothing should take place here!
         */
        try {
            MPI.Init(args);
        } catch (Exception e) {
            e.printStackTrace(System.err);
            try {
                MPI.COMM_WORLD.Abort(2);
            } catch (Exception e2) {
                e2.printStackTrace(System.err);
                System.exit(2);
            }
        }
        /*
         * MPI related variables
         */

        int myRank = -1;
        try {
            myRank = MPI.COMM_WORLD.Rank();
        } catch (Exception e) {
            e.printStackTrace(System.err);
            try {
                MPI.COMM_WORLD.Abort(2);
            } catch (Exception e2) {
                e2.printStackTrace(System.err);
                System.exit(2);
            }
        }


        int noOfProcesses = -1;
        try {
            noOfProcesses = MPI.COMM_WORLD.Size();
        } catch (Exception e) {
            e.printStackTrace(System.err);
            try {
                MPI.COMM_WORLD.Abort(2);
            } catch (Exception e2) {
                e2.printStackTrace(System.err);
                System.exit(2);
            }
        }

        System.out.println("DEBUG: Starting process " + myRank + " of " + noOfProcesses + ".");

        if (noOfProcesses <= 1) {
            // this is not good
            try {
                System.err.println("ERROR: To little processes! " + noOfProcesses);
                MPI.COMM_WORLD.Abort(2);
            } catch (Exception e2) {
                e2.printStackTrace(System.err);
                System.exit(2);
            }
        }

        if (myRank == 0) {
            // MPI master
			/*
             * "Serial" parts: Setup
             */
            // get the thing running part 0a: checking

            double startTime = 0.0;
            try {
                startTime = MPI.Wtime();
            } catch (Exception e) {
                e.printStackTrace(System.err);
            }

            /*
             * When using mpiJava (which is highly encouraged ATM), the arguments look like this:
             * args[0]: The main class, ogolem.core.MainOGOLEM.
             * args[1]: The actual "first" argument for OGOLEM, the input file.
             * args[2+]: Non interesting and MPI related stuff.
             * This actually depends on which MPI backend is used. Regularly it is
             * args[0]: The normal "first" argument for OGOLEM, the input file.
             */


            /*
             * This looks strange and is strange. Apparently the args are different
             * when using MPJexpress.
             * args[0]: the rank (at least I assume so, would fit)
             * args[1]: the link to the MPJe configuration file (inclusive port)
             * args[2]: "niodev"
             * args[3]: the actual args start from here
             */

            /*
             * And the next issue. If we want to use the gcj to compile to native code,
             * the order is (of course) again different
             * args[0] is the input file
             * Does this suck? Yes, it does!
             * And we figured out that at least at the present time the gcj does NOT work.
             */

            //XXX this needs to be adapted for every mpi runtime. too bad.
            // actually, MPI.Init should have a return value of String[] args for MPJe :-)
            final String inputFile = args[3];
            if (inputFile.equalsIgnoreCase("")) {
                System.err.println("Please specify an input file!");
                try {
                    MPI.COMM_WORLD.Abort(2);
                } catch (Exception e2) {
                    e2.printStackTrace(System.err);
                    System.exit(2);
                }
            } else if (inputFile.endsWith(".ogo") != true) {
                System.err.println("Please specify a correct input file!");
                try {
                    MPI.COMM_WORLD.Abort(17);
                } catch (Exception e2) {
                    e2.printStackTrace(System.err);
                    System.exit(17);
                }
            }

            // get the thing running part 0b: set file informations
            final int indexOfOgo = inputFile.indexOf(".ogo");
            final String outputFolder = inputFile.substring(0, indexOfOgo);
            final String outputFile = outputFolder + File.separator + outputFolder + ".out";


            // get the thing running part 1: get configuration
            GlobalConfig globConf = null;
            try {
                globConf = Input.ConfigureMe(inputFile);
            } catch (Exception e) {
                System.err.println("Failure in creating in global configuration!" + e.toString());
                e.printStackTrace(System.err);
                try {
                    MPI.COMM_WORLD.Abort(1);
                } catch (Exception e2) {
                    e2.printStackTrace(System.err);
                    System.exit(1);
                }
            }
            assert(globConf != null);
            globConf.outputFile = outputFile;
            globConf.outputFolder = outputFolder;

            // check the input for sanity
            final Tuple<Boolean, String[]> sanity = InputSanityCheck.isInputSane(globConf);
            if (!sanity.getObject1()) {
                System.err.println("ERROR: Your configuration is insane, I am afraid. Exiting.");
                for (String s : sanity.getObject2()) {
                    System.err.println(s);
                }
                try{
                    MPI.COMM_WORLD.Abort(2);
                } catch(Exception e){
                    e.printStackTrace(System.err);
                    System.exit(2);
                }
            }

            final boolean checkMove = !globConf.restart;
            if(checkMove){
                try{
                    ManipulationPrimitives.moveOldFolders(outputFolder, 9, System.getProperty("user.dir"));
                } catch (Exception e) {
                    System.err.println("ERROR: Failure in initial file/folder operations. " + e.toString());
                    e.printStackTrace(System.err);
                    try {
                        MPI.COMM_WORLD.Abort(21);
                    } catch (Exception e2) {
                        e2.printStackTrace(System.err);
                        System.exit(21);
                    }
                }
            }
            
            double setupTime = 0.0;
            try {
                setupTime = MPI.Wtime();
            } catch (Exception e) {
                e.printStackTrace(System.err);
            }
            
            try{
                final GenericPool<Molecule,Geometry> pool = GenericPool.getInstance(globConf.getGenericPoolConfig(),
                    globConf.getExample());
                final GenericHistory<Molecule,Geometry> history = GenericHistory.getReference(globConf.getGenericHistoryConfig());
        
                final GenericGlobOptJob<Molecule,Geometry> job = new GenericGlobOptJob<>(globConf,
                    pool,history, globConf.getReader(),globConf.getWriter(globConf.outputFolder),
                    globConf.outputFolder,globConf.outputFile);
                GenericMPIOptimization.runAsMaster(job);
            } catch(Exception e){
                System.err.println("Couldn't run master part of the job.");
                e.printStackTrace(System.err);
                try {
                    MPI.COMM_WORLD.Abort(42);
                } catch (Exception e2) {
                    e2.printStackTrace(System.err);
                    System.exit(42);
                }
            }

            double finalTime = 0.0;
            try {
                finalTime = MPI.Wtime();
            } catch (Exception e) {
                e.printStackTrace(System.err);
            }

            final String[] accounting = new String[5];
            accounting[0] = "******************";
            final DecimalFormat format = new DecimalFormat("0.00");
            accounting[1] = "TIMING INFORMATIONS IN SECONDS";
            accounting[2] = "Setup: " + format.format(setupTime - startTime);
            accounting[3] = "Overall: " + format.format(finalTime - startTime);
            accounting[4] = "******************";

            try {
                OutputPrimitives.writeOut(globConf.outputFile, accounting, true);
            } catch (Exception e) {
                System.err.println("Failure in writing to output file." + e.toString());
            }
            
            /*
             * DONE, END, FINISHED!!!
             */
            try {
                MPI.Finalize();
            } catch (Exception e) {
                System.err.println("Couldn't finalize the MPI processes. Aborting.");
                e.printStackTrace(System.err);
                try {
                    MPI.COMM_WORLD.Abort(100);
                } catch (Exception e2) {
                    e2.printStackTrace(System.err);
                    System.exit(100);
                }
            }
            
        }
        /*
         * Actual parallel part
         */
        else {

            final long timeout = 30000; // 30 s timeout

            try {
                GenericMPIOptimization.runAsSlave(timeout);
            } catch (Exception e) {
                System.err.println("Exception in slave process " + myRank);
                e.printStackTrace(System.err);
                try {
                    MPI.COMM_WORLD.Abort(99);
                } catch (Exception e2) {
                    e2.printStackTrace(System.err);
                    System.exit(99);
                }
            }

            /*
             * DONE, END, FINISHED!!!
             */
            try {
                MPI.Finalize();
            } catch (Exception e) {
                System.err.println("Couldn't finalize the MPI processes. Aborting.");
                e.printStackTrace(System.err);
                try {
                    MPI.COMM_WORLD.Abort(1000);
                } catch (Exception e2) {
                    e2.printStackTrace(System.err);
                    System.exit(1000);
                }
            }
        }
    }
}
