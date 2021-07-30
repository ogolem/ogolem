/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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
import java.io.IOException;
import java.text.DecimalFormat;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.mpi.GenericMPIOptimization;
import org.ogolem.generic.mpi.MPIInterface;
import org.ogolem.helpers.Tuple;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.rmi.GenericGlobOptJob;

/**
 * The entry point for any batch-mode, massively parallel usage of the OGOLEM program suite in
 * cluster global optimization mode.
 *
 * @author Johannes Dieterich
 * @version 2020-12-23
 */
public class MainOGOLEM {

  /**
   * The alpha and omega of the OGOLEM program suite. ;-)
   *
   * @param args Arguments, where the first parameter is the input file with {@code .ogo} as an
   *     extension.
   */
  public static void main(String[] args) {
    /*
     * The MPI standard does NOT define, how a MPI-program handles any statements BEFORE
     * the MPI.Init() call.
     * Therefore, out of safety, nothing should take place here!
     */
    MPIInterface mpi = null;
    try {
      mpi = new MPIInterface("ogompi", false);
    } catch (Exception e) {
      e.printStackTrace(System.err);
      System.exit(21);
    }

    final int retInit = mpi.mpiInit(args);
    if (retInit != 0) {
      final int retAbort = mpi.mpiAbort(2);
      ;
      if (retAbort != 0) {
        System.err.println("Failure to MPI abort (init)");
        System.exit(2);
      }
    }

    /*
     * MPI related variables
     */

    int myRank = -1;
    try {
      myRank = mpi.mpiCommRank();
    } catch (Exception e) {
      e.printStackTrace(System.err);
      final int retAbort = mpi.mpiAbort(2);
      if (retAbort != 0) {
        System.err.println("Failure to abort in rank.");
        System.exit(2);
      }
    }

    int noOfProcesses = -1;
    try {
      noOfProcesses = mpi.mpiCommSize();
    } catch (Exception e) {
      e.printStackTrace(System.err);
      final int retAbort = mpi.mpiAbort(2);
      if (retAbort != 0) {
        System.err.println("Failure to abort in size.");
        System.exit(2);
      }
    }

    System.out.println("DEBUG: Starting process " + myRank + " of " + noOfProcesses + ".");

    if (noOfProcesses <= 1) {
      // this is not good
      System.err.println("ERROR: To little processes! " + noOfProcesses);
      final int retAbort = mpi.mpiAbort(2);
      if (retAbort != 0) {
        System.err.println("Failure to abort in #processes.");
        System.exit(2);
      }
    }

    if (myRank == 0) {
      // MPI rank 0 == the queen
      /*
       * "Serial" parts: Setup
       */
      // get the thing running part 0a: checking

      final double startTime = mpi.mpiWtime();

      final String inputFile = args[0];
      if (inputFile.equalsIgnoreCase("")) {
        System.err.println("Please specify an input file!");
        final int retAbort = mpi.mpiAbort(2);
        if (retAbort != 0) {
          System.err.println("Failure to MPI abort (input file)");
          System.exit(2);
        }
      } else if (inputFile.endsWith(".ogo") != true) {
        System.err.println("Please specify a correct input file!");
        final int retAbort = mpi.mpiAbort(17);
        if (retAbort != 0) {
          System.err.println("Failure to MPI abort (correct input)");
          System.exit(17);
        }
      }

      // get the thing running part 0b: set file informations
      final Tuple3D<String, String, String> dirs =
          ManipulationPrimitives.outDirAndBaseName(inputFile);
      final String outputFolder = dirs.getObject2();

      final String outputFile = outputFolder + File.separator + outputFolder + ".out";

      // get the thing running part 1: get configuration
      GlobalConfig globConf = null;
      try {
        globConf = Input.ConfigureMe(inputFile);
      } catch (Exception e) {
        System.err.println("Failure in creating in global configuration!" + e.toString());
        e.printStackTrace(System.err);
        final int retAbort = mpi.mpiAbort(1);
        if (retAbort != 0) {
          System.err.println("Failure to MPI abort (config create)");
          System.exit(1);
        }
      }
      assert (globConf != null);
      globConf.outputFile = outputFile;
      globConf.outputFolder = outputFolder;

      // check the input for sanity
      final Tuple<Boolean, String[]> sanity = InputSanityCheck.isInputSane(globConf);
      if (!sanity.getObject1()) {
        System.err.println("ERROR: Your configuration is insane, I am afraid. Exiting.");
        for (String s : sanity.getObject2()) {
          System.err.println(s);
        }
        final int retAbort = mpi.mpiAbort(2);
        if (retAbort != 0) {
          System.err.println("Failure to MPI abort (config)");
          System.exit(2);
        }
      }

      final boolean checkMove = !globConf.restart;
      if (checkMove) {
        try {
          ManipulationPrimitives.moveOldFolders(outputFolder, 9, System.getProperty("user.dir"));
        } catch (Exception e) {
          System.err.println("ERROR: Failure in initial file/folder operations. " + e.toString());
          e.printStackTrace(System.err);
          final int retAbort = mpi.mpiAbort(21);
          if (retAbort != 0) {
            System.err.println("Failure to MPI abort (initial file/folder)");
            System.exit(21);
          }
        }

        // create folder
        try {
          OutputPrimitives.createAFolder(outputFolder);
        } catch (IOException e) {
          System.err.println("ERROR: Couldn't create output folder!");
          e.printStackTrace(System.err);
          final int retAbort = mpi.mpiAbort(21);
          if (retAbort != 0) {
            System.err.println("Failure to create output folder.");
            System.exit(21);
          }
        }
      }

      final double setupTime = mpi.mpiWtime();

      try {
        final GenericPool<Molecule, Geometry> pool =
            GenericPool.getInstance(globConf.getGenericPoolConfig(), globConf.getExample());
        final GenericHistory<Molecule, Geometry> history =
            GenericHistory.getReference(globConf.getGenericHistoryConfig());

        final GenericGlobOptJob<Molecule, Geometry> job =
            new GenericGlobOptJob<>(
                globConf,
                pool,
                history,
                globConf.getReader(),
                globConf.getWriter(globConf.outputFolder),
                globConf.outputFolder,
                globConf.outputFile);
        GenericMPIOptimization.runAsQueen(mpi, job);
      } catch (Exception e) {
        System.err.println("Couldn't run queen part of the job.");
        e.printStackTrace(System.err);
        final int retAbort = mpi.mpiAbort(42);
        if (retAbort != 0) {
          System.err.println("Failure to MPI abort (queen opt)");
          System.exit(42);
        }
      }

      final double finalTime = mpi.mpiWtime();

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
      final int retFinal = mpi.mpiFinalize();
      if (retFinal != 0) {
        System.err.println("Couldn't finalize the MPI processes. Aborting.");
        final int retAbort = mpi.mpiAbort(1000);
        if (retAbort != 0) {
          System.err.println("Failure to MPI abort (finalize)");
          System.exit(1000);
        }
      }
    }

    /*
     * Actual parallel part
     */
    else {

      final long waittime = 1000; // 1 s wait time after being told "wait"

      try {
        GenericMPIOptimization.runAsDrone(mpi, waittime);
      } catch (Exception e) {
        System.err.println("Exception in drone process " + myRank);
        e.printStackTrace(System.err);
        final int retAbort = mpi.mpiAbort(99);
        if (retAbort != 0) {
          System.err.println("Failure to MPI abort (generic opt)");
          System.exit(99);
        }
      }

      /*
       * DONE, END, FINISHED!!!
       */
      final int retFinal = mpi.mpiFinalize();
      if (retFinal != 0) {
        System.err.println("Couldn't finalize the MPI processes. Aborting.");
        final int retAbort = mpi.mpiAbort(1000);
        if (retAbort != 0) {
          System.err.println("Failure to MPI abort (finalize)");
          System.exit(1000);
        }
      }
    }
  }
}
