/*
Copyright (c) 2014, J. M. Dieterich
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
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Calls ADF (as in: the suite) for geometry optimizations.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class ADFCaller extends AbstractLocOpt {

  private static final long serialVersionUID = (long) 20141027;

  private final boolean DEBUG = true;
  private final String inputStub;
  private final String adfbinDir;

  ADFCaller(final GlobalConfig globConf, final String inputStub) throws Exception {
    super(globConf);
    this.inputStub = inputStub;
    this.adfbinDir = System.getenv("ADFBIN");
    if (adfbinDir == null || adfbinDir.isEmpty()) {
      throw new Exception("ADFBIN environment variable must be populated!");
    }
    // this.DEBUG = (GlobalConfig.DEBUGLEVEL > 0);
  }

  ADFCaller(final ADFCaller orig) {
    super(orig);
    this.inputStub = orig.inputStub;
    this.adfbinDir = orig.adfbinDir;
  }

  @Override
  public ADFCaller copy() {
    return new ADFCaller(this);
  }

  @Override
  protected String myID() {
    return "ADF (the suite) interface";
  }

  @Override
  public String myIDandMethod() {
    return "ADF (the suite) interface using " + this.inputStub;
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      final long lID,
      CartesianCoordinates cartes,
      final boolean[][] baConstraints,
      final boolean isConstricted,
      final BondInfo bonds)
      throws Exception {

    // create a "scratch" directory
    final String dir = "adf_id_" + lID + "_" + cartes.hashCode();
    OutputPrimitives.createAFolder(dir);

    // write the cartesian coordinates
    final String[] xyzData = cartes.createPrintableCartesians();
    OutputPrimitives.writeOut(dir + File.separator + "geom.xyz", xyzData, false);

    // if the stub ends in .adf, it is NOT a built-in and must be linked
    if (inputStub.endsWith(".adf")) {
      OutputPrimitives.createLink(inputStub, dir + File.separator + inputStub);
    }

    // we do the actual calls into ADF in a try/catch as there is more potential for problems and we
    // want to spit more output out
    Process procPrep = null;
    Process procRun = null;
    Process procRep = null;
    String[] xyzEnergyData = null;
    try {

      final Runtime rt = Runtime.getRuntime();
      final File direc = new File(dir);

      // first prepare the input using adfprep
      final String comPrep =
          adfbinDir + File.separator + "adfprep -t " + inputStub + " -m geom.xyz -j adfjob";
      procPrep = rt.exec(comPrep, null, direc);

      // prep gobblers
      final StreamGobbler errorPrep = new StreamGobbler(procPrep.getErrorStream(), "ERROR");
      final StreamGobbler outputPrep = new StreamGobbler(procPrep.getInputStream(), "OUTPUT");
      errorPrep.start();
      outputPrep.start();

      final int codePrep = procPrep.waitFor();
      if (codePrep != 0) {
        if (DEBUG) {
          System.out.println("PREP OUTPUT");
          for (final String s : outputPrep.getData()) {
            System.out.println(s);
          }
          System.out.println("PREP ERROUT");
          for (final String s : errorPrep.getData()) {
            System.out.println(s);
          }
        }
        throw new ConvergenceException("adfprep returns an error code " + codePrep);
      }

      // get the output and write it to a file
      final String[] inpData = outputPrep.getData();
      final String jobFile = dir + File.separator + "jobfile";
      OutputPrimitives.writeOut(jobFile, inpData, false);

      // now execute said jobfile
      ManipulationPrimitives.markExecutable(jobFile);
      final String comRun = "." + File.separator + "jobfile";
      procRun = rt.exec(comRun, null, direc);

      // run gobblers
      final StreamGobbler errorRun = new StreamGobbler(procRun.getErrorStream(), "ERROR");
      final StreamGobbler outputRun = new StreamGobbler(procRun.getInputStream(), "OUTPUT");
      errorRun.start();
      outputRun.start();

      final int codeRun = procRun.waitFor();
      if (codeRun != 0) {
        if (DEBUG) {
          System.out.println("RUN OUTPUT");
          for (final String s : outputRun.getData()) {
            System.out.println(s);
          }
          System.out.println("RUN ERROUT");
          for (final String s : errorRun.getData()) {
            System.out.println(s);
          }
        }
        throw new ConvergenceException("Running the ADF jobfile returns an error code " + codeRun);
      }

      // figure out which results file to take
      String resultsFile = null;
      final File t21File = new File(dir + File.separator + "adfjob.t21");
      final File rkfFile = new File(dir + File.separator + "adfjob.rkf");
      final File runkfFile = new File(dir + File.separator + "adfjob.runkf");
      final File rxkfFile = new File(dir + File.separator + "adfjob.rxkf");
      if (t21File.exists()) {
        resultsFile = "adfjob.t21";
      } else if (rkfFile.exists()) {
        resultsFile = "adfjob.rkf";
      } else if (runkfFile.exists()) {
        resultsFile = "adfjob.runkf";
      } else if (rxkfFile.exists()) {
        resultsFile = "adfjob.rxkf";
      } else {
        throw new Exception(
            "Neither a t21 file, a runkf file, a rxkf nor a rkf file with results found.");
      }

      // get the results, energy and optimized geometry
      final String comRep =
          adfbinDir + File.separator + "adfreport " + resultsFile + " geometry-b energy -plain";
      procRep = rt.exec(comRep, null, direc);

      // run gobblers
      final StreamGobbler errorRep = new StreamGobbler(procRep.getErrorStream(), "ERROR");
      final StreamGobbler outputRep = new StreamGobbler(procRep.getInputStream(), "OUTPUT");
      errorRep.start();
      outputRep.start();

      final int codeRep = procRep.waitFor();
      if (codeRep != 0) {
        if (DEBUG) {
          System.out.println("REP OUTPUT");
          for (final String s : outputRep.getData()) {
            System.out.println(s);
          }
          System.out.println("REP ERROUT");
          for (final String s : errorRep.getData()) {
            System.out.println(s);
          }
        }
        throw new ConvergenceException(
            "Getting results from file " + resultsFile + " returns an error code " + codeRep);
      }

      // get that data
      xyzEnergyData = outputRep.getData();

    } catch (Exception e) {
      System.err.println(e.toString());
      throw new ConvergenceException(
          "ADF calling sequence has a problem (local optimization). For ID " + lID, e);
    } finally {
      if (procPrep != null) {
        procPrep.destroy();
      }
      if (procRun != null) {
        procRun.destroy();
      }
      if (procRep != null) {
        procRep.destroy();
      }
    }

    // create the output cartesian
    final CartesianCoordinates result = cartes.copy();

    // parse the coords first (are in Bohr)
    final double[][] xyz = result.getAllXYZCoord();
    for (int at = 0; at < cartes.getNoOfAtoms(); at++) {
      final String[] data = xyzEnergyData[at].trim().split("\\s+");
      xyz[0][at] = Double.parseDouble(data[1]);
      xyz[1][at] = Double.parseDouble(data[2]);
      xyz[2][at] = Double.parseDouble(data[3]);
    }

    // now try to get the energy
    final String[] eLine = xyzEnergyData[cartes.getNoOfAtoms()].trim().split("\\s+");
    final double e = Double.parseDouble(eLine[0]);
    result.setEnergy(e);

    // delete the directory
    ManipulationPrimitives.remove(dir);

    return result;
  }
}
