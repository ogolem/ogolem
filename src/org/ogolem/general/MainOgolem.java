/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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
package org.ogolem.general;

/**
 * The unified entry point to the OGOLEM framework.
 *
 * @author Johannes Dieterich
 * @version 2020-04-25
 */
public class MainOgolem {

  /**
   * Entry point for everything and everyone. World plus dog. Or no. No dogs allowed. ;-)
   *
   * @param args The arguments. The first defines which program part to use, the rest is part
   *     specific.
   */
  public static void main(final String[] args) {

    if (args == null || args.length == 0) {
      System.err.println("OGOLEM doesn't like to be called without arguments.");
      System.err.println("Try the --help switch.");
      System.exit(1);
    }

    // copy the args array
    final String[] strippedArgs = new String[args.length - 1];
    System.arraycopy(args, 1, strippedArgs, 0, strippedArgs.length);

    // TODO in the long run, add --neural option

    // check which program part to use
    if (args[0].equalsIgnoreCase("--adaptive")) {
      org.ogolem.adaptive.MainOgoAdaptive.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--atom2ogo")) {
      org.ogolem.atom2ogo.MainAtom2Ogo.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--beswitched")) {
      org.ogolem.beswitched.MainBeswitched.main(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--clusters")) {
      org.ogolem.clusters.MainClusters.main(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--core")) {
      org.ogolem.core.MainOGOLEM.main(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--corrfunc")) {
      org.ogolem.corrfunc.MainCorrFunctions.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--dimerizer")) {
      org.ogolem.dimerizer.MainDimerizer.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--dimeranalyzer")) {
      org.ogolem.dimeranalyzer.MainDimerAnalyzer.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--evbqmdff")) {
      org.ogolem.qmdff.MainEVBQMDFF.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--familytree")) {
      org.ogolem.familytree.MainFamilyTree.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--fft")) {
      org.ogolem.fft.MainFFTInterface.exec(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--folder")) {
      org.ogolem.corrfunc.MainCorrelationFolder.execute(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--freqs")) {
      org.ogolem.freqs.MainFreqs.execute(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--lionbench")) {
      org.ogolem.lionbench.MainLion.main(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--ljrefs")) {
      org.ogolem.ljreferences.MainLJRefs.main(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--lhp")) {
      org.ogolem.heat.MainLocalHeat.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--macrobenchmarks")) {
      org.ogolem.macrobenchmarks.MainMacroBenchmarks.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--md")) {
      org.ogolem.md.MainMolecularDynamics.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--microbenchmarks")) {
      org.ogolem.microbenchmarks.MainMicroBenchmarks.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--ogo2atom")) {
      org.ogolem.atom2ogo.MainOgo2Atom.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--parameters")) {
      org.ogolem.parameters.MainParameters.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--powerspec")) {
      org.ogolem.freqs.MainPowerSpec.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--qmdff")) {
      org.ogolem.qmdff.MainQMDFF.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--rmiserver")) {
      org.ogolem.rmi.MainOgolemRMI.main(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--rmiproxy")) {
      org.ogolem.rmi.MainRMIProxy.main(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--rmiclient")) {
      org.ogolem.rmi.MainRMIClient.main(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--rmithreader")) {
      org.ogolem.rmi.MainRMIThreader.main(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--scanalyzer")) {
      org.ogolem.scanalyzer.MainScAnalyzer.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--scanner")) {
      org.ogolem.scanner.MainScanner.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--simplecorrfuncs")) {
      org.ogolem.corrfunc.MainSimpleCorrFunctions.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--switches")) {
      org.ogolem.switches.MainSwitches.main(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--shmem")) {
      org.ogolem.clusters.MainThreadingClusterGlobOpt.execute(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--tests")) {
      org.ogolem.tests.MainTests.run(strippedArgs);
    } else if (args[0].equalsIgnoreCase("--version")) {
      System.out.println(whichVersion());
    } else if (args[0].equalsIgnoreCase("--help") || args[0].equalsIgnoreCase("-h")) {
      System.out.println("");
      System.out.println("OGOLEM is typically called as follows:");
      System.out.println("   java -jar ogolem.jar --SWITCH ARGS");
      System.err.println("");
      System.out.println("where --SWITCH is one of the following:");
      System.out.println("   --adaptive");
      System.out.println("     for global parameter optimization");
      System.out.println("   --aligner");
      System.out.println(
          "     for aligning two sets of cartesian coordinates using a global aligning as presented by M. Eriksson");
      System.out.println("   --atom2ogo");
      System.out.println("     for translating an atom-formatted file to xyz and ogolem bond info");
      System.out.println("   --beswitched");
      System.out.println("     for extracting intermediate switch optimization results");
      System.out.println("   --clusteroverlap");
      System.out.println("     for checking if two cluster structures are identical");
      System.out.println("   --clusterpca");
      System.out.println("     for doing principal component analysis on a set of clusters");
      System.out.println("   --core");
      System.out.println("     for MPI parallelized cluster structure optimization");
      System.out.println("   --corrfunc");
      System.out.println("     to compute correlation functions from MD trajectories");
      System.out.println("   --dimeranalyzer");
      System.out.println("     to analyze the dimerizer output in the crudest form possible");
      System.out.println("   --dimerizer");
      System.out.println("     to do a relaxed grid search for a dimer");
      System.out.println("   --familytree");
      System.out.println("     to create a family tree from the genetic history output");
      System.out.println("   --evbqmdff");
      System.out.println(
          "     to calculate a single energy and gradient for an xyz file with Hartke and Grimme's EVB-QMDFF");
      System.out.println("   --fft");
      System.out.println("     to (\"fast\") fourier transform simple 1D data");
      System.out.println("   --folder");
      System.out.println("     to fold previously computed correlation function");
      System.out.println("   --freqs");
      System.out.println("     for calculating vibrational frequencies");
      System.out.println("   --help");
      System.out.println("     to display this overview");
      System.out.println("   --graphmerger");
      System.out.println(
          "     for merging binary disconnectivity graphs as obtained by the lid subpackage");
      System.out.println("   --lid");
      System.out.println("     for running the lid/threshold algorithm by Schoen and Sibani");
      System.out.println("   --lidgrapher");
      System.out.println("     for preparing a binary disconnectivity graph for plotting");
      System.out.println("   --lionbench");
      System.out.println("     for testing and benchmarking OGOLEM");
      System.out.println("   --ljrefs");
      System.out.println("     for ejecting LJ reference structures (deprecated!)");
      System.out.println("   --lhp");
      System.out.println(
          "     for using local heat pulses as a global optimization (developer-only!)");
      System.out.println("   --md");
      System.out.println("     for running simple molecular dynamics");
      System.out.println("   --macrobenchmarks");
      System.out.println("     for running some macro benchmarks");
      System.out.println("   --microbenchmarks");
      System.out.println("     for running some micro benchmarks");
      System.out.println("   --molecules");
      System.out.println("     to optimize molecules with respect to a or multiple properties");
      System.out.println("   --ogo2atom");
      System.out.println(
          "     for translating a xyz file and ogolem bond info to an Atomdroid atom file");
      System.out.println("   --parameters");
      System.out.println("     for extracting intermediate parametrization results");
      System.out.println("   --powerspec");
      System.out.println("     for calculating a power spectrum from velocity autocorrelation");
      System.out.println("   --qmdff");
      System.out.println(
          "     to calculate a single energy and gradient for an xyz file with Grimme's QMDFF");
      System.out.println("   --rmiserver");
      System.out.println("     for starting an OGOLEM RMI server");
      System.out.println("   --rmiproxy");
      System.out.println("     for starting an OGOLEM RMI proxy");
      System.out.println("   --rmiclient");
      System.out.println("     for starting an OGOLEM RMI client");
      System.out.println("   --rmithreader");
      System.out.println("     for starting an OGOLEM RMI threading client");
      System.out.println("   --scanalyzer");
      System.out.println("     for analyzing the scanning output");
      System.out.println("   --scanner");
      System.out.println("     for scanning along degrees of freedom");
      System.out.println("   --shmem");
      System.out.println("     for thread-based cluster structure optimization");
      System.out.println("   --simplecorrfunc");
      System.out.println(
          "     for calculating simple auto and cross correlation functions (BETA!)");
      System.out.println("   --spectral");
      System.out.println("     for fitting spectra from a set of individual ones (BETA!)");
      System.out.println("   --switches");
      System.out.println("     for thread-based molecular switch optimization");
      System.out.println("   --tests");
      System.out.println("     for testing and development purposes");
      System.out.println("   --version");
      System.out.println("     for the OGOLEM snapshot version");
      System.out.println("");
      System.out.println(
          "and ARGS are the arguments (typical input file(s) and parallelization information)");
      System.out.println("for each part of OGOLEM. Check the manual for details.");
      System.out.println("");
      System.out.println(
          "additionally, some of the packages respond to the help flag, e.g., try --clusterpca help");
      System.out.println("");
    } else {
      System.err.println("Unfortunately OGOLEM does not know what to do with your request.");
      System.err.println("Try the --help switch.");
      System.exit(1);
    }
  }

  private static String whichVersion() {
    return "2016.CURRENT\n Autobuild revision: 1250M";
  }
}
