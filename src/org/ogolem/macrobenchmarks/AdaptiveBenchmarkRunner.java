/*
Copyright (c) 2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.macrobenchmarks;

import java.io.File;
import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.io.InputPrimitives;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Run adaptive optimization benchmarks.
 *
 * @author Johannes Dieterich
 * @version 2020-04-19
 */
class AdaptiveBenchmarkRunner implements BenchmarkRunner {

  private static final Logger LOG = LoggerFactory.getLogger(ClusterBenchmarkRunner.class);

  @Override
  public boolean runBenchmark(String csvLine, int noThreads) throws Exception {

    final String[] csvData = csvLine.trim().split("\\,");
    // format of the line is as follows:
    // 1) the working directory the input resides in
    // 2) the input file, relative to work dir
    // 4) the expected global minimum fitness
    // 5) the file, relative to work dir, containing the reference global minimum
    // 6) threshold for fitness
    // 7) threshold per parameter
    final String workDir = csvData[0].trim();
    final String inputFile = csvData[1].trim();
    final double energy = Double.parseDouble(csvData[2].trim());
    final String refParameters = csvData[3].trim();
    final double fitnessThresh = Double.parseDouble(csvData[4].trim());
    final double paramThresh = Double.parseDouble(csvData[5].trim());

    final File newDir = new File(workDir);
    if (!newDir.exists() || !newDir.isDirectory()) {
      LOG.error(
          "Specified working directory " + workDir + " does not exist or is not a directory.");
      return false;
    }

    LOG.info("Executing benchmark " + inputFile);
    final double startT = System.currentTimeMillis();

    final String[] args = new String[] {"--adaptive", inputFile, "" + noThreads};
    try {
      Helpers.executeJavaProcess(workDir, args);
    } catch (Exception e) {
      throw new RuntimeException("Failure to run benchmark " + csvLine, e);
    }
    final double endT = System.currentTimeMillis();

    // we know the output directory is inputFile w/o .ogo extension
    final String outputDir =
        workDir + File.separator + inputFile.substring(0, inputFile.lastIndexOf(".ogo"));

    // find the rank0 parameters in that output directory
    final File outDir = new File(outputDir);
    final String[] ls = outDir.list(new Helpers.Rank0Filter());
    if (ls == null || ls.length != 1) {
      throw new RuntimeException(
          "No rank0 parameter set for benchmark " + csvLine + " in " + outputDir);
    }

    // let's get two adaptive parameter sets
    final String[] refParams = InputPrimitives.readFileIn(workDir + File.separator + refParameters);
    final AdaptiveParameters refSet = new AdaptiveParameters(refParams, 0);
    final String[] foundParams = InputPrimitives.readFileIn(outputDir + File.separator + ls[0]);
    final AdaptiveParameters foundSet = new AdaptiveParameters(foundParams, 1);

    final double refFit = refSet.getFitness();
    final double foundFit = foundSet.getFitness();
    if (Math.abs(foundFit - refFit) > fitnessThresh) {
      LOG.error("Benchmark fitness wrong. Obtained: " + foundFit + " vs should be " + refFit);
      return false;
    }

    // steps to that individual
    final String ls0NoFront = ls[0].substring("rank0individual".length());
    final long stepsToIndividual = Long.parseLong(ls0NoFront);

    // now the more complex thing - compare parameter values
    final String[] keys = refSet.getAllKeysCopy();
    final String[] oKeys = foundSet.getAllKeysCopy();
    if (keys.length != oKeys.length) {
      LOG.error("Different number of keys for reference and obtained.");
      return false;
    }
    for (final String key : keys) {

      final double[] refParamKey = refSet.getParametersForKey(key);
      final double[] foundParamKey = foundSet.getParametersForKey(key);

      if (refParamKey.length != foundParamKey.length) {
        LOG.error(
            "For key "
                + key
                + " found "
                + foundParamKey.length
                + " should be "
                + refParamKey.length);
        return false;
      }

      for (int i = 0; i < refParamKey.length; i++) {
        if (Math.abs(refParamKey[i] - foundParamKey[i]) > paramThresh) {
          LOG.error(
              "Parameter value " + i + " is " + foundParamKey[i] + " should be " + refParamKey[i]);
          return false;
        }
      }
    }

    LOG.info("Benchmark " + inputFile + " PASSED");
    LOG.info("Benchmark " + inputFile + " took: " + (endT - startT) / 1000 + " s.");
    LOG.info("Benchmark " + inputFile + " took: " + stepsToIndividual + " steps.");

    return true;
  }
}
