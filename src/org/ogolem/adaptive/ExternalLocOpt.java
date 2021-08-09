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
package org.ogolem.adaptive;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.StreamGobbler;
import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.GenericFitnessBackend;
import org.ogolem.generic.GenericLocOpt;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

public class ExternalLocOpt implements GenericLocOpt<Double, AdaptiveParameters> {

  private static final long serialVersionUID = (long) 20200725;

  private final String extCmd;
  private final GenericBackend<Double, AdaptiveParameters> backend;

  ExternalLocOpt(final GenericBackend<Double, AdaptiveParameters> backend) {
    assert (backend != null);
    this.backend = backend;
    final String cmd = System.getenv("OGO_EXTCALLCMD");
    this.extCmd = (cmd == null) ? "ogo_external" : cmd;
  }

  ExternalLocOpt(final ExternalLocOpt orig) {
    this.backend = orig.backend.copy();
    this.extCmd = orig.extCmd;
  }

  @Override
  public ExternalLocOpt copy() {
    return new ExternalLocOpt(this);
  }

  @Override
  public String getMyID() {
    return "external locopt";
  }

  @Override
  public AdaptiveParameters fitness(AdaptiveParameters params, boolean forceOneEval) {

    AdaptiveParameters res = params.copy();
    res.setFitness(FixedValues.NONCONVERGEDENERGY);

    // create a directory and write current parameters out
    final String dirName = "externalcaller-" + params.getID() + "_" + System.currentTimeMillis();
    final String paramsFile = dirName + File.separator + "input.params";
    final String paramsRes = dirName + File.separator + "output.params";
    try {
      OutputPrimitives.createAFolder(dirName);
      final String[] paramsString = params.createFormattedOutput();
      OutputPrimitives.writeOut(paramsFile, paramsString, false);
    } catch (IOException e) {
      e.printStackTrace(System.err);
      return res;
    }

    final List<String> cmdList = new ArrayList<>();
    cmdList.add(extCmd);
    final String mode = (forceOneEval) ? "single" : "locopt";
    cmdList.add(mode);

    final ProcessBuilder pb = new ProcessBuilder(cmdList);
    pb.directory(new File(dirName));

    try {

      final Process proc = pb.start();

      // any error message?
      final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      final int errCode = proc.waitFor();
      if (errCode != 0) {
        throw new Exception(
            "external locopt returns non-zero return value (local optimization). Error code "
                + errCode);
      }
    } catch (Exception e) {
      e.printStackTrace(System.err);
      return res;
    }
    // cmd should(!) have completed normally...

    // read new parameters and their fitness in and return the object
    try {

      final String[] resData = InputPrimitives.readFileIn(paramsRes);
      final AdaptiveParameters inter = new AdaptiveParameters(resData, params.getID());

      // do copy over to make sure everything fits
      final double[] resP = res.getAllParamters();
      final double[] interP = inter.getAllParamters();

      if (resP.length != interP.length) {
        throw new Exception(
            "Externally returned parameters have the wrong length. Should be "
                + resP.length
                + " is "
                + interP.length);
      }

      System.arraycopy(interP, 0, resP, 0, resP.length);
      res.setFitness(inter.getFitness());

    } catch (Exception e) {
      e.printStackTrace(System.err);
      res.setFitness(FixedValues.NONCONVERGEDENERGY);
      return res;
    }

    // clean up
    try {
      ManipulationPrimitives.remove(dirName);
    } catch (Exception e) {
      // not important....
      System.err.println("Exception during external locopt cleanup " + e.getMessage());
    }

    return res;
  }

  @Override
  public GenericFitnessBackend<Double, AdaptiveParameters> getFitnessBackend() {
    return backend;
  }

  @Override
  public GenericBackend<Double, AdaptiveParameters> getBackend() {
    return backend;
  }
}
