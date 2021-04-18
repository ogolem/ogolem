/*
Copyright (c) 2012-2014, J. M. Dieterich
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
package org.ogolem.adaptive;

import java.io.File;
import java.io.FilenameFilter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Constants;
import org.ogolem.core.ConvergenceException;
import org.ogolem.core.StreamGobbler;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Optimizes a Gaussian basis set to the plane-wave results with respect to basis set spillage.
 * Calls proprietary code developed at Princeton University by the Carter group.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class AdaptiveCustomSpillage extends AbstractAdaptivable {

  private static final long serialVersionUID = (long) 20140503;
  private static final boolean DEBUG = false;

  protected List<String> auxData;
  protected final int noChannels;
  protected final int[] noGauss;
  protected final int[] noContr;
  protected long forParamID = -1;
  protected int count = 0;
  protected final DecimalFormat f = new DecimalFormat("####.########");

  AdaptiveCustomSpillage(final int noChannels, final int[] noGauss, final int[] noContr)
      throws Exception {
    this.noChannels = noChannels;
    assert (noChannels == noGauss.length);
    assert (noChannels == noContr.length);
    this.noGauss = noGauss;
    this.noContr = noContr;
    this.auxData = new ArrayList<>();
  }

  AdaptiveCustomSpillage(final AdaptiveCustomSpillage orig) {
    this.noChannels = orig.noChannels;
    this.noGauss = orig.noGauss.clone();
    this.noContr = orig.noContr.clone();
    if (orig.auxData != null) {
      this.auxData = new ArrayList<>();
      orig.auxData.forEach(
          (s) -> {
            this.auxData.add(s);
          });
    } else {
      this.auxData = null;
    }
  }

  @Override
  public AdaptiveCustomSpillage copy() {
    return new AdaptiveCustomSpillage(this);
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      int geomID,
      final BondInfo bonds) {

    if (forParamID != params.getID()) {
      count = 0;
      forParamID = params.getID();
    }

    final String[] out;
    try {
      out = setupAndRunSpillage(cartes, params, geomID, false);
    } catch (Exception e) {
      if (DEBUG) {
        e.printStackTrace(System.err);
      }
      return FixedValues.NONCONVERGEDENERGY;
    }

    // parse the spillage
    double spillage = FixedValues.NONCONVERGEDENERGY;
    for (final String s : out) {
      if (s.trim().startsWith("spillage")) {
        final String[] sa = s.trim().split("\\s+");
        try {
          spillage = Double.parseDouble(sa[2]);
        } catch (Exception e) {
          if (DEBUG) {
            e.printStackTrace(System.err);
          }
          spillage = FixedValues.NONCONVERGEDENERGY;
        }
        break;
      }
    }

    return spillage;
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      int geomID,
      final BondInfo bonds,
      final double[] grad) {

    if (forParamID != params.getID()) {
      count = 0;
      forParamID = params.getID();
    }

    final String[] out;
    try {
      out = setupAndRunSpillage(cartes, params, geomID, true);
    } catch (Exception e) {
      if (DEBUG) {
        e.printStackTrace(System.err);
      }
      return FixedValues.NONCONVERGEDENERGY;
    }

    // parse the spillage
    double spillage = FixedValues.NONCONVERGEDENERGY;
    for (final String s : out) {
      if (s.trim().startsWith("spillage")) {
        final String[] sa = s.trim().split("\\s+");
        try {
          spillage = Double.parseDouble(sa[2]);
        } catch (Exception e) {
          if (DEBUG) {
            e.printStackTrace(System.err);
          }
          spillage = FixedValues.NONCONVERGEDENERGY;
        }
        break;
      }
    }

    // parse the gradient
    for (int i = 0; i < grad.length; i++) {
      grad[i] = 0.0;
    }
    for (int i = 0; i < out.length; i++) {
      if (out[i].trim().startsWith("gradient")) {
        for (int j = 0; j < grad.length; j++) {
          try {
            grad[j] = Double.parseDouble(out[i + 1 + j]);
          } catch (Exception e) {
            if (DEBUG) {
              e.printStackTrace(System.err);
            }
            grad[j] = 0.0;
          }
        }
        break;
      }
    }

    return spillage;
  }

  @Override
  public double[][] minMaxBordersForParams(final AdaptiveParameters params) {

    final double[][] minMax = new double[2][params.getNumberOfParamters()];

    int c = 0;
    for (int channel = 0; channel < noChannels; channel++) {
      for (int gauss = 0; gauss < noGauss[channel]; gauss++) {
        minMax[0][c] = 0.0; // XXX optimize!!!
        minMax[1][c] = 20.0; //
        c++;
      }

      for (int gauss = 0; gauss < noGauss[channel]; gauss++) {
        for (int contr = 0; contr < noContr[channel]; contr++) {
          minMax[0][c] = -5.0; // XXX optimize!!!
          minMax[1][c] = 5.0; //
          c++;
        }
      }
    }

    return minMax;
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      final ArrayList<CartesianCoordinates> refCartes, final String sMethod) {

    int totalVars = 0;
    for (int channel = 0; channel < noChannels; channel++) {
      int noVars = noGauss[channel]; // number of exponents
      // now the matrix
      noVars += noGauss[channel] * noContr[channel];

      totalVars += noVars;
    }

    // create the stub
    final String[] atoms = {"XX"};
    final int[] paramsPerAt = {totalVars};
    final AdaptiveParameters paramStub =
        new AdaptiveParameters(totalVars, -1, atoms, paramsPerAt, sMethod);

    return paramStub;
  }

  protected String[] setupAndRunSpillage(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      int geomID,
      final boolean doGrad)
      throws Exception {

    // create temporary directory
    count++;
    final String dirName =
        "spillage" + params.getID() + "_" + geomID + "_" + params.hashCode() + "_" + count;
    try {
      OutputPrimitives.createAFolder(dirName);
    } catch (Exception e) {
      // may happen if a previous run crashed
      if (DEBUG) {
        e.printStackTrace(System.err);
      }
    }

    // setup the input.dat file
    if (auxData.size() < geomID + 1) {
      // yet unknown input stub
      readAuxFile(geomID);
    }

    final double[][] xyz = cartes.getAllXYZCoord();
    final String[] inpData = new String[xyz[0].length + 4];
    inpData[0] = auxData.get(geomID);
    for (int i = 0; i < xyz[0].length; i++) {
      inpData[i + 1] =
          xyz[0][i] * Constants.BOHRTOANG
              + "   "
              + xyz[1][i] * Constants.BOHRTOANG
              + "   "
              + xyz[2][i] * Constants.BOHRTOANG;
    }
    inpData[inpData.length - 3] = "." + doGrad + "."; // if we want to have a gradient
    inpData[inpData.length - 2] = ".true."; // we want to operate in binarymode (faster I/O)
    inpData[inpData.length - 1] =
        ".false."; // translator mode is not necessary (we already have binaries to start with)

    // setup the basis.info file
    final ArrayList<String> basisDat = new ArrayList<>(params.getNumberOfParamters());
    final double[] p = params.getAllParamters();
    int c = 0;
    for (int channel = 0; channel < noChannels; channel++) {
      final int nog = noGauss[channel];
      final int noc = noContr[channel];
      basisDat.add("" + nog);
      basisDat.add("" + noc);
      // first the exponents
      for (int g = 0; g < nog; g++) {
        basisDat.add(f.format(p[c]));
        c++;
      }
      // then the matrix
      for (int gauss = 0; gauss < nog; gauss++) {
        String s = "";
        for (int contr = 0; contr < noc; contr++) {
          s += f.format(p[c]) + "  ";
          c++;
        }
        basisDat.add(s);
      }
    }

    // write files
    final String sep = System.getProperty("file.separator");
    OutputPrimitives.writeOut(dirName + sep + "basis.info", basisDat, false);
    OutputPrimitives.writeOut(dirName + sep + "input.dat", inpData, false);

    // link the correct wfn.* files to the directory
    final File folder = new File(System.getProperty("user.dir"));
    final String prefix = "geom" + geomID + "_wfn";
    final File[] listOfFiles = folder.listFiles(new WfnFilter(prefix));
    for (final File wfn : listOfFiles) {
      // strip off the prefix
      final String wfnName = wfn.getName();
      final String realName = wfnName.substring(prefix.length() - 3);
      final String link = dirName + sep + realName;
      OutputPrimitives.createLink(wfnName, link);
    }

    // run the spillage code
    Process proc = null;
    String[] out;
    try {

      String spillCmd = System.getenv("OGO_SPILLAGECMD");
      if (spillCmd == null) {
        spillCmd = "driver";
      }

      Runtime rt = Runtime.getRuntime();
      File dir = new File(dirName);
      proc = rt.exec(spillCmd, null, dir);

      // any error message?
      StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      final int exitValue = proc.waitFor();
      if (DEBUG) {
        System.out.println("DEBUG: Output");
        for (final String s : outputGobbler.getData()) {
          System.out.println("DEBUG: " + s);
        }
        System.out.println("DEBUG: Error");
        for (final String s : errorGobbler.getData()) {
          System.out.println("DEBUG: err" + s);
        }
      }
      if (exitValue != 0) {
        throw new ConvergenceException(
            "WARNING: SPILLAGE returns non-zero return value (local optimization): "
                + exitValue
                + ".");
      } // else: Spillage should(!) have completed normally...

    } catch (Exception e) {
      e.printStackTrace(System.err);
      throw new ConvergenceException("SPILLAGE has a problem (local optimization).", e);
    } finally {
      if (proc != null) {
        proc.destroy();
      }
    }

    out = InputPrimitives.readFileIn(dirName + sep + "out");

    // remove directory
    ManipulationPrimitives.remove(dirName);

    if (DEBUG) {
      for (final String s : out) {
        System.out.println("DEBUG: " + s);
      }
    }

    return out;
  }

  protected static class WfnFilter implements FilenameFilter {

    private final String prefix;

    WfnFilter(final String prefix) {
      this.prefix = prefix;
    }

    @Override
    public boolean accept(File dir, String name) {
      return name.startsWith(prefix);
    }
  }

  protected void readAuxFile(final int geomID) throws Exception {

    final String[] tmp = Input.ReadFile("spillage-aux" + geomID + ".ogo");

    String s = "";
    final String lineEnd = System.getProperty("line.separator");
    for (int i = 0; i < tmp.length; i++) {
      s += tmp[i];
      if (i < tmp.length - 1) {
        s += lineEnd;
      }
    }
    auxData.add(geomID, s);
  }
}
