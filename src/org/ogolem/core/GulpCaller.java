/*
Copyright (c) 2012-2013, J. M. Dieterich
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

import org.ogolem.io.InputPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Calls Gulp for local optimizations of periodic structures. Needs an input stub file and assumes
 * this file to specify fractional coordinates.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class GulpCaller extends AbstractLocOpt {

  private static final long serialVersionUID = (long) 20130104;
  private final String[] auxData;
  private final double[][] cell;
  private final double[] anker;
  private final String gulpCmd;

  GulpCaller(final GlobalConfig globConf, final String auxFile, final String ank) throws Exception {
    super(globConf);
    final String[] tmp = InputPrimitives.readFileIn(auxFile);
    final String endl = System.getProperty("line.separator");
    String dat = "";
    for (final String s : tmp) {
      dat += s + endl;
    }
    this.auxData = dat.split("OGOCOORD");
    if (!globConf.isPeriodic) {
      throw new Exception("No periodic optimization chosen with Gulp. Not possible.");
    }
    this.cell = globConf.periodicCell.clone();
    this.anker = new double[3];
    final String[] sa = ank.trim().split("\\,");
    for (int i = 0; i < 3; i++) {
      anker[i] = Double.parseDouble(sa[i].trim());
    }
    final String tmp2 = System.getenv("OGO_GULPCMD");
    if (tmp == null) {
      this.gulpCmd = "gulp";
    } else {
      this.gulpCmd = tmp2;
    }
  }

  GulpCaller(final GulpCaller orig) {
    super(orig);
    this.auxData = orig.auxData.clone();
    this.cell = orig.cell.clone();
    this.anker = orig.anker.clone();
    this.gulpCmd = orig.gulpCmd;
  }

  @Override
  public GulpCaller copy() {
    return new GulpCaller(this);
  }

  @Override
  public String myID() {
    return "Gulp";
  }

  @Override
  public String myIDandMethod() {
    return "Gulp (custom input)";
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      long lID,
      CartesianCoordinates cartes,
      final boolean[][] baConstraints,
      final boolean isConstricted,
      final BondInfo bonds)
      throws InitIOException, ConvergenceException, CastException, Exception {

    // translate the coordinates to their periodic equivalent
    final double[][] per = PeriodicCoordTranslator.coordsToPeriodic(cartes, cell, anker, true);
    assert (per[0].length == auxData.length / 3 + 2);

    // merge the input file
    String inp = auxData[0];
    int count = 1;
    for (int at = 0; at < per[0].length; at++) {
      for (int c = 0; c < 3; c++) {
        inp += per[c][at] + auxData[count];
        count++;
      }
    }

    // write input there
    final String inputFile = "gulp" + lID + ".dat";
    try {
      OutputPrimitives.writeOut(inputFile, inp, false);
    } catch (Exception e) {
      throw new InitIOException("Exception when writing gulp input for id " + lID, e);
    }

    // run gulp
    final String[] output = runGulp(gulpCmd, inputFile);

    // parse output
    double e = Double.MAX_VALUE;
    for (int i = 0; i < output.length; i++) {
      final String line = output[i].trim();
      if (line.startsWith("Final energy =")) {
        final String[] sa = line.split("\\s+");
        e = Double.parseDouble(sa[3]) * Constants.EVTOHARTREE;
      } else if (line.startsWith("Final fractional coordinates of atoms")) {
        // the actual coordinates start six lines later
        for (int at = 0; at < per[0].length; at++) {
          final String[] atline = output[at + i + 6].trim().split("\\s+");
          per[0][at] = Double.parseDouble(atline[3]);
          per[1][at] = Double.parseDouble(atline[4]);
          per[2][at] = Double.parseDouble(atline[5]);
        }
      }
    }

    // delete input file
    ManipulationPrimitives.remove(inputFile);

    // translate periodic coordinates to cartesians
    final double[][] xyzOpt = PeriodicCoordTranslator.coordsFromPeriodic(per, cell);

    // setup cartes object
    final CartesianCoordinates cartOpt = new CartesianCoordinates(cartes);
    cartOpt.setAllXYZ(xyzOpt);
    cartOpt.setEnergy(e);

    return cartOpt;
  }

  private static String[] runGulp(final String gulpCmd, final String input) throws Exception {

    final Runtime rt = Runtime.getRuntime();
    final Process proc = rt.exec(gulpCmd + " " + input);

    // any error message?
    final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

    // any output?
    final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

    // kick them off
    errorGobbler.start();
    outputGobbler.start();

    // any error???
    if (proc.waitFor() != 0) {
      throw new Exception("Gulp returns non-zero return value (local optimization).");
    } else {
      // gulp should(!) have completed normally...
      return outputGobbler.getData();
    }
  }
}
