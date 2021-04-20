/*
Copyright (c) 2010-2015, J. M. Dieterich
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
package org.ogolem.freqs;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import org.ogolem.corrfunc.CorrelationData;
import org.ogolem.helpers.Fortune;
import org.ogolem.io.InputPrimitives;

/**
 * Entry point to calculate a power spectrum from the FFT transform of an velocity autocorrelation
 * function. Please note that this is *NOT* working in atomic units as the rest of ogolem but
 * instead directly works in wavenumbers. Implemented using Tinker's basic algorithm.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class MainPowerSpec {

  /**
   * Prefactor to convert cm-1 for the Fourier transform. Last multiplier is the speed of light in
   * cm per picosecond.
   */
  private static final double PREFAC = Math.PI * 2 * 2.99792458E-2;

  /**
   * Entry point to calculate the power spectrum.
   *
   * @param args First argument: the velocity autocorrelation file. Second: the step length in fs
   *     (!). Last (optional): maximum wavelength for analysis in wavenumbers (default: 5000).
   */
  public static void run(final String[] args) {

    if (args[0].equalsIgnoreCase("help")) {
      System.out.println(
          "This calculates power spectra from correlation functions. It needs as parameters:");
      System.out.println(" * the correlation file");
      System.out.println(" * the step length in fs");
      System.out.println("Optionally, you can specify:");
      System.out.println(
          " * the maximum wavelength for the analysis in wavenumbers (default: 5000)");
      System.out.println(" * an output file to use. Default: System.out");
      System.exit(0);
    }

    final String corrFile = args[0];
    final double stepLength = Double.parseDouble(args[1]) / 1000.0;
    final int maxWave = (args.length >= 3) ? Integer.parseInt(args[2]) : 5000;
    PrintStream ps = System.out;
    if (args.length >= 4) {
      final String outFile = args[3];
      final File f = new File(outFile);
      try {
        final boolean fileNotExists = f.createNewFile();
        if (!fileNotExists) {
          System.out.println("INFO: File " + outFile + "already exists.");
        }
        ps = new PrintStream(f, Charset.forName("UTF-8"));
      } catch (Exception e) {
        System.err.println("Failure to open PrintStream for output file " + outFile);
        e.printStackTrace(System.err);
      }
    }

    // let's read in the previously prepared data
    CorrelationData corrData = null;
    try {
      final String[] data = InputPrimitives.readFileIn(corrFile);
      corrData = new CorrelationData(data);
    } catch (Exception e) {
      System.err.println(
          "ERROR: Couldn't read and/or parse velocity autocorrelation file " + corrFile);
      e.printStackTrace(System.err);
      System.exit(57);
    }

    System.err.println(
        "FIXME: Also (if the correlation function should not be unitless) assumes A/ps or something for the norm/average of the correlation function");

    // Fourier transforming it for every individual wavelength
    final double[] intensities = new double[maxWave];
    for (int wave = 1; wave < maxWave; wave++) {
      final double freq = PREFAC * wave;
      double d = 0.0;
      for (final CorrelationData.CorrelationDataPoint pt : corrData) {
        final double time = stepLength * pt.k;
        d += pt.norm * Math.cos(freq * time);
      }
      intensities[wave] = 1000 * stepLength * d;
    }

    // done, print results and say goodbye
    final DecimalFormat form = new DecimalFormat("##########.######");
    ps.println("# -----------------------------------------------");
    ps.println("# Power spectrum from the velocity autocorrelation");
    ps.println("# Source file " + corrFile);
    ps.println("# Frequency [cm-1]            Intensity");
    ps.println("# -----------------------------------------------");
    for (int wave = 1; wave < intensities.length; wave++) {
      ps.println(" " + wave + " \t " + form.format(intensities[wave]));
    }
    ps.println("# " + Fortune.randomFortune());
  }
}
