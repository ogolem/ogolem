/*
Copyright (c) 2010-2014, J. M. Dieterich
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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.io.OutputPrimitives;

/**
 * Value object holding the frequencies.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public final class Frequencies implements Serializable {

  private static final long serialVersionUID = (long) 20130411;

  /**
   * Taken from http://users.mccammon.ucsd.edu/~blu/Research-Handbook/physical-constant.html,
   * confirmed from http://web.utk.edu/~rcompton/constants, conversion factor to be multiplied with
   * Eh to yield wavenumbers. Also wiki features this value now. 219474.6313705 This on the other
   * hand, after long fiddling with BXH is the CORRECT conversion factor for this purpose. The a.u.
   * time and the speed of light in cm instead m. Sounds funny, is funny.
   */
  private static final double AUFREQSTOWAVENUMBERS = 1 / 2.41888428E-17 / 2.99792458E10;

  private final List<Frequency> frequencies;
  private final boolean anyImag;

  public Frequencies(
      double[] imagFreqs, double[] realFreqs, double[] realIntensities, double[][] vecs) {

    this.frequencies = new ArrayList<>();
    assert (imagFreqs.length == realFreqs.length);
    boolean aImag = false;
    for (int i = 0; i < realFreqs.length; i++) {

      final double inten = (realIntensities == null) ? 1.0 : realIntensities[i];
      final double[] eigenvec = new double[vecs.length];
      for (int j = 0; j < vecs.length; j++) {
        eigenvec[j] = vecs[j][i];
      }
      final boolean b = (imagFreqs[i] != 0.0);
      if (b) {
        aImag = true;
      }
      final Frequency freq = new Frequency(realFreqs[i], inten, eigenvec, b);
      frequencies.add(freq);
    }
    this.anyImag = aImag;

    Collections.sort(frequencies);
  }

  public List<Frequency> getAllFrequncies() {
    return frequencies;
  }

  public boolean anyImaginaryFrequency() {
    return anyImag;
  }

  public String[] printableFrequencies() {

    final String[] printout = new String[frequencies.size() + 2];
    printout[0] = "# Frequencies";
    printout[1] = "# counter\t wavelength [a.u.]\t wavelength [cm-1]\t imaginary?\t intensity";

    final String freqFormat = "%8.2f";
    final String freqFormatAU = "%17.15f";
    final String countFormat = "%9d";
    final String intsFormat = "%8.2f";

    int counter = 0;
    for (int i = 0; i < frequencies.size(); i++) {

      final Frequency freq = frequencies.get(i);

      printout[counter + 2] =
          String.format(Locale.US, countFormat, counter)
              + "\t"
              + String.format(Locale.US, freqFormatAU, freq.freq)
              + "\t"
              + String.format(Locale.US, freqFormat, freq.freq * AUFREQSTOWAVENUMBERS);
      if (freq.isImaginary) {
        printout[counter + 2] += "\t\t Y \t\t";
      } else {
        printout[counter + 2] += "\t\t N \t\t";
      }
      printout[counter + 2] += String.format(Locale.US, intsFormat, 1.0);

      counter++;
    }

    return printout;
  }

  /**
   * Prints out all frequencies. Really untested at the moment!!!
   *
   * @param filePrefix
   * @param cartes
   */
  void printAllFrequencies(final String filePrefix, final CartesianCoordinates cartes) {

    final String freqFormat = "%9.8f";
    final double convFac = 1.0;
    final int noAts = cartes.getNoOfAtoms();

    for (int i = 0; i < frequencies.size(); i++) {

      final Frequency freq = frequencies.get(i);
      final String file = filePrefix + "_" + i + ".xyz";
      final String[] print = cartes.createPrintableCartesians();

      // add the frequency we are talking about
      print[1] += "\t Frequency: " + freq.freq * AUFREQSTOWAVENUMBERS;

      // PLEASE NOTE: eigen vectors are FORTRAN-style ordered
      // add the vector information
      for (int j = 0; j < noAts; j++) {
        print[j + 2] +=
            "\t"
                + String.format(Locale.US, freqFormat, convFac * freq.eigenvec[j])
                + "\t"
                + String.format(Locale.US, freqFormat, convFac * freq.eigenvec[j + noAts])
                + "\t"
                + String.format(Locale.US, freqFormat, convFac * freq.eigenvec[j + 2 * noAts]);
      }

      try {
        OutputPrimitives.writeOut(file, print, true);
      } catch (Exception e) {
        e.printStackTrace(System.err);
      }
    }
  }

  public static class Frequency implements Comparable<Frequency>, Serializable {

    private static final long serialVersionUID = (long) 20141121;

    private final double intensity;
    private final double freq;
    private final double[] eigenvec;
    private final boolean isImaginary;

    Frequency(
        final double freq, final double inten, final double[] eigenvec, final boolean isImaginary) {
      this.intensity = inten;
      this.freq = freq;
      this.isImaginary = isImaginary;
      this.eigenvec = eigenvec;
    }

    @Override
    public int compareTo(final Frequency f) {

      final double d1 = this.getFreq();
      final double d2 = f.getFreq();

      return Double.compare(d1, d2);
    }

    /** @return the intensity */
    public double getIntensity() {
      return intensity;
    }

    /** @return the frequency */
    public double getFreq() {
      return freq;
    }

    /** @return the eigen vector of this frequency */
    public double[] getEigenvec() {
      return eigenvec;
    }

    /** @return wether it is imaginary */
    public boolean isIsImaginary() {
      return isImaginary;
    }
  }
}
