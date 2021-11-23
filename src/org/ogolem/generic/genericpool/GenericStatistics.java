/*
Copyright (c) 2010-2014, J. M. Dieterich
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
package org.ogolem.generic.genericpool;

import java.io.File;
import java.io.RandomAccessFile;
import java.io.Serializable;
import java.nio.charset.Charset;
import java.text.SimpleDateFormat;
import java.util.Calendar;

/**
 * This is a really simple stub for a statistics object.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public final class GenericStatistics implements Serializable {

  private static final long serialVersionUID = (long) 20150803;

  private static final String DATE_FORMAT_NOW = "yyyy-MM-dd HH:mm:ss";

  private final String logFile;
  private final long noBestCountsToFlush;

  private boolean beSilent = false;
  private long lastBestID = -1000l;
  private double lastBestFitness = Double.MAX_VALUE;
  private long noBestCount = 0l; // counter since the last time a new best one was registered
  private boolean lastLineWasStatus = false;
  private long lastLinePointer = 0l;

  /**
   * Constructs the statistics object.
   *
   * @param log The path to the log file.
   * @param countsToFlush how many counts that are NOT new best individuals we wait before printing
   *     a line.
   */
  public GenericStatistics(final String log, final long countsToFlush) {
    this.logFile = log;
    this.noBestCountsToFlush = countsToFlush;

    // touch the log file
    final File f = new File(log);
    if (!f.exists()) {
      try {
        final boolean success = f.createNewFile();
        if (!success) {
          throw new RuntimeException("File reports no success with file creation.");
        }
      } catch (Exception e) {
        throw new RuntimeException("Logfile " + log + " could not be created!");
      }
    }
  }

  /**
   * Constructs the statistics object. Option to not touch the file.
   *
   * @param log The path to the log file.
   * @param countsToFlush how many counts that are NOT new best individuals we wait before printing
   *     a line.
   * @param touchFile whether or not to touch the log file as trial
   */
  public GenericStatistics(final String log, final long countsToFlush, final boolean touchFile) {
    this.logFile = log;
    this.noBestCountsToFlush = countsToFlush;

    if (touchFile) {
      // touch the log file
      final File f = new File(log);
      if (!f.exists()) {
        try {
          final boolean success = f.createNewFile();
          if (!success) {
            throw new RuntimeException("File reports no success with file creation.");
          }
        } catch (Exception e) {
          throw new RuntimeException("Logfile " + log + " could not be created!");
        }
      }
    }
  }

  /**
   * If the individual is the new best individual, a line is printed to the log file with the
   * current date and time.
   *
   * @param individual the ID of the individual
   * @param position the position in the pool
   * @param fitness the fitness of this individual
   */
  public void individualAddedToPool(
      final long individual, final int position, final double fitness) {

    if (beSilent) {
      return;
    }

    noBestCount++;

    if (position == 0) {
      // print a a line to file
      final Calendar cal = Calendar.getInstance();
      final SimpleDateFormat sdf = new SimpleDateFormat(DATE_FORMAT_NOW);
      final String date = sdf.format(cal.getTime());
      final String line =
          "["
              + date
              + "]"
              + "   new best individual "
              + individual
              + "\t fitness "
              + fitness
              + System.lineSeparator();

      appendToFile(line, true);
      lastBestID = individual;
      lastBestFitness = fitness;
      noBestCount = 0;
      lastLineWasStatus = false;
    }

    if (noBestCount >= noBestCountsToFlush) {
      // print a different line
      final String line = assembleNoNewLine(individual);

      appendToFile(line, false);
      noBestCount = 0;
      lastLineWasStatus = true;
    }
  }

  /**
   * Registers that an attempt was made to add an individual to the pool but it never made it for
   * whatever reason.
   *
   * @param individual the ID of this (the last) individual
   */
  public void registerIndividualNotAdded(final long individual) {

    if (beSilent) {
      return;
    }

    noBestCount++;

    if (noBestCount >= noBestCountsToFlush) {
      // print a different line
      final String line = assembleNoNewLine(individual);

      appendToFile(line, false);
      noBestCount = 0;
      lastLineWasStatus = true;
    }
  }

  private String assembleNoNewLine(final long individual) {

    final Calendar cal = Calendar.getInstance();
    final SimpleDateFormat sdf = new SimpleDateFormat(DATE_FORMAT_NOW);
    final String date = sdf.format(cal.getTime());
    final String line =
        "["
            + date
            + "]"
            + "   best "
            + lastBestID
            + " current "
            + individual
            + "      fitness "
            + lastBestFitness
            + System.lineSeparator();

    return line;
  }

  private void appendToFile(final String line, final boolean isNew) {

    try (final RandomAccessFile f = new RandomAccessFile(logFile, "rw")) {
      if (lastLineWasStatus) {
        // remove the last line
        f.setLength(lastLinePointer);
      }

      // append
      f.seek(lastLinePointer);
      f.write(line.getBytes(Charset.forName("UTF-8")));

      if (isNew) {
        lastLinePointer = f.getFilePointer();
      }
    } catch (Exception e) {
      System.err.println("WARNING: Couldn't write to log file. " + e.toString());
    }
  }

  void silence() {
    this.beSilent = true;
  }
}
