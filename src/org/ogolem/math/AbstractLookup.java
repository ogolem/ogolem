/*
Copyright (c) 2012-2014, J. M. Dieterich
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
package org.ogolem.math;

import java.io.Serializable;
import org.ogolem.generic.Copyable;

/**
 * An abstract lookup table for functions.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public abstract class AbstractLookup implements Serializable, Copyable {

  private static final long serialVersionUID = (long) 20120612;

  public static final int[] POSSIBLE_SIZES = {
    1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152
  };

  protected final int ENTRIES;
  protected final int ENTRIESDECR;
  protected double[] table;
  protected final double st;
  protected final double en;
  protected final double dis;
  protected final double disEntries;

  public AbstractLookup(final int entries, final double start, final double end) {
    ENTRIES = entries;
    ENTRIESDECR = ENTRIES - 1;
    table = new double[ENTRIES];
    st = start;
    en = end;
    dis = end - start;
    disEntries = ENTRIESDECR / dis;
    final double incr = dis / (ENTRIESDECR);
    double x = start;
    for (int i = 0; i < ENTRIES; i++) {
      table[i] = func(x);
      x += incr;
    }
  }

  public AbstractLookup(
      final int entries, final double start, final double end, final boolean beLazy) {
    ENTRIES = entries;
    ENTRIESDECR = ENTRIES - 1;
    st = start;
    en = end;
    dis = end - start;
    disEntries = ENTRIESDECR / dis;
    if (!beLazy) {
      table = new double[ENTRIES];
      final double incr = dis / (ENTRIESDECR);
      double x = start;
      for (int i = 0; i < ENTRIES; i++) {
        table[i] = func(x);
        x += incr;
      }
    }
  }

  public AbstractLookup(final AbstractLookup orig) {
    this.ENTRIES = orig.ENTRIES;
    this.ENTRIESDECR = orig.ENTRIESDECR;
    this.st = orig.st;
    this.en = orig.en;
    this.dis = orig.dis;
    this.disEntries = orig.disEntries;
    this.table = orig.table.clone();
  }

  @Override
  public abstract AbstractLookup copy();

  /**
   * Always gives the canonical function result back.
   *
   * @param x Must be in the interval [end,start]
   * @return the canonical func(x).
   */
  public final double canonical(final double x) {
    return func(x);
  }

  /**
   * A lookup function for the exp including linear interpolation.
   *
   * @param x Must be in the interval [end,start]
   * @return The looked up and interpolated func(x). Canonical func(x) if x is outside the interval.
   */
  public double funcInter(final double x) {

    if (table == null) {
      // lazy init, do now
      table = new double[ENTRIES];
      final double incr = dis / (ENTRIESDECR);
      double xxx = st;
      for (int i = 0; i < ENTRIES; i++) {
        table[i] = func(xxx);
        xxx += incr;
      }
    }

    if (x >= en || x < st) return func(x);

    final double poi = (x - st) * disEntries;
    final int point = (int) poi;
    final double rest = poi - point;
    // jump in
    final double uncorr = table[point];
    final double uncorr2 = (point == ENTRIESDECR) ? uncorr : table[point + 1];
    return (uncorr + rest * (uncorr2 - uncorr));
  }

  /**
   * A lookup function for the exp without interpolation.
   *
   * @param x Must be in the interval [start,end]
   * @return The looked up func(x). Canonical func(x) if x is outside the interval.
   */
  protected double funcNonInter(final double x) {

    if (table == null) {
      // lazy init, do now
      table = new double[ENTRIES];
      final double incr = dis / (ENTRIESDECR);
      double xxx = st;
      for (int i = 0; i < ENTRIES; i++) {
        table[i] = func(xxx);
        xxx += incr;
      }
    }

    if (x >= en || x < st) return func(x);

    final double poi = (x - st) * disEntries;
    final int point = (int) poi;
    // jump in
    return table[point];
  }

  /**
   * All implementations must override this.
   *
   * @param x key
   * @return value
   */
  protected abstract double func(final double x);
}
