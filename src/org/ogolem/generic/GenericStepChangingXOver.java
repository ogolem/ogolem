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
package org.ogolem.generic;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.helpers.Tuple;

/**
 * A mutation which changes character based on the step we are dealing with.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class GenericStepChangingXOver<E, T extends Optimizable<E>>
    implements GenericCrossover<E, T> {

  private static final long serialVersionUID = (long) 20140727;
  private final long totalSteps;
  private final List<GenericCrossover<E, T>> xovers;
  private final long[] offsets;
  private GenericCrossover<E, T> ref;

  public GenericStepChangingXOver(
      final long totSteps,
      final List<GenericCrossover<E, T>> xovers,
      final List<Double> percComplete) {
    this.totalSteps = totSteps;
    if (xovers.size() != percComplete.size()) {
      throw new RuntimeException("Mismatch in length of percentages and mutation operators!");
    }

    this.xovers = xovers;
    this.offsets = new long[xovers.size()];
    long old = 0l;
    for (int i = 0; i < xovers.size(); i++) {
      long size = (long) Math.ceil(totalSteps * percComplete.get(i));
      System.out.println(" " + percComplete.get(i) + "\t" + xovers.get(i).getMyID());
      offsets[i] = old;
      old += size;
    }

    if (old < totalSteps - 10) {
      throw new RuntimeException(
          "Too little percents specified in stepwise changing xover. " + old + "\t" + totalSteps);
    }
    if (old > totalSteps + 10) {
      throw new RuntimeException(
          "Too many percents specified in stepwise changing xover." + old + "\t" + totalSteps);
    }
  }

  public GenericStepChangingXOver(final GenericStepChangingXOver<E, T> orig) {
    this.totalSteps = orig.totalSteps;
    this.offsets = orig.offsets.clone();
    this.xovers = new ArrayList<>();
    for (int i = 0; i < orig.xovers.size(); i++) {
      this.xovers.add(orig.xovers.get(i).copy());
    }
  }

  @Override
  public GenericStepChangingXOver<E, T> copy() {
    return new GenericStepChangingXOver<>(this);
  }

  @Override
  public String getMyID() {
    String s = "STEPWISE CHANGING CROSSOVER:\n";
    for (int i = 0; i < xovers.size(); i++) {
      s += "\t\t step complete offset: " + offsets[i] + "\t " + xovers.get(i).getMyID();
    }

    return s;
  }

  @Override
  public Tuple<T, T> crossover(final T mother, final T father, final long futureID) {

    assert (offsets.length == xovers.size());

    for (int i = 0; i < offsets.length; i++) {
      if (offsets[i] <= futureID) {
        final Tuple<T, T> tup = xovers.get(i).crossover(mother, father, futureID);
        ref = xovers.get(i);
        return tup;
      }
    }

    throw new RuntimeException("Apparently doing more steps than we initially anticipated?!");
  }

  @Override
  public short hasPriority() {
    return ref.hasPriority();
  }
}
