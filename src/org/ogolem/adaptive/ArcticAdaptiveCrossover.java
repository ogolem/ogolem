/*
Copyright (c) 2019, M. Dittner, B. Hartke
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * Neither the name of the copyright holder nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

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

import java.util.List;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * An n-point, generic genotype, real-number-based _intermediate_ crossover operator. A multiple
 * mixing portugal variant, i.e. not a _discrete_, but _intermediate_ recombination (known as
 * frequently used variant in evolution strategies). Aimed purpose: after several runs (farming) and
 * adaption of parameter boundaries, this can be used to map down to best mixed solution. Of course,
 * when using giant boundaries and no seeding, usually a mean value of father and mother is not
 * supposed to be a good solution. Frequently used as "hill climbing" crossover step that is able to
 * sample new values instead of new combinations of values.
 *
 * <p>See: T.Weise's book "Global Optimization Algorithms -- Theory and Application", 3rd ed., pp.
 * 337, "(Weighted) Average Crossover, especially item 2. a random weight more at the center between
 * two alleles.
 *
 * @author Mark Dittner
 * @version 2020-12-29
 */
class ArcticAdaptiveCrossover implements GenericCrossover<Double, AdaptiveParameters> {

  private static final long serialVersionUID = (long) 20200425;
  private final int noMixing;
  private final boolean bRandomAvg;
  private final Lottery rd;
  private final int iOrder;

  public ArcticAdaptiveCrossover(final int noMixing, final boolean bRandomAvg, final int order) {
    this.noMixing = noMixing;
    this.bRandomAvg = bRandomAvg;
    this.rd = Lottery.getInstance();
    this.iOrder = order;
    assert (iOrder >= 1);
  }

  @Override
  public ArcticAdaptiveCrossover copy() {
    return new ArcticAdaptiveCrossover(noMixing, bRandomAvg, iOrder);
  }

  @Override
  public String getMyID() {
    return "arcticMixing\n\tn-Point (random: "
        + bRandomAvg
        + ") mixing of AdaptiveParameters, points: "
        + noMixing;
  }

  @Override
  public Tuple<AdaptiveParameters, AdaptiveParameters> crossover(
      AdaptiveParameters mother, AdaptiveParameters father, long futureID) {

    final AdaptiveParameters child1 = mother.copy();
    final AdaptiveParameters child2 = father.copy();
    final Double[] genomeChild1 = child1.getGenomeCopy();
    final Double[] genomeChild2 = child2.getGenomeCopy();

    if (genomeChild1 == null || genomeChild2 == null) {
      return new Tuple<>(null, null);
    }

    assert (genomeChild1.length == genomeChild2.length);

    if (noMixing
        > genomeChild1
            .length) { // not '>=', as one needs not crossings, but array elements for mixing
      System.err.println(
          "ERROR: OGOLEM does NOT believe that you want to use this amount of "
              + "mixing recombination points! Sorry. ;-)");
      // MxD: TODO: Actually such checks should be made during a sanity check before the GA or
      // throw a specialized (to be defined/handled) Exception for communicating this
      // fundamentally wrong input setting: combination of a small problem size with a crossover
      // that is wrongly configured.
      return new Tuple<>(null, null);
    }

    final List<Integer> mixPoints = RandomUtils.listOfPoints(noMixing, genomeChild1.length);

    if (!bRandomAvg) {
      // exactly in between both (arithmetic mean), both children are the same
      for (int locus : mixPoints) {
        final double mix = (genomeChild1[locus] + genomeChild2[locus]) / 2.0;
        genomeChild1[locus] = mix;
        genomeChild2[locus] = mix;
      }
    } else {
      // random weight between 0 and 1 with somewhat perferred centering in the middle
      double dRandom = rd.nextDouble();
      final boolean whichDirection = rd.nextBoolean();
      for (int j = 0; j < iOrder - 1; j++) {
        // decides on how much it is "in between" for exploring volume of crossover hypercube!
        dRandom *= dRandom;
      }
      final double weight = whichDirection ? 0.5 + 0.5 * dRandom : 0.5 - 0.5 * dRandom;
      for (int locus : mixPoints) {
        final double mix1 = (1 - weight) * genomeChild1[locus] + weight * genomeChild2[locus];
        final double mix2 = (1 - weight) * genomeChild2[locus] + weight * genomeChild1[locus];
        genomeChild1[locus] = mix1;
        genomeChild2[locus] = mix2;
      }
    }
    // TODO: possible to define another subversion of this operator with _one_ complete crossover:
    // One weight
    //  for the full sub-genotype, not a new random double for each separate gene

    // put back in
    child1.setGenome(genomeChild1);
    child2.setGenome(genomeChild2);

    return new Tuple<>(child1, child2);
  }

  @Override
  public short hasPriority() {
    return -1;
  }
}
