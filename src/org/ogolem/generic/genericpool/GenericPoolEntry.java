/*
Copyright (c) 2012, J. M. Dieterich
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
package org.ogolem.generic.genericpool;

import java.io.Serializable;
import org.ogolem.generic.Copyable;
import org.ogolem.generic.Optimizable;

/**
 * An individual for the generic pool
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class GenericPoolEntry<E, T extends Optimizable<E>> implements Serializable, Copyable {

  private static final long serialVersionUID = (long) 20120215;

  private double fitness = Double.MAX_VALUE;
  private T individual = null;
  private Niche niche = null;

  GenericPoolEntry() {}

  GenericPoolEntry(final T individual, final double fitness) {
    this.fitness = fitness;
    this.individual = individual;
    this.niche = null;
  }

  GenericPoolEntry(final T individual, final double fitness, final Niche niche) {
    this.fitness = fitness;
    this.individual = individual;
    this.niche = niche;
  }

  @SuppressWarnings("unchecked")
  GenericPoolEntry(final GenericPoolEntry<E, T> orig) {
    this.fitness = orig.fitness;
    this.individual = (orig.individual == null) ? null : (T) orig.individual.copy();
    this.niche = (orig.niche == null) ? null : orig.niche.copy();
  }

  @Override
  public GenericPoolEntry<E, T> copy() {
    return new GenericPoolEntry<>(this);
  }

  public void setFitness(final double fit) {
    fitness = fit;
  }

  public double getFitness() {
    return fitness;
  }

  public void setIndividual(final T individuum) {
    individual = individuum;
  }

  public T getIndividual() {
    return individual;
  }

  public void setNiche(final Niche n) {
    niche = n;
  }

  public Niche getNiche() {
    return niche;
  }
}
