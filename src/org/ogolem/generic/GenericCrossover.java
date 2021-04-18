/*
Copyright (c) 2013-2020, J. M. Dieterich and B. Hartke
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

import java.io.Serializable;
import org.ogolem.helpers.Tuple;

/**
 * Interface for a generic crossover.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public interface GenericCrossover<E, T extends Optimizable<E>> extends Serializable, Copyable {

  @Override
  GenericCrossover<E, T> copy();

  String getMyID();

  /**
   * Generic crossover. Returns (in best case) two children T, if something went wrong, one of both
   * may be null.
   *
   * @param mother the T to be used as a mother in the crossover
   * @param father the T to be used as a father in the crossover
   * @param futureID the ID of the children in the future
   * @return a tuple of two children, content may be null
   */
  public Tuple<T, T> crossover(final T mother, final T father, final long futureID);

  /**
   * If this crossover w.r.t. to the last carried out crossover operation has a priority which of
   * the two children to use.
   *
   * @return -1: no priority, 0: child 0, 1: child 1, all other values are not allowed
   */
  public short hasPriority();
}
