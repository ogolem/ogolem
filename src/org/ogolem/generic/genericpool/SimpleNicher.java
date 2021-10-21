/*
Copyright (c) 2012-2013, J. M. Dieterich
              2016-2020, J. M. Dieterich and B. Hartke
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

import org.ogolem.generic.Optimizable;
import org.ogolem.helpers.Tuple;

/**
 * An abstract implementation of the nicher.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class SimpleNicher<E, T extends Optimizable<E>> extends AbstractNicher<E, T> {

  private static final long serialVersionUID = (long) 20130403;
  private final int maxIndividualsPerNiche;

  public SimpleNicher(final int maxIndividualsPerNiche) {
    super();
    this.maxIndividualsPerNiche = maxIndividualsPerNiche;
  }

  @Override
  public boolean cleanUp(final GenericPool<E, T> pool) {

    final int poolSize = pool.getCurrentPoolSize();

    // check for most populated niche
    int maxPop = 0;
    Niche maxPopNiche = null;
    for (final Tuple<Niche, Integer> tup : nichePopulation) {
      if (tup.getObject2() > maxPop) {
        maxPop = tup.getObject2();
        maxPopNiche = tup.getObject1().copy();
      }
    }

    if (maxPop <= maxIndividualsPerNiche || maxPopNiche == null) {
      return false;
    }

    // delete the worst individual (and really NEVER delete the best one!)
    boolean deleted = false;
    PoolLoop:
    for (int j = poolSize - 1; j >= 0; j--) {
      final GenericPoolEntry<E, T> entry = pool.getEntryAtPosition(j);
      if (entry.niche().comp(maxPopNiche)) {
        if (j == 0) {
          System.err.println(
              "ERROR: Trying to delete the best individual. Exiting, first some hopefully helpful output...");
          System.err.println(" maxPop " + maxPop);
          System.err.println(" maxIndividualsPerNiche " + maxIndividualsPerNiche);
          System.err.println(" maxPopNiche " + maxPopNiche.getID());
          for (int x = 0; x < pool.getCurrentPoolSize(); x++) {
            final GenericPoolEntry<E, T> en = pool.getEntryAtPosition(x);
            System.err.println(" Pool entry " + x + " " + en.niche().getID());
          }
          System.err.println("");
          nichePopulation.forEach(
              (tup) -> {
                System.err.println(
                    " Niche pop: " + tup.getObject1().getID() + " " + tup.getObject2());
              });
          System.err.println(" Serializing pool for post-mortem analysis...");
          pool.serializeMe();
          System.err.println("");
          System.err.println("This ain't Roy's reaper algorithm. Therefore: Bye. ;-)");
          System.exit(42);
        }
        if (DEBUG) {
          System.out.println(
              "DEBUG: Deleting individual at position " + j + " with ID " + maxPopNiche.getID());
        }
        pool.removeIndividualAtPos(
            j); // Please note that this function DOES delete the niche from the nicher, hence, we
        // do NOT need to update stats
        deleted = true;
        break PoolLoop;
      }
    }

    if (!deleted) {
      System.err.println(
          "ERROR: Although there is a niche w/ population exceeding maxPop, the nicher didn't delete anything. Contact author(s).");
      System.err.println("maxPop " + maxPop);
      System.err.println("maxPopNiche " + maxPopNiche.getID());
      return deleted;
    }

    return true;
  }
}
