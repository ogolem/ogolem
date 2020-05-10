/**
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

import org.ogolem.helpers.Tuple;

/**
 * Adaptor to translate turn a mutation into a crossover.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class GenericMutationAsXOver<E, T extends Optimizable<E>> implements GenericCrossover<E,T> {

    private static final long serialVersionUID = (long) 20200429;
    private final GenericMutation<E,T> mutation;
    
    public GenericMutationAsXOver(final GenericMutation<E,T> mut){
        this.mutation = mut;
    }
    
    public GenericMutationAsXOver(final GenericMutationAsXOver<E,T> orig){
        this.mutation = orig.mutation.clone();
    }
    
    @Override
    public GenericMutationAsXOver<E,T> clone() {
        return new GenericMutationAsXOver<>(this);
    }

    @Override
    public String getMyID() {
        return "MUTATION AS CROSSOVER ADAPTOR\n\tmutation: " + mutation.getMyID();
    }

    @SuppressWarnings("unchecked")
    @Override
    public Tuple<T, T> crossover(final T mother, final T father, final long futureID) {
        
        final T work1 = (T) mother.copy();
        work1.setID(futureID);
        final T work2 = (T) father.copy();
        work2.setID(futureID);
        
        final T child1 = mutation.mutate(work1);
        final T child2 = mutation.mutate(work2);
        
        return new Tuple<>(child1,child2);
    }

    @Override
    public short hasPriority() {
        // odds are: no
        return -1;
    }
}
