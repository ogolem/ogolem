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
 * A no crossover, essentially a double copy machine.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class NoCrossover<E, T extends Optimizable<E>> implements GenericCrossover<E,T> {

    private static final long serialVersionUID = (long) 20200429;
    
    @Override
    public GenericCrossover<E, T> clone() {
        return new NoCrossover<>();
    }

    @Override
    public String getMyID() {
        return "NO CROSSOVER";
    }

    @Override
    @SuppressWarnings("unchecked")
    public Tuple<T, T> crossover(final T mother, final T father, final long futureID) {
        
        final T child1 = (T) mother.copy();
        child1.setID(futureID);
        final T child2 = (T) father.copy();
        child2.setID(futureID);
        return new Tuple<>(child1,child2);
    }

    @Override
    public short hasPriority() {
        return -1;
    }
}