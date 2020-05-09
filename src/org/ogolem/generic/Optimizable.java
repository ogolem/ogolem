/**
Copyright (c) 2010, J. M. Dieterich and B. Hartke
              2012-2013, J. M. Dieterich
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

import java.io.Serializable;

/**
 * Defined to unify all possible optimizable problems.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public interface Optimizable<T extends Object> extends Serializable, Copyable {

    Optimizable<T> copy();
    
    boolean isDiscreteProblem();
    
    long getID();
    
    void setID(final long id);

    long getFatherID();

    void setFatherID(final long fatherID);
    
    long getMotherID();

    void setMotherID(final long motherID);
    
    double getFitness();
    
    void setFitness(final double fitness);
    
    /**
     * May return NULL!
     * @return Genom as a double or null if this is not implemented.
     */
    double[] getGenomeAsDouble();
    
    /**
     * @return this genom as a T array. Must be returning something valid.
     */
    T[] getGenomeCopy();
    
    void setGenome(final T[] genome);
}
