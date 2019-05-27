/**
Copyright (c) 2013, J. M. Dieterich
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

/**
 * A simple generic darwin implementation.
 * @author Johannes Dieterich
 * @version 2013-11-30
 */
public class SimpleGenericDarwin<E,T extends Optimizable<E>> extends GenericAbstractDarwin<E,T>{
    
    private static final long serialVersionUID = (long) 20131130;
    
    public SimpleGenericDarwin(final GenericCrossover<E,T> cross, final GenericMutation<E,T> mut,
            final GenericSanityCheck<E,T> sanity, final GenericFitnessFunction<E,T> fitness,
            final IndividualWriter<T> writer, final double crossPoss, final double mutPoss,
            final boolean printBeforeFitness, final int noOfTries){
        super(cross, mut, sanity, fitness, writer, crossPoss, mutPoss,
            printBeforeFitness, noOfTries);
    }
    
    public SimpleGenericDarwin(final SimpleGenericDarwin<E,T> orig){
        super(orig);
    }
    
    @Override
    public SimpleGenericDarwin<E,T> clone(){
        return new SimpleGenericDarwin<>(this);
    }
    
    @Override
    public String getMyID(){
        return "simple darwin using: "
                + "\n\txover    " + xover.getMyID()
                + "\n\tmutation " + mutation.getMyID()
                + "\n\tfitness  " + fitness.getMyID();
    }
    
    @Override
    public void postXOver(final T individual1, final T individual2, final long futureID){}
    
    @Override
    public void postMutation(final T individual){}
    
    @Override
    public void runAfterEachTry(){}
}
