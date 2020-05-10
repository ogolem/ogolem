/**
Copyright (c) 2013, J. M. Dieterich and B. Hartke
              2013, J. M. Dieterich and B. Hartke
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

import java.util.List;
import org.ogolem.random.RandomUtils;
import org.ogolem.helpers.Tuple;

/**
 * An n-point, generic genotype, real-number based crossover operator.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class GenericPortugalCrossover <E,T extends Optimizable<E>> implements GenericCrossover<E,T>{
    
    private static final long serialVersionUID = (long) 20200429;
    private final int noCrosses;
    
    public GenericPortugalCrossover(final int noCrosses){
        this.noCrosses = noCrosses;
    }
    
    @Override
    public GenericPortugalCrossover<E,T> clone(){
        return new GenericPortugalCrossover<>(noCrosses);
    }
    
    @Override
    public String getMyID(){
        return "portugal " + noCrosses;
    }
    
    @Override
    public Tuple<T,T> crossover(final T mother, final T father, final long futureID){
        
        @SuppressWarnings("unchecked")
        final T child1 = (T) mother.copy();
        @SuppressWarnings("unchecked")
        final T child2 = (T) father.copy();
        final E[] genomeChild1 = child1.getGenomeCopy();
        final E[] genomeChild2 = child2.getGenomeCopy();
        
        if(genomeChild1 == null || genomeChild2 == null){
            return new Tuple<>(null,null);
        }
        
        assert(genomeChild1.length == genomeChild2.length);
        
        if(noCrosses >= genomeChild1.length){
            System.err.println("ERROR: OGOLEM does NOT believe that you want to use this amount of crossing points! Sorry. ;-)");
            return new Tuple<>(null,null);
        }
        
        final List<Integer> crossPts = RandomUtils.rndListOfPoints(noCrosses, genomeChild1.length);
                
        // swap it
        boolean mothFirst = true;
        for(int i = 0; i < genomeChild1.length; i++){
            
            if(!mothFirst){
                final E mo = genomeChild1[i];
                genomeChild1[i] = genomeChild2[i];
                genomeChild2[i] = mo;
            } // else: no swapping needed, mother data is in child1 already

            // check whether we need to swap
            if(crossPts.contains(i)){
                mothFirst = !mothFirst;
            }
        }
        
        // put back in
        child1.setGenome(genomeChild1);
        child2.setGenome(genomeChild2);
        
        return new Tuple<>(child1,child2);
    }    

    @Override
    public short hasPriority() {
        return -1;
    }
}
