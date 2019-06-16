/**
Copyright (c) 2012-2013, J. M. Dieterich
              2016, J. M. Dieterich and B. Hartke
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

import java.util.ArrayList;
import java.util.List;
import org.ogolem.generic.Optimizable;
import org.ogolem.helpers.Tuple;

/**
 * An abstract implementation of the nicher.
 * @author Johannes Dieterich
 * @version 2016-01-14
 */
abstract class AbstractNicher<E,T extends Optimizable<E>>  implements Nicher<E,T> {

    private static final long serialVersionUID = (long) 20130403;
    protected static final boolean DEBUG = false;
    protected final List<Tuple<Niche,Integer>> nichePopulation;
    
    public AbstractNicher(){
        this.nichePopulation = new ArrayList<>();
    }
    
    @Override
    public void report(final Niche added){
        
        for(final Tuple<Niche,Integer> tup : nichePopulation){
            if(tup.getObject1().comp(added)){
                final int curr = tup.getObject2();
                tup.setObject2(curr + 1);
                return;
            }
        }
        
        // apparently unknown so far
        final Tuple<Niche, Integer> tup = new Tuple<>(added.clone(),1);
        nichePopulation.add(tup);
    }
    
    @Override
    public void delete(final Niche deleted){
        
        for(int i = 0; i < nichePopulation.size(); i++){
            final Tuple<Niche,Integer> tup = nichePopulation.get(i);
            if(tup.getObject1().comp(deleted)){
                final int curr = tup.getObject2();
                if(curr == 1){
                    // after the decrement this niche would vanish
                    nichePopulation.remove(i);
                    return;
                } else {
                    tup.setObject2(curr - 1);
                    return;
                }
            }
        }
        
        // apparently unknown so far, which shouldn't happen
        System.err.println("ERROR: Trying to delete a niche unknown to the nicher. This is a bug, notify the author(s).");
        System.err.println("Niche: " + deleted.getID());
        nichePopulation.forEach((tup) -> {
            System.err.println(" Niche pop: " + tup.getObject1().getID() + " " + tup.getObject2());
        });
        System.exit(84);
    }
        
    @Override
    public List<Tuple<Niche,Integer>> calculateNichePopulation(final boolean print, final GenericPool<E,T> pool){

        assert(pool != null);
        assert(pool.getNicheOfIndividualAtPos(0) != null);
        
        final List<Tuple<Niche,Integer>> niches = new ArrayList<>();
        for(final Tuple<Niche,Integer> tup : nichePopulation){
            final Tuple<Niche,Integer> copy = new Tuple<>(tup.getObject1().clone(),tup.getObject2());
            niches.add(copy);
        }
        
        if(DEBUG){
            final List<Niche> alNiches = new ArrayList<>();
            final List<Integer> alPops = new ArrayList<>();
        
            alNiches.add(pool.getNicheOfIndividualAtPos(0).clone());
            alPops.add(1);

            PopLoop:
            for(int i = 1; i < pool.getCurrentPoolSize(); i++){
                final Niche niche = pool.getNicheOfIndividualAtPos(i);
                for(int j = 0; j < alNiches.size(); j++){
                    final Niche compNiche = alNiches.get(j);
                    if(compNiche.comp(niche)){
                        int iPop = alPops.get(j);
                        iPop++;
                        alPops.set(j, iPop);
                        continue PopLoop;
                    }

                    if(j == alNiches.size() -1){
                        // new niche
                        alNiches.add(niche.clone());
                        alPops.add(1);
                        continue PopLoop;
                    }
                }
            }
            
            // compare niche population with this fresh thing
            TupLoop: for(final Tuple<Niche,Integer> tup : niches){
                final Niche n = tup.getObject1();
                for(int i = 0; i < alNiches.size(); i++){
                    if(n.comp(alNiches.get(i))){
                        // let's see...
                        final int pop = tup.getObject2();
                        final int freshPop = alPops.get(i);
                        if(pop != freshPop){
                            System.err.println("ERROR: Something wrong in niche housekeeping. For niche " + n.getID() + " fresh number " + freshPop + " vs. " + pop + ".");
                            nichePopulation.forEach((tup2) -> {
                                System.err.println(" Niche pop: " + tup2.getObject1().getID() + " " + tup2.getObject2());
                            });
                            System.exit(126);
                        }
                        continue TupLoop;
                    }
                    if(i == alNiches.size()-1){
                        System.err.println("ERROR: Niche " + n.getID() + " vanished from niches!");
                        System.exit(168);
                    }
                }
            }
        }

        if(print){
            System.out.println("*********************************************");
            for(int i = 0; i < niches.size(); i++){
                System.out.println("Niche " + i + " has a population of " + niches.get(i).getObject2() + ". Identifier: " + niches.get(i).getObject1().getID());
            }
            System.out.println("*********************************************");
        }

        return niches;
    }
}
