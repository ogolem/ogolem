/**
Copyright (c) 2014, J. M. Dieterich
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

import java.util.ArrayList;
import java.util.List;

/**
 * A mutation which changes character based on the step we are dealing with.
 * @author Johannes Dieterich
 * @version 2014-07-27
 */
public class GenericStepChangingMut<E, T extends Optimizable<E>> implements GenericMutation<E,T>{

    private static final long serialVersionUID = (long) 20140727;
    private final long totalSteps;
    private final List<GenericMutation<E,T>> mutations;
    private final long[] offsets;
    
    public GenericStepChangingMut(final long totSteps, final List<GenericMutation<E,T>> muts,
            final List<Double> percComplete){
        this.totalSteps = totSteps;
        if(muts.size() != percComplete.size()){throw new RuntimeException("Mismatch in length of percentages and mutation operators!");}
        
        this.mutations = muts;
        this.offsets = new long[muts.size()];
        long old = 0l;
        for(int i = 0; i < muts.size(); i++){
            long size = (long) Math.ceil(totalSteps*percComplete.get(i));
            offsets[i] = old;
            old += size;
        }
        
        if(old < totalSteps-10){throw new RuntimeException("Too little percents specified in stepwise changing mutation. " + old + "\t" + totalSteps);}
        if(old > totalSteps+10){throw new RuntimeException("Too many percents specified in stepwise changing mutation." + old + "\t" + totalSteps);}
    }
    
    public GenericStepChangingMut(final GenericStepChangingMut<E,T> orig){
        this.totalSteps = orig.totalSteps;
        this.offsets = orig.offsets.clone();
        this.mutations = new ArrayList<>();
        for(int i = 0; i < orig.mutations.size(); i++){
            this.mutations.add(orig.mutations.get(i).clone());
        }
    }
    
    @Override
    public GenericStepChangingMut<E, T> clone() {
        return new GenericStepChangingMut<>(this);
    }

    @Override
    public String getMyID() {
        String s = "STEPWISE CHANGING MUTATION:\n";
        for(int i = 0; i < mutations.size(); i++){
            s += "\t\t step complete offset: " + offsets[i] + "\t "+ mutations.get(i).getMyID();
        }
        
        return s;
    }

    @Override
    public T mutate(final T orig) {
   
        assert(offsets.length == mutations.size());
        
        final long lID = orig.getID();
        for(int i = 0; i < offsets.length; i++){
            if(offsets[i] <= lID){
                return mutations.get(i).mutate(orig);
            }
        }
        
        throw new RuntimeException("Apparently doing more steps than we initially anticipated?!");
    }
}
