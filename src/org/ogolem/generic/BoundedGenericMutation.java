/**
Copyright (c) 2013, J. M. Dieterich and B. Hartke
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
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * A bounded generic mutation.
 * @author Johannes Dieterich
 * @version 2013-12-04
 */
public class BoundedGenericMutation<T extends Optimizable<Double>> implements GenericMutation<Double,T>{
    
    public static final int SINGLEMUT = 0;
    public static final int MULTIMUT = 1;
    
    private static final long serialVersionUID = (long) 20131201;
    private final Lottery r;
    private final int mode;
    private final double[] upper;
    private final double[] lower;
    
    public BoundedGenericMutation(final int mode, final double[] low, final double[] up){
        this.mode = mode;
        this.upper = up;
        this.lower = low;
        this.r = Lottery.getInstance();
    }
    
    public BoundedGenericMutation(final BoundedGenericMutation<T> orig){
        this.r = Lottery.getInstance();
        this.mode = orig.mode;
        this.upper = orig.upper.clone();
        this.lower = orig.lower.clone();
    }

    @Override
    public BoundedGenericMutation<T> clone() {
        return new BoundedGenericMutation<>(this);
    }

    @Override
    public String getMyID() {
        return "generic bounded mutation\n\tmode " + mode;
    }

    @Override
    public T mutate(final T orig) {
        
        @SuppressWarnings("unchecked")
        final T mut = (T) orig.clone();
        
        final Double[] params = mut.getGenomeCopy();
        final int dims = params.length;
        
        if(mode == SINGLEMUT){
            // find one spot
            final int loc = r.nextInt(dims);
            final double d = r.nextDouble();
            assert(upper[loc] >= lower[loc]);
            params[loc] = d*(upper[loc]-lower[loc])+lower[loc];
        } else if(mode == MULTIMUT){
            // find a list
            final int no = r.nextInt(dims);
            final List<Integer> spots = RandomUtils.listOfPoints(no, dims);
            for(final int loc : spots){
                assert(upper[loc] >= lower[loc]);
                final double d = r.nextDouble();
                params[loc] = d*(upper[loc]-lower[loc])+lower[loc];
            }
        } else{
            throw new RuntimeException("Mode " + mode + " unknown in bounded generic mutation.");
        }
        
        mut.setGenome(params);
        
        return mut;
    }
}
