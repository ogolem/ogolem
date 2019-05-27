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
package org.ogolem.core;

import org.ogolem.generic.GenericMutation;

/**
 * A Monte Carlo mutation operator for a Molecule.
 * @author Johannes Dieterich
 * @version 2014-03-27
 */
public class MonteCarloMoleculeMutation implements GenericMutation<Double,Molecule>{
    
    public static final int SINGLEMOVEMODE = 0;
    public static final int ALLMOVEMODE = 1;
    public static final int PARTIALMOVEMODE = 2;
    
    private static final long serialVersionUID = (long) 20140327;
    private final int mode;
    private final double maxMove;
    
    public MonteCarloMoleculeMutation(final int mode, final double maxMove){
        assert(maxMove > 0.0);
        assert(mode == 0 || mode == 1 || mode == 2);
        this.maxMove = maxMove;
        this.mode = mode;
    }
    
    public MonteCarloMoleculeMutation(final MonteCarloMoleculeMutation orig){
        this.maxMove = orig.maxMove;
        this.mode = orig.mode;
    }
    
    @Override
    public MonteCarloMoleculeMutation clone() {
        return new MonteCarloMoleculeMutation(this);
    }

    @Override
    public String getMyID() {
        return "Monte Carlo mutation: \n\tmode: " + mode + "\n\tmaxmove: " + maxMove;
    }

    @Override
    public Molecule mutate(final Molecule orig) {
        return MonteCarloMutation.mutate(orig, mode, maxMove);
    }
}
