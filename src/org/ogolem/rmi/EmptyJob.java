/**
Copyright (c) 2010-2014, J. M. Dieterich
              2015-2016, J. M. Dieterich and B. Hartke
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
package org.ogolem.rmi;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.rmi.RMICodes.JOBSTATE;

/**
 * An empty job as a fallback.
 * @author Johannes Dieterich
 * @version 2016-01-25
 */
final class EmptyJob implements Job<Boolean>{
    
    private static final long serialVersionUID = (long) 20141101;

    @Override
    public Task<Boolean> nextTask(){
        return null;
    }

    @Override
    public RMICodes.JOBSTATE submitResult(Result<Boolean> result){
        return JOBSTATE.FINISH;
    }

    @Override
    public boolean jobWaiting(){
        return false;
    }

    @Override
    public boolean jobFinished(){
        return true;
    }    

    @Override
    public List<Task<Boolean>> nextInitTasks(final int maxTasks, final int noProxies) {
        return new ArrayList<>();
    }

    @Override
    public List<Boolean> mergePools(final int noOfAssocResults, final List<Boolean> clientPool, final int maxStructsBack, final long lastStart) {
        return clientPool;
    }

    @Override
    public Tuple3D<Boolean, Long, Integer> nextGlobOptChunk(final int maxTasks) {
        return new Tuple3D<>(false,0l,0);
    }
}