/**
Copyright (c) 2010-2011, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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

import java.io.Serializable;

/**
 * A generic equipped class for computation results.
 * @author Johannes Dieterich
 * @version 2015-08-02
 */
public final class Result<T> implements Serializable{

    private static final long serialVersionUID = (long) 20110329;

    private final T result;
    private final boolean wasOK;
    private final int clientID;
    private final long[] family;

    Result(final T res, final boolean ok, final int clientID, final long[] family){
        assert(family != null);
        this.result = res;
        this.wasOK = ok;
        this.clientID = clientID;
        this.family = family;
    }

    Result(final T res, final boolean ok, final int clientID){
        this.result = res;
        this.wasOK = ok;
        this.clientID = clientID;
        this.family = new long[]{-1,-1,-1};
    }

    T getResult(){
        return result;
    }

    boolean wasOK(){
        return wasOK;
    }

    int getClientID(){
        return clientID;
    }

    long[] getFamily(){
        assert(family != null);
        return family;
    }
}
