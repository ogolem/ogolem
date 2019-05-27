/**
Copyright (c) 2010     , J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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
package org.ogolem.helpers;

import java.io.Serializable;
import java.lang.reflect.Method;

/**
 * A 3D Tupel.
 * @author Johannes Dieterich
 * @version 2014-07-10
 */
public class Tuple3D<E,T,V> extends Tuple<E,T> implements Serializable, Cloneable {

    private static final long serialVersionUID = (long) 20110314;

    protected V object3;

    public Tuple3D(E ob1, T ob2, V ob3) {
        super(ob1, ob2);
        object3 = ob3;
    }
    
    @SuppressWarnings("unchecked")
    @Override
    public Tuple3D<E,T,V> clone() throws CloneNotSupportedException {
        
        try {
            final Class<?> cE = object1.getClass();
	    final Method cloneE = cE.getDeclaredMethod("clone");
	    final Object clE = cloneE.invoke(object1, new Object[]{});
            
            final Class<?> cT = object2.getClass();
	    final Method cloneT = cT.getDeclaredMethod("clone");
	    final Object clT = cloneT.invoke(object2, new Object[]{});
            
            final Class<?> cV = object3.getClass();
	    final Method cloneV = cV.getDeclaredMethod("clone");
	    final Object clV = cloneV.invoke(object3, new Object[]{});
            
            return new Tuple3D<>((E) clE, (T) clT, (V) clV);
        } catch (Throwable t){
            throw new CloneNotSupportedException(t.toString());
        }
    }

    public V getObject3(){
        return object3;
    }
    
    public void setObject3(V v){
        object3 = v;
    }
}
