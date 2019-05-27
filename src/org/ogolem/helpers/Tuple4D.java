/**
Copyright (c) 2015, J. M. Dieterich and B. Hartke
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
 * A 4D tuple.
 * @author Johannes Dieterich
 * @version 2015-07-25
 */
public class Tuple4D<E,T,V,W> extends Tuple3D<E,T,V> implements Serializable, Cloneable {
    
    private static final long serialVersionUID = (long) 20140727;

    private W object4;

    public Tuple4D(final E ob1, final T ob2, final V ob3, final W ob4) {
        super(ob1, ob2, ob3);
        object4 = ob4;
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
                        
            final Class<?> cW = object4.getClass();
	    final Method cloneW = cW.getDeclaredMethod("clone");
	    final Object clW = cloneW.invoke(object4, new Object[]{});
            
            return new Tuple4D<>((E) clE, (T) clT, (V) clV, (W) clW);
        } catch (Throwable t){
            throw new CloneNotSupportedException(t.toString());
        }
    }

    public W getObject4(){
        return object4;
    }
    
    public void setObject4(W w){
        object4 = w;
    }
}
