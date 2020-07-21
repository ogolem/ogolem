/**
Copyright (c) 2013, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.generic.threading;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.ogolem.generic.Copyable;
import org.ogolem.helpers.Tuple;

/**
 * A generic object cache. Thread-safe (well, that is the purpose of it).
 * @author Johannes Dieterich
 * @version 2020-04-29
 * @param <T> must be an implementation of Copyable
 */
public class ObjectCache<T extends Copyable> {
    
    private static final int MAXTRIES = 1000;
    private final T ref;
    private final List<Tuple<Boolean,T>> cacheContent;
    
    public ObjectCache(final int noOfThreads, final T ref) throws Exception {
        this.cacheContent = Collections.synchronizedList(new ArrayList<>(noOfThreads));//new ArrayList<Tuple<Boolean,T>>(noOfThreads);//Collections.synchronizedList(new ArrayList<Tuple<Boolean,T>>(noOfThreads));
        this.ref = ref;
        for(int i = 0; i < 2*noOfThreads; i++){ // a little bit extra does not hurt. :-)
            cacheContent.add(new Tuple<>(false, dirtyClone(ref)));
        }
    }
    
    public T getOriginalEntry(){
        return ref;
    }
    
    /**
     * Returns a thread private cached object.
     * @return a tuple containing a flag that must be set to FALSE after being done
     */
    public synchronized Tuple<Boolean,T> getUnusedEntry(){
        
        int c = 0;
        while(c < MAXTRIES){
            for(final Tuple<Boolean,T> tup : cacheContent){
                if(tup.getObject1()) {continue;}
                tup.setObject1(true);
                return tup;
            }
            c++;
        }
        
        System.err.println("ERROR: No unused entry in cache.");
        System.err.println("ERROR:  The likeliest cause for this problem is an exception being caused in the parallel region. Please run with a higher DebugLevel to obtain stack traces.");
        throw new RuntimeException("No unused entry in cache.");
    }
    
    @SuppressWarnings("unchecked")
    private T dirtyClone(final T orig) throws Exception {
        
    	final T rtn = (T) orig.copy();

        return rtn;
    }
    
    /*@SuppressWarnings("unchecked")
    private T dirtyDeepCopy(final T orig) throws Exception {
        
        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
        final ObjectOutputStream oos = new ObjectOutputStream(baos);
        oos.writeObject(orig);
        final ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
        final ObjectInputStream ois = new ObjectInputStream(bais);
        final T deepCopy = (T) ois.readObject();
        
        return deepCopy;
    }*/
}
