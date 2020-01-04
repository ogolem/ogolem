/**
Copyright (c) 2017, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.ogolem.adaptive.genericfitness.GenericReferencePoint;
import org.ogolem.adaptive.genericfitness.ReferenceInputData;
import org.ogolem.properties.Property;

/**
 * A library of reference points. Indexing included. The implementation is fugly.
 * The only way that one could make it more crazy is by introducing a librarian
 * that is a monk... errr APE.
 * @author Johannes Dieterich
 * @version 2017-12-17
 */
@SuppressWarnings({"unchecked", "rawtypes"})
class ReferencePointLibrary implements Serializable {
    
    private static final long serialVersionUID = (long) 20171217;
    private static final boolean DEBUG = false;
    
    private final Map<String,Integer> nameToListPlace;
    private final List<List<GenericReferencePoint>> allReferencePoints;
    
    ReferencePointLibrary(){
        this.nameToListPlace = new HashMap<>();
        this.allReferencePoints = new ArrayList<>();
    }
    
    <T extends Property, V extends ReferenceInputData<T>> void addReferencePoint(final GenericReferencePoint<T,V> point, final String name){
        final Integer listPlace = nameToListPlace.get(name);
        if(listPlace == null){
            // not yet known
            final List<GenericReferencePoint> list = new ArrayList<>();
            list.add(point);
            allReferencePoints.add(list);
            nameToListPlace.put(name, allReferencePoints.size()-1);
            if(DEBUG){System.out.println("DEBUG: Adding " + name + " as new in position " + (allReferencePoints.size()-1));}
        } else {
            // known
            final List<GenericReferencePoint> listForProp = allReferencePoints.get(listPlace);
            listForProp.add(point);
            if(DEBUG){System.out.println("DEBUG: Adding " + name + " to existing in position " + listPlace);}
        }
    }
    
    <T extends Property, V extends ReferenceInputData<T>> List<GenericReferencePoint<T,V>> retrieveReferencePointsForName(final String name) throws Error {
        
        final Integer listPlace = nameToListPlace.get(name);
        if(listPlace == null){
            final List<GenericReferencePoint<T,V>> listForPropGeneric = new ArrayList<>();
            if(DEBUG){System.out.println("DEBUG: Returning empty for " + name);}
            return listForPropGeneric;
        } else {
            final List listForProp = allReferencePoints.get(listPlace);
            final List<GenericReferencePoint<T,V>> listForPropGeneric = (List<GenericReferencePoint<T,V>>) listForProp; // we are kinda certain this works
            if(DEBUG){System.out.println("DEBUG: Returning " + name + " position " + listPlace);}
            return listForPropGeneric;
        }
    }
    
    List<List<GenericReferencePoint>> retrieveAllReferencePoints(){
        return allReferencePoints;
    }
    
    Set<String> getAllKeysAdded(){
        return nameToListPlace.keySet();
    }
}
