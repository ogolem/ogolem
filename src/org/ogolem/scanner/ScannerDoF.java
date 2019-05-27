/**
Copyright (c) 2013, J. M. Dieterich
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
package org.ogolem.scanner;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * A degree of freedom. Essentially, a pointer to one entry in the a z-Matrix
 * with an associated set of bounds and an increment.
 * @author Johannes Dieterich
 * @version 2013-09-22
 */
class ScannerDoF implements Iterable<Double>{
   
    private static final boolean DEBUG = false;
    private final int molecule;
    private final int atom;
    private final int dof;
    private final double start;
    private final double end;
    private final double incr;
    
    ScannerDoF(final int mol, final int at, final int dof, final double start, final double end,
            final double incr){
        this.molecule = mol;
        this.atom = at;
        this.dof = dof;
        this.start = start;
        this.end = end;
        this.incr = incr;
        assert(molecule >= 0);
        assert(atom >= 0);
        assert(end >= start);
        assert(incr > 0.0);
        assert(dof >= 0 && dof <= 2);
    }
    
    @Override
    public Iterator<Double> iterator(){
        return new ScannerDoFIterator();
    }

    public int getAtom() {
        return atom;
    }

    public int getMolecule() {
        return molecule;
    }

    public int getDoF() {
        return dof;
    }
    
    class ScannerDoFIterator implements Iterator<Double> {
        
        private double curr;
        ScannerDoFIterator(){
            curr = start;
        }
        
        @Override
        public void remove(){
            throw new UnsupportedOperationException("remove() operation not available in ScannerDoF iterator.");
        }
        
        @Override
        public boolean hasNext(){
            if(DEBUG){System.out.println("DEBUG: For molecule " + molecule + " atom "
                    + atom + " and dof " + dof + " we have more? " + (curr <= end)
                    + " curr " + curr + " incr " + incr + " end " + end);}
            return (curr <= end);
        }
        
        @Override
        public Double next(){
            if(curr > end){throw new NoSuchElementException("No more element smaller than end");}
            if(DEBUG){System.out.println("DEBUG: About to return " + curr + " for mol " + molecule + " atom " + atom + " dof " + dof);}
            final double old = curr;
            curr += incr;
            
            return old;
        }
    }
}
