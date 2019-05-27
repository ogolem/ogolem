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
package org.ogolem.fft;

import org.apache.commons.math3.complex.Complex;

/**
 * Pads data for butterfly FFTs. 
 * @author Johannes Dieterich
 * @version 2013-12-23
 */
public class Padder {
    
    // disallow instantiation
    private Padder(){}
    
    /**
     * Finds the power of two padding for a given array length.
     * @param length the original length
     * @return the padded length, a power of two
     */
    public static int findPowTwoPad(final int length){
        
        assert(length > 0);
        final int origLength = length;
        int padLength = 2;
        while(padLength < origLength){
            padLength *= 2;
        }
        
        return padLength;
    }
    
    /**
     * Pads a given input array to a power of two array of at least original size (padded with zeros)
     * @param orig the original data as a double array
     * @return the zero-padded array of at least original size
     */
    public static double[] powTwoPad(final double[] orig){
        
        assert(orig != null);
        assert(orig.length > 0);
        
        final int origLength = orig.length;
        final int padLength = findPowTwoPad(origLength);
        
        // we now know the padded length, allocate a new array and copy stuff
        final double[] padded = new double[padLength];
        System.arraycopy(orig, 0, padded, 0, origLength);
        
        return padded;
    }
    
    /**
     * Pads a APACHE COMMONS Complex type as an array to a power of two.
     * @param orig the original complex array.
     * @return a padded array containing copies and zero'd objects.
     */
    public static Complex[] powTwoPad(final Complex[] orig){
        
        assert(orig != null);
        assert(orig.length > 0);
        
        final int origLength = orig.length;
        final int padLength = findPowTwoPad(origLength);
        
        // we now know the padded length, allocate a new array and copy stuff
        final Complex[] padded = new Complex[padLength];
        for(int i = 0; i < origLength; i++){
            padded[i] = new Complex(orig[i].getReal(), orig[i].getImaginary());
        }
        
        for(int i = origLength; i < padLength; i++){
            padded[i] = new Complex(0.0); // pad with zeros
        }
        
        return padded;
    }
}
