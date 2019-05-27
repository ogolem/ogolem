/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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
package org.ogolem.math;

/**
 * A collection of fast (and therefore inaccurate!) functions.
 * @author Johannes Dieterich
 * @version 2013-03-27
 */
public class FastFunctions {
    
    // avoid instantiation
    private FastFunctions(){}
    
    /**
     * From: http://martin.ankerl.com/2007/02/11/optimized-exponential-functions-for-java/
     * See: Nical N. Schraudolph "A Fast, Compact Approximation of the exponential function"
     * Has an 11 bit accuracy, is apparently 5 times faster than Math.exp().
     * Not exact enough for most applications (including gradients).
     * @param x
     * @return An approximation to exp(x)
     */
    public static double fastExp(final double x) {
        final long tmp = (long) (1512775 * x + 1072632447);
        return Double.longBitsToDouble(tmp << 32);
    }
    
    public static double fastExp2(final double val) {
        final long tmp = (long) (1512775 * val) + (1072693248 - 60801);
        return Double.longBitsToDouble(tmp << 32);
    }

    public static double fastCorrExp(final double val) {
        System.err.println("you shouldn't use fastCorrExp!");
        if(true) {return Math.exp(val);}
        final long tmp = (long) (1512775 * val) + 1072693248;
        final long mantissa = tmp & 0x000FFFFF;
        long error = mantissa >> 7;   // remove chance of overflow
        error = (error - mantissa * mantissa) / 186; // subtract mantissa^2 * 64
                                                     // 64 / 186 = 1/2.90625
        return Double.longBitsToDouble((tmp - error) << 32);
}
    
    /**
     * Approximated exp(x) through lim_{n\rightarrow\infty}(1+\dfrac{1}{n})^n
     * with n=128. Bernoulli approximation.
     * @param x
     * @return the approximated exp(x) 
     */
    public static double lim128Exp(final double x) {
        double y = 1d + x / 128d;
        y *= y; y *= y; y *= y; y *= y;
        y *= y; y *= y; y *= y;
        return y;
    }
    
    /**
     * Approximated exp(x) through lim_{n\rightarrow\infty}(1+\dfrac{1}{n})^n
     * with n=256. Bernoulli approximation.
     * @param x
     * @return the approximated exp(x) 
     */
    public static double lim256Exp(final double x) {
        double y = 1d + x / 256d;
        y *= y; y *= y; y *= y; y *= y;
        y *= y; y *= y; y *= y; y *= y;
        return y;
    }
    
    /**
     * Approximated exp(x) through lim_{n\rightarrow\infty}(1+\dfrac{1}{n})^n
     * with n=512. Bernoulli approximation.
     * @param x
     * @return the approximated exp(x) 
     */
    public static double lim512Exp(final double x) {
        double y = 1d + x / 512d;
        y *= y; y *= y; y *= y; y *= y;
        y *= y; y *= y; y *= y; y *= y;
        y *= y;
        return y;
    }

    /**
     * Approximated exp(x) through lim_{n\rightarrow\infty}(1+\dfrac{1}{n})^n
     * with n=1024. Bernoulli approximation.
     * @param x
     * @return the approximated exp(x)
     */
    public static double lim1024Exp(final double x) {
        double y = 1d + x / 1024d;
        y *= y; y *= y; y *= y; y *= y;
        y *= y; y *= y; y *= y; y *= y;
        y *= y; y *= y;
        return y;
    }
    
       /**
     * Approximated exp(x) through lim_{n\rightarrow\infty}(1+\dfrac{1}{n})^n
     * with n=2048. Bernoulli approximation.
     * @param x
     * @return the approximated exp(x) 
     */
    public static double lim2048Exp(final double x) {
        double y = 1d + x / 2048d;
        y *= y; y *= y; y *= y; y *= y;
        y *= y; y *= y; y *= y; y *= y;
        y *= y; y *= y; y *= y;
        return y;
    }
    
       /**
     * Approximated exp(x) through lim_{n\rightarrow\infty}(1+\dfrac{1}{n})^n
     * with n=4096. Bernoulli approximation.
     * @param x
     * @return the approximated exp(x)
     */
    public static double lim4096Exp(final double x) {
        double y = 1d + x / 4096d;
        y *= y; y *= y; y *= y; y *= y;
        y *= y; y *= y; y *= y; y *= y;
        y *= y; y *= y; y *= y; y *= y; 
        return y;
    }
    
    /**
     * Returns the square of a double.
     * @param x the value to be squared
     * @return the square
     */
    public static double sq(final double x){
        return x*x;
    }
    
    /**
     * Fast pow() call for integer exponents doing some fancy operations (found
     * on the internet).
     * @param base the base
     * @param exp the exponent
     * @return base to the power of exp
     */
    public static double pow (final double base, final int exp){
        assert(exp >= 0); // the below is not capable to check for oddness for negative numbers according to findbugs.
        return exp==0?1:sq(pow(base,exp/2))*(exp%2==1?base:1);
    }
}
