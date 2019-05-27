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
package org.ogolem.math;

import static org.junit.Assert.*;
import org.junit.Test;

/**
 *
 * @author Johannes Dieterich
 * @version 2013-03-27
 */
public class FastFunctionsTest {
    
    public static final double LOWACC = 1E-3;
    public static final double MEDACC = 5E-4;
    public static final double MEDPLUSACC = 1E-4;
    public static final double HIGHACC = 1E-8;
    
    public FastFunctionsTest() {
    }
    
    /**
     * Test of fastExp method, of class FastFunctions.
     */
    @Test
    public void testFastExp() {
        System.out.println("fastExp");
        double x = -20.0;
        while(x < 10.0){
            final double expResult = Math.exp(x);
            final double result = FastFunctions.fastExp(x);
            assertEquals(expResult, result, Math.max(MEDACC,Math.abs(0.05*expResult))); // maximum of 5E-4 or 2%
            x += 0.05;
        }
    }

    /**
     * Test of fastExp2 method, of class FastFunctions.
     */
    @Test
    public void testFastExp2() {
        System.out.println("fastExp2");
        double x = -20.0;
        while(x < 10.0){
            final double expResult = Math.exp(x);
            final double result = FastFunctions.fastExp2(x);
            assertEquals(expResult, result, Math.max(MEDACC,Math.abs(0.05*expResult))); // maximum of 5E-4 or 2%
            x += 0.05;
        }
    }

    /**
     * Test of fastCorrExp method, of class FastFunctions.
     */
    // currently, fastCorrExp() is broken and returns Math.exp()
    /*@Test
    public void testFastCorrExp() {
        System.out.println("fastCorrExp");
        double x = 0.0;
        while(x < 10.0){
            final double expResult = Math.exp(x);
            final double result = FastFunctions.fastCorrExp(x);
            assertEquals(expResult, result, Math.max(MEDPLUSACC,Math.abs(0.01*expResult))); // maximum of 1E-4 or 1%
            System.out.println("fine: " + x + " " + expResult);
            x += 0.05;
        }
    }

    /**
     * Test of lim128Exp method, of class FastFunctions.
     */
    @Test
    public void testLim128Exp() {
        System.out.println("lim128Exp");
        double x = -10.0;
        while(x < 5.0){
            final double expResult = Math.exp(x);
            final double result = FastFunctions.lim128Exp(x);
            assertEquals(expResult, result, Math.max(LOWACC,Math.abs(0.1*expResult))); // maximum of 1E-3 or 10%
            x += 0.05;
        }
    }

    /**
     * Test of lim256Exp method, of class FastFunctions.
     */
    @Test
    public void testLim256Exp() {
        System.out.println("lim256Exp");
        double x = -20.0;
        while(x < 7.0){
            final double expResult = Math.exp(x);
            final double result = FastFunctions.lim256Exp(x);
            assertEquals(expResult, result, Math.max(LOWACC,Math.abs(0.1*expResult))); // maximum of 1E-3 or 10%
            x += 0.05;
        }
    }

    /**
     * Test of lim512Exp method, of class FastFunctions.
     */
    @Test
    public void testLim512Exp() {
        System.out.println("lim512Exp");
        double x = -20.0;
        while(x < 10.0){
            final double expResult = Math.exp(x);
            final double result = FastFunctions.lim512Exp(x);
            assertEquals(expResult, result, Math.max(LOWACC,Math.abs(0.1*expResult))); // maximum of 1E-3 or 10%
            x += 0.05;
        }
    }

    /**
     * Test of lim1024Exp method, of class FastFunctions.
     */
    @Test
    public void testLim1024Exp() {
        System.out.println("lim1024Exp");
        double x = -20.0;
        while(x < 10.0){
            final double expResult = Math.exp(x);
            final double result = FastFunctions.lim1024Exp(x);
            assertEquals(expResult, result, Math.max(MEDACC,Math.abs(0.05*expResult))); // maximum of 5E-4 or 5%
            x += 0.05;
        }
    }

    /**
     * Test of lim2048Exp method, of class FastFunctions.
     */
    @Test
    public void testLim2048Exp() {
        System.out.println("lim2048Exp");
        double x = -20.0;
        while(x < 10.0){
            final double expResult = Math.exp(x);
            final double result = FastFunctions.lim2048Exp(x);
            assertEquals(expResult, result, Math.max(MEDACC,Math.abs(0.03*expResult))); // maximum of 5E-4 or 3%
            x += 0.05;
        }
    }

    /**
     * Test of lim4096Exp method, of class FastFunctions.
     */
    @Test
    public void testLim4096Exp() {
        System.out.println("lim4096Exp");
        double x = -20.0;
        while(x < 10.0){
            final double expResult = Math.exp(x);
            final double result = FastFunctions.lim4096Exp(x);
            assertEquals(expResult, result, Math.max(MEDPLUSACC,Math.abs(0.02*expResult))); // maximum of 1E-4 or 2%
            x += 0.05;
        }
    }

    /**
     * Test of sq method, of class FastFunctions.
     */
    @Test
    public void testSq() {
        System.out.println("sq");
        assertEquals(0.0, FastFunctions.sq(0.0), HIGHACC);
        assertEquals(4.0, FastFunctions.sq(2.0), HIGHACC);
        assertEquals(9.0, FastFunctions.sq(3.0), HIGHACC);
        assertEquals(2.5*2.5, FastFunctions.sq(2.5), HIGHACC);
        assertEquals(1.0, FastFunctions.sq(-1.0), HIGHACC);
        assertEquals(4.0, FastFunctions.sq(-2.0), HIGHACC);
        assertEquals(0.01, FastFunctions.sq(0.1), HIGHACC);
    }

    /**
     * Test of pow method, of class FastFunctions.
     */
    @Test
    public void testPow() {
        System.out.println("pow");
        assertEquals(1.0, FastFunctions.pow(0.0, 0), HIGHACC);
        assertEquals(1.0, FastFunctions.pow(1.0, 0), HIGHACC);
        assertEquals(1.0, FastFunctions.pow(1.0, 1), HIGHACC);
        assertEquals(4.0, FastFunctions.pow(2.0, 2), HIGHACC);
        assertEquals(27.0, FastFunctions.pow(3.0, 3), HIGHACC);
    }
}
