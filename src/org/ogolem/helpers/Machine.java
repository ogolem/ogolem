/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2011, J. M. Dieterich
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

/**
 * Calculates various machine related features, such as the machine precision.
 * See, e.g.:
 * http://archive.devx.com/java/free/book.asp
 * for a good introduction into this topic.
 * @author Johannes Dieterich
 * @version 2011-01-14
 */
public class Machine {

    /** Radix used by floating-point numbers. */
    private final static int iRadix = computeRadix();
    /** Largest positive value which, when added to 1.0, yields 0 */
    private final static double dMachinePrecision = calcMachinePrecision();
    /** Typical meaningful precision for numerical calculations. */
    private final static double dNumericalPrecision = Math.sqrt(dMachinePrecision);

    public static double calcMachinePrecision() {
        double dFloatingRadix = getRadix();
        double dInverseRadix = 1.0d / dFloatingRadix;
        double dMachPrecision = 1.0d;
        double dTmp = 1.0d + dMachPrecision;
        while (dTmp - 1.0d != 0.0d) {
            dMachPrecision *= dInverseRadix;
            dTmp = 1.0d + dMachPrecision;
        }
        return dMachPrecision;
    }

    public static int getRadix() {
      return iRadix;
    }

    public static double getMachinePrec(){
        return dMachinePrecision;
    }

    public static double getNumericalPrec(){
        return dNumericalPrecision;
    }

    private static int computeRadix() {
        int iCalcRadix = 0;
        double dA = 1.0d;
        double dTmp1, dTmp2;
        do {
            dA += dA;
            dTmp1 = dA + 1.0d;
            dTmp2 = dTmp1 - dA;
        } while (dTmp2 - 1.0d != 0.0d);
        double dB = 1.0d;
        while (iCalcRadix == 0) {
            dB += dB;
            dTmp1 = dA + dA;
            iCalcRadix = (int) (dTmp1 - dA);
        }
        return iCalcRadix;
    }
}

