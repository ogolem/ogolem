/**
Copyright (c) 2014, J. M. Dieterich
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
package org.ogolem.core;

import java.util.Random;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Johannes Dieterich
 * @version 2014-12-30
 */
public class CoordTranslationTest {
    
    private static final double NUMACC = 1e-10;
    private final Random r;
    
    public CoordTranslationTest(){
        this.r = new Random();
    }
    
    /**
     * Test of sanitizePhi method, of class CoordTranslation.
     */
    @Test
    public void testSanitizePhi() {
        System.out.println("sanitizePhi");
        
        final int TRIALS = 1000;
        final double lower = -Math.PI;
        final double upper = Math.PI;
        final double period = 2*Math.PI;
        final int maxPer = 10;
        
        for(int i = 0; i < TRIALS; i++){
            final double ang1 = lower + r.nextDouble()*period;
            final double san1 = CoordTranslation.sanitizePhi(ang1);
            assertEquals(ang1,san1,NUMACC);
            
            // now one period on top
            final double ang2 = ang1+period;
            final double san2 = CoordTranslation.sanitizePhi(ang2);
            assertEquals(ang1,san2,NUMACC);
            
            // remove one period
            final double ang3 = ang1-period;
            final double san3 = CoordTranslation.sanitizePhi(ang3);
            assertEquals(ang1,san3,NUMACC);
            
            // add rnd periods
            final int pers1 = r.nextInt(maxPer);
            final double ang4 = ang1+pers1*period;
            final double san4 = CoordTranslation.sanitizePhi(ang4);
            assertEquals(ang1,san4,NUMACC);
            
            
            // remove rnd periods
            final int pers2 = r.nextInt(maxPer);
            final double ang5 = ang1-pers2*period;
            final double san5 = CoordTranslation.sanitizePhi(ang5);
            assertEquals(ang1,san5,NUMACC);
            
        }        
    }

    /**
     * Test of sanitizeOmega method, of class CoordTranslation.
     */
    @Test
    public void testSanitizeOmega() {
        System.out.println("sanitizeOmega");
        final int TRIALS = 1000;
        final double lower = -0.5*Math.PI;
        final double upper = 0.5*Math.PI;
        final double period = Math.PI;
        final int maxPer = 10;
        
        for(int i = 0; i < TRIALS; i++){
            final double ang1 = lower + r.nextDouble()*period;
            final double san1 = CoordTranslation.sanitizeOmega(ang1);
            assertEquals(ang1,san1,NUMACC);
            
            // now one period on top
            final double ang2 = ang1+period;
            final double san2 = CoordTranslation.sanitizeOmega(ang2);
            assertEquals(ang1,san2,NUMACC);
            
            // remove one period
            final double ang3 = ang1-period;
            final double san3 = CoordTranslation.sanitizeOmega(ang3);
            assertEquals(ang1,san3,NUMACC);
            
            // add rnd periods
            final int pers1 = r.nextInt(maxPer);
            final double ang4 = ang1+pers1*period;
            final double san4 = CoordTranslation.sanitizeOmega(ang4);
            assertEquals(ang1,san4,NUMACC);
            
            
            // remove rnd periods
            final int pers2 = r.nextInt(maxPer);
            final double ang5 = ang1-pers2*period;
            final double san5 = CoordTranslation.sanitizeOmega(ang5);
            assertEquals(ang1,san5,NUMACC);
            
        }        
    }

    /**
     * Test of sanitizePsi method, of class CoordTranslation.
     */
    @Test
    public void testSanitizePsi() {
        System.out.println("sanitizePsi");
        
        final int TRIALS = 1000;
        final double lower = -Math.PI;
        final double upper = Math.PI;
        final double period = 2*Math.PI;
        final int maxPer = 10;
        
        for(int i = 0; i < TRIALS; i++){
            final double ang1 = lower + r.nextDouble()*period;
            final double san1 = CoordTranslation.sanitizePsi(ang1);
            assertEquals(ang1,san1,NUMACC);
            
            // now one period on top
            final double ang2 = ang1+period;
            final double san2 = CoordTranslation.sanitizePsi(ang2);
            assertEquals(ang1,san2,NUMACC);
            
            // remove one period
            final double ang3 = ang1-period;
            final double san3 = CoordTranslation.sanitizePsi(ang3);
            assertEquals(ang1,san3,NUMACC);
            
            // add rnd periods
            final int pers1 = r.nextInt(maxPer);
            final double ang4 = ang1+pers1*period;
            final double san4 = CoordTranslation.sanitizePsi(ang4);
            assertEquals(ang1,san4,NUMACC);
            
            
            // remove rnd periods
            final int pers2 = r.nextInt(maxPer);
            final double ang5 = ang1-pers2*period;
            final double san5 = CoordTranslation.sanitizePsi(ang5);
            assertEquals(ang1,san5,NUMACC);
            
        }
    }
}
