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
package org.ogolem.helpers;

import java.util.List;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 * Tests the random utilities.
 * @author Johannes Dieterich
 * @version 2013-04-22
 */
public class RandomUtilsTest {
    
    public RandomUtilsTest(){}
    
    /**
     * Test of gaussDouble method, of class RandomUtils
     */
    @Test
    public void testGaussDouble(){
        System.out.println("gaussDouble");
        for(int it = 0; it < 100; it++){
            final double d = RandomUtils.gaussDouble(-5.0, 5.0, 2.0);
            if(d > 5.0 || d < -5.0){fail("gaussDouble outside valid range I.");}
        }
        
        for(int it = 0; it < 100; it++){
            final double d = RandomUtils.gaussDouble(-5.0, 10.0, 4.0);
            if(d > 10.0 || d < -5.0){fail("gaussDouble outside valid range II.");}
        }
        
        for(int it = 0; it < 100; it++){
            final double d = RandomUtils.gaussDouble(0.0, 5.0, 0.4);
            if(d > 5.0 || d < 0.0){fail("gaussDouble outside valid range III.");}
        }
        
        for(int it = 0; it < 100; it++){
            final double d = RandomUtils.gaussDouble(-5.0, -1.0, 1.4);
            if(d > -1.0 || d < -5.0){fail("gaussDouble outside valid range IV.");}
        }
    }
    
    @Test
    public void testHalfgaussDouble(){
        System.out.println("halfgaussDouble");
        for(int it = 0; it < 100; it++){
            final double d = RandomUtils.halfgaussDouble(-5.0, 5.0, 2.0);
            if(d > 5.0 || d < -5.0){fail("halfgaussDouble outside valid range I.");}
        }
        
        for(int it = 0; it < 100; it++){
            final double d = RandomUtils.halfgaussDouble(-5.0, 10.0, 4.0);
            if(d > 10.0 || d < -5.0){fail("halfgaussDouble outside valid range II.");}
        }
        
        for(int it = 0; it < 100; it++){
            final double d = RandomUtils.halfgaussDouble(0.0, 5.0, 0.4);
            if(d > 5.0 || d < 0.0){fail("halfgaussDouble outside valid range III.");}
        }
        
        for(int it = 0; it < 100; it++){
            final double d = RandomUtils.halfgaussDouble(-5.0, -1.0, 1.4);
            if(d > -1.0 || d < -5.0){fail("halfgaussDouble outside valid range IV.");}
        }
    }
    
    @Test
    public void testListOfPoints(){
        System.out.println("listOfPoints");
        for(int it = 0; it < 100; it++){
            final int low = 0;
            final int high = 10;
            final int no = 3;
            final List<Integer> points = RandomUtils.listOfPoints(no, low, high);
            if(points.size() != no){
                fail("Wrong number of points returned I. " + points.size());
            }
            for(final int i : points){
                if(i < low || i >= high){
                    fail("Wrong point inside test I." + i);
                }
            }
        }
        
        for(int it = 0; it < 100; it++){
            final int low = -10;
            final int high = -2;
            final int no = 4;
            final List<Integer> points = RandomUtils.listOfPoints(no, low, high);
            if(points.size() != no){
                fail("Wrong number of points returned I. " + points.size());
            }
            for(final int i : points){
                if(i < low || i >= high){
                    fail("Wrong point inside test I." + i);
                }
            }
        }
        
        for(int it = 0; it < 100; it++){
            final int low = -5;
            final int high = 10;
            final int no = 7;
            final List<Integer> points = RandomUtils.listOfPoints(no, low, high);
            if(points.size() != no){
                fail("Wrong number of points returned I. " + points.size());
            }
            for(final int i : points){
                if(i < low || i >= high){
                    fail("Wrong point inside test I." + i);
                }
            }
        }
        
        for(int it = 0; it < 100; it++){
            final int low = 5;
            final int high = 10;
            final int no = 1;
            final List<Integer> points = RandomUtils.listOfPoints(no, low, high);
            if(points.size() != no){
                fail("Wrong number of points returned I. " + points.size());
            }
            for(final int i : points){
                if(i < low || i >= high){
                    fail("Wrong point inside test I." + i);
                }
            }
        }
    }
    
    @Test
    public void testVector(){
        System.out.println("testVector");
        for(int it = 1; it < 100; it++){
            final double[] x = new double[it];
            RandomUtils.randomVector(x);
            
            // get the norm
            double d = 0.0;
            for(int i = 0; i < it; i++){
                d += x[i]*x[i];
            }
            final double norm = Math.sqrt(d);
            if(Math.abs(norm-1.0) > 1E-8){
                fail("Norm of vector not one " + norm + " length " + it);
            }
        }
    }
    
    @Test
    public void testEulers(){
        System.out.println("testEulers");
        final double[] x = new double[3];
        for(int it = 1; it < 1000; it++){
            RandomUtils.randomEulers(x);
            
            if(x[0] > Math.PI || x[0] < -Math.PI
                    || x[1] > 0.5*Math.PI || x[1] < -0.5*Math.PI
                    || x[2] > Math.PI || x[2] < -Math.PI){
                fail("Something is wrong with Euler angles " + x[0] + " " + x[1] + " " + x[2]);
            }
        }
    }
}
