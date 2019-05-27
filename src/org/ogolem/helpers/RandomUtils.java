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

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.ogolem.core.CoordTranslation;

/**
 * Some helper functions for random numbers.
 * @author Johannes Dieterich
 * @version 2013-08-02
 */
public class RandomUtils {
    
    public static final long TRIESTOEMERGENCY=10000l;
    private static final boolean DEBUG = false;
    private static final Random r = new Random();
    
    // disallow instantiation
    private RandomUtils(){}
        
    /**
     * Generates a Gaussian-distributed random number in between bounds.
     * @param low lower bound for the number
     * @param high upper bound for the number
     * @return a random number
     */
    public static double gaussDouble(final double low, final double high){
        return gaussDouble(low,high,1.0);
    }
    
    /**
     * Generates a Gaussian-distributed random number in between bounds with a
     * defined standard deviation.
     * @param low lower bound for the number
     * @param high upper bound for the number
     * @param stdDev the standard deviation of the Gaussian distribution
     * @return a random number
     */
    public static double gaussDouble(final double low, final double high, final double stdDev){
        
        assert(high != low);
        assert(high > low);
        assert(stdDev > 0.0);
        
        final double mid = (high-low)/2+low;
        final double std = Math.abs(stdDev*(high-low)/2);
        
        long c = 0;
        do {
            final double d = r.nextGaussian()*std+mid;
            c++;
            if(d >= low && d <= high){return d;}
        } while(c < TRIESTOEMERGENCY);
                
        System.err.println("WARNING: Emergency occured in gaussDouble(). Returning " + mid + ".");
        return mid;
    }
    
    /**
     * Generates a Gaussian-distributed random number in between bounds.
     * @param low lower bound for the number
     * @param high upper bound for the number
     * @return a random number
     */
    public static double halfgaussDouble(final double low, final double high){
        return halfgaussDouble(low,high,1.0);
    }
    
    /**
     * Generates a Gaussian-distributed random number in between bounds.
     * @param low lower bound for the number
     * @param high upper bound for the number
     * @param stdDev the standard deviation of the Gaussian distribution
     * @return a random number
     */
    public static double halfgaussDouble(final double low, final double high, final double stdDev){
        
        assert(high != low);
        assert(high > low);
        assert(stdDev > 0.0);
        
        final double std = Math.abs(stdDev*(high-low));
        
        long c = 0;
        do {
            final double d = low+Math.abs(r.nextGaussian())*std;
            c++;
            if(d < high){return d;}
        } while(c < TRIESTOEMERGENCY);
        
        System.err.println("WARNING: Emergency occured in halfgaussDouble(). Returning " + low + ".");
        return low;
    }
    
    public static List<Integer> rndListOfPoints(final int noPoints, final int high){
        return rndListOfPoints(noPoints,0,high);
    }
    
    public static List<Integer> rndListOfPoints(final int noPoints, final int low, final int high){
        
        final List<Integer> list = listOfPoints(noPoints,low,high);
        // shake it
        for(int i = 0; i < noPoints; i++){
            final int mv = r.nextInt(noPoints);
            final int curr = list.remove(i);
            list.add(mv, curr);
        }
        
        return list;
    }
    
    /**
     * Returns a list of random points in between 0 (inclusive) and high (exclusive).
     * @param noPoints number of random points to draw
     * @param high the upper border (exclusive)
     * @return a list of points (in ascending order)
     */
    public static List<Integer> listOfPoints(final int noPoints, final int high){
        final List<Integer> list = new ArrayList<>(noPoints);
        listOfPoints(noPoints, 0, high,list);
        
        return list;
    }
    
    /**
     * Returns a list of random points in between low (inclusive) and high (exclusive).
     * @param noPoints number of random points to draw
     * @param low the lower border (inclusive)
     * @param high the upper border (exclusive)
     * @return a list of points (in ascending order)
     */
    public static List<Integer> listOfPoints(final int noPoints, final int low, final int high){
        final List<Integer> list = new ArrayList<>(noPoints);
        listOfPoints(noPoints, low, high,list);
        
        return list;
    }
    
    /**
     * Returns a list of random points in between low (inclusive) and high (exclusive).
     * @param noPoints number of random points to draw
     * @param low the lower border (inclusive)
     * @param high the upper border (exclusive)
     * @param points a list object, supposed to eliminate the overhead of object allocation
     */
    public static void listOfPoints(final int noPoints, final int low, final int high, final List<Integer> points){
        
        assert(noPoints > 0);
        assert(high > low);
        final int diff = high-low;
        assert(diff >= noPoints);
        
        points.clear();
        for(int i = low; i < high; i++){
            points.add(i);
        }
        
        // remove some
        while(points.size() > noPoints){
            final int which = r.nextInt(points.size());
            points.remove(which);
        }
        
        assert(points.size() == noPoints);
    }
    
    /**
     * Returns a shuffled list of integers from [0 to (exclusive) no]. The sole
     * idea of this routine is to be fast and somewhat random, we DO NOT claim
     * these numbers to be in any way properly randomized.
     * @param no the number of integers to shuffle
     * @param list Our list. Changed on exit.
     */
    public static void randomList(final int no, final List<Integer> list){
        
        assert(list != null);
        assert(no > 0);
       
        list.clear(); // just to make sure
        list.add(0);
        
        for(int i = 1; i < no; i++){
            final int pos = r.nextInt(list.size()+1);
            list.add(pos,i);
        }
    }
    
    /**
     * Returns a vector with a norm of 1.0.
     * @param x Our vector. Changed on exit.
     */
    public static void randomVector(final double[] x){
        randomVector(x, 1.0);
    }
    
    /**
     * Returns a vector with a defined norm.
     * @param x Our vector. Changed on exit.
     * @param norm The norm.
     */
    public static void randomVector(final double[] x, final double norm){
        
        assert(x != null);
        assert(x.length > 0);
        
        double d = 0.0;
        for(int i = 0; i < x.length; i++){
            final double rd = (r.nextBoolean()) ? - r.nextDouble() : r.nextDouble();
            x[i] = rd;
            d += rd*rd;
        }
        
        final double n = norm/Math.sqrt(d);
        for(int i = 0; i < x.length; i++){
            x[i] *= n;
        }
        
        assert(Math.abs(norm - norm(x)) < 1E-8);
    }
    
    private static double norm(final double[] x){
        
        double norm = 0.0;
        for(int i = 0; i < x.length; i++){
            norm += x[i]*x[i];
        }
        norm = Math.sqrt(norm);
        
        return norm;
    }
    
    /**
     * Returns a set of random Euler angles in order phi/omega/psi.
     * @param x The array in which the Euler angles will be stored on exit
     */
    public static void randomEulers(final double[] x){
        
        assert(x.length == 3);
        
        // phi
        x[0] = (r.nextBoolean()) ? Math.PI*r.nextDouble() : -Math.PI*r.nextDouble();
        // omega
        x[1] = (r.nextBoolean()) ? 0.5*Math.PI*r.nextDouble(): -0.5*Math.PI*r.nextDouble();
        // psi
        x[2] = (r.nextBoolean()) ? Math.PI*r.nextDouble() : -Math.PI*r.nextDouble();
    }
    
    /**
     * Increments a set of Euler angles randomly.
     * @param eulers The starting euler angles, changed on exit.
     * @param strength The strength of the increment as a fraction of the total interval of the eulers
     */
    public static void randomEulerIncrements(final double[] eulers, final double strength){
        
        assert(eulers.length == 3);
        assert(strength >= 0.0);
        assert(strength <= 1.0);
        
        final double incrPhi = strength*r.nextDouble()*2*Math.PI;
        final double incrOmega = strength*r.nextDouble()*Math.PI;
        final double incrPsi = strength*r.nextDouble()*2*Math.PI;
        
        // now add them on top, sanitizing if needed
        eulers[0] = (r.nextBoolean()) ? eulers[0] + incrPhi : eulers[0] - incrPhi;
        eulers[1] = (r.nextBoolean()) ? eulers[1] + incrOmega : eulers[1] - incrOmega;
        eulers[2] = (r.nextBoolean()) ? eulers[2] + incrPsi : eulers[2] - incrPsi;
        
        /*while(!(eulers[0] <= Math.PI && eulers[0] >= -Math.PI)){
            if(eulers[0] > Math.PI){assert(eulers[0]-Math.PI > 0.0); eulers[0] = -Math.PI+(eulers[0]-Math.PI);}
            else if(eulers[0] < -Math.PI){assert(eulers[0]-Math.PI < 0.0); eulers[0] = Math.PI+(eulers[0]-Math.PI);}
        }

        while(!(eulers[1] <= 0.5*Math.PI && eulers[1] >= -0.5*Math.PI)){
            if(eulers[1] > 0.5*Math.PI){assert(eulers[1]-0.5*Math.PI > 0.0); eulers[1] = -0.5*Math.PI+(eulers[1]-0.5*Math.PI);}
            else if(eulers[1] < -0.5*Math.PI){assert(eulers[1]-0.5*Math.PI < 0.0); eulers[1] = 0.5*Math.PI+(eulers[1]-0.5*Math.PI);}
        }
        
        while(!(eulers[2] <= Math.PI && eulers[2] >= -Math.PI)){
            if(eulers[2] > Math.PI){assert(eulers[2]-Math.PI > 0.0); eulers[2] = -Math.PI+(eulers[2]-Math.PI);}
            else if(eulers[2] < -Math.PI){assert(eulers[2]-Math.PI < 0.0); eulers[2] = Math.PI+(eulers[2]-Math.PI);}
        }*/
        
        eulers[0] = CoordTranslation.sanitizePhi(eulers[0]);
        eulers[1] = CoordTranslation.sanitizeOmega(eulers[1]);
        eulers[2] = CoordTranslation.sanitizePsi(eulers[2]);

        assert(eulers[0] <= Math.PI && eulers[0] >= -Math.PI);
        assert(eulers[1] <= 0.5*Math.PI && eulers[1] >= -0.5*Math.PI);
        assert(eulers[2] <= Math.PI && eulers[2] >= -Math.PI);
    }
    
    /**
     * Rotates a set of coordinates randomly.
     * @param coords coodinates to be rotated
     * @return rotated coordinates, may be the same object if no rotation has taken place
     */
    public static double[][] randomRotation(final double[][] coords){
        
        assert(coords.length == 3);
        assert(coords[0].length == coords[1].length);
        assert(coords[0].length == coords[2].length);
        
        final double[] eulers = new double[3];
        randomEulers(eulers);
        
        return CoordTranslation.rotateXYZ(coords, eulers);
    }
    
    /**
     * Gets a random point in a sphere of radius r.
     * @param radius the radius of the sphere
     * @param xyz the cartesian coordinates of the point in the sphere
     */
    public static void randomPointInSphere(final double radius, final double[] xyz){
        
        assert(xyz.length == 3);
        assert(radius > 0.0);
        
        final double[] spherical = new double[3];
        randomSphericals(spherical,radius);
        
        CoordTranslation.sphericalToCartesianCoord(spherical, xyz);
    }
    
    /**
     * Gives random spherical coordinates.
     * @param spherical The spherical coordinates, overriden on exit. Order: r/phi/omega
     * @param radius The radius of the sphere.
     */
    public static void randomSphericals(final double[] spherical, final double radius){
        
        assert(spherical.length == 3);
        spherical[0] = radius*r.nextDouble();
        spherical[1] = Math.PI * r.nextDouble();
        spherical[2] = 2*Math.PI*r.nextDouble();
    }
    
    public static void incrementDihedrals(final int[] whichDihedrals, final double[] dihedrals, final double maxStrength){
        
        assert(whichDihedrals.length <= dihedrals.length);
        assert(maxStrength >= 0.0);
        
        // manipulate the dihedrals
        for(final int which : whichDihedrals){
            final double incr = 2*Math.PI*maxStrength*r.nextDouble();
            // now add them on top, sanitizing if needed
            dihedrals[which] = (r.nextBoolean()) ? dihedrals[which] + incr : dihedrals[which] - incr;
            
            if(dihedrals[which] > Math.PI) {dihedrals[which] = dihedrals[which]-2*Math.PI;}
            else if(dihedrals[which] < -Math.PI) {dihedrals[which] = 2*Math.PI+dihedrals[which];}
            
            if(DEBUG){System.out.println("DEBUG: After sanitizing: new dihedral " + which +  " is " + dihedrals[which] + " or " + Math.toDegrees(dihedrals[which]));}
            
            assert(dihedrals[which] < Math.PI && dihedrals[which] > -Math.PI);
        }
    }
}
