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

/**
 * A really limited (but growing?) set of statistic primitives.
 * @author Johannes Dieterich
 * @version 2013-09-23
 */
public class StatisticUtils {
    
    private static final boolean DEBUG = true;
    private StatisticUtils(){}
    
    /**
     * Computes the arithmetic mean and the standard deviation (Bessel corrected) of a set of points.
     * @param dat
     * @return tupel of the mean and the std. dev. StdDev is null if too little points. Mean is null if no points.
     */
    public static Tuple<Double,Double> meanAndStdDev(final List<Double> dat){
        
        if(dat.isEmpty()){return new Tuple<>(null,null);}
        
        double mean = 0.0;
        for(final double d : dat){
            mean += d;
        }
        
        mean /= dat.size();
        if(dat.size() < 3){return new Tuple<>(mean,null);}
        
        double stdDev = 0;
        for(final double d : dat){
            final double diff = d-mean;
            stdDev += diff*diff;
        }
        
        stdDev /= (dat.size()-1);
        stdDev = Math.sqrt(stdDev);
        
        return new Tuple<>(mean,stdDev);
    }
    
    /**
     * Computes the arithmetic mean and the standard deviation (Bessel corrected) of a set of points.
     * @param dat
     * @return tupel of the mean and the std. dev. StdDev is null if too little points. Mean is null if no points.
     */
    public static Tuple<Double,Double> meanAndStdDev(final double[] dat){
        
        if(dat.length == 0){return new Tuple<>(null,null);}
        
        double mean = 0.0;
        for(final double d : dat){
            mean += d;
        }
        
        mean /= dat.length;
        if(dat.length < 3){return new Tuple<>(mean,null);}
        
        double stdDev = 0;
        for(final double d : dat){
            final double diff = d-mean;
            stdDev += diff*diff;
        }
        
        stdDev /= (dat.length-1);
        stdDev = Math.sqrt(stdDev);
        
        return new Tuple<>(mean,stdDev);
    }
    
    /**
     * Bins a set of data.
     * @param data the data as an array of doubles
     * @param noBins the number of bins
     * @return a tupel containg the binning information as a integer array of length noBins and the binning offsets as an double array of length noBins.
     */
    public static Tuple<int[],double[]> binData(final double[] data, final int noBins){
        
        assert(noBins > 0);
        assert(data != null);
        assert(data.length > 0);
        
        final int[] bins = new int[noBins];
        final double[] binOffs = new double[noBins];
        
        // find max and min
        double minD = Double.MAX_VALUE;
        double maxD = Double.NEGATIVE_INFINITY;
        
        for(final double d : data){
            minD = Math.min(minD,d);
            maxD = Math.max(maxD,d);
        }
        
        // compute the offsets
        final double incr = (maxD-minD)/(noBins);
        for(int i = 0; i < noBins; i++){
            binOffs[i] = incr*i;
        }
        
        // now bin
        for(final double d : data){
            int bin = (int) Math.floor((d-minD)/incr);
            if(bin >= noBins) {bin = noBins-1;}
            else if(bin < 0){bin = 0;}
            bins[bin] = bins[bin]+1;
        }
        
        return new Tuple<>(bins,binOffs);
    }
    
    /**
     * Bins a set of data.
     * @param data the data as an array of doubles
     * @param noBins the number of bins
     * @return a tupel containg the binning information as a integer array of length noBins and the binning offsets as an double array of length noBins.
     */
    public static Tuple<int[],double[]> binData(final List<Double> data, final int noBins){
        
        assert(noBins > 0);
        assert(data != null);
        assert(!data.isEmpty());
        
        final int[] bins = new int[noBins];
        final double[] binOffs = new double[noBins];
        
        // find max and min
        double minD = Double.MAX_VALUE;
        double maxD = Double.NEGATIVE_INFINITY;
        
        for(final double d : data){
            minD = Math.min(minD,d);
            maxD = Math.max(maxD,d);
        }
        
        // compute the offsets
        final double incr = (maxD-minD)/(noBins);
        for(int i = 0; i < noBins; i++){
            binOffs[i] = incr*i;
        }
        
        // now bin
        for(final double d : data){
            int bin = (int) Math.floor((d-minD)/incr);
            if(bin >= noBins) {bin = noBins-1;}
            else if(bin < 0){bin = 0;}
            bins[bin] = bins[bin]+1;
        }
        
        return new Tuple<>(bins,binOffs);
    }
}
