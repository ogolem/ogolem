/**
Copyright (c) 2016-2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.random;

import java.io.Serializable;
import java.util.Random;

/**
 * A wrapper around a random number generator (typically java.util.Random) to make random
 * numbers... less random.
 * @author Johannes Dieterich
 * @version 2020-04-25
 */
public final class Lottery implements Serializable {
    
    private static final long serialVersionUID = (long) 20200425;
    
    private static Lottery instance = null;
    
    private RNGenerator rng = null;
    
    private Lottery(final RNGenerator randNG){
        // no instantiation from outside!
        rng = randNG;
    }
    
    public static synchronized void setGenerator(final RNGenerator rng){
        instance = new Lottery(rng);
    }
    
    public static synchronized Lottery getInstance(){
        
        if(instance == null){
            System.err.println("Lottery got no random number generator! Initializing...");
            
            final Random r = new Random();
            final long seed = r.nextLong();
            final RNGenerator rngen= new StandardRNG(seed);
            System.err.println("Intialized: " + rngen.getInformation());
            
            instance = new Lottery(rngen);
            
            System.err.println("How did you get here (for developers!):");
            
            final StackTraceElement[] trace = new Throwable().getStackTrace();
            for(final StackTraceElement el : trace){
                System.err.println(el.toString());
            }
        }
        
        return instance;
    }
    
    public String getInformation(){
        return rng.getInformation();
    }
    
    public boolean nextBoolean(){
        return rng.nextBoolean();
    }
    
    public double nextDouble(){
        return rng.nextDouble();
    }
    
    public float nextFloat(){
        return rng.nextFloat();
    }
    
    public double nextGaussian(){
        return rng.nextGaussian();
    }
    
    public int nextInt(){
        return rng.nextInt();
    }
    
    public int nextInt(final int n){
        return rng.nextInt(n);
    }
    
    public long nextLong(){
        return rng.nextLong();
    }
}
