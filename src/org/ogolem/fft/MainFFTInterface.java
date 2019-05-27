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

import org.ogolem.io.InputPrimitives;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.TransformType;
import org.ogolem.helpers.Fortune;

/**
 * A VERY SIMPLE way of doing FFT transforms of data stored on disk. Why very simple?
 * We employ apache commons, 1D, in-core, padded-to-pow-of-two, so our code is really
 * just a wrapper and nothing else. Basically just for people (like me) that cannot
 * FFT functions in their head (like BXH).
 * @author Johannes Dieterich
 * @version 2013-12-23
 */
public class MainFFTInterface {
    
    public static void exec(final String[] args){
        
        if(args[0].equalsIgnoreCase("help")){
            System.out.println("This is the basic FFT transform interface of ogolem. Needed input:");
            System.out.println(" * the file name of the data");
            System.out.println(" * complex yes/no?");
            System.out.println(" * forward or backward?");
            System.exit(0);
        }
        
        final String fileN = args[0];
        final boolean isCompl = Boolean.parseBoolean(args[1]);
        final boolean fwdOrBwd = Boolean.parseBoolean(args[2]);
        
        if(isCompl){
            // we are operating on COMPLEX data
            System.out.println("# FFT transformed from COMPLEX data and file " + fileN);
            Complex[] unpadded = null;
            try {
                final String[] dataCont = InputPrimitives.readFileIn(fileN);
                unpadded = new Complex[dataCont.length];
                int c = 0;
                for (final String dat : dataCont) {
                    if (dat.trim().startsWith("#") || dat.trim().startsWith("//")) {
                        continue;
                    }
                    final String[] sa = dat.trim().split("\\s+");
                    final double real = Double.parseDouble(sa[0]);
                    final double imag = Double.parseDouble(sa[1]);
                    unpadded[c] = new Complex(real,imag);
                    c++;
                }
            } catch (Exception e) {
                System.err.println("ERROR: Failure in reading or parsing the FFT data.");
                e.printStackTrace(System.err);
                System.exit(1);
            }

            // pad the stuff
            final Complex[] padded = Padder.powTwoPad(unpadded);

            // FFT :-)
            final FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
            final Complex[] trafo = fft.transform(padded, (fwdOrBwd) ? TransformType.FORWARD : TransformType.INVERSE);

            System.out.println("# FFT transformed, padded results coming (first col real, second col imag part)");
            for (final Complex c : trafo) {
                System.out.println("" + c.getReal() + "   " + c.getImaginary());
            }
            System.out.println("# " + Fortune.randomFortune());
            
        } else{
            // we are operating on REAL data
            System.out.println("# FFT transformed from REAL data and file " + fileN);
            double[] unpadded = null;
            try {
                final String[] dataCont = InputPrimitives.readFileIn(fileN);
                unpadded = new double[dataCont.length];
                int c = 0;
                for (final String dat : dataCont) {
                    if (dat.trim().startsWith("#") || dat.trim().startsWith("//")) {
                        continue;
                    }
                    unpadded[c] = Double.parseDouble(dat.trim());
                    c++;
                }
            } catch (Exception e) {
                System.err.println("ERROR: Failure in reading or parsing the FFT data.");
                e.printStackTrace(System.err);
                System.exit(1);
            }

            // pad the stuff
            final double[] padded = Padder.powTwoPad(unpadded);

            // FFT :-)
            final FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
            final Complex[] trafo = fft.transform(padded, (fwdOrBwd) ? TransformType.FORWARD : TransformType.INVERSE);

            System.out.println("# FFT transformed, padded results coming (first col real, second col imag part)");
            for (final Complex c : trafo) {
                System.out.println("" + c.getReal() + "   " + c.getImaginary());
            }
            System.out.println("# " + Fortune.randomFortune());
        }
    }
}
