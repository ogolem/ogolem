/**
Copyright (c) 2015, J. M. Dieterich
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

import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericLocOpt;

/**
 * Builds fitness functions.
 * @author Johannes Dieterich
 * @version 2015-04-27
 */
public class FitnessFunctionFactory {
    
    public static GenericFitnessFunction<Molecule,Geometry> build(final GlobalConfig globConf, final GenericLocOpt<Molecule,Geometry> backend, final String fitnessFunctionConfig) throws Exception {
        
        if(fitnessFunctionConfig.startsWith("energy")){
            return backend;
        } else if(fitnessFunctionConfig.startsWith("spectral:")){

            int corrFrames = -1;
            int corrFramesPerBlock = -1;
            int corrGap = -1;
            String refSpecFile = null;
            int maxWaveNumbers = -1;
            double fullThresh = -1.0;
            boolean toBaseline = true;
            boolean doLocOpt = true;
            int start = 0;
            int end = 5000;
            int dipoleStepsToSnap = 1;

            final String sub = fitnessFunctionConfig.substring(9).trim();
            final String[] tokens = sub.split("\\;");
            for(final String token : tokens){
                if(token.startsWith("corrframes=")){
                    corrFrames = Integer.parseInt(token.substring("corrframes=".length()).trim());
                } else if(token.startsWith("corrframesperblock=")){
                    corrFramesPerBlock = Integer.parseInt(token.substring("corrframesperblock=".length()).trim());
                } else if(token.startsWith("corrgap=")){
                    corrGap = Integer.parseInt(token.substring("corrgap=".length()).trim());
                } else if(token.startsWith("maxwavenumbers=")){
                    maxWaveNumbers = Integer.parseInt(token.substring("maxwavenumbers=".length()).trim());
                } else if(token.startsWith("fullthresh=")){
                    fullThresh = Double.parseDouble(token.substring("fullthresh=".length()).trim());
                } else if(token.startsWith("refspecfile=")){
                     refSpecFile = token.substring("refspecfile=".length()).trim();
                } else if(token.startsWith("intervalstart=")){
                    start = Integer.parseInt(token.substring("intervalstart=".length()).trim());
                } else if(token.startsWith("intervalend=")){
                    end = Integer.parseInt(token.substring("intervalend=".length()).trim());
                } else if(token.equalsIgnoreCase("nottobaseline")){
                    toBaseline = false;
                } else if(token.equalsIgnoreCase("nocoefflocopt")){
                    doLocOpt = false;
                } else if(token.startsWith("dipolestepstosnap=")){
                    dipoleStepsToSnap = Integer.parseInt(token.substring("dipolestepstosnap=".length()).trim());
                } else {
                    throw new RuntimeException("Unknown token " + token + " for spectral optimization target.");
                }
            }

            if(corrFrames <= 0){throw new RuntimeException("spectral fitness function: corrframes must be positive.");}
            if(corrFramesPerBlock <= 0){throw new RuntimeException("spectral fitness function: corrframesperblock must be positive.");}
            if(corrGap <= 0){throw new RuntimeException("spectral fitness function: corrgap must be positive.");}
            if(refSpecFile == null){throw new RuntimeException("spectral fitness function: refspec must be defined ");}
            if(maxWaveNumbers <= 0){throw new RuntimeException("spectral fitness function: maxWaveNumbers must be positive.");}
            if(fullThresh <= 0.0){throw new RuntimeException("spectral fitness function: fullThresh must be positive.");}
            if(start > end){throw new RuntimeException("spectral fitness function: start must be smaller than end");}
            if(start < 0){throw new RuntimeException("spectral fitness function: start must be at least 0 wavenumbers");}
            if(end > maxWaveNumbers){throw new RuntimeException("spectral fitness function: end must be smaller than maxwavenumbers");}
            if(dipoleStepsToSnap < 1){throw new RuntimeException("spectral fitness function: dipolestepstosnap must be positive");}
            
            return new SingleGeomSpectralFitnessFunction(backend, corrFrames, corrFramesPerBlock, corrGap,
                refSpecFile, maxWaveNumbers, fullThresh, toBaseline, doLocOpt, start, end,
                dipoleStepsToSnap, globConf.mdConf, globConf.blowFacBondDetect, globConf.blowFacDissocDetect);
        }
                            
        throw new RuntimeException("Unknown optimization target " + fitnessFunctionConfig);
    }
}
