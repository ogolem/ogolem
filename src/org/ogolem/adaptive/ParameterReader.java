/**
Copyright (c) 2014, J. M. Dieterich
              2017-2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive;

import org.ogolem.generic.IndividualReader;
import org.ogolem.io.InputPrimitives;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A reader for parameters.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class ParameterReader implements IndividualReader<AdaptiveParameters>{

    private static final long serialVersionUID = (long) 20200429;
    private static final Logger LOG = LoggerFactory.getLogger(ParameterReader.class);
    
    @Override
    public IndividualReader<AdaptiveParameters> copy() {
        return new ParameterReader();
    }

    @Override
    public void populateIndividualFromFile(final AdaptiveParameters individual, final String file) throws Exception {
        
        System.out.println("INFO: Please note, the reading from a seed for the parameters may not always work for purely technical reasons."
                + " If this should now break and you are sure that the parameter set is identical: ping the author(s).");
        
        final String[] sa = InputPrimitives.readFileIn(file);
        
        final AdaptiveParameters pTmp = new AdaptiveParameters(sa,-1);
        
        // do some checks (extend if necessary)
        if(pTmp.getNumberOfParamters() != individual.getNumberOfParamters()){
            throw new RuntimeException("Mismatch in number of parameters: " + pTmp.getNumberOfParamters() + " vs " + individual.getNumberOfParamters());
        }
        
        final String[] keysOrig = individual.getAllKeysCopy();
        final String[] keysTmp = pTmp.getAllKeysCopy();
        if(keysOrig.length != keysTmp.length){
            throw new RuntimeException("Not the same number of keys in read in and original parameter set " + keysOrig.length + " vs "
                + keysTmp.length);
        }
        
        for(final String key : keysOrig){
            final int origStart = individual.getStartPointForKey(key);
            final int tmpStart = pTmp.getStartPointForKey(key);
            if(origStart != tmpStart){
                throw new RuntimeException("Not the right offset in read in and original parameter set "
                    + origStart + " vs " + tmpStart + " for key " + key);
            }
        }
        
        final double[] pVals = pTmp.getAllParamters();
        individual.copyParameters(pVals);
        individual.setFitness(pTmp.getFitness());
        
        final String[] outP = individual.createFormattedOutput();
        String sCom = "";
        for(final String s : outP){
            sCom += s + "\n";
        }
        LOG.info("Read seed from file " + file);
        LOG.info(sCom); // dump it in the info stream
    }
}
