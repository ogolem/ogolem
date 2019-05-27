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
package org.ogolem.adaptive;

import org.ogolem.generic.IndividualWriter;
import org.ogolem.io.OutputPrimitives;

/**
 * Writes one geometry out.
 * @author Johannes Dieterich
 * @version 2014-05-02
 */
public class ParameterWriter implements IndividualWriter<AdaptiveParameters> {
    
    private static final long serialVersionUID = (long) 20140502;
    private final String outputFolder;
    
    ParameterWriter(final String outputFolder){
        this.outputFolder = outputFolder;
    }
    
    @Override
    public ParameterWriter clone(){
        return new ParameterWriter(outputFolder);
    }
    
    @Override
    public void writeIndividual(final AdaptiveParameters p){
        // write xyz out
        try {
            final String[] sa = p.createFormattedOutput();
            OutputPrimitives.writeOut(outputFolder, sa, true);
        } catch (Exception e) {
            System.err.println("WARNING: Couldn't write parameters " + p.getID() + ". " + e.toString());
            e.printStackTrace(System.err);
        }
    }

    @Override
    public void writeIndividual(AdaptiveParameters p, String toFile) {
        // write xyz out
        try {
            final String[] sa = p.createFormattedOutput();
            OutputPrimitives.writeOut(toFile, sa, true);
        } catch (Exception e) {
            System.err.println("WARNING: Couldn't write parameters " + p.getID() + ". " + e.toString());
            e.printStackTrace(System.err);
        }
    }
}
