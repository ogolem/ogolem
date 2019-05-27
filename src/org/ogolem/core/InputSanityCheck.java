/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

import org.ogolem.helpers.Tuple;

/**
 * Checks whether the input is sane.
 * @author Johannes Dieterich
 * @version 2014-05-02
 */
public final class InputSanityCheck {

    // do not allow instantiation
    private InputSanityCheck(){
    }

    // XXX could probably be extended ;-)
    // XXX spit out all warnings at the same time

    public static Tuple<Boolean,String[]> isInputSane(final GlobalConfig globConf){
        
        /*
         * first check for some general settings
         */
        if (globConf.maxCellSize[0] <= 0 || globConf.maxCellSize[1] <= 0 || globConf.maxCellSize[2] <= 0) {
            String[] message = {"Your cell is in at least one dimension too small. Fix it please."};

            return new Tuple<>(false, message);
        }

        if(globConf.poolSize <= 1){
            String[] message = {"Your pool size is too small. Fix it please."};

            return new Tuple<>(false, message);
        }

        if(globConf.blowFacDissocDetect <= 0.0){
            String[] message = {"Your blow factor for the dissociation detection is too small. Fix it please."};

            return new Tuple<>(false, message);
        }

        /*
         * check whether the global optimization and local optimization objects are set
         */
        if(globConf.opter == null){
            String[] message = {"No global optimization engine available."};

            return new Tuple<>(false, message);
        }

        if(globConf.refNewton == null){
            String[] message = {"No local optimization engine available."};

            return new Tuple<>(false, message);
        }

        /*
         * now check the geometry
         */
        final GeometryConfig gc = globConf.geoConf;
        if(gc == null || gc.geomMCs.isEmpty()){
            String[] message = {"No geometry definition existing."};
            
            return new Tuple<>(false,message);
        }

        for (final MoleculeConfig mc : gc.geomMCs) {
            if(mc.sID.equalsIgnoreCase("N/A") || mc.sID.isEmpty()){
                String[] message = {"Molecules String ID not properly set. Please inform developer(s)."};
                
                return new Tuple<>(false, message);
            }
            
            if(mc.constricted && mc.flexy) {
                String[] message = {"At the present moment, we do not support molecules which are both flexible and constricted.",
                "Please contact the author(s) if you really want to use this."};

                return new Tuple<>(false, message);
            }
        }
        
        /*
         * if we end up here, everything is fine
         */
        final Tuple<Boolean, String[]> result = new Tuple<>(true, new String[0]);

        return result;
    }
}
