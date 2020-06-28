/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
              2015-2016, J. M. Dieterich and B. Hartke
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

import java.util.List;

/**
 * Creates the fitting environment object based on the data provided.
 * @author Johannes Dieterich
 * @version 2016-02-29
 */
final class EnvironmentFactory {

    static enum KIND{SIMPLE};
    
    /**
     * Creates the object.
     * @param whichEnvType the kind of environment
     * @return Will return null, if it doesn't understand the specified environment
     * type.
     */
    static Environment createEnvironment(final KIND whichEnvType,
            final CartesianCoordinates envCartes, final boolean flexyEnvironment,
            final List<Integer> listSecAtoms, final double blowFacInitEnvBonds,
            final CollisionDetection.CDTYPE whichCollDetect, final DissociationDetection.DDTYPE whichDissDetect, final double blowColl,
            final double blowDiss, Atom[] referencePoints, AllowedSpace space){

        Environment env = null;

        switch(whichEnvType){
            /*
             * SIMPLE ENVIRONMENT
             */
            case SIMPLE:
                final BondInfo bonds = createEnvBonding(envCartes, blowFacInitEnvBonds);
                env = new SimpleEnvironment(envCartes, flexyEnvironment, bonds, listSecAtoms,
                        whichCollDetect,blowColl, referencePoints, space);
                break;
            /*
             * DEFAULT
             */
            default:
                System.err.println("ERROR: We are in an unknown enviroment type. " +
                        "Contact the author(s).");
                break;
        }
        
        return env;
    }

    private static BondInfo createEnvBonding(final CartesianCoordinates env,
            final double blowFac){
        
        return CoordTranslation.checkForBonds(env, blowFac);
    }
}
