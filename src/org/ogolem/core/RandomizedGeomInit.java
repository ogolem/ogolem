/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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

/**
 * Does a standard randomized geometry initialization.
 * @author Johannes Dieterich
 * @version 2015-04-28
 */
final class RandomizedGeomInit implements GeometryInitialization{

    private static final long serialVersionUID = (long) 20140401;
    private static final boolean DEBUG = false;
    private final boolean useDissDetect;

    RandomizedGeomInit(final boolean useAnyDissDetect){
        this.useDissDetect = useAnyDissDetect;
    }

    @Override
    public Geometry initTheGeometry(Geometry geom, final CollisionDetection.CDTYPE whichCollisionDetection,
            final DissociationDetection.DDTYPE whichDissociationDetection, final double[] cellSize, 
            final double dBlowDissDetect, final double dBlowCollDetect, final float explDoFRatio,
            final boolean molecularCD){

        int iEmergencyCounter = 0;
        boolean bFitsToEnv = false;
        while(!bFitsToEnv && iEmergencyCounter < FixedValues.MAXTOEMERGENCY){
            // we just transfer through to Geometry's own initializing routine
            geom.setRandomExtCoordMolecule(whichCollisionDetection, whichDissociationDetection,
                    cellSize, dBlowDissDetect, dBlowCollDetect, useDissDetect, explDoFRatio,
                    molecularCD);

            if(geom.containsEnvironment()){
                // init that one
                geom.initializeEnvironment();

                // check whether they fit together
                bFitsToEnv = geom.doesFitWithEnvironment();
            } else{
                // it'll always fit ;-)
                bFitsToEnv = true;
            }
            iEmergencyCounter++;
            // in case of an "emergency", the last geom will be returned and will fail to optimize, too bad ;-)
        }
        
        if(DEBUG){
            final String[] sa = geom.makePrintableAbsoluteCoord(true);
            System.out.println("DEBUG: Geometry " + geom.getID() + " after init.");
            for(final String s : sa){
                System.out.println(s);
            }
        }

        return geom;
    }
}
