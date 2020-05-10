/**
Copyright (c) 2014, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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
import org.ogolem.generic.GenericInitializer;

/**
 * An adaptor from the old interface to the generic one. Hopefully temporary...
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class GeometryInitializationToGenericAdaptor implements GenericInitializer<Molecule, Geometry> {
    
    private static final long serialVersionUID = (long) 20200429;
    private final GenericFitnessFunction<Molecule,Geometry> fitness;
    private final GeometryInitialization init;
    private final CollisionDetection.CDTYPE whichCollDetect;
    private final DissociationDetection.DDTYPE whichDissDetect;
    private final double[] cellSize;
    private final double blowDiss;
    private final double blowColl;
    private final float explDoFRatio;
    private final boolean molecularCD;

    public GeometryInitializationToGenericAdaptor(final GeometryInitialization init, final CollisionDetection.CDTYPE whichCollDetect,
            final DissociationDetection.DDTYPE whichDissDetect, final double[] cellSize, final double blowDiss, final double blowColl,
            final float explDoFRatio, final boolean molecularCD, final GenericFitnessFunction<Molecule,Geometry> fitness){
        this.blowColl = blowColl;
        this.blowDiss = blowDiss;
        this.cellSize = cellSize;
        this.explDoFRatio = explDoFRatio;
        this.init = init;
        this.molecularCD = molecularCD;
        this.whichCollDetect = whichCollDetect;
        this.whichDissDetect = whichDissDetect;
        this.fitness = fitness;
    }
    
    GeometryInitializationToGenericAdaptor(final GeometryInitializationToGenericAdaptor orig){
        this.blowColl = orig.blowColl;
        this.blowDiss = orig.blowDiss;
        this.cellSize = orig.cellSize.clone();
        this.explDoFRatio = orig.explDoFRatio;
        this.init = orig.init; // does not support clone...
        this.molecularCD = orig.molecularCD;
        this.whichCollDetect = orig.whichCollDetect;
        this.whichDissDetect = orig.whichDissDetect;
        this.fitness = orig.fitness.copy();
    }
    
    @Override
    public GenericInitializer<Molecule, Geometry> copy() {
        return new GeometryInitializationToGenericAdaptor(this);
    }

    @Override
    public Geometry initialize(final Geometry ref, final long futureID) {
        
        final Geometry initGeom = init.initTheGeometry(ref, whichCollDetect,
                whichDissDetect, cellSize, blowDiss, blowColl, explDoFRatio,
                molecularCD);
        initGeom.setID(futureID);
        
        final Geometry opt = fitness.fitness(initGeom, false);
        
        return opt;               
    }
    
}
