/**
Copyright (c) 2014, J. M. Dieterich
              2016, J. M. Dieterich and B. Hartke
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

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.ogolem.generic.GenericMutation;
import org.ogolem.helpers.RandomUtils;

/**
 * A orientation-based mutation.
 * @author Johannes Dieterich
 * @version 2016-03-22
 */
public class FoehrGeometryMutation implements GenericMutation<Molecule,Geometry> {
    
    private static final long serialVersionUID = (long) 20150217;

    public static enum MODUS {INCREMENT,FULLYRANDOM};
    public static enum HOWMANY {SINGLE, MULTIPLE, MULTIPLEGAUSS, ALL};
    
    private static final double INCRMUTSTRENGTH = 0.1; // XXX hard coded, but this seems to be a sane choice (10% of euler interval)
    
    private final Random r;
    private final MODUS mode;
    private final HOWMANY noMutations;
    
    FoehrGeometryMutation(final MODUS mode, final HOWMANY noOfMuts){
        this.mode = mode;
        this.noMutations = noOfMuts;
        this.r = new Random();
    }
    
    FoehrGeometryMutation(final FoehrGeometryMutation orig){
        this.mode = orig.mode;
        this.noMutations = orig.noMutations;
        this.r = new Random();
    }
    
    @Override
    public FoehrGeometryMutation clone() {
        return new FoehrGeometryMutation(this);
    }

    @Override
    public String getMyID() {        
        return "FOEHR GEOMETRY MUTATION:\n\t mutation mode: " + mode.toString();
    }

    @Override
    public Geometry mutate(final Geometry orig) {
        
        final int noMols = orig.getNumberOfIndieParticles();
        
        final List<Integer> mutations = new ArrayList<>();
        switch(noMutations){
            case SINGLE:
                final int which = r.nextInt(noMols);
                mutations.add(which);
                break;
            case MULTIPLE:
                final int howMany = r.nextInt(noMols);
                RandomUtils.listOfPoints(howMany, 0, noMols, mutations);
                break;
            case MULTIPLEGAUSS:
                final int howMany2 = (int) Math.round(RandomUtils.halfgaussDouble(0, 1)*noMols);
                RandomUtils.listOfPoints(howMany2, 0, noMols, mutations);
                break;
            case ALL:
                for(int i = 0; i < noMols; i++){mutations.add(i);}
                break;
        }
        
        
        final Geometry mutated = orig.clone();
        for(final int mutMol : mutations){
            
            final Molecule mol = mutated.getMoleculeAtPosition(mutMol);
            final double[] eulers = mol.getOrientation();
            
            switch(mode){
                case INCREMENT:
                    RandomUtils.randomEulerIncrements(eulers, INCRMUTSTRENGTH);
                    break;
                case FULLYRANDOM:
                    RandomUtils.randomEulers(eulers);
                    break;
            }
        }
        
        return mutated;
    }
}
