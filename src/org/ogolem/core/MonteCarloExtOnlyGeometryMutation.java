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
package org.ogolem.core;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.ogolem.generic.GenericMutation;
import org.ogolem.helpers.RandomUtils;

/**
 * A Monte Carlo mutation operator for a Geometry operating only in the external
 * coordinates (i.e., good only for non-flexible small molecules).
 * @author Johannes Dieterich
 * @version 2014-07-13
 */
public class MonteCarloExtOnlyGeometryMutation implements GenericMutation<Molecule,Geometry>{
    
    public static final int SINGLEMOVEMODE = 0;
    public static final int ALLMOVEMODE = 1;
    public static final int PARTIALMOVEMODE = 2;
    
    private static final long serialVersionUID = (long) 20140713;
    private final Random r = new Random();
    private final int mode;
    private final double maxMoveCOM;
    private final double maxMoveEuler;
    
    public MonteCarloExtOnlyGeometryMutation(final int mode, final double maxMoveCOM,
            final double maxMoveEuler){
        assert(maxMoveCOM >= 0.0);
        assert(maxMoveEuler >= 0.0);
        assert(maxMoveEuler <= 1.0);
        assert(mode == 0 || mode == 1 || mode == 2 || mode == 10 || mode == 11 || mode == 12);
        this.maxMoveCOM = maxMoveCOM;
        this.maxMoveEuler = maxMoveEuler;
        this.mode = mode;
    }
    
    public MonteCarloExtOnlyGeometryMutation(final MonteCarloExtOnlyGeometryMutation orig){
        this.maxMoveCOM = orig.maxMoveCOM;
        this.maxMoveEuler = orig.maxMoveEuler;
        this.mode = orig.mode;
    }

    @Override
    public MonteCarloExtOnlyGeometryMutation clone() {
        return new MonteCarloExtOnlyGeometryMutation(this);
    }

    @Override
    public String getMyID() {
        return "Monte Carlo (ext only) mutation: \n\tmode: " + mode
                + "\n\tmaxmove (COM): " + maxMoveCOM
                + "\n\tmaxmove (Euler): " + maxMoveEuler;
    }

    @Override
    public Geometry mutate(final Geometry orig) {
        
        final int noMols = orig.getNumberOfIndieParticles();
        
        List<Integer> moveables;
        if(mode == SINGLEMOVEMODE){
            final int which = r.nextInt(noMols);
            moveables = new ArrayList<>();
            moveables.add(which);
        } else if(mode == ALLMOVEMODE){
            moveables = new ArrayList<>();
            for(int i = 0; i < noMols; i++){
                moveables.add(i);
            }
        } else if(mode == PARTIALMOVEMODE){
            final int noPoints = r.nextInt(noMols);
            moveables = RandomUtils.listOfPoints(noPoints, noMols);
        } else {
            throw new RuntimeException("Illegal move mode " + mode + " in extonly MC mutation.");
        }
        
        
        final Geometry work = new Geometry(orig);
        final double[] comMove = new double[3];
        moveables.forEach((moveMol) -> {
            final double comM = maxMoveCOM*r.nextDouble();
            RandomUtils.randomVector(comMove, comM);
            
            final double eulerM = maxMoveEuler*r.nextDouble();
            
            final Molecule m = work.getMoleculeAtPosition(moveMol);
            final double[] origCOM = m.getExternalCenterOfMass();
            final double sign1 = (r.nextBoolean()) ? -1.0 : 1.0;
            origCOM[0] += sign1*comMove[0];
            origCOM[1] += sign1*comMove[1];
            origCOM[2] += sign1*comMove[2];
    
            RandomUtils.randomEulerIncrements(m.getOrientation(), eulerM);
        });
        
        return work;
    }
}
