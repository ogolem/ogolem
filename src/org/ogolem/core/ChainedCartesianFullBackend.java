/**
Copyright (c) 2016-2017, J. M. Dieterich and B. Hartke
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

/**
 * Chains multiple cartesian full backends after each other.
 * @author Johannes Dieterich
 * @version 2017-03-03
 */
class ChainedCartesianFullBackend implements CartesianFullBackend {

    private static final long serialVersionUID = (long) 20160122;

    private final List<CartesianFullBackend> backends;
    
    ChainedCartesianFullBackend(final List<CartesianFullBackend> backends){
        
        assert(backends != null);
        assert(backends.size() > 0);
        
        this.backends = backends;
    }
    
    private ChainedCartesianFullBackend(final ChainedCartesianFullBackend orig){
        
        this.backends = new ArrayList<>();
        for(final CartesianFullBackend back : orig.backends){
            backends.add(back.clone());
        }
    }
    
    @Override
    public ChainedCartesianFullBackend clone() {
        return new ChainedCartesianFullBackend(this);
    }

    @Override
    public String getMethodID() {
        
        String s = "chained:\n";
        for(final CartesianFullBackend back : backends){
            s += "\t\t" + back.getMethodID() + "\n";
        }
        
        return s;
    }

    @Override
    public void gradientCalculation(final long lID, final int iIteration, final double[] xyz1D, final String[] saAtomTypes,
            final short[] atomNos, final int[] atsPerMol, final double[] energyparts, final int iNoOfAtoms, final float[] faCharges,
            final short[] spins, final BondInfo bonds, final Gradient grad, final boolean hasRigidEnv) {
        
        final Gradient workGrad = new Gradient(3,iNoOfAtoms);
        
        // just add up
        double e = 0.0;
        grad.zeroGradient();
        final double[][] gradMat = grad.getTotalGradient();
        for (final CartesianFullBackend back : backends){
            back.gradientCalculation(lID, iIteration, xyz1D, saAtomTypes, atomNos, atsPerMol, energyparts, iNoOfAtoms, faCharges, spins, bonds, workGrad, hasRigidEnv);
            
            final double[][] workGradMat = workGrad.getTotalGradient();
            for(int coord = 0; coord < 3; coord++){
                for(int at = 0; at < iNoOfAtoms; at++){
                    gradMat[coord][at] += workGradMat[coord][at];
                }
            }
            
            e += workGrad.getFunctionValue();
        }
        
        grad.setFunctionValue(e);
    }

    @Override
    public double energyCalculation(final long lID, final int iIteration, final double[] xyz1D, final String[] saAtomTypes,
            final short[] atomNos, final int[] atsPerMol, final double[] energyparts, final int iNoOfAtoms, final float[] faCharges,
            final short[] spins, final BondInfo bonds, final boolean hasRigidEnv) {
        
        // just add up
        double e = 0.0;
        for (final CartesianFullBackend back : backends){
            e += back.energyCalculation(lID, iIteration, xyz1D, saAtomTypes, atomNos, atsPerMol, energyparts, iNoOfAtoms, faCharges, spins, bonds, hasRigidEnv);
        }
        
        return e;
    }
}
