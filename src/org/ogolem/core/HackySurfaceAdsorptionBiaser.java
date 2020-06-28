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

/**
 * A very hacky way to enforce/bias towards a cluster actually adsorbing on top the surface
 * and not below it. XX/"0" atoms are understood as the lowest most layer and, if any
 * other atom is below this in z component, large energy and gradient contributions
 * are added for this atom to push it up.
 * This version uses a step potential for eIncr>0 && gradIncr=0,
 * or a linear ramp potential for eIncr=0 && gradIncr>0, or no potential otherwise.
 * @author Johannes Dieterich
 * @version 2017-03-03
 */
class HackySurfaceAdsorptionBiaser implements CartesianFullBackend {

    private static final long serialVersionUID = (long) 20161019;
    
    private final double eIncr;
    private final double gradIncr;
    
    HackySurfaceAdsorptionBiaser(final double eIncrement, final double gradIncrement){
        
        assert(eIncrement >= 0.0);
        assert(gradIncrement >= 0.0);
        
        this.eIncr = eIncrement;
        this.gradIncr = gradIncrement;
    }
    
    private HackySurfaceAdsorptionBiaser(final HackySurfaceAdsorptionBiaser orig){
        this.eIncr = orig.eIncr;
        this.gradIncr = orig.gradIncr;
    }
    
    @Override
    public HackySurfaceAdsorptionBiaser clone() {
        return new HackySurfaceAdsorptionBiaser(this);
    }

    @Override
    public String getMethodID() {
        return "hacky top surface adsorption biaser";
    }

    @Override
    public void gradientCalculation(final long lID, final int iIteration, final double[] xyz1D, final String[] saAtomTypes,
            final short[] atomNos, final int[] atsPerMol, final double[] energyparts, final int iNoOfAtoms, final float[] faCharges,
            final short[] spins, final BondInfo bonds, final Gradient grad, final boolean hasRigidEnv) {
        
        grad.zeroGradient();
        
        int beginEnv = 0;
        for(int i = 0; i < atsPerMol.length - 1; i++){
            beginEnv += atsPerMol[i];
        }
        
        final double[][] gradMat = grad.getTotalGradient();
        double eAddition = 0.0;
        for(int at1 = beginEnv; at1 < iNoOfAtoms; at1++){
            
            final short at1No = atomNos[at1];
            if(at1No != 0){continue;}
            
            final double z1 = xyz1D[2*iNoOfAtoms + at1];
            
            for(int at2 = 0; at2 < beginEnv; at2++){
                
                final short at2No = atomNos[at2];
                if(at2No == 0){continue;} // not between dummies
                
                final double z2 = xyz1D[2*iNoOfAtoms + at2];
                
                if(z2 < z1){
                    if(eIncr > 0.0 && gradIncr == 0.0){ // step potential (which has zero gradient)
                        eAddition += eIncr;
                    }else if(eIncr == 0.0 && gradIncr > 0.0){ // linear ramp potential
                        eAddition += gradIncr*(z1-z2);
                        gradMat[2][at2] -= gradIncr;                        
                    }
                    // else: no additional potential & gradient
                }
            }
        }
        
        grad.setTotalEnergy(eAddition);
    }

    @Override
    public double energyCalculation(final long lID, final int iIteration, final double[] xyz1D, final String[] saAtomTypes,
            final short[] atomNos, final int[] atsPerMol, final double[] energyparts, final int iNoOfAtoms, final float[] faCharges,
            final short[] spins, final BondInfo bonds, final boolean hasRigidEnv) {
        
        int beginEnv = 0;
        for(int i = 0; i < atsPerMol.length - 1; i++){
            beginEnv += atsPerMol[i];
        }
        
        double eAddition = 0.0;
        for(int at1 = beginEnv; at1 < iNoOfAtoms; at1++){
            
            final short at1No = atomNos[at1];
            if(at1No != 0){continue;}
            
            final double z1 = xyz1D[2*iNoOfAtoms + at1];
            
            for(int at2 = 0; at2 < beginEnv; at2++){
                
                final short at2No = atomNos[at2];
                if(at2No == 0){continue;} // not between dummies
                
                final double z2 = xyz1D[2*iNoOfAtoms + at2];
                
                if(z2 < z1){
                    if(eIncr > 0.0 && gradIncr == 0.0){ // step potential
                        eAddition += eIncr;
                    }else if(eIncr == 0.0 && gradIncr > 0.0){
                        eAddition += gradIncr*(z1-z2);                    
                    }
                    // else: no additional potential
                }
            }
        }
        
        return eAddition;
    }
}
