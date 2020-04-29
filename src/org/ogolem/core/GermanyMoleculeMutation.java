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

import java.util.List;
import org.ogolem.generic.GenericMutation;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * A molecule specific genotype mutation.
 * @author Johannes Dieterich
 * @version 2014-03-27
 */
public class GermanyMoleculeMutation implements GenericMutation<Double,Molecule> {
    
    public static final double DEFAULTBONDFACTOR = 0.25;
    
    private static final long serialVersionUID = (long) 20140327;
    private final Lottery random = Lottery.getInstance();
    private final boolean beAggressive;
    private final double bondFactor;
    
    GermanyMoleculeMutation(final boolean beAggressive, final double bondFactor){
        assert(bondFactor > 0.0);
        this.beAggressive = beAggressive;
        this.bondFactor = bondFactor;
    }
    
    GermanyMoleculeMutation(final GermanyMoleculeMutation orig){
        this.beAggressive = orig.beAggressive;
        this.bondFactor = orig.bondFactor;
    }

    @Override
    public GermanyMoleculeMutation clone() {
        return new GermanyMoleculeMutation(this);
    }

    @Override
    public String getMyID() {
        return "GERMANY\nmolecule genotype mutation\n\t aggressive? " + beAggressive
                + "\n\t bondFactor: " + bondFactor;
    }

    @Override
    public Molecule mutate(final Molecule orig) {
        if(beAggressive){
            return genotypeMutation(orig);
        } else{
            return mildGenotypeMutation(orig);
        }
    }
    
    private Molecule genotypeMutation(final Molecule start){
        
        //TODO this is too aggressive and does not work!

        final Molecule end = new Molecule(start);

        final ZMatrix zmat = start.getZMatrix();
        final boolean[][] dofs = start.getDegreesOfFreedom();

        int noDoFs = 0;
        for (final boolean[] dof : dofs) {
            for (int j = 0; j < dofs[0].length; j++) {
                if (dof[j]) {
                    noDoFs++;
                }
            }
        }

        if (noDoFs == 0) {
            return end;
        }

        final int muts = random.nextInt(noDoFs);
        final List<Integer> which = RandomUtils.listOfPoints(muts, noDoFs);
        
        //XXX this is a little inefficient. implement better.
        for(int mut : which){
            // translate back
            int c = 0;
            int kind = 2;
            int atom = 0;
            for (int i = 0; i < dofs.length; i++) {
                for (int j = 0; j < dofs[0].length; j++) {
                    if (dofs[i][j]) {
                        c++;
                    }
                    if (c == mut) {
                        kind = j;
                        atom = i;
                        break;
                    }
                }
            }

            // mutate
            final double d = random.nextDouble();
            double val = 0.0;
            if(kind == 0){
                //XXX this might need some more thinking and tinkering...
                val = zmat.getABondLength(atom) * d;
            } else if(kind == 1){
                val = Math.PI * d;
            } else if(kind == 2){
                final boolean b = random.nextBoolean();
                val = (b) ? Math.PI * d : -Math.PI * d;
            }

            zmat.setAValue(kind, atom, val);
        }

        end.setZMatrix(zmat);
        final CartesianCoordinates refCartes = zmat.translateToCartesianAndAlign(end.getCartesians());
        end.setReferenceCartesian(refCartes);

        return end;
    }

    private Molecule mildGenotypeMutation(final Molecule start){

        final Molecule endMolecule = new Molecule(start);

        // mutate
        final ZMatrix zmat = start.getZMatrix();

        final boolean[][] dofs = start.getDegreesOfFreedom();

        int howManyDoFs = 0;
        for (final boolean[] dofRow : dofs) {
            for (int j = 0; j < dofs[0].length; j++) {
                if (dofRow[j]) {
                    howManyDoFs++;
                }
            }
        }

        if (howManyDoFs == 0) {
            /*
             * we do nothing since this case can happen. the user
             * might, e.g., specify a molecule as being flexible
             * but forget to say which degree of freedom is used.
             */
            System.out.println("INFO: Molecular mutation started but no DoFs specified?!");
            return endMolecule;
        }
        

        // first drag a random telling which coord gets manipulated
        final int whichCoord = random.nextInt(howManyDoFs) + 1; // the random should be unequal zero

        // now translate back to what this means
        int tmp = 0;
        int whichKind = 2;
        int onWhichAtom = 0;
        for (int i = 0; i < dofs.length; i++) {
            for (int j = 0; j < dofs[0].length; j++) {
                if (dofs[i][j]) {
                    tmp++;
                }
                if (tmp == whichCoord) {
                    whichKind = j;
                    onWhichAtom = i;
                    break;
                }
            }
        }

        if (whichKind == 0) {
            /*
             * maxBondStretch
             */
            final double bond = zmat.getABondLength(onWhichAtom);
            final double factor = random.nextDouble();
            final boolean b = random.nextBoolean();
            
            final double newBond = (b) ? bond + factor*bondFactor : bond - factor*bondFactor;

            // put in
            zmat.setABondLength(onWhichAtom, newBond);
            
        } else {
            // some kind of angle, we should also be able make the angle *bigger*
            final boolean b = random.nextBoolean();
            final double angleFactor = (b) ? random.nextDouble() : random.nextDouble()+1;
            
            if(whichKind == 1){
                final double angle = zmat.getABondAngle(onWhichAtom);
                double newAngle = angle*angleFactor;
                if(newAngle < 0.0){
                    while(newAngle < 0.0){
                        newAngle += Math.PI;
                    }
                } else if(newAngle > Math.PI){
                    while(newAngle > Math.PI){
                        newAngle -= Math.PI;
                    }
                }
                zmat.setABondAngle(onWhichAtom, newAngle);
            } else {
                final double angle = zmat.getADihedral(onWhichAtom);
                double newAngle = angle*angleFactor;
                if(newAngle < -Math.PI){
                    while(newAngle < -Math.PI){
                        newAngle = Math.PI - (Math.PI+newAngle); // substract what is too negative from the upper bound
                    }
                } else if(newAngle > Math.PI){
                    while(newAngle < -Math.PI){
                        newAngle = -Math.PI + (newAngle-Math.PI); // add what is too positive to the lower bound
                    }
                }
                zmat.setADihedral(onWhichAtom, newAngle);
            }
        }
                
        endMolecule.setZMatrix(zmat);
        final CartesianCoordinates refCartes = zmat.translateToCartesianAndAlign(endMolecule.getCartesians());
        endMolecule.setReferenceCartesian(refCartes);

        return endMolecule;
    }
}
