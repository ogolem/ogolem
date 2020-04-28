/**
Copyright (c) 2012-2014, J. M. Dieterich
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

import java.io.Serializable;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * Implements a mutation operator which in essence is only a MC step.
 * @author Johannes Dieterich
 * @version 2016-12-18
 */
public class MonteCarloMutation implements Serializable {
    
    private static final long serialVersionUID = (long) 20140327;
    
    private static final Lottery random = Lottery.getInstance();
    
    public static Molecule mutate(final Molecule mol, final int mode, final double maxMove){
        
        final CartesianCoordinates cartes = mol.getCartesians();
        final double[][] xyz = cartes.getAllXYZCoord();
        final boolean[][] constr = mol.getConstraints();
        
        if(mode == 0){
            moveOne(xyz,constr,maxMove);
        } else if(mode == 1){
            moveAll(xyz,constr,maxMove);
        } else if(mode == 2){
            moveSome(xyz,constr,maxMove);
        } else{
            System.err.println("ERROR: Unknown case in molecular MonteCarloMutation. Contact author(s)!");
        }
        
        final Molecule molNew = new Molecule (cartes, mol.getMolPosition(), mol.getSID(),
                mol.getFlexy(), mol.getDegreesOfFreedom(), mol.isConstricted(), mol.getConstraints());
        
        return molNew;
    }
    
    public static Geometry mutate(final Geometry g, final int mode, final double maxMove){
        
        if(mode < 10){
            final CartesianCoordinates cartes = new CartesianCoordinates(g.getCartesians());
            final double[][] xyz = cartes.getAllXYZCoord();
            final boolean[][] constr = g.getAllConstraintsXYZ(false);
        
            if(mode == 0){
                moveOne(xyz,constr,maxMove);
            } else if(mode == 1){
                moveAll(xyz,constr,maxMove);
            } else if(mode == 2){
                moveSome(xyz,constr,maxMove);
            } else{
                System.err.println("ERROR: Unknown case in geometric MonteCarloMutation. Contact author(s)!");
            }
            
            final Geometry gNew = new Geometry(cartes, g.getID(), g.getNumberOfIndieParticles(),
                cartes.getAllAtomsPerMol(), g.getAllFlexies(), g.getExplicitDoFs(), g.getAllConstraints(false),
                g.getAllConstraintsXYZ(false), g.getSIDs(), g.getBondInfo().clone());

            return gNew;
        } else {
            
            final Geometry gNew = new Geometry(g);
            final boolean[] constr = g.getAllConstraints(false);
            
            if(mode == 10){
                moveOne(gNew,constr,maxMove);
            } else if(mode == 11){
                moveAll(gNew,constr,maxMove);
            } else if(mode == 12){
                moveSome(gNew,constr,maxMove);
            } else{
                System.err.println("ERROR: Unknown case in geometric MonteCarloMutation2. Contact author(s)!");
            }
            
            return gNew;
        }
    }
    
    private static void moveOne(final double[][] xyz, final boolean[][] constr, final double maxMove){
        
        final double[] r = new double[3];
        
        int which;
        boolean cont = true;
        do{
            which = random.nextInt(xyz[0].length);
            if(!constr[0][which] || !constr[1][which] || !constr[2][which]) cont = false;
        } while(cont);
        
        // move
        randomizer(xyz,constr,maxMove,which,r);
    }
    
    private static void moveAll(final double[][] xyz, final boolean[][] constr, final double maxMove){
        
        final double[] r = new double[3];
        
        for(int i = 0; i < xyz[0].length; i++){
            if(!constr[0][i] || !constr[1][i] || !constr[2][i]){
                randomizer(xyz,constr,maxMove,i,r);
            }
        }
    }
    
    private static void moveSome(final double[][] xyz, final boolean[][] constr, final double maxMove){
        
        final double target = random.nextDouble();
        final double[] r = new double[3];
        
        for(int i = 0; i < xyz[0].length; i++){
            if((!constr[0][i] || !constr[1][i] || !constr[2][i]) && (target > random.nextDouble())){
                randomizer(xyz,constr,maxMove,i,r);
            }
        }
    }
    
    private static void moveOne(final Geometry g, final boolean[] constr, final double maxMove){
        
        int which;
        boolean cont = true;
        do{
            which = random.nextInt(g.getNumberOfIndieParticles());
            if(!constr[which]) cont = false;
        } while(cont);
        
        // move
        final Molecule newMol = mutate(g.getMoleculeAtPosition(which), 2, maxMove);
        g.setMoleculeAtPosition(which, newMol);
    }
    
    private static void moveAll(final Geometry g, final boolean[] constr, final double maxMove){
        
        for(int i = 0; i < g.getNumberOfIndieParticles(); i++){
            if(!constr[i]){
                final Molecule newMol = mutate(g.getMoleculeAtPosition(i), 2, maxMove);
                g.setMoleculeAtPosition(i, newMol);
            }
       }
    }
    
    private static void moveSome(final Geometry g, final boolean[] constr, final double maxMove){

        final double target = random.nextDouble();
        for(int i = 0; i < g.getNumberOfIndieParticles(); i++){
            if(!constr[i] && (target > random.nextDouble())){
                final Molecule newMol = mutate(g.getMoleculeAtPosition(i), 2, maxMove);
                g.setMoleculeAtPosition(i, newMol);
            }
       }
    }
    
    private static void randomizer(final double[][] xyz, final boolean[][] constr, final double maxMove, final int which, final double[] r){
        RandomUtils.randomVector(r);
        final boolean b = random.nextBoolean();
        final double scaling = (b) ? maxMove * random.nextDouble() : -maxMove * random.nextDouble();
        if(!constr[0][which]) xyz[0][which] += scaling*r[0];
        if(!constr[1][which]) xyz[1][which] += scaling*r[1];
        if(!constr[2][which]) xyz[2][which] += scaling*r[2];
    }
}
