/**
Copyright (c) 2012-2014, J. M. Dieterich
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

import static java.lang.Math.PI;
import static java.lang.Math.sqrt;
import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;

/**
 * Niche computation based on hydrogen bonds. Works only for hydrogen bonds
 * based on water molecules where atoms are defined as OHH in the input.
 * Ignores other building blocks/molecules.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class HydrogenBondNicheComp implements NicheComputer<Molecule,Geometry>{
    
    private static final long serialVersionUID = (long) 20200429;
    private static final boolean DEBUG = false;
    private final double maxDistOH;
    private final double maxDistOO;
    private final double maxAngle;
    
    public HydrogenBondNicheComp(final double maxDistOH, final double maxDistOO, final double maxAngle){
        this.maxDistOO = maxDistOO;
        this.maxDistOH = maxDistOH;
        this.maxAngle = maxAngle*PI/180.0;
    }
    
    private HydrogenBondNicheComp(final HydrogenBondNicheComp orig){
        this.maxAngle = orig.maxAngle;
        this.maxDistOH = orig.maxDistOH;
        this.maxDistOO = orig.maxDistOO;
    }
    
    @Override
    public HydrogenBondNicheComp copy(){
        return new HydrogenBondNicheComp(this);
    }
    
    @Override
    public Niche computeNiche(final Geometry g){

        // get cartesian coordinates
        final CartesianCoordinates c = g.getCartesians();
        final short[] atomNos = c.getAllAtomNumbers();
        final double[][] xyz = c.getAllXYZCoord();
        final int[] noAtsPerMol = c.getAllAtomsPerMol();
        
        int noHBonds = 0;
                
        int atOff1 = 0;
        for(int mol1 = 0; mol1 < c.getNoOfMolecules()-1; mol1++){
            if(noAtsPerMol[mol1] != 3 || !allAtsCorrect(atomNos, atOff1)){
                atOff1 += noAtsPerMol[mol1];
                continue;
            }
            
            int atOff2 = atOff1;
            for(int mol2 = mol1+1; mol2 < c.getNoOfMolecules(); mol2++){
                atOff2 += noAtsPerMol[mol2];
                if(noAtsPerMol[mol2] != 3 || !allAtsCorrect(atomNos, atOff2)){
                    //atOff2 += noAtsPerMol[mol2];
                    continue;
                }
                
                // check o-o distance (fastest criterion)
                final double ooX = xyz[0][atOff1]-xyz[0][atOff2];
                final double ooY = xyz[1][atOff1]-xyz[1][atOff2];
                final double ooZ = xyz[2][atOff1]-xyz[2][atOff2];          
                final double distOO = sqrt(ooX*ooX + ooY*ooY + ooZ*ooZ);
                if(DEBUG) System.out.println("DEBUG: dist " + distOO + " between " + atOff1 + " and " + atOff2);
                if(distOO > maxDistOO) continue;
                
                // check o-h distances and angles h-o-o (multiple permutations)
                final double ohO1X = xyz[0][atOff1]-xyz[0][atOff2+1];
                final double ohO1Y = xyz[1][atOff1]-xyz[1][atOff2+1];
                final double ohO1Z = xyz[2][atOff1]-xyz[2][atOff2+1];
                final double distO1 = sqrt(ohO1X*ohO1X + ohO1Y*ohO1Y + ohO1Z*ohO1Z);
                if(distO1 <= maxDistOH){
                    final double angle = CoordTranslation.calcAngle(xyz, atOff2+1, atOff2, atOff1);
                    if (angle <= maxAngle) {
                        noHBonds++;
                        continue;
                    }
                }
                
                final double ohO2X = xyz[0][atOff1]-xyz[0][atOff2+2];
                final double ohO2Y = xyz[1][atOff1]-xyz[1][atOff2+2];
                final double ohO2Z = xyz[2][atOff1]-xyz[2][atOff2+2];
                final double distO2 = sqrt(ohO2X*ohO2X + ohO2Y*ohO2Y + ohO2Z*ohO2Z);
                if(distO2 <= maxDistOH){
                    final double angle = CoordTranslation.calcAngle(xyz, atOff2+2, atOff2, atOff1);
                    if (angle <= maxAngle) {
                        noHBonds++;
                        continue;
                    }
                }
                
                final double oh1OX = xyz[0][atOff1+1]-xyz[0][atOff2];
                final double oh1OY = xyz[1][atOff1+1]-xyz[1][atOff2];
                final double oh1OZ = xyz[2][atOff1+1]-xyz[2][atOff2];
                final double dist1O = sqrt(oh1OX*oh1OX + oh1OY*oh1OY + oh1OZ*oh1OZ);
                if(dist1O <= maxDistOH){
                    final double angle = CoordTranslation.calcAngle(xyz, atOff1+1, atOff1, atOff2);
                    if (angle <= maxAngle) {
                        noHBonds++;
                        continue;
                    }
                }
                
                final double oh2OX = xyz[0][atOff1+2]-xyz[0][atOff2];
                final double oh2OY = xyz[1][atOff1+2]-xyz[1][atOff2];
                final double oh2OZ = xyz[2][atOff1+2]-xyz[2][atOff2];
                final double dist2O = sqrt(oh2OX*oh2OX + oh2OY*oh2OY + oh2OZ*oh2OZ);
                if(dist2O <= maxDistOH){
                    final double angle = CoordTranslation.calcAngle(xyz, atOff1+2, atOff1, atOff2);
                    if (angle <= maxAngle) {
                        noHBonds++;
                        continue;
                    }
                }
                
                //atOff2 += noAtsPerMol[mol2];
            }
            atOff1 += noAtsPerMol[mol1];
        }
        
        return new Niche("hbonds" + noHBonds);
    }
    
    private boolean allAtsCorrect(final short[] atoms, final int offset){
        if(atoms.length < offset+3) return false; // safty net
        return (atoms[offset] == 8 && atoms[offset+1] == 1 && atoms[offset+2] == 1);
    }
}
