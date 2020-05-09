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

import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;

/**
 * Calculates the niche based on the percentage of near 90 degree angles in the
 * clusters. Routine similar to the one implemented in Phenix and described in
 * B. Hartke, Z. Phys. Chem., 214, 1251 (2000).
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class WaterAngleNicheComp implements NicheComputer<Molecule,Geometry> {

    private static final long serialVersionUID = (long) 20200429;
    private final static double ANGLOW = Math.toRadians(80.0);
    private final static double ANGHIGH = Math.toRadians(100.0);
    private final static double PLOW = 50.0;
    private final static double PCUT = 75.0; // was: 50.0
    private final static double BONDSTSQ = 3.3*Constants.ANGTOBOHR;
    private final static double BONDST = BONDSTSQ*BONDSTSQ;
    private int nicheWidth = 1;

    WaterAngleNicheComp(){};
    WaterAngleNicheComp(final int nicheWidth){
        this.nicheWidth = nicheWidth;
    }
    
    private WaterAngleNicheComp(final WaterAngleNicheComp orig){
        this.nicheWidth = orig.nicheWidth;
    }
    
    @Override
    public WaterAngleNicheComp copy(){
        return new WaterAngleNicheComp(this);
    }
    
    /**
     * Makes a statistic of the bond angles for the cluster geometry in x,y,z,
     * but only for actually bonded triples. From this, a certain measure nproz
     * for the "clustering" of near-90-degree bond angles (and hence for the 
     * amount of "cubicity" in the geometry) is derived, see comments below.
     * @param g the geometry
     * @return the niche this geometry can be associated to 
     */
    @Override
    public Niche computeNiche(final Geometry g){
        
        // get cartesian coordinates (the combined arrays come in handy)
        final CartesianCoordinates c = g.getCartesians();
        final short[] atomNos = c.getAllAtomNumbers();
        final int[] noAtsPerMol = c.getAllAtomsPerMol();
        final int noMols = c.getNoOfMolecules();
        
        // tmp arrays
        final int[] nang = new int[noMols];
        final int[] nrect = new int[noMols];
        final double[] perc = new double[noMols];
        final int[] offs = new int[noMols];
        
        for(int i = 1; i < noMols; i++){
            offs[i] = offs[i-1]+noAtsPerMol[i-1];
        }
        
        int count = 0;
        int totcount = 0;
        for(int i = 0; i < noMols; i++){
            
            if(noAtsPerMol[i] != 3 || !allAtsCorrect(atomNos, offs[i])) continue;
            final double[] comI = g.getCOM(i);
            
            for(int j = 0; j < noMols; j++){
                if(i == j) continue;
                if(noAtsPerMol[j] != 3 || !allAtsCorrect(atomNos, offs[j])) continue;
                final double[] comJ = g.getCOM(j);
                
                final double dX = comI[0]-comJ[0];
                final double dY = comI[1]-comJ[1];
                final double dZ = comI[2]-comJ[2];
                final double dist = dX*dX+dY*dY+dZ*dZ;
                if(dist > BONDST) continue;
                
                for(int k = j+1; k < noMols; k++){
                    if(k == i) continue;
                    if(noAtsPerMol[k] != 3 || !allAtsCorrect(atomNos, offs[k])) continue;
                    final double[] comK = g.getCOM(k);
                
                    final double dX2 = comI[0]-comK[0];
                    final double dY2 = comI[1]-comK[1];
                    final double dZ2 = comI[2]-comK[2];
                    final double dist2 = dX2*dX2+dY2*dY2+dZ2*dZ2;
                    if(dist2 > BONDST) continue;
                    
                    // triple i,j,k is unique and actually bound; determine bond angle
                    final double tal=dX*dX2+dY*dY2+dZ*dZ2;
                    final double alpha=Math.acos(tal/Math.sqrt(dist)/Math.sqrt(dist2));
                    totcount++;
                    nang[i]++;
                    if (alpha >= ANGLOW && alpha <= ANGHIGH){
                        count++;
                        nrect[i]++;
                    }
                }
                
            }
        }
        
        // calculate percentage of near-90deg angles at each atom
        for(int i = 0; i < noMols; i++){
            perc[i] = nrect[i]/(double)nang[i]*100.0;
            // zero them to make them reusable
            nang[i] = 0;
            nrect[i] = 0;
        }
        
        // calculate total number of neighbors of each atom, and number of neighbors
        // for which both(!) atoms have a value above plow for the above percentage
        for(int i = 0; i < noMols-1; i++){
            final double[] comI = g.getCOM(i);
            for(int j = i+1; j < noMols; j++){
                final double[] comJ = g.getCOM(j);
                final double dX = comI[0]-comJ[0];
                final double dY = comI[1]-comJ[1];
                final double dZ = comI[2]-comJ[2];
                final double dist = dX*dX+dY*dY+dZ*dZ;
                if(dist > BONDST) continue;
                
                // i and j are bonded (interacting) neighbors
                nang[i]++;
                nang[j]++;
                if(perc[i] >= PLOW && perc[j] >= PLOW){
                    nrect[i]++;
                    nrect[j]++;
                }
            }
        }
        
        // calculate percentage of above-limit neighbors for each atom
        count = 0;
        for(int i = 0; i < noMols; i++){
            final double proz = nrect[i]/(double) nang[i]*100.0;
            if(proz >= PCUT) count++;
        }
        
        // finally, summarize this info into the percentage of atoms for which
        // the value proz is above pcut
        final int nproz = (int) Math.round(count/(double)noMols*100.0);

        //bxh: adding in a modification that makes the niches wider than just 1
        final int magic = nproz/nicheWidth; // INTENTIONAL integer division
        final int loEnd = magic*nicheWidth;
        final int hiEnd = loEnd + nicheWidth -1;
        
//        return new Niche("waterangle" + nproz);
        return new Niche("waterangle" + loEnd + "-" + hiEnd);
    }
    
    private boolean allAtsCorrect(final short[] atoms, final int offset){
        if(atoms.length < offset+3) return false; // safty net
        return (atoms[offset] == 8 && atoms[offset+1] == 1 && atoms[offset+2] == 1);
    }
}
