/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2015, J. M. Dieterich
              2019-2020, J. M. Dieterich and B. Hartke
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
import static org.ogolem.core.Constants.ANGTOBOHR;
import static org.ogolem.core.Constants.BOHRTOANG;
import org.ogolem.math.GenericLookup;

/**
 * Provides the TIP3P force field for an simple (and very biased) description of water clusters.
 * Numerical data from http://www1.lsbu.ac.uk/water/water_models.html and wiki
 * @author Johannes Dieterich
 * @version 2020-02-08
 */
public class TIP3PForceField implements RigidBodyBackend {

    // the ID
    private static final long serialVersionUID = (long) 20200208;
    private static final boolean DEBUG = false;
    
    // the LJ cutoff in bohr
    private static final double CUTOFFLJSQ = 10*ANGTOBOHR*10*ANGTOBOHR;
    
    // partial charge on H's
    private final double CHARGEH;
    
    // partial charges on O's
    private final double CHARGEO;
    
    private final double CHARGEOO;
    private final double CHARGEOH;
    private final double CHARGEHH;
    
    // unit is 10^-3 kcal*A^{12}/mol before the translation
    private final double LJAOXYGEN;

    //unit is kcal * A^{6}/mol before the translation
    private final double LJBOXYGEN;
    
    private final double rOH;
    private final double angHOH;

    private static final double ANTICOLDEND = 0.0; // 4.0 bohr
    private final TIPnPHelpers.AntiColdFusionFunction anticold;
    private final GenericLookup anticoldMed;    
    
    public TIP3PForceField(){
        // use standard TIP3P parameters
        final TIPnPParameters.TIP3PParameters params = new TIPnPParameters.StandardTIP3PParameters();
        this.CHARGEH = params.getChargeH();
        this.CHARGEO = params.getChargeO();
        this.CHARGEHH = CHARGEH*CHARGEH;
        this.CHARGEOH = CHARGEO*CHARGEH;
        this.CHARGEOO = CHARGEO*CHARGEO;
        this.LJAOXYGEN = params.getLJA();
        this.LJBOXYGEN = params.getLJB();
        this.rOH = params.getOH();
        this.angHOH = params.getHOH();
        // create Bernd's anti cold fusion term
        this.anticold = new TIPnPHelpers.AntiColdFusionFunction();
        this.anticoldMed = new GenericLookup(anticold,GenericLookup.POSSIBLE_SIZES[1],0.0,ANTICOLDEND); // we want a very compact table
    }

    public TIP3PForceField(final TIPnPParameters.TIP3PParameters params){
        this.CHARGEH = params.getChargeH();
        this.CHARGEO = params.getChargeO();
        this.CHARGEHH = CHARGEH*CHARGEH;
        this.CHARGEOH = CHARGEO*CHARGEH;
        this.CHARGEOO = CHARGEO*CHARGEO;
        this.LJAOXYGEN = params.getLJA();
        this.LJBOXYGEN = params.getLJB();
        this.rOH = params.getOH();
        this.angHOH = params.getHOH();
        // create Bernd's anti cold fusion term
        this.anticold = new TIPnPHelpers.AntiColdFusionFunction();
        this.anticoldMed = new GenericLookup(anticold,GenericLookup.POSSIBLE_SIZES[1],0.0,ANTICOLDEND); // we want a very compact table
    }
    
    TIP3PForceField(final TIP3PForceField orig){
        this.CHARGEH = orig.CHARGEH;
        this.CHARGEHH = orig.CHARGEHH;
        this.CHARGEO = orig.CHARGEO;
        this.CHARGEOH = orig.CHARGEOH;
        this.CHARGEOO = orig.CHARGEOO;
        this.LJAOXYGEN = orig.LJAOXYGEN;
        this.LJBOXYGEN = orig.LJBOXYGEN;
        this.rOH = orig.rOH;
        this.angHOH = orig.angHOH;
        // create Bernd's anti cold fusion term
        this.anticold = new TIPnPHelpers.AntiColdFusionFunction();
        this.anticoldMed = new GenericLookup(anticold,GenericLookup.POSSIBLE_SIZES[1],0.0,ANTICOLDEND); // we want a very compact table
    }
    
    @Override
    public TIP3PForceField clone(){
        return new TIP3PForceField(this);
    }

    @Override
    public String getMethodID(){
        return "TIP3P (dirty implementation)";
    }

    /**
     * Checks if the geometry is actually suitable for TIP3P
     * @param ref the reference geometry
     * @return true, if all molecules are waters in atom order O/H/H
     */
    @Override
    public boolean suitableForThisBackend(final Geometry ref){
        
        assert(ref != null);
        assert(ref.getNumberOfIndieParticles() > 0);
        
        // first figure out whether the the first three could be a water
        final short[] atomNos = ref.getMoleculeAtPosition(0).getAtomNumbers();
        if(!(atomNos.length == 3 && atomNos[0] == 8
                && atomNos[1] == 1
                && atomNos[2] == 1)){
            System.err.println("ERROR: We figured out that this is no suitable water (in O/H/H order) cluster!");
            return false;
        }

        // now the rest
        for(int mol = 0; mol < ref.getNumberOfIndieParticles(); mol++){
            
            final short[] atomNosMol = ref.getMoleculeAtPosition(mol).getAtomNumbers();
            // check whether all molecules are having the right order
            if (!(atomNosMol.length == 3 && atomNosMol[0] == 8 && atomNosMol[1] == 1 && atomNosMol[2] == 1)) {
                // something is f***** up
                System.err.println("ERROR: Some parts of your cluster aren't water (in O/H/H order).");
                return false;
            }
        }
        
        return true;
    }
    
    @Override
    public CartesianCoordinates adjustCartesians(final Geometry ref, final CartesianCoordinates unadjusted, final int molID){
        return unadjusted;
    }
    
    @Override
    public void rigidify(final Molecule mol){
        TIPnPHelpers.sanitizeWater(mol, rOH, angHOH);
    }

    @Override
    public double energy(final Geometry ref, final List<CartesianCoordinates> cartes,
            final double[][] coms, final int counter){
        
        double energy = 0.0;
        for(int mol1 = 0; mol1 < cartes.size(); mol1++){
            
            // adjust the first water molecules coordinates
            final CartesianCoordinates c1 = cartes.get(mol1);
            final double[] com1 = coms[mol1];
            final double[][] xyz1O = c1.getAllXYZCoord();
            
            final double o1X = xyz1O[0][0] + com1[0];
            final double o1Y = xyz1O[1][0] + com1[1];
            final double o1Z = xyz1O[2][0] + com1[2];
            final double h11X = xyz1O[0][1] + com1[0];
            final double h11Y = xyz1O[1][1] + com1[1];
            final double h11Z = xyz1O[2][1] + com1[2];
            final double h12X = xyz1O[0][2] + com1[0];
            final double h12Y = xyz1O[1][2] + com1[1];
            final double h12Z = xyz1O[2][2] + com1[2];
            
            for(int mol2 = mol1+1; mol2 < cartes.size(); mol2++){
                
                // adjust the second water molecules coordinates
                final double[] com2 = coms[mol2];
                final CartesianCoordinates c2 = cartes.get(mol2);
                final double[][] xyz2O = c2.getAllXYZCoord();
                
                final double o2X = xyz2O[0][0] + com2[0];
                final double o2Y = xyz2O[1][0] + com2[1];
                final double o2Z = xyz2O[2][0] + com2[2];
                final double h21X = xyz2O[0][1] + com2[0];
                final double h21Y = xyz2O[1][1] + com2[1];
                final double h21Z = xyz2O[2][1] + com2[2];
                final double h22X = xyz2O[0][2] + com2[0];
                final double h22Y = xyz2O[1][2] + com2[1];
                final double h22Z = xyz2O[2][2] + com2[2];
                                
                // XXX cutoff for the electrostatics
                // for very large clusters, it may pay off to just compute the diff between the COMs and based on that cut off at e.g. 40 Angstrom
                
                // get some distance diffs
                final double diffOOx = o1X - o2X;
                final double diffOOy = o1Y - o2Y;
                final double diffOOz = o1Z - o2Z;
                final double distOOSq = diffOOx*diffOOx + diffOOy*diffOOy + diffOOz*diffOOz;
                final double distOO = Math.sqrt(distOOSq);
                if(distOO <= ANTICOLDEND){energy += anticoldMed.inter(distOO);}
                
                final double diffH1H1x = h11X - h21X;
                final double diffH1H1y = h11Y - h21Y;
                final double diffH1H1z = h11Z - h21Z;
                final double distH1H1Sq = diffH1H1x*diffH1H1x + diffH1H1y*diffH1H1y + diffH1H1z*diffH1H1z;
                final double distH1H1 = Math.sqrt(distH1H1Sq);
                if(distH1H1 <= ANTICOLDEND){energy += anticoldMed.inter(distH1H1);}
                
                final double diffH2H2x = h12X - h22X;
                final double diffH2H2y = h12Y - h22Y;
                final double diffH2H2z = h12Z - h22Z;
                final double distH2H2Sq = diffH2H2x*diffH2H2x + diffH2H2y*diffH2H2y + diffH2H2z*diffH2H2z;
                final double distH2H2 = Math.sqrt(distH2H2Sq);
                if(distH2H2 <= ANTICOLDEND){energy += anticoldMed.inter(distH2H2);}
                
                final double diffH1H2x = h11X - h22X;
                final double diffH1H2y = h11Y - h22Y;
                final double diffH1H2z = h11Z - h22Z;
                final double distH1H2Sq = diffH1H2x*diffH1H2x + diffH1H2y*diffH1H2y + diffH1H2z*diffH1H2z;
                final double distH1H2 = Math.sqrt(distH1H2Sq);
                if(distH1H2 <= ANTICOLDEND){energy += anticoldMed.inter(distH1H2);}
                
                final double diffH2H1x = h12X - h21X;
                final double diffH2H1y = h12Y - h21Y;
                final double diffH2H1z = h12Z - h21Z;
                final double distH2H1Sq = diffH2H1x*diffH2H1x + diffH2H1y*diffH2H1y + diffH2H1z*diffH2H1z;
                final double distH2H1 = Math.sqrt(distH2H1Sq);
                if(distH2H1 <= ANTICOLDEND){energy += anticoldMed.inter(distH2H1);}
                
                final double diffOH1x = o1X - h21X;
                final double diffOH1y = o1Y - h21Y;
                final double diffOH1z = o1Z - h21Z;
                final double distOH1Sq = diffOH1x*diffOH1x + diffOH1y*diffOH1y + diffOH1z*diffOH1z;
                final double distOH1 = Math.sqrt(distOH1Sq);
                if(distOH1 <= ANTICOLDEND){energy += anticoldMed.inter(distOH1);}
                
                final double diffOH2x = o1X - h22X;
                final double diffOH2y = o1Y - h22Y;
                final double diffOH2z = o1Z - h22Z;
                final double distOH2Sq = diffOH2x*diffOH2x + diffOH2y*diffOH2y + diffOH2z*diffOH2z;
                final double distOH2 = Math.sqrt(distOH2Sq);
                if(distOH2 <= ANTICOLDEND){energy += anticoldMed.inter(distOH2);}
                
                final double diffH1Ox = h11X - o2X;
                final double diffH1Oy = h11Y - o2Y;
                final double diffH1Oz = h11Z - o2Z;
                final double distH1OSq = diffH1Ox*diffH1Ox + diffH1Oy*diffH1Oy + diffH1Oz*diffH1Oz;
                final double distH1O = Math.sqrt(distH1OSq);
                if(distH1O <= ANTICOLDEND){energy += anticoldMed.inter(distH1O);}
                
                final double diffH2Ox = h12X - o2X;
                final double diffH2Oy = h12Y - o2Y;
                final double diffH2Oz = h12Z - o2Z;
                final double distH2OSq = diffH2Ox*diffH2Ox + diffH2Oy*diffH2Oy + diffH2Oz*diffH2Oz;
                final double distH2O = Math.sqrt(distH2OSq);
                if(distH2O <= ANTICOLDEND){energy += anticoldMed.inter(distH2O);}
                
                // now compute the actual interaction energies
                // first: the electrostatics, we assume that we never can cut that one off
                energy += (CHARGEOO)/ distOO + (CHARGEOH)/ distH1O + (CHARGEOH)/ distH2O
                        + (CHARGEOH)/ distOH1 + (CHARGEOH)/ distOH2 + (CHARGEHH)/ distH1H1
                        + (CHARGEHH)/ distH2H2 + + (CHARGEHH)/ distH1H2 + (CHARGEHH)/ distH2H1;
                
                // second: the LJ potential, we assume we CAN cut this off based on the COM-COM distance of the O's                
                if(distOOSq >= CUTOFFLJSQ){continue;}
                
                final double distOO6 = distOOSq*distOOSq*distOOSq;
                final double distOO12 = distOO6*distOO6;
                
                energy += LJAOXYGEN/distOO12 - LJBOXYGEN/distOO6;
            }
        }
        
        return energy;
    }

    @Override
    public double gradient(final Geometry ref, final List<CartesianCoordinates> cartes,
            final double[][] coms, final List<double[][]> gradient, final int counter){
        
        assert(gradient != null);
        assert(gradient.size() == cartes.size());

        for(int i = 0; i < gradient.size(); i++){
            final double[][] gradPart = gradient.get(i);
            for(int coord = 0; coord < 3; coord++){
                for(int atom = 0; atom < gradPart[0].length; atom++){
                    gradPart[coord][atom] = 0.0;
                }
            }
        }
        
        double energy = 0.0;
        for(int mol1 = 0; mol1 < cartes.size(); mol1++){
            
            // adjust the first water molecules coordinates
            final CartesianCoordinates c1 = cartes.get(mol1);
            final double[] com1 = coms[mol1];
            final double[][] xyz1O = c1.getAllXYZCoord();
            
            final double o1X = xyz1O[0][0] + com1[0];
            final double o1Y = xyz1O[1][0] + com1[1];
            final double o1Z = xyz1O[2][0] + com1[2];
            final double h11X = xyz1O[0][1] + com1[0];
            final double h11Y = xyz1O[1][1] + com1[1];
            final double h11Z = xyz1O[2][1] + com1[2];
            final double h12X = xyz1O[0][2] + com1[0];
            final double h12Y = xyz1O[1][2] + com1[1];
            final double h12Z = xyz1O[2][2] + com1[2];
                        
            final double[][] gradMol1 = gradient.get(mol1);
            
            for(int mol2 = mol1+1; mol2 < cartes.size(); mol2++){
                
                // adjust the second water molecules coordinates
                final double[] com2 = coms[mol2];
                final CartesianCoordinates c2 = cartes.get(mol2);
                final double[][] xyz2O = c2.getAllXYZCoord();
                
                final double o2X = xyz2O[0][0] + com2[0];
                final double o2Y = xyz2O[1][0] + com2[1];
                final double o2Z = xyz2O[2][0] + com2[2];
                final double h21X = xyz2O[0][1] + com2[0];
                final double h21Y = xyz2O[1][1] + com2[1];
                final double h21Z = xyz2O[2][1] + com2[2];
                final double h22X = xyz2O[0][2] + com2[0];
                final double h22Y = xyz2O[1][2] + com2[1];
                final double h22Z = xyz2O[2][2] + com2[2];
                
                final double[][] gradMol2 = gradient.get(mol2);
                                
                // XXX cutoff for the electrostatics
                // for very large clusters, it may pay off to just compute the diff between the COMs and based on that cut off at e.g. 40 Angstrom
                
                // get some distance diffs
                final double diffOOx = o1X - o2X;
                final double diffOOy = o1Y - o2Y;
                final double diffOOz = o1Z - o2Z;
                final double distOOSq = diffOOx*diffOOx + diffOOy*diffOOy + diffOOz*diffOOz;
                final double distOO = Math.sqrt(distOOSq);
                if(distOO <= ANTICOLDEND){energy += anticoldMed.inter(distOO);}
                
                final double diffH1H1x = h11X - h21X;
                final double diffH1H1y = h11Y - h21Y;
                final double diffH1H1z = h11Z - h21Z;
                final double distH1H1Sq = diffH1H1x*diffH1H1x + diffH1H1y*diffH1H1y + diffH1H1z*diffH1H1z;
                final double distH1H1 = Math.sqrt(distH1H1Sq);
                if(distH1H1 <= ANTICOLDEND){energy += anticoldMed.inter(distH1H1);}
                
                final double diffH2H2x = h12X - h22X;
                final double diffH2H2y = h12Y - h22Y;
                final double diffH2H2z = h12Z - h22Z;
                final double distH2H2Sq = diffH2H2x*diffH2H2x + diffH2H2y*diffH2H2y + diffH2H2z*diffH2H2z;
                final double distH2H2 = Math.sqrt(distH2H2Sq);
                if(distH2H2 <= ANTICOLDEND){energy += anticoldMed.inter(distH2H2);}
                
                final double diffH1H2x = h11X - h22X;
                final double diffH1H2y = h11Y - h22Y;
                final double diffH1H2z = h11Z - h22Z;
                final double distH1H2Sq = diffH1H2x*diffH1H2x + diffH1H2y*diffH1H2y + diffH1H2z*diffH1H2z;
                final double distH1H2 = Math.sqrt(distH1H2Sq);
                if(distH1H2 <= ANTICOLDEND){energy += anticoldMed.inter(distH1H2);}
                
                final double diffH2H1x = h12X - h21X;
                final double diffH2H1y = h12Y - h21Y;
                final double diffH2H1z = h12Z - h21Z;
                final double distH2H1Sq = diffH2H1x*diffH2H1x + diffH2H1y*diffH2H1y + diffH2H1z*diffH2H1z;
                final double distH2H1 = Math.sqrt(distH2H1Sq);
                if(distH2H1 <= ANTICOLDEND){energy += anticoldMed.inter(distH2H1);}
                
                final double diffOH1x = o1X - h21X;
                final double diffOH1y = o1Y - h21Y;
                final double diffOH1z = o1Z - h21Z;
                final double distOH1Sq = diffOH1x*diffOH1x + diffOH1y*diffOH1y + diffOH1z*diffOH1z;
                final double distOH1 = Math.sqrt(distOH1Sq);
                if(distOH1 <= ANTICOLDEND){energy += anticoldMed.inter(distOH1);}
                
                final double diffOH2x = o1X - h22X;
                final double diffOH2y = o1Y - h22Y;
                final double diffOH2z = o1Z - h22Z;
                final double distOH2Sq = diffOH2x*diffOH2x + diffOH2y*diffOH2y + diffOH2z*diffOH2z;
                final double distOH2 = Math.sqrt(distOH2Sq);
                if(distOH2 <= ANTICOLDEND){energy += anticoldMed.inter(distOH2);}
                
                final double diffH1Ox = h11X - o2X;
                final double diffH1Oy = h11Y - o2Y;
                final double diffH1Oz = h11Z - o2Z;
                final double distH1OSq = diffH1Ox*diffH1Ox + diffH1Oy*diffH1Oy + diffH1Oz*diffH1Oz;
                final double distH1O = Math.sqrt(distH1OSq);
                if(distH1O <= ANTICOLDEND){energy += anticoldMed.inter(distH1O);}
                
                final double diffH2Ox = h12X - o2X;
                final double diffH2Oy = h12Y - o2Y;
                final double diffH2Oz = h12Z - o2Z;
                final double distH2OSq = diffH2Ox*diffH2Ox + diffH2Oy*diffH2Oy + diffH2Oz*diffH2Oz;
                final double distH2O = Math.sqrt(distH2OSq);
                if(distH2O <= ANTICOLDEND){energy += anticoldMed.inter(distH2O);}
                
                // now compute the actual interaction energies
                // first: the electrostatics, we assume that we never can cut that one off
                energy += (CHARGEOO)/ distOO + (CHARGEOH)/ distH1O + (CHARGEOH)/ distH2O
                        + (CHARGEOH)/ distOH1 + (CHARGEOH)/ distOH2 + (CHARGEHH)/ distH1H1
                        + (CHARGEHH)/ distH2H2 + + (CHARGEHH)/ distH1H2 + (CHARGEHH)/ distH2H1;
                
                // gradient for the electrostatics                
                final double gradEOO   = -(CHARGEOO)/ distOOSq;
                final double gradEH1O  = -(CHARGEOH)/ distH1OSq;
                final double gradEH2O  = -(CHARGEOH)/ distH2OSq;
                final double gradEOH1  = -(CHARGEOH)/ distOH1Sq;
                final double gradEOH2  = -(CHARGEOH)/ distOH2Sq;
                final double gradEH1H1 = -(CHARGEHH)/ distH1H1Sq;
                final double gradEH2H2 = -(CHARGEHH)/ distH2H2Sq;
                final double gradEH1H2 = -(CHARGEHH)/ distH1H2Sq;
                final double gradEH2H1 = -(CHARGEHH)/ distH2H1Sq;
                
                // divide the distances in all three dimensions by the total distance to cover for the coordinate system
                
                final double h1OG = gradEH1O/distH1O;
                final double h1H1G = gradEH1H1/distH1H1;
                final double h1H2G = gradEH1H2/distH1H2;
                
                final double h1OxG = h1OG*diffH1Ox;
                final double h1H1xG = h1H1G*diffH1H1x;
                final double h1H2xG = h1H2G*diffH1H2x;
                
                gradMol1[0][1] += h1OxG;
                gradMol2[0][0] -= h1OxG;
                gradMol1[0][1] += h1H1xG;
                gradMol2[0][1] -= h1H1xG;
                gradMol1[0][1] += h1H2xG;
                gradMol2[0][2] -= h1H2xG;
                
                final double h2OG = gradEH2O/distH2O;
                final double h2H2G = gradEH2H2/distH2H2;
                final double h2H1G = gradEH2H1/distH2H1;
                
                final double h2OxG = h2OG*diffH2Ox;
                final double h2H2xG = h2H2G*diffH2H2x;
                final double h2H1xG = h2H1G*diffH2H1x;
                
                gradMol1[0][2] += h2OxG;
                gradMol2[0][0] -= h2OxG;
                gradMol1[0][2] += h2H2xG;
                gradMol2[0][2] -= h2H2xG;
                gradMol1[0][2] += h2H1xG;
                gradMol2[0][1] -= h2H1xG;
                
                final double oOG = gradEOO/distOO;
                final double oH1G = gradEOH1/distOH1;
                final double oH2G = gradEOH2/distOH2;
                
                final double oOxG = oOG*diffOOx;
                final double oH1xG = oH1G*diffOH1x;
                final double oH2xG = oH2G*diffOH2x;
                
                gradMol1[0][0] += oOxG;
                gradMol2[0][0] -= oOxG;
                gradMol1[0][0] += oH1xG;
                gradMol2[0][1] -= oH1xG;
                gradMol1[0][0] += oH2xG;
                gradMol2[0][2] -= oH2xG;
                
                
                final double h1OyG = h1OG*diffH1Oy;
                final double h1H1yG = h1H1G*diffH1H1y;
                final double h1H2yG = h1H2G*diffH1H2y;
                
                gradMol1[1][1] += h1OyG;
                gradMol2[1][0] -= h1OyG;
                gradMol1[1][1] += h1H1yG;
                gradMol2[1][1] -= h1H1yG;
                gradMol1[1][1] += h1H2yG;
                gradMol2[1][2] -= h1H2yG;
                
                final double h2OyG = h2OG*diffH2Oy;
                final double h2H2yG = h2H2G*diffH2H2y;
                final double h2H1yG = h2H1G*diffH2H1y;
                
                gradMol1[1][2] += h2OyG;
                gradMol2[1][0] -= h2OyG;
                gradMol1[1][2] += h2H2yG;
                gradMol2[1][2] -= h2H2yG;
                gradMol1[1][2] += h2H1yG;
                gradMol2[1][1] -= h2H1yG;
                
                final double oOyG = oOG*diffOOy;
                final double oH1yG = oH1G*diffOH1y;
                final double oH2yG = oH2G*diffOH2y;
                
                gradMol1[1][0] += oOyG;
                gradMol2[1][0] -= oOyG;
                gradMol1[1][0] += oH1yG;
                gradMol2[1][1] -= oH1yG;
                gradMol1[1][0] += oH2yG;
                gradMol2[1][2] -= oH2yG;
                
                final double h1OzG = h1OG*diffH1Oz;
                final double h1H1zG = h1H1G*diffH1H1z;
                final double h1H2zG = h1H2G*diffH1H2z;
                
                gradMol1[2][1] += h1OzG;
                gradMol2[2][0] -= h1OzG;
                gradMol1[2][1] += h1H1zG;
                gradMol2[2][1] -= h1H1zG;
                gradMol1[2][1] += h1H2zG;
                gradMol2[2][2] -= h1H2zG;
                
                final double h2OzG = h2OG*diffH2Oz;
                final double h2H2zG = h2H2G*diffH2H2z;
                final double h2H1zG = h2H1G*diffH2H1z;
                
                gradMol1[2][2] += h2OzG;
                gradMol2[2][0] -= h2OzG;
                gradMol1[2][2] += h2H2zG;
                gradMol2[2][2] -= h2H2zG;
                gradMol1[2][2] += h2H1zG;
                gradMol2[2][1] -= h2H1zG;
                
                final double oOzG = oOG*diffOOz;
                final double oH1zG = oH1G*diffOH1z;
                final double oH2zG = oH2G*diffOH2z;
                
                gradMol1[2][0] += oOzG;
                gradMol2[2][0] -= oOzG;
                gradMol1[2][0] += oH1zG;
                gradMol2[2][1] -= oH1zG;
                gradMol1[2][0] += oH2zG;
                gradMol2[2][2] -= oH2zG;
                
                // second: the LJ potential, we assume we CAN cut this off based on the COM-COM distance of the O's                
                if(distOOSq >= CUTOFFLJSQ){continue;}
                
                final double distOO6 = distOOSq*distOOSq*distOOSq;
                final double distOO12 = distOO6*distOO6;
                
                final double ljA = LJAOXYGEN/distOO12;
                final double ljB = LJBOXYGEN/distOO6;
                final double lj = ljA - ljB;
                energy += lj;
                
                if(DEBUG){
                    System.out.println("DEBUG: LJ contrib " + lj + " " + distOO*BOHRTOANG);
                }
                
                // LJ gradient terms
                final double distOOInv = 1.0/distOO;
                final double gradLJ = -12*ljA*distOOInv + 6*ljB*distOOInv;
                final double gradLJx = gradLJ*diffOOx*distOOInv;
                final double gradLJy = gradLJ*diffOOy*distOOInv;
                final double gradLJz = gradLJ*diffOOz*distOOInv;
                                
                gradMol1[0][0] += gradLJx;
                gradMol2[0][0] -= gradLJx;
                gradMol1[1][0] += gradLJy;
                gradMol2[1][0] -= gradLJy;
                gradMol1[2][0] += gradLJz;
                gradMol2[2][0] -= gradLJz;
            }
        }
        
        return energy;
    }
}
