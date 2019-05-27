/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2015, J. M. Dieterich
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
import org.apache.commons.math3.util.FastMath;
import static org.ogolem.core.Constants.ANGTOBOHR;
import static org.ogolem.core.Constants.BOHRTOANG;
import org.ogolem.math.GenericLookup;
import org.ogolem.math.TrivialLinearAlgebra;

/**
 * Provides the TIP4P force field for an simple description of water clusters.
 * Numerical data from http://www1.lsbu.ac.uk/water/water_models.html and wiki
 * @author Johannes Dieterich
 * @version 2015-01-08
 */
public class TIP4PForceField implements RigidBodyBackend {

    // the ID
    private static final long serialVersionUID = (long) 20140108;
    private static final boolean DEBUG = false;
    
    // the LJ cutoff in bohr
    private static final double CUTOFFLJSQ = 10*ANGTOBOHR*10*ANGTOBOHR;
    
    // partial charge on H's
    private final double CHARGEH;
    
    // partial charges on M's
    private final double CHARGEM;
    
    // distance of O to M
    private final double DISTOM;
    
    private final double CHARGEMM;
    private final double CHARGEMH;
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
    
    private final double[][] xyz1 = new double[3][4]; // account for the M-site!
    private final double[][] xyz2 = new double[3][4];

    public TIP4PForceField(){
        // use standard TIP4P parameters
        final TIPnPParameters.TIP4PParameters params = new TIPnPParameters.StandardTIP4PParameters();
        this.CHARGEH = params.getChargeH();
        this.CHARGEM = params.getChargeM();
        this.CHARGEHH = CHARGEH*CHARGEH;
        this.CHARGEMH = CHARGEH*CHARGEM;
        this.CHARGEMM = CHARGEM*CHARGEM;
        this.LJAOXYGEN = params.getLJA();
        this.LJBOXYGEN = params.getLJB();
        this.DISTOM = params.getOMDist();
        this.rOH = params.getOH();
        this.angHOH = params.getHOH();
        // create Bernd's anti cold fusion term
        this.anticold = new TIPnPHelpers.AntiColdFusionFunction();
        this.anticoldMed = new GenericLookup(anticold,GenericLookup.POSSIBLE_SIZES[1],0.0,ANTICOLDEND); // we want a very compact table
    }
    
    public TIP4PForceField(final TIPnPParameters.TIP4PParameters params){
        this.CHARGEH = params.getChargeH();
        this.CHARGEM = params.getChargeM();
        this.CHARGEHH = CHARGEH*CHARGEH;
        this.CHARGEMH = CHARGEH*CHARGEM;
        this.CHARGEMM = CHARGEM*CHARGEM;
        this.LJAOXYGEN = params.getLJA();
        this.LJBOXYGEN = params.getLJB();
        this.DISTOM = params.getOMDist();
        this.rOH = params.getOH();
        this.angHOH = params.getHOH();
        // create Bernd's anti cold fusion term
        this.anticold = new TIPnPHelpers.AntiColdFusionFunction();
        this.anticoldMed = new GenericLookup(anticold,GenericLookup.POSSIBLE_SIZES[1],0.0,ANTICOLDEND); // we want a very compact table
    }
    
    TIP4PForceField(final TIP4PForceField orig){
        this.CHARGEH = orig.CHARGEH;
        this.CHARGEHH = orig.CHARGEHH;
        this.CHARGEM = orig.CHARGEM;
        this.CHARGEMM = orig.CHARGEMM;
        this.CHARGEMH = orig.CHARGEMH;
        this.DISTOM = orig.DISTOM;
        this.LJAOXYGEN = orig.LJAOXYGEN;
        this.LJBOXYGEN = orig.LJBOXYGEN;
        this.rOH = orig.rOH;
        this.angHOH = orig.angHOH;
        // create Bernd's anti cold fusion term
        this.anticold = new TIPnPHelpers.AntiColdFusionFunction();
        this.anticoldMed = new GenericLookup(anticold,GenericLookup.POSSIBLE_SIZES[1],0.0,ANTICOLDEND); // we want a very compact table
    }
    
    @Override
    public TIP4PForceField clone(){
        return new TIP4PForceField(this);
    }

    @Override
    public String getMethodID(){
        return "TIP4P (dirty implementation)";
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
        if(atomNos[0] == 8
                && atomNos[1] == 1
                && atomNos[2] == 1){
        } else{
            System.err.println("ERROR: We figured out that this is no suitable water (in O/H/H order) cluster!");
            return false;
        }

        // now the rest
        for(int mol = 0; mol < ref.getNumberOfIndieParticles(); mol++){
            
            final short[] atomNosMol = ref.getMoleculeAtPosition(mol).getAtomNumbers();
            // check whether all molecules are having the right order
            if (atomNosMol[0] == 8 && atomNosMol[1] == 1 && atomNosMol[2] == 1) {
                // good
            } else {
                // something is f***** up
                System.err.println("ERROR: Some parts of your cluster aren't water (in O/H/H order).");
                return false;
            }
        }
        
        return true;
    }
    
    @Override
    public CartesianCoordinates adjustCartesians(final Geometry ref, final CartesianCoordinates unadjusted, final int molID){
        
        // here, we just add the M-site to the CartesianCoordinates as to avoid recalculating it all the time
        // everything else: does basically not matter
        final CartesianCoordinates adjCartes = new CartesianCoordinates(4,1,new int[]{4});
        final String[] atoms = adjCartes.getAllAtomTypes();
        atoms[0] = "O";
        atoms[1] = "H";
        atoms[2] = "H";
        atoms[3] = "XX";
        final short[] nos = new short[4];
        nos[0] = 8;
        nos[1] = 1;
        nos[2] = 1;
        nos[3] = 0;
        adjCartes.setAtomNumbers(nos);
        final double[][] xyzAdj = adjCartes.getAllXYZCoord();
        final double[][] xyzOrg = unadjusted.getAllXYZCoord();
        System.arraycopy(xyzOrg[0], 0, xyzAdj[0], 0, 3);
        System.arraycopy(xyzOrg[1], 0, xyzAdj[1], 0, 3);
        System.arraycopy(xyzOrg[2], 0, xyzAdj[2], 0, 3);
        final double[] scr1 = new double[3];
        final double[] scr2 = new double[3];
        final double[] mPos = new double[3];
        calcMPosition(xyzOrg, scr1, scr2, mPos);
        xyzAdj[0][3] = mPos[0];
        xyzAdj[1][3] = mPos[1];
        xyzAdj[2][3] = mPos[2];
        
        return adjCartes;
    }
    
    @Override
    public void rigidify(final Molecule mol){
        TIPnPHelpers.sanitizeWater(mol, rOH, angHOH);
    }

    @Override
    public double energy(final Geometry ref, final List<CartesianCoordinates> cartes,
            final double[][] coms, final int counter){
                
        double energy = 0.0;
        for(int mol1 = 0; mol1 < cartes.size()-1; mol1++){
            
            // adjust the first water molecules coordinates
            final CartesianCoordinates c1 = cartes.get(mol1);
            final double[] com1 = coms[mol1];
            final double[][] xyz1O = c1.getAllXYZCoord();
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 4; j++){
                    xyz1[i][j] = xyz1O[i][j] + com1[i];
                }
            }
            
            if(DEBUG){
                final boolean ok = checkCartesForMPosition(c1,DISTOM);
                if(!ok){
                    throw new RuntimeException("Something is off with the M position in molecule " + mol1);
                }
            }

            
            if(DEBUG){
                System.out.println("No " + mol1);
                System.out.println("" + 4);
                System.out.println("");
                System.out.println("O " + xyz1[0][0]*BOHRTOANG + " "  + xyz1[1][0]*BOHRTOANG + " "  + xyz1[2][0]*BOHRTOANG);
                System.out.println("H " + xyz1[0][1]*BOHRTOANG + " "  + xyz1[1][1]*BOHRTOANG + " "  + xyz1[2][1]*BOHRTOANG);
                System.out.println("H " + xyz1[0][2]*BOHRTOANG + " "  + xyz1[1][2]*BOHRTOANG + " "  + xyz1[2][2]*BOHRTOANG);
                System.out.println("XX " + xyz1[0][3]*BOHRTOANG + " "  + xyz1[1][3]*BOHRTOANG + " "  + xyz1[2][3]*BOHRTOANG);
            }
            
            for(int mol2 = mol1+1; mol2 < cartes.size(); mol2++){
                
                // adjust the second water molecules coordinates
                final double[] com2 = coms[mol2];
                final CartesianCoordinates c2 = cartes.get(mol2);
                final double[][] xyz2O = c2.getAllXYZCoord();
                for(int i = 0; i < 3; i++){
                    for(int j = 0; j < 4; j++){
                        xyz2[i][j] = xyz2O[i][j] + com2[i];
                    }
                }
                
                if(DEBUG){
                    final boolean ok = checkCartesForMPosition(c2,DISTOM);
                    if(!ok){
                        throw new RuntimeException("Something is off with the M position in molecule " + mol2);
                    }
                }
                                
                if(DEBUG){
                    System.out.println("No " + mol2);
                    System.out.println("" + 4);
                    System.out.println("");
                    System.out.println("O " + xyz2[0][0]*BOHRTOANG + " "  + xyz2[1][0]*BOHRTOANG + " "  + xyz2[2][0]*BOHRTOANG);
                    System.out.println("H " + xyz2[0][1]*BOHRTOANG + " "  + xyz2[1][1]*BOHRTOANG + " "  + xyz2[2][1]*BOHRTOANG);
                    System.out.println("H " + xyz2[0][2]*BOHRTOANG + " "  + xyz2[1][2]*BOHRTOANG + " "  + xyz2[2][2]*BOHRTOANG);
                    System.out.println("XX " + xyz2[0][3]*BOHRTOANG + " "  + xyz2[1][3]*BOHRTOANG + " "  + xyz2[2][3]*BOHRTOANG);
                }
                                
                // XXX cutoff for the electrostatics
                // for very large clusters, it may pay off to just compute the diff between the COMs and based on that cut off at e.g. 40 Angstrom
                
                // get some distance diffs
                final double diffMMx = xyz1[0][3] - xyz2[0][3];
                final double diffMMy = xyz1[1][3] - xyz2[1][3];
                final double diffMMz = xyz1[2][3] - xyz2[2][3];
                final double distMMSq = diffMMx*diffMMx + diffMMy*diffMMy + diffMMz*diffMMz;
                final double distMM = FastMath.sqrt(distMMSq);
                if(distMM <= ANTICOLDEND){energy += anticoldMed.inter(distMM);}
                
                final double diffOOx = xyz1[0][0] - xyz2[0][0];
                final double diffOOy = xyz1[1][0] - xyz2[1][0];
                final double diffOOz = xyz1[2][0] - xyz2[2][0];
                final double distOOSq = diffOOx*diffOOx + diffOOy*diffOOy + diffOOz*diffOOz;
                final double distOO = FastMath.sqrt(distOOSq);
                if(distOO <= ANTICOLDEND){energy += anticoldMed.inter(distOO);}
                
                final double diffH1H1x = xyz1[0][1] - xyz2[0][1];
                final double diffH1H1y = xyz1[1][1] - xyz2[1][1];
                final double diffH1H1z = xyz1[2][1] - xyz2[2][1];
                final double distH1H1Sq = diffH1H1x*diffH1H1x + diffH1H1y*diffH1H1y + diffH1H1z*diffH1H1z;
                final double distH1H1 = FastMath.sqrt(distH1H1Sq);
                if(distH1H1 <= ANTICOLDEND){energy += anticoldMed.inter(distH1H1);}
                
                final double diffH2H2x = xyz1[0][2] - xyz2[0][2];
                final double diffH2H2y = xyz1[1][2] - xyz2[1][2];
                final double diffH2H2z = xyz1[2][2] - xyz2[2][2];
                final double distH2H2Sq = diffH2H2x*diffH2H2x + diffH2H2y*diffH2H2y + diffH2H2z*diffH2H2z;
                final double distH2H2 = FastMath.sqrt(distH2H2Sq);
                if(distH2H2 <= ANTICOLDEND){energy += anticoldMed.inter(distH2H2);}
                
                final double diffH1H2x = xyz1[0][1] - xyz2[0][2];
                final double diffH1H2y = xyz1[1][1] - xyz2[1][2];
                final double diffH1H2z = xyz1[2][1] - xyz2[2][2];
                final double distH1H2Sq = diffH1H2x*diffH1H2x + diffH1H2y*diffH1H2y + diffH1H2z*diffH1H2z;
                final double distH1H2 = FastMath.sqrt(distH1H2Sq);
                if(distH1H2 <= ANTICOLDEND){energy += anticoldMed.inter(distH1H2);}
                
                final double diffH2H1x = xyz1[0][2] - xyz2[0][1];
                final double diffH2H1y = xyz1[1][2] - xyz2[1][1];
                final double diffH2H1z = xyz1[2][2] - xyz2[2][1];
                final double distH2H1Sq = diffH2H1x*diffH2H1x + diffH2H1y*diffH2H1y + diffH2H1z*diffH2H1z;
                final double distH2H1 = FastMath.sqrt(distH2H1Sq);
                if(distH2H1 <= ANTICOLDEND){energy += anticoldMed.inter(distH2H1);}
                
                final double diffMH1x = xyz1[0][3] - xyz2[0][1];
                final double diffMH1y = xyz1[1][3] - xyz2[1][1];
                final double diffMH1z = xyz1[2][3] - xyz2[2][1];
                final double distMH1Sq = diffMH1x*diffMH1x + diffMH1y*diffMH1y + diffMH1z*diffMH1z;
                final double distMH1 = FastMath.sqrt(distMH1Sq);
                if(distMH1 <= ANTICOLDEND){energy += anticoldMed.inter(distMH1);}
                
                final double diffMH2x = xyz1[0][3] - xyz2[0][2];
                final double diffMH2y = xyz1[1][3] - xyz2[1][2];
                final double diffMH2z = xyz1[2][3] - xyz2[2][2];
                final double distMH2Sq = diffMH2x*diffMH2x + diffMH2y*diffMH2y + diffMH2z*diffMH2z;
                final double distMH2 = FastMath.sqrt(distMH2Sq);
                if(distMH2 <= ANTICOLDEND){energy += anticoldMed.inter(distMH2);}
                
                final double diffH1Mx = xyz1[0][1] - xyz2[0][3];
                final double diffH1My = xyz1[1][1] - xyz2[1][3];
                final double diffH1Mz = xyz1[2][1] - xyz2[2][3];
                final double distH1MSq = diffH1Mx*diffH1Mx + diffH1My*diffH1My + diffH1Mz*diffH1Mz;
                final double distH1M = FastMath.sqrt(distH1MSq);
                if(distH1M <= ANTICOLDEND){energy += anticoldMed.inter(distH1M);}
                
                final double diffH2Mx = xyz1[0][2] - xyz2[0][3];
                final double diffH2My = xyz1[1][2] - xyz2[1][3];
                final double diffH2Mz = xyz1[2][2] - xyz2[2][3];
                final double distH2MSq = diffH2Mx*diffH2Mx + diffH2My*diffH2My + diffH2Mz*diffH2Mz;
                final double distH2M = FastMath.sqrt(distH2MSq);
                if(distH2M <= ANTICOLDEND){energy += anticoldMed.inter(distH2M);}
                
                
                // now compute the actual interaction energies
                // first: the electrostatics, we assume that we never can cut that one off
                final double elect = (CHARGEMM)/ distMM + (CHARGEMH)/ distH1M + (CHARGEMH)/ distH2M
                        + (CHARGEMH)/ distMH1 + (CHARGEMH)/ distMH2 + (CHARGEHH)/ distH1H1
                        + (CHARGEHH)/ distH2H2 + + (CHARGEHH)/ distH1H2 + (CHARGEHH)/ distH2H1;
                
                if(DEBUG){
                    System.out.println("DEBUG: distMM " + (CHARGEMM)/ distMM + " " + distMM*BOHRTOANG);
                    System.out.println("DEBUG: distH1M " + (CHARGEMH)/ distH1M + " " + distH1M*BOHRTOANG);
                    System.out.println("DEBUG: distH2M " + (CHARGEMH)/ distH2M + " " + distH2M*BOHRTOANG);
                    System.out.println("DEBUG: distMH1 " + (CHARGEMH)/ distMH1 + " " + distMH1*BOHRTOANG);
                    System.out.println("DEBUG: distMH2 " + (CHARGEMH)/ distMH2 + " " + distMH2*BOHRTOANG);
                    System.out.println("DEBUG: distH1H1 " + (CHARGEHH)/ distH1H1 + " " + distH1H1*BOHRTOANG);
                    System.out.println("DEBUG: distH2H2 " + (CHARGEHH)/ distH2H2 + " " + distH2H2*BOHRTOANG);
                    System.out.println("DEBUG: distH1H2 " + (CHARGEHH)/ distH1H2 + " " + distH1H2*BOHRTOANG);
                    System.out.println("DEBUG: distH2H1 " + (CHARGEHH)/ distH2H1 + " " + distH2H1*BOHRTOANG);
                    System.out.println("DEBUG: elec " + elect);
                }
                energy += elect;
                
                // second: the LJ potential, we assume we CAN cut this off based on the COM-COM distance of the O's                
                if(distOOSq >= CUTOFFLJSQ){continue;}
                
                final double distOO6 = distOOSq*distOOSq*distOOSq;
                final double distOO12 = distOO6*distOO6;
                
                final double lj = LJAOXYGEN/distOO12 - LJBOXYGEN/distOO6;
                energy += lj;
                
                if(DEBUG){
                    System.out.println("DEBUG: LJ contrib " + lj + " " + distOO*BOHRTOANG);
                }
            }
        }
        
        if(DEBUG){
            System.out.println("DEBUG: At iter " + counter + " energy is " + energy);
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
        for(int mol1 = 0; mol1 < cartes.size()-1; mol1++){
            
            // adjust the first water molecules coordinates
            final CartesianCoordinates c1 = cartes.get(mol1);
            final double[] com1 = coms[mol1];
            final double[][] xyz1O = c1.getAllXYZCoord();
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 4; j++){
                    xyz1[i][j] = xyz1O[i][j] + com1[i];
                }
            }
            
            if(DEBUG){
                final boolean ok = checkCartesForMPosition(c1,DISTOM);
                if(!ok){
                    throw new RuntimeException("Something is off with the M position in molecule " + mol1);
                }
            }
            
            final double[][] gradMol1 = gradient.get(mol1);
            
            if(DEBUG){
                System.out.println("No " + mol1);
                System.out.println("" + 4);
                System.out.println("");
                System.out.println("O " + xyz1[0][0]*BOHRTOANG + " "  + xyz1[1][0]*BOHRTOANG + " "  + xyz1[2][0]*BOHRTOANG);
                System.out.println("H " + xyz1[0][1]*BOHRTOANG + " "  + xyz1[1][1]*BOHRTOANG + " "  + xyz1[2][1]*BOHRTOANG);
                System.out.println("H " + xyz1[0][2]*BOHRTOANG + " "  + xyz1[1][2]*BOHRTOANG + " "  + xyz1[2][2]*BOHRTOANG);
                System.out.println("XX " + xyz1[0][3]*BOHRTOANG + " "  + xyz1[1][3]*BOHRTOANG + " "  + xyz1[2][3]*BOHRTOANG);
            }
            
            for(int mol2 = mol1+1; mol2 < cartes.size(); mol2++){
                
                // adjust the second water molecules coordinates
                final double[] com2 = coms[mol2];
                final CartesianCoordinates c2 = cartes.get(mol2);
                final double[][] xyz2O = c2.getAllXYZCoord();
                for(int i = 0; i < 3; i++){
                    for(int j = 0; j < 4; j++){
                        xyz2[i][j] = xyz2O[i][j] + com2[i];
                    }
                }
                
                if(DEBUG){
                    final boolean ok = checkCartesForMPosition(c2,DISTOM);
                    if(!ok){
                        throw new RuntimeException("Something is off with the M position in molecule " + mol2);
                    }
                }
                
                final double[][] gradMol2 = gradient.get(mol2);
                
                if(DEBUG){
                    System.out.println("No " + mol2);
                    System.out.println("" + 4);
                    System.out.println("");
                    System.out.println("O " + xyz2[0][0]*BOHRTOANG + " "  + xyz2[1][0]*BOHRTOANG + " "  + xyz2[2][0]*BOHRTOANG);
                    System.out.println("H " + xyz2[0][1]*BOHRTOANG + " "  + xyz2[1][1]*BOHRTOANG + " "  + xyz2[2][1]*BOHRTOANG);
                    System.out.println("H " + xyz2[0][2]*BOHRTOANG + " "  + xyz2[1][2]*BOHRTOANG + " "  + xyz2[2][2]*BOHRTOANG);
                    System.out.println("XX " + xyz2[0][3]*BOHRTOANG + " "  + xyz2[1][3]*BOHRTOANG + " "  + xyz2[2][3]*BOHRTOANG);
                }
                                
                // XXX cutoff for the electrostatics
                // for very large clusters, it may pay off to just compute the diff between the COMs and based on that cut off at e.g. 40 Angstrom
                                
                // get some distance diffs
                final double diffMMx = xyz1[0][3] - xyz2[0][3];
                final double diffMMy = xyz1[1][3] - xyz2[1][3];
                final double diffMMz = xyz1[2][3] - xyz2[2][3];
                final double distMMSq = diffMMx*diffMMx + diffMMy*diffMMy + diffMMz*diffMMz;
                final double distMM = FastMath.sqrt(distMMSq);
                if(distMM <= ANTICOLDEND){energy += anticoldMed.inter(distMM);}
                
                final double diffOOx = xyz1[0][0] - xyz2[0][0];
                final double diffOOy = xyz1[1][0] - xyz2[1][0];
                final double diffOOz = xyz1[2][0] - xyz2[2][0];
                final double distOOSq = diffOOx*diffOOx + diffOOy*diffOOy + diffOOz*diffOOz;
                final double distOO = FastMath.sqrt(distOOSq);
                if(distOO <= ANTICOLDEND){energy += anticoldMed.inter(distOO);}
                
                final double diffH1H1x = xyz1[0][1] - xyz2[0][1];
                final double diffH1H1y = xyz1[1][1] - xyz2[1][1];
                final double diffH1H1z = xyz1[2][1] - xyz2[2][1];
                final double distH1H1Sq = diffH1H1x*diffH1H1x + diffH1H1y*diffH1H1y + diffH1H1z*diffH1H1z;
                final double distH1H1 = FastMath.sqrt(distH1H1Sq);
                if(distH1H1 <= ANTICOLDEND){energy += anticoldMed.inter(distH1H1);}
                
                final double diffH2H2x = xyz1[0][2] - xyz2[0][2];
                final double diffH2H2y = xyz1[1][2] - xyz2[1][2];
                final double diffH2H2z = xyz1[2][2] - xyz2[2][2];
                final double distH2H2Sq = diffH2H2x*diffH2H2x + diffH2H2y*diffH2H2y + diffH2H2z*diffH2H2z;
                final double distH2H2 = FastMath.sqrt(distH2H2Sq);
                if(distH2H2 <= ANTICOLDEND){energy += anticoldMed.inter(distH2H2);}
                
                final double diffH1H2x = xyz1[0][1] - xyz2[0][2];
                final double diffH1H2y = xyz1[1][1] - xyz2[1][2];
                final double diffH1H2z = xyz1[2][1] - xyz2[2][2];
                final double distH1H2Sq = diffH1H2x*diffH1H2x + diffH1H2y*diffH1H2y + diffH1H2z*diffH1H2z;
                final double distH1H2 = FastMath.sqrt(distH1H2Sq);
                if(distH1H2 <= ANTICOLDEND){energy += anticoldMed.inter(distH1H2);}
                
                final double diffH2H1x = xyz1[0][2] - xyz2[0][1];
                final double diffH2H1y = xyz1[1][2] - xyz2[1][1];
                final double diffH2H1z = xyz1[2][2] - xyz2[2][1];
                final double distH2H1Sq = diffH2H1x*diffH2H1x + diffH2H1y*diffH2H1y + diffH2H1z*diffH2H1z;
                final double distH2H1 = FastMath.sqrt(distH2H1Sq);
                if(distH2H1 <= ANTICOLDEND){energy += anticoldMed.inter(distH2H1);}
                
                final double diffMH1x = xyz1[0][3] - xyz2[0][1];
                final double diffMH1y = xyz1[1][3] - xyz2[1][1];
                final double diffMH1z = xyz1[2][3] - xyz2[2][1];
                final double distMH1Sq = diffMH1x*diffMH1x + diffMH1y*diffMH1y + diffMH1z*diffMH1z;
                final double distMH1 = FastMath.sqrt(distMH1Sq);
                if(distMH1 <= ANTICOLDEND){energy += anticoldMed.inter(distMH1);}
                
                final double diffMH2x = xyz1[0][3] - xyz2[0][2];
                final double diffMH2y = xyz1[1][3] - xyz2[1][2];
                final double diffMH2z = xyz1[2][3] - xyz2[2][2];
                final double distMH2Sq = diffMH2x*diffMH2x + diffMH2y*diffMH2y + diffMH2z*diffMH2z;
                final double distMH2 = FastMath.sqrt(distMH2Sq);
                if(distMH2 <= ANTICOLDEND){energy += anticoldMed.inter(distMH2);}
                
                final double diffH1Mx = xyz1[0][1] - xyz2[0][3];
                final double diffH1My = xyz1[1][1] - xyz2[1][3];
                final double diffH1Mz = xyz1[2][1] - xyz2[2][3];
                final double distH1MSq = diffH1Mx*diffH1Mx + diffH1My*diffH1My + diffH1Mz*diffH1Mz;
                final double distH1M = FastMath.sqrt(distH1MSq);
                if(distH1M <= ANTICOLDEND){energy += anticoldMed.inter(distH1M);}
                
                final double diffH2Mx = xyz1[0][2] - xyz2[0][3];
                final double diffH2My = xyz1[1][2] - xyz2[1][3];
                final double diffH2Mz = xyz1[2][2] - xyz2[2][3];
                final double distH2MSq = diffH2Mx*diffH2Mx + diffH2My*diffH2My + diffH2Mz*diffH2Mz;
                final double distH2M = FastMath.sqrt(distH2MSq);
                if(distH2M <= ANTICOLDEND){energy += anticoldMed.inter(distH2M);}
                
                
                // now compute the actual interaction energies
                // first: the electrostatics, we assume that we never can cut that one off
                final double elect = (CHARGEMM)/ distMM + (CHARGEMH)/ distH1M + (CHARGEMH)/ distH2M
                        + (CHARGEMH)/ distMH1 + (CHARGEMH)/ distMH2 + (CHARGEHH)/ distH1H1
                        + (CHARGEHH)/ distH2H2 + (CHARGEHH)/ distH1H2 + (CHARGEHH)/ distH2H1;
                
                                if(DEBUG){
                    System.out.println("DEBUG: distMM " + (CHARGEMM)/ distMM + " " + distMM*BOHRTOANG);
                    System.out.println("DEBUG: distH1M " + (CHARGEMH)/ distH1M + " " + distH1M*BOHRTOANG);
                    System.out.println("DEBUG: distH2M " + (CHARGEMH)/ distH2M + " " + distH2M*BOHRTOANG);
                    System.out.println("DEBUG: distMH1 " + (CHARGEMH)/ distMH1 + " " + distMH1*BOHRTOANG);
                    System.out.println("DEBUG: distMH2 " + (CHARGEMH)/ distMH2 + " " + distMH2*BOHRTOANG);
                    System.out.println("DEBUG: distH1H1 " + (CHARGEHH)/ distH1H1 + " " + distH1H1*BOHRTOANG);
                    System.out.println("DEBUG: distH2H2 " + (CHARGEHH)/ distH2H2 + " " + distH2H2*BOHRTOANG);
                    System.out.println("DEBUG: distH1H2 " + (CHARGEHH)/ distH1H2 + " " + distH1H2*BOHRTOANG);
                    System.out.println("DEBUG: distH2H1 " + (CHARGEHH)/ distH2H1 + " " + distH2H1*BOHRTOANG);
                    System.out.println("DEBUG: elec " + elect);
                }
                energy += elect;
                
                // gradient for the electrostatics                
                final double gradEMM   = -(CHARGEMM)/ distMMSq;
                final double gradEH1M  = -(CHARGEMH)/ distH1MSq;
                final double gradEH2M  = -(CHARGEMH)/ distH2MSq;
                final double gradEMH1  = -(CHARGEMH)/ distMH1Sq;
                final double gradEMH2  = -(CHARGEMH)/ distMH2Sq;
                final double gradEH1H1 = -(CHARGEHH)/ distH1H1Sq;
                final double gradEH2H2 = -(CHARGEHH)/ distH2H2Sq;
                final double gradEH1H2 = -(CHARGEHH)/ distH1H2Sq;
                final double gradEH2H1 = -(CHARGEHH)/ distH2H1Sq;
                
                // divide the distances in all three dimensions by the total distance to cover for the coordinate system
                                
                final double h1MG = gradEH1M/distH1M;
                final double h1H1G = gradEH1H1/distH1H1;
                final double h1H2G = gradEH1H2/distH1H2;
                
                final double h1MxG = h1MG*diffH1Mx;
                final double h1H1xG = h1H1G*diffH1H1x;
                final double h1H2xG = h1H2G*diffH1H2x;
                
                gradMol1[0][1] += h1MxG;
                gradMol2[0][3] -= h1MxG;
                gradMol1[0][1] += h1H1xG;
                gradMol2[0][1] -= h1H1xG;
                gradMol1[0][1] += h1H2xG;
                gradMol2[0][2] -= h1H2xG;
                
                final double h2MG = gradEH2M/distH2M;
                final double h2H2G = gradEH2H2/distH2H2;
                final double h2H1G = gradEH2H1/distH2H1;
                
                final double h2MxG = h2MG*diffH2Mx;
                final double h2H2xG = h2H2G*diffH2H2x;
                final double h2H1xG = h2H1G*diffH2H1x;
                
                gradMol1[0][2] += h2MxG;
                gradMol2[0][3] -= h2MxG;
                gradMol1[0][2] += h2H2xG;
                gradMol2[0][2] -= h2H2xG;
                gradMol1[0][2] += h2H1xG;
                gradMol2[0][1] -= h2H1xG;
                
                final double mMG = gradEMM/distMM;
                final double mH1G = gradEMH1/distMH1;
                final double mH2G = gradEMH2/distMH2;
                
                final double mMxG = mMG*diffMMx;
                final double mH1xG = mH1G*diffMH1x;
                final double mH2xG = mH2G*diffMH2x;
                
                gradMol1[0][3] += mMxG;
                gradMol2[0][3] -= mMxG;
                gradMol1[0][3] += mH1xG;
                gradMol2[0][1] -= mH1xG;
                gradMol1[0][3] += mH2xG;
                gradMol2[0][2] -= mH2xG;
                
                
                final double h1MyG = h1MG*diffH1My;
                final double h1H1yG = h1H1G*diffH1H1y;
                final double h1H2yG = h1H2G*diffH1H2y;
                
                gradMol1[1][1] += h1MyG;
                gradMol2[1][3] -= h1MyG;
                gradMol1[1][1] += h1H1yG;
                gradMol2[1][1] -= h1H1yG;
                gradMol1[1][1] += h1H2yG;
                gradMol2[1][2] -= h1H2yG;
                
                final double h2MyG = h2MG*diffH2My;
                final double h2H2yG = h2H2G*diffH2H2y;
                final double h2H1yG = h2H1G*diffH2H1y;
                
                gradMol1[1][2] += h2MyG;
                gradMol2[1][3] -= h2MyG;
                gradMol1[1][2] += h2H2yG;
                gradMol2[1][2] -= h2H2yG;
                gradMol1[1][2] += h2H1yG;
                gradMol2[1][1] -= h2H1yG;
                
                final double mMyG = mMG*diffMMy;
                final double mH1yG = mH1G*diffMH1y;
                final double mH2yG = mH2G*diffMH2y;
                
                gradMol1[1][3] += mMyG;
                gradMol2[1][3] -= mMyG;
                gradMol1[1][3] += mH1yG;
                gradMol2[1][1] -= mH1yG;
                gradMol1[1][3] += mH2yG;
                gradMol2[1][2] -= mH2yG;
                
                final double h1MzG = h1MG*diffH1Mz;
                final double h1H1zG = h1H1G*diffH1H1z;
                final double h1H2zG = h1H2G*diffH1H2z;
                
                gradMol1[2][1] += h1MzG;
                gradMol2[2][3] -= h1MzG;
                gradMol1[2][1] += h1H1zG;
                gradMol2[2][1] -= h1H1zG;
                gradMol1[2][1] += h1H2zG;
                gradMol2[2][2] -= h1H2zG;
                
                final double h2MzG = h2MG*diffH2Mz;
                final double h2H2zG = h2H2G*diffH2H2z;
                final double h2H1zG = h2H1G*diffH2H1z;
                
                gradMol1[2][2] += h2MzG;
                gradMol2[2][3] -= h2MzG;
                gradMol1[2][2] += h2H2zG;
                gradMol2[2][2] -= h2H2zG;
                gradMol1[2][2] += h2H1zG;
                gradMol2[2][1] -= h2H1zG;
                
                final double mMzG = mMG*diffMMz;
                final double mH1zG = mH1G*diffMH1z;
                final double mH2zG = mH2G*diffMH2z;
                
                gradMol1[2][3] += mMzG;
                gradMol2[2][3] -= mMzG;
                gradMol1[2][3] += mH1zG;
                gradMol2[2][1] -= mH1zG;
                gradMol1[2][3] += mH2zG;
                gradMol2[2][2] -= mH2zG;
                
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
        
        if(DEBUG){
            System.out.println("DEBUG: At iter " + counter + " energy in gradient is " + energy);
        }
        
        return energy;
    }
    
    private void calcMPosition(final double[][] xyz, final double[] scr1, final double[] scr2, final double[] m){
        
        // calculate O-H1 and O-H2 vectors
        scr1[0] = xyz[0][1] - xyz[0][0];
        scr1[1] = xyz[1][1] - xyz[1][0];
        scr1[2] = xyz[2][1] - xyz[2][0];
        
        scr2[0] = xyz[0][2] - xyz[0][0];
        scr2[1] = xyz[1][2] - xyz[1][0];
        scr2[2] = xyz[2][2] - xyz[2][0];
        
        final double normoh1 = 1./FastMath.sqrt(scr1[0]*scr1[0] + scr1[1]*scr1[1] + scr1[2]*scr1[2]);
        final double normoh2 = 1./FastMath.sqrt(scr2[0]*scr2[0] + scr2[1]*scr2[1] + scr2[2]*scr2[2]);
        
        // norm the vectors
        scr1[0] *= normoh1;
        scr1[1] *= normoh1;
        scr1[2] *= normoh1;
        
        scr2[0] *= normoh2;
        scr2[1] *= normoh2;
        scr2[2] *= normoh2;
        
        // calculate M position as the angle halfing vector, multiplied by the distance to the O and from the O position
        m[0] = (scr1[0]+scr2[0]);
        m[1] = (scr1[1]+scr2[1]);
        m[2] = (scr1[2]+scr2[2]);
        
        final double normM = DISTOM/FastMath.sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
        
        m[0] = normM*m[0] + xyz[0][0];
        m[1] = normM*m[1] + xyz[1][0];
        m[2] = normM*m[2] + xyz[2][0];
    }
    
        
    private static boolean checkCartesForMPosition(final CartesianCoordinates c, final double distMO){
        
        // checks if the M position is still right (as in: is it the correct distance away from the O and is it in 
        // the plane with the Hs and the O
        
        // we DO NOT check ATM the angles H1OM and MOH2
        final double MYNUMACC = 1e-6;
        
        final double[][] xyz = c.getAllXYZCoord();
        final double distOMx = xyz[0][0]-xyz[0][3];
        final double distOMy = xyz[1][0]-xyz[1][3];
        final double distOMz = xyz[2][0]-xyz[2][3];
        final double thisDistOM = FastMath.sqrt(distOMx*distOMx + distOMy*distOMy + distOMz*distOMz);
        if(Math.abs(distMO-thisDistOM) > MYNUMACC){
            assert(Math.abs(distMO-thisDistOM) > MYNUMACC);
            return false;
        }
        
        // calc the normal of the OHH plane
        final double[] scr1 = new double[3];
        final double[] scr2 = new double[3];
        
        for(int i = 0; i < 3; i++){
            scr1[i] = xyz[i][1] - xyz[i][0];
            scr2[i] = xyz[i][2] - xyz[i][0];
        }
        
        final double[] norm = new double[3];
        TrivialLinearAlgebra.crossProduct(scr1, scr2, norm);
        
        // calc the MO diff vector
        for(int i = 0; i < 3; i++){
            scr1[i] = xyz[i][3] - xyz[i][0];
        }
        
        // calc the dotproduct
        final double dot = TrivialLinearAlgebra.dotProduct(scr1, norm);
        
        // if the dot product is below MYNUMACC, all is well
        if(Math.abs(dot) > MYNUMACC){
            assert(Math.abs(dot) > MYNUMACC);
            return false;
        }
        
        if(DEBUG){
            System.out.println("DEBUG: all is well in M positioning.");
        }
        
        return true;
    }

}
