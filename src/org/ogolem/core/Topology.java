/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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
import java.util.ArrayList;
import java.util.List;
import org.ogolem.math.AcosLookup;

/**
 * A value object representing a topology. Not everything needs to be filled in
 * and it is highly recommended to extend this when needed.
 * @author Johannes Dieterich
 * @version 2016-02-24
 */
public class Topology implements Cloneable, Serializable {

    private static final long serialVersionUID = (long) 20160224;
    private BondInfo bonds;
    private String[] atomNames;
    private short[] atomicNumbers;
    private float[] charges;
    private short[] spins;
    private double[][] positions;
    // calculated on-the-fly
    private double[][] distances;
    private double[][][] angles;
    private double[][][] bondedAngles;
    private double[][][][] dihedrals;
    // calculated on-the-fly (upon request)
    private boolean calcDiffs = false;
    private double[][][] diffPos;
    // just if applicable
    private List<int[]> contributions13;
    private int[][] contri13;
    private List<int[]> contributions14;
    private int[][] contri14;
    
    // important for computations
    private final int wAcos;
    private final AcosLookup acosLook;
    
    // some internal scratch space
    private final double[] scr1 = new double[3];
    private final double[] scr2 = new double[3];
    private final double[] scr3 = new double[3];

    public Topology(CartesianCoordinates cartes, BondInfo bonds){
        this(cartes,bonds,false);
    }
    
    public Topology(CartesianCoordinates cartes, BondInfo bonds, int wAcos){
        this(cartes,bonds,false, wAcos);
    }
    
    public Topology(CartesianCoordinates cartes, BondInfo bonds, boolean calcAllDiffs){
        this(cartes,bonds,false,0);
    }
        
    public Topology(CartesianCoordinates cartes, BondInfo bonds, boolean calcAllDiffs, int whichAcos){
        this.calcDiffs = calcAllDiffs;
        this.bonds = bonds;
        this.atomNames = cartes.getAllAtomTypes();
        this.atomicNumbers = cartes.getAllAtomNumbers();
        this.charges = cartes.getAllCharges();
        this.spins = cartes.getAllSpins();
        
        /*
         * our lookup table (potential overhead!)
         */
        switch(whichAcos){
            case 0: acosLook = null; wAcos = 0; break;
            case 1: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[0]); wAcos = 1; break;
            case 2: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[1]); wAcos = 1; break;
            case 3: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[2]); wAcos = 1; break;
            case 4: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[3]); wAcos = 1; break;
            case 5: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[4]); wAcos = 1; break;
            case 11: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[0]); wAcos = 2; break;
            case 12: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[1]); wAcos = 2; break;
            case 13: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[2]); wAcos = 2; break;
            case 14: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[3]); wAcos = 2; break;
            case 15: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[4]); wAcos = 2; break;
            case 20: acosLook = null; wAcos = 3; break;
            default: acosLook = null; wAcos = 0; break;
       }
        

        // calculate distances
        final int iNoOfAtoms = cartes.getNoOfAtoms();
        final double[][] daXYZ = cartes.getAllXYZCoord();

        this.positions = daXYZ;

        this.distances = new double[iNoOfAtoms][iNoOfAtoms];

        if(calcDiffs) {this.diffPos = new double[iNoOfAtoms][iNoOfAtoms][3];}
        if(calcDiffs){
            for(int i = 0; i < iNoOfAtoms -1; i++){
                for(int j = i+1; j < iNoOfAtoms; j++){
                    final double dX = daXYZ[0][i] - daXYZ[0][j];
                    final double dY = daXYZ[1][i] - daXYZ[1][j];
                    final double dZ = daXYZ[2][i] - daXYZ[2][j];
                    
                    distances[i][j] = Math.sqrt(dX*dX + dY*dY + dZ*dZ);
                    diffPos[i][j][0] = dX;
                    diffPos[i][j][1] = dY;
                    diffPos[i][j][2] = dZ;
                    
                    // it's kind of symmetric
                    diffPos[j][i][0] = -dX;
                    diffPos[j][i][1] = -dY;
                    diffPos[j][i][2] = -dZ;
                    distances[j][i] = distances[i][j];
                }
            }
        } else{
            for(int i = 0; i < iNoOfAtoms -1; i++){
                for(int j = i+1; j < iNoOfAtoms; j++){
                    distances[i][j] = CoordTranslation.distance(daXYZ, i, j);
                    
                    // it's symmetric
                    distances[j][i] = distances[i][j];
                }
            }
        }
    }
    
    public Topology(String[] atoms, double[][] xyz, BondInfo bonds,
            float[] charges, short[] spins, short[] numbers,
            List<int[]> contr13, List<int[]> contr14){
        this(atoms, xyz, bonds, charges, spins, numbers, contr13, contr14, false);
    }
    
    public Topology(String[] atoms, double[][] xyz, BondInfo bonds,
            float[] charges, short[] spins, short[] numbers,
            List<int[]> contr13, List<int[]> contr14, int whichAcos){
        this(atoms, xyz, bonds, charges, spins, numbers, contr13, contr14, false, whichAcos);
    }

    public Topology(String[] atoms, double[][] xyz, BondInfo bonds,
            float[] charges, short[] spins, short[] numbers,
            List<int[]> contr13, List<int[]> contr14,
            boolean calcAllDiffs){
        this(atoms, xyz, bonds, charges, spins, numbers, contr13, contr14, false, 0);
    }
        
    public Topology(String[] atoms, double[][] xyz, BondInfo bonds,
            float[] charges, short[] spins, short[] numbers,
            List<int[]> contr13, List<int[]> contr14,
            boolean calcAllDiffs, int whichAcos){
        this.calcDiffs = calcAllDiffs;
        this.bonds = bonds;
        this.atomNames = atoms;
        this.atomicNumbers = numbers;
        this.positions = xyz;
        this.charges = charges;
        this.spins = spins;
        
        /*
         * our lookup table (potential overhead!)
         */
        switch(whichAcos){
            case 0: acosLook = null; wAcos = 0; break;
            case 1: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[0]); wAcos = 1; break;
            case 2: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[1]); wAcos = 1; break;
            case 3: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[2]); wAcos = 1; break;
            case 4: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[3]); wAcos = 1; break;
            case 5: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[4]); wAcos = 1; break;
            case 11: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[0]); wAcos = 2; break;
            case 12: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[1]); wAcos = 2; break;
            case 13: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[2]); wAcos = 2; break;
            case 14: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[3]); wAcos = 2; break;
            case 15: acosLook = new AcosLookup(AcosLookup.POSSIBLE_SIZES[4]); wAcos = 2; break;
            default: acosLook = null; wAcos = 0; break;
       }

        this.contributions13 = contr13;
        if(contr13 != null){
            this.contri13 = new int[contr13.size()][3];
            for(int i = 0; i < contr13.size(); i++) {contri13[i] = contr13.get(i);}
        }

        this.contributions14 = contr14;
        if(contr14 != null){
            this.contri14 = new int[contr14.size()][4];
            for(int i = 0; i < contr14.size(); i++) {contri14[i] = contr14.get(i);}
        }

        // calculate distances
        final int iNoOfAtoms = atoms.length;
        final double[][] daXYZ = xyz;

        this.distances = new double[iNoOfAtoms][iNoOfAtoms];
        if(calcDiffs) {this.diffPos = new double[iNoOfAtoms][iNoOfAtoms][3];}

        if(calcDiffs){
            for(int i = 0; i < iNoOfAtoms -1; i++){
                for(int j = i+1; j < iNoOfAtoms; j++){
                    final double dX = daXYZ[0][i] - daXYZ[0][j];
                    final double dY = daXYZ[1][i] - daXYZ[1][j];
                    final double dZ = daXYZ[2][i] - daXYZ[2][j];
                    
                    distances[i][j] = Math.sqrt(dX*dX + dY*dY + dZ*dZ);
                    diffPos[i][j][0] = dX;
                    diffPos[i][j][1] = dY;
                    diffPos[i][j][2] = dZ;
                    
                    // it's kind of symmetric
                    diffPos[j][i][0] = -dX;
                    diffPos[j][i][1] = -dY;
                    diffPos[j][i][2] = -dZ;
                    distances[j][i] = distances[i][j];
                }
            }
        } else{
            for(int i = 0; i < iNoOfAtoms -1; i++){
                for(int j = i+1; j < iNoOfAtoms; j++){
                    distances[i][j] = CoordTranslation.distance(daXYZ, i, j);
                    
                    // it's symmetric
                    distances[j][i] = distances[i][j];
                }
            }
        }
    }        
        
    public Topology(final Topology orig, final boolean shallow){
        /*
         * we need to be careful not to rise any NullPointerExceptions!
         * Thou shalt not follow the NULL pointer, for chaos and madness await
         * thee at its end.
         * ;-)
         * Some things (bonds and contributions) will always only be shallow copies!
         */
         if(orig.bonds != null) {this.bonds = orig.bonds.clone();}
         if(orig.contributions13 != null) {this.contributions13 = new ArrayList<>(orig.contributions13);}
         if(orig.contri13 != null) {this.contri13 = orig.contri13.clone();}
         if(orig.contri14 != null) {this.contri14 = orig.contri14.clone();}
         if(orig.contributions14 != null) {this.contributions14 = new ArrayList<>(orig.contributions14);}
         if(orig.charges != null) {this.charges = orig.charges.clone();}
         if(orig.spins != null) {this.spins = orig.spins.clone();}
         if(orig.atomNames != null) {this.atomNames = orig.atomNames.clone();}
         if(orig.atomicNumbers != null) {this.atomicNumbers = orig.atomicNumbers.clone();}
         if(orig.positions != null){
             // deep copy
             this.positions = new double[orig.positions.length][orig.positions[0].length];
             for(int i = 0; i < orig.positions.length; i++){
                 System.arraycopy(orig.positions[i], 0, positions[i], 0, orig.positions[0].length);
             }
         }
         if(orig.distances != null){
             if(!shallow){
                 // deep copy
                 this.distances = new double[orig.distances.length][orig.distances[0].length];
                 for(int i = 0; i < orig.distances.length; i++){
                     System.arraycopy(orig.distances[i], 0, distances[i], 0, orig.distances[0].length);
                }
             } else{
                 this.distances = orig.distances.clone();
             }
         }
         this.calcDiffs = orig.calcDiffs;
         if(orig.diffPos != null){
             if(!shallow){
                 // deep copy
                 this.diffPos = new double[orig.diffPos.length][orig.diffPos.length][3];
                 for(int i = 0; i < orig.diffPos.length; i++){
                     for(int j = 0; j < orig.diffPos.length; j++){
                        System.arraycopy(orig.diffPos[i][j], 0, this.diffPos[i][j], 0, 3);
                    }
                }
             } else{
                 this.diffPos = orig.diffPos.clone();
             }
         }
         
         this.wAcos = orig.wAcos;
         if(orig.acosLook != null) {this.acosLook = orig.acosLook.clone();} //XXX do we need a clone w/ deep copy?!
         else {this.acosLook = null;}
         
         if(orig.angles != null){
             if(!shallow){
                 this.angles = new double[orig.angles.length][orig.angles[0].length][orig.angles[0][0].length];
                 for(int i = 0; i < angles.length; i++){
                     for(int j = 0; j < angles[j].length; j++){
                        System.arraycopy(orig.angles[i][j], 0, angles[i][j], 0, angles[i][j].length);
                     }
                 }
             } else{
                 this.angles = orig.angles.clone(); // XXX think about this...
             }
         }
         
         if(orig.bondedAngles != null){
             if(!shallow){
                 //TODO implement
             } else{
                 this.bondedAngles = orig.bondedAngles.clone(); //XXX think about this!
             }
         }
         
         /*
          * We leave the rest null'd on purpose since we assume that the cost to
          * deep copy it is not much less than computing it again (if needed).
          * Also it is unlikely to be != null.
          */
    }
    
    public Topology(Topology orig){
        this(orig, false);
    }
    
    /**
     * 
     * @return a shallow copy of this topology
     */
    @Override
    public Topology clone(){
        return new Topology(this, true);
    }
    
    public int getNumberOfAtoms(){
        return atomNames.length;
    }

    public void setPositions(final double[][] positions){
        this.positions = positions;
        // make the differences invalid
        this.positions = null;
        this.diffPos = null;
    }

    public double[][] getPositions(){
        return positions;
    }

    public void setCharges(final float[] charges){
        this.charges = charges;
    }

    public float[] getCharges(){
        return charges;
    }

    public void setSpins(final short[] spins){
        this.spins = spins;
    }

    public short[] getSpins(){
        return spins;
    }

    public void setBonds(final BondInfo bonding){
        bonds = bonding;
    }

    public BondInfo getBonds(){
        return bonds;
    }

    public void setAtomNames(final String[] atomNames){
        this.atomNames = atomNames;
    }

    public String[] getAtomNames(){
        return atomNames;
    }
    
    public void setAtomicNumbers(final short[] nos){
        atomicNumbers = nos;
    }
    
    public short[] getAtomicNumbers(){
        return atomicNumbers;
    }

    public void setDistances(final double[][] distances){
        this.distances = distances;
        // make the position differences invalid (what they are now)
        this.diffPos = null;
    }

    public double[][] getDistances(){
        //XXX think about calculating them only upon request...
        return distances;
    }
    
    public double[][][] getPosDiffs(){
        
        if(diffPos != null) {return diffPos;}
        
        this.diffPos = new double [distances.length][distances.length][3];    
        for (int i = 0; i < distances.length - 1; i++) {
            for (int j = i + 1; j < distances.length; j++) {
                final double dX = positions[0][i] - positions[0][j];
                final double dY = positions[1][i] - positions[1][j];
                final double dZ = positions[2][i] - positions[2][j];

                diffPos[i][j][0] = dX;
                diffPos[i][j][1] = dY;
                diffPos[i][j][2] = dZ;

                // it's kind of symmetric
                diffPos[j][i][0] = -dX;
                diffPos[j][i][1] = -dY;
                diffPos[j][i][2] = -dZ;
            }
        }

        return diffPos;
    }

    public void setAngles(final double[][][] angles){
        this.angles = angles;
    }
    
    public double getAngle(final int i, final int j, final int k){
        
        assert(i != j);
        assert(j != k);
        assert(i != k);
        
        // only compute if necessary
        if(angles != null) {return angles[i][j][k];}
        //XXX acos()?!
        else {return CoordTranslation.calcAngle(positions, i, j, k, scr1, scr2);}
    }
    
    /**
     * Computes all internal angles in a triangle IJK.
     * @param i
     * @param j
     * @param k
     * @param angles will contain all internal angles in order IJK/JKI/KIJ.
     */
    public void getTriangle(final int i, final int j, final int k, final double[] angles){
        
        assert(i != j);
        assert(j != k);
        assert(i != k);
        assert(angles.length >= 3);
        
        // only compute if necessary
        if(this.angles != null){
            angles[0] = this.angles[i][j][k];
            angles[1] = this.angles[j][k][i];
            angles[2] = this.angles[k][i][j];
        } else{
            // compute (but try to be smart)
            if(diffPos != null){
                final double[] dIJ = diffPos[i][j];
                final double[] dIK = diffPos[i][k];
                final double[] dKJ = diffPos[k][j];
                final double lIJ = Math.sqrt(dIJ[0]*dIJ[0] + dIJ[1]*dIJ[1] + dIJ[2]*dIJ[2]);
                final double lIK = Math.sqrt(dIK[0]*dIK[0] + dIK[1]*dIK[1] + dIK[2]*dIK[2]);
                final double lKJ = Math.sqrt(dKJ[0]*dKJ[0] + dKJ[1]*dKJ[1] + dKJ[2]*dKJ[2]);
                final double dotIJwKJ = org.ogolem.math.TrivialLinearAlgebra.dotProduct(dIJ, dKJ);
                final double dotJKwIK = - org.ogolem.math.TrivialLinearAlgebra.dotProduct(dKJ, dIK);
                final double normIJ_KJ = dotIJwKJ/(lIJ*lKJ);
                final double normJK_IK = dotJKwIK/(lKJ*lIK);
                angles[0] = acos(normIJ_KJ);
                angles[1] = acos(normJK_IK);
            } else{
                for(int z = 0; z < 3; z++){
                    scr1[z] = positions[z][i] - positions[z][j];
                    scr2[z] = positions[z][i] - positions[z][k];
                    scr3[z] = positions[z][k] - positions[z][j];
                }
                final double lIJ = Math.sqrt(scr1[0] * scr1[0] + scr1[1] * scr1[1] + scr1[2] * scr1[2]);
                final double lIK = Math.sqrt(scr2[0] * scr2[0] + scr2[1] * scr2[1] + scr2[2] * scr2[2]);
                final double lKJ = Math.sqrt(scr3[0] * scr3[0] + scr3[1] * scr3[1] + scr3[2] * scr3[2]);
                final double dotIJwKJ = org.ogolem.math.TrivialLinearAlgebra.dotProduct(scr1, scr3);
                final double dotJKwIK = -org.ogolem.math.TrivialLinearAlgebra.dotProduct(scr3, scr2);
                final double normIJ_KJ = dotIJwKJ / (lIJ * lKJ);
                final double normJK_IK = dotJKwIK / (lKJ * lIK);
                angles[0] = acos(normIJ_KJ);
                angles[1] = acos(normJK_IK);
            }

            angles[2] = Math.PI - angles[0] -angles[1];
        }
    }

    /**
     * Should only be used when knowing what the memory footprint of this method
     * is.
     * @return angles A field of all angles in the system.
     */
    public double[][][] getAngles(){
        
        if(angles == null) {calcAllAngles();}
        return angles;
    }
    
    public void calcAllAngles(){
        
        final int noOfAngles = distances.length;
        angles = new double[noOfAngles][noOfAngles][noOfAngles];
        for(int i = 0; i < noOfAngles; i++){
            for(int j = 0; j < noOfAngles; j++){
                if(i == j) {continue;}
                for(int k = 0; k < noOfAngles; k++){
                    if(i == k || j == k) {continue;}
                    //XXX symmetry!
                    angles[i][j][k] = CoordTranslation.calcAngle(positions, i, j, k, scr1, scr2);
                }
            }
        }
    }

    public double[][][] getBondedAngles(){
        //XXX this method is unused and makes little sense, IMHO.
        if(bondedAngles != null) {return bondedAngles;}
        
        final int noOfAngles = distances.length;
        bondedAngles = new double[noOfAngles][noOfAngles][noOfAngles];
        
        for(int i = 0; i < noOfAngles; i++){
            for(int j = 0; j < noOfAngles; j++){
                if(i == j || !bonds.hasBond(i, j)) {continue;}
                for(int k = 0; k < noOfAngles; k++){
                    if(i == k || j == k || !bonds.hasBond(j, k)) {continue;}
                    //XXX symmetry!
                    bondedAngles[i][j][k] = CoordTranslation.calcAngle(positions, i, j, k, scr1, scr2);
                }
            }
        }
        
        return bondedAngles;
    }

    public void setDihedrals(double[][][][] dihedrals){
        this.dihedrals = dihedrals;
    }

    /**
     * Should only be used when knowing what the memory footprint of this method
     * is.
     * @return dihedrals A field of all dihedrals in the system.
     */
    public double[][][][] getDihedrals(){

        if(dihedrals != null) {return dihedrals;}

        // initialize dihedrals
        int iNoOfAngles = distances.length;
        dihedrals = new double[iNoOfAngles][iNoOfAngles][iNoOfAngles][iNoOfAngles];

        for(int i = 0; i < iNoOfAngles; i++){

            final double[] daPos1 = {positions[0][i],positions[1][i], positions[2][i]};

            for(int j = 0; j < iNoOfAngles; j++){

                final double[] daPos2 = {positions[0][j],positions[1][j], positions[2][j]};

                for(int k = 0; k < iNoOfAngles; k++){

                    final double[] daPos3 = {positions[0][k],positions[1][k], positions[2][k]};
                    
                    for(int l = 0; l < iNoOfAngles; l++){

                        final double[] daPos4 = {positions[0][l],positions[1][l], positions[2][l]};

                        if ((i == j) || (j == k) || (i == k) || (i == l) || (j == l) || (k == l) ){
                            dihedrals[i][j][k][l] = Double.NaN;
                            continue;
                        }

                        // compute dihedral
                        //XXX is also not acos aware... if it matters? who knows...
                        dihedrals[i][j][k][l] = CoordTranslation.calcDihedral(daPos1, daPos2, daPos3, daPos4);
                    }
                }
            }
        }

        return dihedrals;
    }

    public boolean does13exist(int i, int j){
        
        for(final int[] c13 : contributions13){
            if((i == c13[0] && j == c13[2]) || (j == c13[0] && i == c13[2])) {return true;}
        }
        
        return false;
    }

    public boolean does14exist(int i, int j){
        
        for(final int[] c14 : contributions14){
            if((i == c14[0] && j == c14[3]) || (j == c14[0] && i == c14[3])) {return true;}
        }
        
        return false;
    }

    public List<int[]> get13Contributions(){
        return contributions13;
    }

    public List<int[]> get14Contributions(){
        return contributions14;
    }
    
    public int[][] get13ContributionsField(){
        return contri13;
    }

    public int[][] get14ContributionsField(){
        return contri14;
    }
    
    public String[] createPrintableCartesians(){
        final int noAtoms = positions[0].length;
        final String[] sa = new String[noAtoms+2];
        sa[0] = "" + noAtoms;
        sa[1] = "";
        for(int i = 0; i < noAtoms; i++){
            sa[i+2] = atomNames[i] + "\t"
                    + positions[0][i]*Constants.BOHRTOANG + "\t"
                    + positions[1][i]*Constants.BOHRTOANG + "\t"
                    + positions[2][i]*Constants.BOHRTOANG;
        }
        
        return sa;
    }
    
    private double acos(final double x){
        switch(wAcos){
            case 0: return Math.acos(x);
            case 1: return acosLook.acosInter(x);
            case 2: return acosLook.acosNonInter(x);
            case 3: return org.apache.commons.math3.util.FastMath.acos(x);
            default: return Math.acos(x);
        }
    }
}
