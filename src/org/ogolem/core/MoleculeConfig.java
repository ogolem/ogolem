/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

/**
 * The configuration object for a Molecules class providing default values.
 * @author Johannes Dieterich
 * @version 2014-09-06
 */
public class MoleculeConfig implements Serializable, Cloneable {
	
    //the serial version ID
    private static final long serialVersionUID = (long) 20130925;
	
    // String ID
    String sID = "N/A";
    // the energy of this particular molcule
    double molecularEnergy = Double.MAX_VALUE;
    // Integer ID
    int iID = 0;
    // boolean semi-stiff or flexible molecule
    boolean flexy = false;
    // boolean stiff parts or not?
    public boolean constricted = false;
    // external center of mass as array of doubles
    double[] externalCOM = new double[3];
    // orientation
    double[] externalOrient = new double[3];
    // number of atoms in molecule
    public int noOfAtoms = 1;
    // the cartesian coordinates
    double[][] refXYZ;
    // the atom types
    String[] atomTypes;
    // the atom numbers
    short[] atomNumbers;
    // the z matrix
    ZMatrix zmat = null;
    // to support charges
    public float[] charges;
    // to support spins
    public short[] spins;
    // represents the degrees of freedom in the molecule
    public boolean[][] degreesOfFreedom = null;
    // represents the constrains in the molecules (based on XYZ!!!)
    public boolean[][] constraints = null;
    
    MoleculeConfig(final boolean dontAlloc){
        if(!dontAlloc) {
            this.externalCOM = new double[3];
            this.externalOrient = new double[3];
        }
    }

    MoleculeConfig(final MoleculeConfig orig){
        this.constricted = orig.constricted;
        this.flexy = orig.flexy;
        if(orig.constraints != null){
            this.constraints = orig.constraints.clone();
        }
        if(orig.degreesOfFreedom != null){
            this.degreesOfFreedom = orig.degreesOfFreedom.clone();
        }
        this.molecularEnergy = orig.molecularEnergy;
        this.externalCOM = orig.externalCOM.clone();
        this.externalOrient = orig.externalOrient.clone();
        this.charges = orig.charges.clone();
        this.iID = orig.iID;
        this.noOfAtoms = orig.noOfAtoms;
        this.spins = orig.spins.clone();
        this.atomTypes = orig.atomTypes.clone();
        this.atomNumbers = (orig.atomNumbers == null) ? null : orig.atomNumbers.clone();
        this.refXYZ = new double[3][this.noOfAtoms];
        System.arraycopy(orig.refXYZ[0], 0, this.refXYZ[0], 0, this.noOfAtoms);
        System.arraycopy(orig.refXYZ[1], 0, this.refXYZ[1], 0, this.noOfAtoms);
        System.arraycopy(orig.refXYZ[2], 0, this.refXYZ[2], 0, this.noOfAtoms);
        this.sID = orig.sID;
        this.zmat = (orig.zmat == null) ? null : new ZMatrix(orig.zmat);
    }
    
    MoleculeConfig(final MoleculeConfig orig, boolean deep){
        this.constricted = orig.constricted;
        this.flexy = orig.flexy;
        if(orig.constraints != null){
            if(!deep) {this.constraints = orig.constraints.clone();}
            else{
                this.constraints = new boolean[orig.constraints.length][];
                for(int i = 0; i < orig.constraints.length; i++) {this.constraints[i] = orig.constraints[i].clone();}
            }
        }
        if(orig.degreesOfFreedom != null){
            if(!deep) {this.degreesOfFreedom = orig.degreesOfFreedom.clone();}
            else{
                this.degreesOfFreedom = new boolean[orig.degreesOfFreedom.length][];
                for(int i = 0 ; i < orig.degreesOfFreedom.length; i++) {this.degreesOfFreedom[i] = orig.degreesOfFreedom[i].clone();}
            }
        }
        this.molecularEnergy = orig.molecularEnergy;
        this.externalCOM = orig.externalCOM.clone();
        this.externalOrient = orig.externalOrient.clone();
        this.charges = orig.charges.clone();
        this.iID = orig.iID;
        this.noOfAtoms = orig.noOfAtoms;
        this.spins = orig.spins;
        this.atomTypes = orig.atomTypes.clone();
        this.atomNumbers = (orig.atomNumbers == null) ? null : orig.atomNumbers.clone();
        this.refXYZ = new double[3][this.noOfAtoms];
        System.arraycopy(orig.refXYZ[0], 0, this.refXYZ[0], 0, this.noOfAtoms);
        System.arraycopy(orig.refXYZ[1], 0, this.refXYZ[1], 0, this.noOfAtoms);
        System.arraycopy(orig.refXYZ[2], 0, this.refXYZ[2], 0, this.noOfAtoms);
        this.sID = orig.sID;
        this.zmat = (orig.zmat == null) ? null : new ZMatrix(orig.zmat);
    }
    
    @Override
    public MoleculeConfig clone(){
        return new MoleculeConfig(this, true);
    }
    
    public CartesianCoordinates toCartesians(){
        return CoordTranslation.moleculeToCartesian(new Molecule(this),false);
    }

    private void writeObject(ObjectOutputStream oos) throws IOException {
        oos.writeInt(iID);
        oos.writeObject(sID);
        oos.writeInt(noOfAtoms);
        oos.writeDouble(molecularEnergy);
        oos.writeBoolean(flexy);
        oos.writeBoolean(constricted);
        oos.writeObject(externalCOM);
        oos.writeObject(externalOrient);
        oos.writeObject(zmat);
        oos.writeObject(charges);
        oos.writeObject(spins);
        oos.writeObject(degreesOfFreedom);
        oos.writeObject(constraints);
        oos.writeObject(atomTypes);
        oos.writeObject(atomNumbers);
        oos.writeObject(refXYZ);
    }

    private void readObject(ObjectInputStream ois) throws IOException,
            ClassNotFoundException, ClassCastException {

        iID = ois.readInt();
        sID = (String) ois.readObject();
        noOfAtoms = ois.readInt();
        molecularEnergy = ois.readDouble();
        flexy = ois.readBoolean();
        constricted = ois.readBoolean();
        externalCOM = (double[]) ois.readObject();
        externalOrient = (double[]) ois.readObject();
        zmat = (ZMatrix) ois.readObject();
        charges = (float[]) ois.readObject();
        spins = (short[]) ois.readObject();
        degreesOfFreedom = (boolean[][]) ois.readObject();
        constraints = (boolean[][]) ois.readObject();
        atomTypes = (String[]) ois.readObject();
        atomNumbers = (short[]) ois.readObject();
        refXYZ = (double[][]) ois.readObject();
    }
}
