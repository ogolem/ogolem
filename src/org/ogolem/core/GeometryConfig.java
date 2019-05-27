/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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
import java.util.ArrayList;
import java.util.List;

/**
 * This is the configuration object for the Geometry class and provides default
 * values.
 * @author Johannes Dieterich
 * @version 2013-09-26
 */
public final class GeometryConfig implements Serializable, Cloneable {

    // the serial version ID
    private static final long serialVersionUID = (long) 20130926;

    // long ID
    long lID;
    // ID of "father" geometry
    long fatherID;
    // ID of "mother" geometry
    long motherID;
    // Number of independent particles
    public int noOfParticles;
    // Fitness of this Geometry
    double fitness;
    // ArrayList of MoleculeConfigs
    public List<MoleculeConfig> geomMCs;
    // the environment
    public Environment env = null;
    // the bonds
    public BondInfo bonds;
    
    public GeometryConfig(){      
    }

    @Override
    public GeometryConfig clone(){
        
        final GeometryConfig gc = new GeometryConfig();
        gc.lID = this.lID;
        gc.motherID = this.motherID;
        gc.fatherID = this.fatherID;
        assert(this.noOfParticles > 0);
        gc.noOfParticles = this.noOfParticles;
        gc.fitness = this.fitness;
        gc.env = (this.env == null) ? null : this.env.clone();
        gc.geomMCs = new ArrayList<>();
        assert(geomMCs.size() == this.noOfParticles);
        for(final MoleculeConfig mc : this.geomMCs){
            gc.geomMCs.add(mc.clone());
        }
        assert(this.bonds != null);
        gc.bonds = this.bonds.clone();
        
        return gc;
    }

    private void writeObject(ObjectOutputStream oos) throws IOException {
        oos.writeLong(lID);
        oos.writeLong(fatherID);
        oos.writeLong(motherID);
        oos.writeInt(noOfParticles);
        oos.writeDouble(fitness);
        oos.writeObject(geomMCs);
        oos.writeObject(env);
        oos.writeObject(bonds);
    }

    @SuppressWarnings(value = "unchecked")
    private void readObject(ObjectInputStream ois) throws IOException,
            ClassNotFoundException, ClassCastException {
        lID = ois.readLong();
        fatherID = ois.readLong();
        motherID = ois.readLong();
        noOfParticles = ois.readInt();
        fitness = ois.readDouble();
        geomMCs = (ArrayList<MoleculeConfig>) ois.readObject();
        env = (Environment) ois.readObject();
        bonds = (BondInfo) ois.readObject();
    }
}
