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
 * This configures the Atom class and provides default values.
 * @author Johannes Dieterich
 * @version 2014-09-06
 */
class AtomConfig implements Serializable{
	
    // the serial version ID
    private static final long serialVersionUID = (long) 20121201;

    // String ID
    String sID = "N/A";
    // atom number
    short atomNo = -1;
    // int ID
    int iID = 0;
    // part of energy caused by this atom as double
    double energypart = 0.0;
    // the position in an arbitrary coordinate system (cartesian though)
    double[] position = {0.0, 0.0, 0.0};

    AtomConfig(){
        // everything as is
    }

    AtomConfig(AtomConfig orig){
        this.energypart = orig.energypart;
        this.iID = orig.iID;
        this.atomNo = orig.atomNo;
        this.sID = orig.sID;
        this.position = orig.position.clone();
    }

    private void writeObject(ObjectOutputStream oos) throws IOException{
         oos.writeDouble(energypart);
         oos.writeInt(iID);
         oos.writeInt(atomNo);
         oos.writeObject(sID);
         oos.writeObject(position);
     }

    private void readObject(ObjectInputStream ois) throws IOException,
            ClassNotFoundException, ClassCastException {

        energypart = ois.readDouble();
        iID = ois.readInt();
        atomNo = ois.readShort();
        sID = (String) ois.readObject();
        position = (double[]) ois.readObject();
    }
		
}
