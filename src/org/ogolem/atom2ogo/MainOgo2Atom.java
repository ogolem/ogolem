/**
Copyright (c) 2013, J. M. Dieterich
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
package org.ogolem.atom2ogo;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Constants;
import org.ogolem.core.Input;
import org.ogolem.core.SimpleBondInfo;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Converts an xyz with OGOLEM-style bond info into an Atomdroid atom-formatted file.
 * Eventually, once we support more formats (will we ever?) one could picture
 * renaming this package to org.ogolem.translator or something.
 * @author Johannes Dieterich
 * @version 2013-10-04
 */
public class MainOgo2Atom {
    
    public static void run(final String[] args){
        
        if(args[0].equalsIgnoreCase("help")){
            System.out.println("This translates xyz+bond info to an Atomdroid file.");
            System.out.println("We need as arguments:");
            System.out.println(" * the path to a .xyz file");
            System.out.println(" * the path to a bond list");
            System.out.println(" * the path to the new .atom file");
            return;
        }
        
        final String xyzFile = args[0];
        final String bondsFile = args[1];
        final String atomFile = args[2];
        
        // read the cartesians and bonds in
        CartesianCoordinates cartes = null;
        BondInfo info = null;
        try{
            cartes = Input.readCartesFromFile(xyzFile);
            info = Input.readBondsListIn(bondsFile, cartes.getNoOfAtoms());
        } catch(Exception e){
            System.err.println("ERROR: Couldn't read coordinate and/or bond info file.");
            e.printStackTrace(System.err);
            System.exit(1);
        }
        
        // stub our nice atomdroid format together
        final List<String> atomData = new ArrayList<>(2*cartes.getNoOfAtoms());
        final String[] coords = cartes.createPrintableCartesians();
        // first the coordinates
        atomData.add(coords[0]);
        for(int i = 0; i < cartes.getNoOfAtoms(); i++){
            final String line = coords[i+2];
            atomData.add(line + "  0.00"); // no partial charge
        }
        
        // then the bonds
        for(int at1 = 0; at1 < cartes.getNoOfAtoms()-1; at1++){
            for(int at2 = at1+1; at2 < cartes.getNoOfAtoms(); at2++){
                final short b = info.bondType(at1, at2);
                if(b == BondInfo.NOBOND || b == BondInfo.VDW) {continue;} // no bond or vdW "bond"
                else if(b == BondInfo.UNCERTAIN){
                    // atomdroid will get a single bond, whatever...
                    final String line = at1 + " " + at2 + " " + BondInfo.SINGLE;
                    atomData.add(line);
                } else {
                    final String line = at1 + " " + at2 + " " + b;
                    atomData.add(line);
                }
            }
        }
        
        // plot atom file
        try{
            OutputPrimitives.writeOut(atomFile, atomData, true);
        } catch(Exception e){
            System.err.println("ERROR: Couldn't write translated atom file.");
            e.printStackTrace(System.err);
            System.exit(2);
        }
    }
}
