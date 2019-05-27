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

import java.util.List;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Constants;
import org.ogolem.core.SimpleBondInfo;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Converts Atomdroids atom-formatted files into proper xyz and a set of
 * bond information. Eventually, once we support more formats (will we ever)
 * one could picture renaming this package to org.ogolem.translator or something.
 * @author Johannes Dieterich
 * @version 2013-10-03
 */
public class MainAtom2Ogo {
    
    public static void run(final String[] args){
        
        if(args[0].equalsIgnoreCase("help")){
            System.out.println("This translates Atomdroid formatted files to ogolem bond info.");
            System.out.println("We need as arguments:");
            System.out.println(" * the path to a .atom file");
            System.out.println(" * the path to the new .xyz file");
            return;
        }
        
        final String atomFile = args[0];
        final String xyzFile = args[1];
        
        // read atom file
        String[] atomData = null;
        try{
            atomData = InputPrimitives.readFileIn(atomFile);
        } catch(Exception e){
            System.err.println("ERROR: Failed to read atom file. " + atomFile);
            e.printStackTrace(System.err);
            System.exit(1);
        }
        
        // now parse info and translate on the fly (not that difficult)
        final int noAtoms = Integer.parseInt(atomData[0].trim());
        final CartesianCoordinates cartes = new CartesianCoordinates(noAtoms,1,new int[]{noAtoms});
        final BondInfo info = new SimpleBondInfo(noAtoms);
        final String[] atoms = cartes.getAllAtomTypes();
        final double[][] xyz = cartes.getAllXYZCoord();
        for(int i = 0; i < noAtoms; i++){
            final String[] line = atomData[i+1].trim().split("\\s+");
            atoms[i] = line[0];
            xyz[0][i] = Double.parseDouble(line[1])*Constants.ANGTOBOHR;
            xyz[1][i] = Double.parseDouble(line[2])*Constants.ANGTOBOHR;
            xyz[2][i] = Double.parseDouble(line[3])*Constants.ANGTOBOHR;
        }
        
        // parse bonds
        for(int i = noAtoms+1; i < atomData.length; i++){
            final String[] line = atomData[i].trim().split("\\s+");
            final int at1 = Integer.parseInt(line[0]);
            final int at2 = Integer.parseInt(line[1]);
            final short bondType = Short.parseShort(line[2]);
            info.setBond(at1, at2, bondType);
        }
        
        final String[] xyzData = cartes.createPrintableCartesians();
        xyzData[1] = "ogolem translated atom file";
        
        // plot xyz
        try{
            OutputPrimitives.writeOut(xyzFile, xyzData, true);
        } catch(Exception e){
            System.err.println("ERROR: Couldn't write translated xyz file.");
            e.printStackTrace(System.err);
            System.exit(2);
        }
        
        // plot the bonds to output stream
        final List<String> bondInfo = info.translateToInput();
        for(final String s : bondInfo){
            System.out.println(s);
        }
    }
}
