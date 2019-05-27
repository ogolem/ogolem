/**
Copyright (c) 2010-2014, J. M. Dieterich
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
package org.ogolem.freqs;

import org.ogolem.core.CartesianFullBackend;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Geometry;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Newton;

/**
 * Calculate frequencies for cartesian coordinates.
 * @author Johannes Dieterich
 * @version 2014-11-21
 */
public class MainFreqs {
    
    /**
     * Entry point to frequency calculations.
     * @param args First argument is the path to an ogo file, then to the
     * xyz file, second one the method choice, third one the parameter filename,
     * fourth one (if present) enabling full
     * output for visualization of each frequency (experimental!).
     */
    public static void execute(final String[] args){
        
        if(args == null || args[0].equalsIgnoreCase("help")){
            System.out.println("This is the harmonic frequencies subpackage. Mandatory input arguments:");
            System.out.println(" * the configuration file (a .ogo one)");
            System.out.println(" * the geometry of to be analyzed as a .xyz.");
            System.out.println("Optional input arguments:");
            System.out.println(" * a boolean whether or not full visualization (i.e., the actual vibrations) is wished. Default: false.");
            System.out.println(" * prefix for the vibrational files. Default: freqs.");
            System.exit(0);
        }
        
        // setup cartesian coordinates
        CartesianCoordinates cartes = null;
        Geometry geom = null;
        GlobalConfig conf = null;
        try{
            conf = org.ogolem.core.Input.ConfigureMe(args[0]);
            geom = new Geometry(conf.geoConfCopy());
            cartes = geom.getCartesians();
            
            // update with the real coordinates
            final CartesianCoordinates c2 = org.ogolem.core.Input.readCartesFromFile(args[1]);
            final double[][] xyzReal = c2.getAllXYZCoordsCopy();
            cartes.setAllXYZ(xyzReal);
            
        } catch(Exception e){
            System.err.println("ERROR: Couldn't read cartesian coordinates file. Aborting.");
            e.printStackTrace(System.err);
            System.exit(1);
        }
        
        // setup method: we know that any FrequencyMethod must be backend, so we can try there
        FrequencyMethod method = null;
        try{
            final Newton newton = conf.getRefNewton();
            if(newton instanceof FrequencyMethod){
                method = (FrequencyMethod) newton;
            } else {
                final CartesianFullBackend back = newton.getBackend();
                if(back instanceof FrequencyMethod){
                    method = (FrequencyMethod) back;
                } else {
                    throw new Exception("Neither local optimization nor backend are instances of FrequencyMethod.");
                }
            }
        } catch(Exception e){
            System.err.println("ERROR: Couldn't setup frequency method. Aborting.");
            e.printStackTrace(System.err);
            System.exit(2);
        }
        
        // full visualization?
        boolean fullVisual = false;
        String prefix = "freqs";
        try {
            fullVisual = Boolean.parseBoolean(args[2]);
            prefix = args[3];
        } catch (Exception e) {
            System.err.println("INFO: Couldn't parse if full frequency visualization is wanted. Not doing it." + e.toString());
        }

        // calculate frequencies for this combination of coordinates and method
        final Frequencies freqs = HarmonicFrequencyCalculator.calculateFrequencies(cartes, method, geom.getBondInfo());

        // printout frequencies
        final String[] freqsOut = freqs.printableFrequencies();

        System.out.println("Frequency calculation finished successfully. Information coming...");
        for(final String s : freqsOut){
            System.out.println(s);
        }
        
        if(fullVisual){
            System.out.println("INFO: Printing out all frequencies for visualization with prefix " + prefix);
            freqs.printAllFrequencies(prefix, cartes);
        }
   }
}
