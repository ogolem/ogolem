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
package org.ogolem.scanner;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.Constants;

/**
 * The configuration of the scanner.
 * @author Johannes Dieterich
 * @version 2013-09-22
 */
class ScannerConfig {
    
    ScannerConfig(final String[] input) throws Exception {
        
        this.scanDoFs = new ArrayList<>();
        
        for(int i = 0; i < input.length; i++){
            final String line = input[i].trim();
            if(line.equalsIgnoreCase("<SCANCONF>")){
                i++;
                while(!input[i].trim().equalsIgnoreCase("</SCANCONF>")){
                    final String l = input[i].trim();
                    if(l.startsWith("!") || l.startsWith("//") || l.startsWith("#")){
                        i++; continue;
                    } else if(l.startsWith("QuenchScanStructures")){
                        this.optimize = true;
                    } else if(l.startsWith("EnergyScanStructures")){
                        this.singlePoints = true;
                    } else {
                        // a DoF specification
                        final String[] sa = l.split("\\;");
                        final int mol = Integer.parseInt(sa[0].trim());
                        final int atom = Integer.parseInt(sa[1].trim());
                        final int dof = Integer.parseInt(sa[2].trim());
                        if(dof == 0){
                            // bond, in angstrom
                            final double start = Double.parseDouble(sa[3].trim())*Constants.ANGTOBOHR;
                            final double end = Double.parseDouble(sa[4].trim())*Constants.ANGTOBOHR;
                            final double incr = Double.parseDouble(sa[5].trim())*Constants.ANGTOBOHR;
                            final ScannerDoF thisDof = new ScannerDoF(mol,atom,dof,start,end,incr);
                            this.scanDoFs.add(thisDof);
                        } else {
                            // angle or torsion, in degrees
                            final double start = Math.toRadians(Double.parseDouble(sa[3].trim()));
                            final double end = Math.toRadians(Double.parseDouble(sa[4].trim()));
                            final double incr = Math.toRadians(Double.parseDouble(sa[5].trim()));
                            final ScannerDoF thisDof = new ScannerDoF(mol,atom,dof,start,end,incr);
                            this.scanDoFs.add(thisDof);
                        }
                    }
                    i++;
                }
                return;
            }
        }
    }
    
    /**
     * If we want to quench the scanned structures as well.
     */
    boolean optimize = false;
    
    /**
     * If we want to compute the energy of all the single point structures.
     */
    boolean singlePoints = false;
    
    /**
     * A list of all the degrees of freedom to be scanned over.
     */
    List<ScannerDoF> scanDoFs;
}
