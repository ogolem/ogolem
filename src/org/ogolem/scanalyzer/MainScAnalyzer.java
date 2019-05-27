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
package org.ogolem.scanalyzer;

import java.io.File;
import org.ogolem.core.FixedValues;
import org.ogolem.io.InquiryPrimitives;

/**
 * Does some small statistics on the output of a scanner run.
 * @author Johannes Dieterich
 * @version 2013-09-23
 */
public class MainScAnalyzer {
    
    public static void run(String args[]) {
    
        if(args[0].equalsIgnoreCase("help")){
            System.out.println("Please specify the following:");
            System.out.println(" * a directory containing all the scanned structures optionally followed by");
            System.out.println(" * usenounopt (if we do not want to analyze the not optimized structures) and/or");
            System.out.println(" * useopt (if you want to analyze the optimized structures)");
            System.out.println(" * nobins=XXX where XXX is the number of bins (default: 10)");
            System.out.println(" * cutoff=XXX where XXX is the cutoff in Hartee (default: NONCONVERGEDENEGY/1005.0)");
            return;
        }
        
        //TODO one could add a custom flag and allow arbitrary suffixes to be specified?
        
        boolean nonOptStructures = true;
        boolean optStructures = false;
        int noBins = 10;
        final String structFolder = args[0];
        double cutoff = FixedValues.NONCONVERGEDENERGY;
        for(int i = 1; i < args.length; i++){
            if(args[i].equalsIgnoreCase("usenounopt")){
                nonOptStructures = false;
            } else if(args[i].equalsIgnoreCase("useopt")){
                optStructures = true;
            } else if(args[i].startsWith("nobins=")){
                final String s = args[i].substring(7);
                noBins = Integer.parseInt(s);
            } else if(args[i].startsWith("cutoff=")){
                final String s = args[i].substring(7);
                cutoff = Double.parseDouble(s);
            } else{
                System.err.println("ERROR: Unknown argument " + args[i]);
                System.exit(1);
            }
        }
        
        final File f = new File(structFolder);
        if(!(f.exists() && f.isDirectory())){
            System.err.println("ERROR: Are you sure folder " + structFolder + " exists?");
            System.exit(2);
        }
        
        // find a list of all the nonoptstructure files and analyze (if wanted)
        if (nonOptStructures) {
            try {
                final String[] nonopt = InquiryPrimitives.fileListWithSuffix("-nonopt.xyz", structFolder);
                System.out.println("Analyzing non-optimized structure energies:");
                Analyzer.analyze(structFolder,nonopt,noBins,cutoff);
            } catch (Exception e) {
                System.err.println("WARNING: Some problem during analyzation of the non-opt structures");
                e.printStackTrace(System.err);
            }
        }
        
        // find a list of all the optstructure files and analyze (if wanted)
        if (optStructures) {
            try {
                final String[] nonopt = InquiryPrimitives.fileListWithSuffix("-opt.xyz", structFolder);
                System.out.println("Analyzing optimized structure energies:");
                Analyzer.analyze(structFolder,nonopt,noBins,cutoff);
            } catch (Exception e) {
                System.err.println("WARNING: Some problem during analyzation of the opt structures");
                e.printStackTrace(System.err);
            }
        }
    }
}
