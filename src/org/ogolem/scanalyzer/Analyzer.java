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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.Constants;
import org.ogolem.helpers.StatisticUtils;
import org.ogolem.helpers.Tuple;

/**
 * Analyzes a set of structures for typical energies etc pp.
 * @author Johannes Dieterich
 * @version 2013-09-23
 */
class Analyzer {
    private static final boolean DEBUG = false;
    private Analyzer(){}
    
    static void analyze(final String folderName, final String[] fileList, final int noBins, final double cutoff) throws Exception {
        
        final List<Double> energies = new ArrayList<>(fileList.length);
        
        // parse all the energies, kind of tedious...
        for(int i = 0; i < fileList.length; i++){
            final String fileName = folderName + File.separator + fileList[i];
            if(DEBUG){System.out.println("DEBUG: Reading in " + fileName);}
            try(final BufferedReader buffreader = new BufferedReader(new FileReader(fileName))){
                buffreader.readLine(); // number of atoms, who cares
                final String line = buffreader.readLine();
                if(line == null) {throw new RuntimeException("Line should be non-null, is null.");}
                final String[] sa = line.trim().split("\\s+");
                final double e = Double.parseDouble(sa[1])*Constants.KJTOHARTREE;
                if(DEBUG){System.out.println("DEBUG: read in energy " + e + " in kj/mol " + e*Constants.HARTREETOKJ);}
                if(e >= cutoff){continue;}
                energies.add(e);
            } catch (IOException e) {
                throw e;
            }
            
        }
        
        // now do some basic analyzation of them
        final Tuple<Double,Double> tup = StatisticUtils.meanAndStdDev(energies);
        System.out.println("  mean:    " + tup.getObject1()*Constants.HARTREETOKJ + " kJ/mol");
        System.out.println("  std.dev: " + tup.getObject2()*Constants.HARTREETOKJ + " kJ/mol");
        
        // lets see if we could bin?
        final Tuple<int[],double[]> tup2 = StatisticUtils.binData(energies, noBins);
        final int[] binned = tup2.getObject1();
        final double[] binOffs = tup2.getObject2();
        for(int i = 0; i < noBins; i++){
            System.out.println("  bin " + i + " start " + binOffs[i] + " population " + binned[i] + " percents " + (binned[i]*100.0/energies.size()));
        }
    }
}
