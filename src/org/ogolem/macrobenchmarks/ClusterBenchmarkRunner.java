/**
Copyright (c) 2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.macrobenchmarks;

import java.io.File;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Constants;
import org.ogolem.core.GeometryConfig;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.Input;
import org.ogolem.core.MoleculeConfig;
import org.ogolem.io.InputPrimitives;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Run cluster optimization benchmarks.
 * @author Johannes Dieterich
 * @version 2020-05-24
 */
class ClusterBenchmarkRunner implements BenchmarkRunner {

    private static final Logger LOG = LoggerFactory.getLogger(ClusterBenchmarkRunner.class);
    
    @Override
    public boolean runBenchmark(String csvLine, int noThreads) throws Exception {

        final String[] csvData = csvLine.trim().split("\\,");
        // format of the line is as follows:
        // 1) the working directory the input resides in
        // 2) the input file, relative to work dir
        // 4) the expected global minimum energy
        // 5) the file, relative to work dir, containing the reference global minimum
        // 6) threshold for energy
        // 7) threshold for cluster structure
        final String workDir = csvData[0].trim();
        final String inputFile = csvData[1].trim();
        final double energy = Double.parseDouble(csvData[2].trim());
        final String refGeometry = csvData[3].trim();
        final double energyThresh = Double.parseDouble(csvData[4].trim());
        final double structThresh = Double.parseDouble(csvData[5].trim());
        
        final File newDir = new File(workDir);
        if(!newDir.exists() || !newDir.isDirectory()){
            LOG.error("Specified working directory " + workDir + " does not exist or is not a directory.");
            return false;
        }
        
        // let's parse us a GlobalConfig
        final GlobalConfig globConf = Input.ConfigureMe(workDir + File.separator + inputFile, true);
        final GeometryConfig gc = globConf.geoConfCopy();
        
        LOG.info("Executing benchmark " + inputFile);
        final double startT = System.currentTimeMillis();
        
        final String[] args = new String[]{"--shmem", inputFile, "" + noThreads};
        try {
            Helpers.executeJavaProcess(workDir, args);
        } catch(Exception e){
            throw new RuntimeException("Failure to run benchmark " + csvLine, e);
        }
        final double endT = System.currentTimeMillis();
        
        // we know the output directory is inputFile w/o .ogo extension
        final String outputDir = workDir + File.separator + inputFile.substring(0, inputFile.lastIndexOf(".ogo"));
        
        // find the rank0 geometry in that output directory
        final File outDir = new File(outputDir);
        final String[] ls = outDir.list(new Helpers.Rank0Filter());
        if(ls == null){
            final String[] files = outDir.list();
            System.out.println("ERROR: didn't find a rank0 - found files:");
            for(final String file : files){
                System.err.println(file);
            }
            throw new RuntimeException("Null'd filtered output directory for rank0 geometry for benchmark " + csvLine + " in " + outputDir);
        } else if (ls.length == 0) {
            throw new RuntimeException("No rank0 geometry for benchmark " + csvLine + " in " + outputDir);
        } else if (ls.length > 1) {
            throw new RuntimeException("Multiple rank0 geometry for benchmark " + csvLine + " in " + outputDir);
        }
        
        final String rank0 = outputDir + File.separator + ls[0];
        final String[] content = InputPrimitives.readFileIn(rank0);
        final String[] eLine = content[1].trim().split("\\s+");
        final double foundE = Double.parseDouble(eLine[1].trim())*Constants.KJTOHARTREE;
        if(Math.abs(foundE - energy) > energyThresh){
            LOG.error("Benchmark energy wrong. Obtained: " + foundE + " vs should be " + energy);
            return false;
        }

        // steps to that individual
        final String ls0NoFront = ls[0].substring("rank0individual".length());
        final String ls0NoBack = ls0NoFront.substring(0, ls0NoFront.length()-".xyz".length());
        final long stepsToIndividual = Long.parseLong(ls0NoBack);
        
        // now the more complex thing - compare structures
        
        int totAtoms = 0;
        for(final MoleculeConfig mc : gc.geomMCs){
            totAtoms += mc.noOfAtoms;
        }
        
        final int[] atsPerMol = new int[gc.noOfParticles];
        final short[] spins = new short[totAtoms];
        final float[] charges = new float[totAtoms];
        final String[] sids = new String[gc.noOfParticles];
        int off = 0;
        for(int i = 0; i < gc.geomMCs.size(); i++){
            final MoleculeConfig mc = gc.geomMCs.get(i);
            sids[i] = "mol" + i;
            atsPerMol[i] = mc.noOfAtoms;
            for(int j = 0; j < mc.noOfAtoms; j++){
                spins[j+off] = mc.spins[j];
                charges[j+off] = mc.charges[j];
            }
            off += mc.noOfAtoms;
        }
        
        final CartesianCoordinates refCartes = Input.readCartesFromFile(workDir + File.separator + refGeometry,
                gc.noOfParticles, atsPerMol, spins, charges);
        
        final CartesianCoordinates resultCartes = Input.readCartesFromFile(rank0,
                gc.noOfParticles, atsPerMol.clone(), spins.clone(), charges.clone());
        
        //TODO XXX
        // needs either lidheat branch (Hundt) or the graph thingy.
        LOG.info("No support yet for comparing structures.");
        
        LOG.info("Benchmark " + inputFile + " PASSED");
        LOG.info("Benchmark " + inputFile + " took: " + (endT-startT)/1000 + " s.");
        LOG.info("Benchmark " + inputFile + " took: " + stepsToIndividual + " steps.");
        
        return true;
    }
}
