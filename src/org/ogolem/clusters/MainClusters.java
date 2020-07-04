/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.clusters;

import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import org.ogolem.core.CartesianCoordinates;
// get some physical constants at hand
import static org.ogolem.core.Constants.*;
import org.ogolem.core.Geometry;
import org.ogolem.core.GeometryWriter;
import org.ogolem.core.Input;
import org.ogolem.core.Molecule;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.GenericPoolEntry;
import org.ogolem.generic.genericpool.Niche;
import org.ogolem.io.OutputPrimitives;

/**
 * A posteriori analysis of ogolem runs.
 * @author Johannes Dieterich
 * @author Bernd Hartke
 * @version 2020-06-07
 */
public class MainClusters{

    /**
     * The actual entry point to analysis.
     * -i inputfile (default: pool.bin)
     * -o outputfile (default: cluster-analysis.out)
     * -dissblow blow factor for the dissociation detection (default: 1.8)
     * -threshsame threshold for the moments of inertia, defining when they are "equal" (default: 1E-1)
     * -threshdiff threshold for the moments of inertia, defining when they are "different" (default: 2.5E-1)
     * -noofbins the number of bins (default: 100)
     * -noenergies disables scanning of energies (default: on)
     * -noinertias disables moments of inertia analysis (default: on)
     * -nocomdiffs disables difference of COMs analysis (default: on)
     * -nodissdetect disables dissociation detection (default: on)
     * -ljstrains (default: off)
     * -getstructs reads all structures out of a pool (default: off)
     * -binstructs enable binning by amount of structures instead of energy (default: off)
     * @param args
     */
    public static void main(final String[] args) {
        
        if(args == null || args.length == 0 || args[0].equalsIgnoreCase("help")){
            System.out.println("This is the cluster analysis subpackage. Possible options: ");
            System.out.println(" * -i POOLFILE, to use POOLFILE as the binary pool. Deault: pool.bin");
            System.out.println(" * -o OUTFILE, to use OUTFILE as the output file. Default: cluster-analysis.out");
            System.out.println(" * -dissblow XX.X, to set the DD blow factor to XX.X. Default: 1.8");
            System.out.println(" * -threshame XX.X, to set the moments of inertia threshold for identical moments to XX.X. Default: 1E-1");
            System.out.println(" * -threshdiff XX.X, to set the moments of inertia threshold for different moments to XX.X. Default: 2.5E-1");
            System.out.println(" * -noenergies, disables energy analysis");
            System.out.println(" * -noinertias, disables moments of intertia analysis");
            System.out.println(" * -nodissdectect, disables the DD");
            System.out.println(" * -ljstrains, enable LJ strain calculation");
            System.out.println(" * -getstructs, to get all structures in the binary pool");
            System.out.println(" * -binstructs, to bin by number of structures instead of energy");
            System.out.println(" * -getniches, to see the distribution of pool individuals into niches, in OUTFILE");
            System.exit(0);
        }

        // first read in what we are supposed to do and what not, depending on the args
        String sOutputFile = "clusters-analysis.out";
        String sOutBinsShape = "clusters-shapebins.dat";
        String sOutBinsDiss = "clusters-dissbins.dat";
        String sInputFile = "pool.bin";
        String sStructsFolder = "structs";

        boolean bEnergies = true;
        boolean bMomentsOfInertia = true;
        boolean bCOMDiffs = true;
        boolean bDissociations = true;
        boolean bBinStructs = false;
        boolean bGetStructs = false;
        boolean bLJStrains = false;
        boolean bGetNiches = false;
        
        double dBlowDissoc = 2.5;
        double dThreshSame = 1E-1;
        double dThreshDiff = 2.5E-1;

        int iNoOfBins = 100;

        for(int i = 0; i <args.length; i++){
            if(args[i].equalsIgnoreCase("-o")){
                // the next entry in the array contains the outputfile
                if(args[i+1].startsWith("-")){
                    System.err.println("ERROR: You need to specify an output file after the -o switch. Continuing.");
                } else{
                    sOutputFile = args[i+1];
                }
                i = i + 1;
                continue;
            } else if(args[i].equalsIgnoreCase("-i")){
                // the next entry in the array contains the inputfile
                if(args[i+1].startsWith("-")){
                    System.err.println("ERROR: You need to specify an input file after the -i switch. Continuing.");
                } else{
                    sInputFile = args[i+1];
                }
                i = i + 1;
                continue;
            } else if(args[i].equalsIgnoreCase("-dissblow")){
                // the next entry in the array contains the blow factor
                if(args[i+1].startsWith("-")){
                    System.err.println("ERROR: You need to specify a blow factor after " +
                            "the -dissblow switch. Continuing.");
                } else{
                    try{
                        dBlowDissoc = Double.parseDouble(args[i+1]);
                    } catch(Exception e){
                        System.err.println("ERROR: Failure to parse the blow factor for " +
                                "the dissociation detection. Continuing with default.");
                    }
                }
                i = i + 1;
                continue;
            } else if(args[i].equalsIgnoreCase("-threshsame")){
                // the next entry in the array contains the blow factor
                if(args[i+1].startsWith("-")){
                    System.err.println("ERROR: You need to specify a factor after " +
                            "the -threshsame switch. Continuing.");
                } else{
                    try{
                        dThreshSame = Double.parseDouble(args[i+1]);
                    } catch(Exception e){
                        System.err.println("ERROR: Failure to parse the equal threshhold factor for " +
                                "the moments of inertia. Continuing with default.");
                    }
                }
                i = i + 1;
                continue;
            } else if(args[i].equalsIgnoreCase("-threshdiff")){
                // the next entry in the array contains the blow factor
                if(args[i+1].startsWith("-")){
                    System.err.println("ERROR: You need to specify a factor after " +
                            "the -threshdiff switch. Continuing.");
                } else{
                    try{
                        dThreshDiff = Double.parseDouble(args[i+1]);
                    } catch(Exception e){
                        System.err.println("ERROR: Failure to parse the difference threshhold factor for " +
                                "the moments of inertia. Continuing with default.");
                    }
                }
                i = i + 1;
                continue;
            } else if(args[i].equalsIgnoreCase("-noofbins")){
                // the next entry in the array contains the number of bins
                if(args[i+1].startsWith("-")){
                    System.err.println("ERROR: You need to specify the number of bins after " +
                            "the -noofbins switch. Continuing.");
                } else{
                    try{
                        iNoOfBins = Integer.parseInt(args[i+1]);
                    } catch(Exception e){
                        System.err.println("ERROR: Failure to parse the number of bins. " +
                                "Continuing with default.");
                    }
                }
                i = i + 1;
                continue;
            } else if(args[i].equalsIgnoreCase("-noenergies")){
                bEnergies = false;
                continue;
            } else if(args[i].equalsIgnoreCase("-noinertias")){
                bMomentsOfInertia = false;
                continue;
            } else if(args[i].equalsIgnoreCase("-nocomdiffs")){
                bCOMDiffs = false;
                continue;
            } else if(args[i].equalsIgnoreCase("-nodissdetect")){
                bDissociations = false;
                continue;
            } else if(args[i].equalsIgnoreCase("-binstructs")){
                bBinStructs = true;
                continue;
            } else if(args[i].equalsIgnoreCase("-getstructs")){
                bGetStructs = true;
                continue;
            } else if(args[i].equalsIgnoreCase("-ljstrains")){
                bLJStrains = true;
                continue;
            } else if(args[i].equalsIgnoreCase("-getniches")){
                bGetNiches = true;
                continue;
            } else{
                // if we end up here, there was a non valid option
                System.err.println("ERROR: Unrecognized option " + args[i] + " Continuing.");
            }
        }

        /*
         * configuring the program is done, start the analysis by reading the pool in
         */
        GenericPool<Molecule,Geometry> pool = null;
        try{
            String sDir = System.getProperty("user.dir");
            pool = Input.ReadSerializedPool(sDir, sInputFile);
        } catch(Exception e){
            // this is crucial
            System.err.println("ERROR: Failure to read the serialized pool in. Exiting ");
            e.printStackTrace(System.err);
            System.exit(42);
        }

        // for the binning, we need a couple of things: first all energies
        final double[] daFitnesses = pool.getAllFitnesses();

        final int[] iaWhichBinForGeo = new int[daFitnesses.length];
        final double[] daBins = new double[iNoOfBins];

        if (!bBinStructs) {
            // calculate the binsize and store the binoffset
            double dBinOffset = daFitnesses[0];
            double dBinSize = (daFitnesses[daFitnesses.length - 1] - daFitnesses[0]) / iNoOfBins;
            dBinSize = Math.abs(dBinSize);

            double dTemp = dBinOffset;
            for (int i = 0; i < iNoOfBins; i++) {
                daBins[i] = dTemp * HARTREETOKJ;
                dTemp += dBinSize;
            }

            // which bin for which struct
            for(int i = 0; i < daFitnesses.length; i++){
                double dWhichBin = Math.abs(daFitnesses[i] - dBinOffset) / dBinSize;
                dWhichBin = Math.floor(dWhichBin);
                int iWhichBin = (int) dWhichBin;
                if(iWhichBin >= iNoOfBins){
                    iWhichBin = iNoOfBins - 1;
                }
                iaWhichBinForGeo[i] = iWhichBin;
            }

        } else {
            final int iBinSize = daFitnesses.length / iNoOfBins;

            if((iBinSize * iNoOfBins) != daFitnesses.length){
                System.err.println("ERROR: There is a mismatch between the amount of structures and " +
                        "the amount of bins. Please change that. Exiting.");
                System.exit(112);
            }

            // populate bin information
            for (int i = 0; i < iNoOfBins; i++) {
                daBins[i] = (double) iBinSize * i;
            }

            // which bin for which struct
            int iTempCountBins = 0;
            int iTempCountStruct = -1;
            for(int i = 0; i < daFitnesses.length; i++){
                if(iTempCountStruct < iBinSize - 1){
                    // this is the correct bin
                    iaWhichBinForGeo[i] = iTempCountBins;
                    iTempCountStruct++;
                } else {
                    // next bin
                    iTempCountBins++;
                    iTempCountStruct = 0;
                    iaWhichBinForGeo[i] = iTempCountBins;
                }
            }
        }

        final LinkedList<String> llOutput = new LinkedList<>();

        final int iNumberOfStructs = pool.getCurrentPoolSize();

        final DecimalFormat formPercentages = new DecimalFormat("0.00");

        /*
         * energies?
         */
        if(bEnergies){
            DecimalFormat formEnergies = new DecimalFormat("0.0000");

            llOutput.add("Scanning of energies");
            llOutput.add("\t rank\t energy [kJ/mol]");

            CartesianCoordinates cartesTemp;

            for(int i = 0; i < iNumberOfStructs; i++){
                cartesTemp = pool.getIndividualAtPosition(i).getCartesians();
                double dEnergy = cartesTemp.getEnergy() * HARTREETOKJ;

                llOutput.add("\t" + i + "\t" + formEnergies.format(dEnergy));
            }

            llOutput.add("");
        }

        /*
         * should we get the structures?
         */
        if(bGetStructs){
            llOutput.add("GETTING ALL STRUCTURES FROM THE POOL");
            llOutput.add("putting data to folder " + sStructsFolder);

            try {
                final File fold = new File(sStructsFolder);
                if(!fold.exists()){
                    final boolean suc = fold.mkdir();
                    if(!suc){throw new RuntimeException("Folder " + sStructsFolder + " could not be created");}
                } else if(fold.exists() && !fold.isDirectory()){
                    throw new RuntimeException("Folder " + sStructsFolder + "is not a folder!");
                }
                
                final GeometryWriter writer = new GeometryWriter(sStructsFolder);
                int rank = 0;
                for(final GenericPoolEntry<Molecule,Geometry> entry : pool){
                    final Geometry g = entry.getIndividual();
                    final String toFile = sStructsFolder + File.separator + "rank" + rank + "geometry" + g.getID();
                    writer.writeIndividual(g, toFile);
                    rank++;
                }
            } catch (Exception e) {
                System.err.println("Error in final writing!" + e.toString());
                System.exit(1);
            }
        }

        /*
         * moments of intertia?
         */
        if(bMomentsOfInertia){

            DecimalFormat formInertias = new DecimalFormat("0.000");

            llOutput.add("ANALYZATION OF MOMENTS OF INERTIA");
            llOutput.add("\t rank\t moment 2\t moment 3\t shape");

            CartesianCoordinates cartesTemp;

            int[] iaShapeTypes = new int[4];

            int[][] iaBinningInfo = new int[iNoOfBins][4];

            for(int i = 0; i < iNumberOfStructs; i++){
                cartesTemp = pool.getIndividualAtPosition(i).getCartesians();
                double[] daTempMoments = ClusterAnalyzation.calcMomentsOfInertia(cartesTemp);

                final ClusterAnalyzation.Shape shape = ClusterAnalyzation.shapeOfCluster(daTempMoments, dThreshSame, dThreshDiff);

                String sShape = "";
                int iShape = 0;
                if(shape == ClusterAnalyzation.Shape.SPHERICAL){
                    sShape = "spherical";
                    iaShapeTypes[0]++;
                    iShape = 0;
                } else if(shape == ClusterAnalyzation.Shape.PROLATE){
                    sShape = "prolate";
                    iaShapeTypes[1]++;
                    iShape = 1;
                } else if(shape == ClusterAnalyzation.Shape.OBLATE){
                    sShape = "oblate";
                    iaShapeTypes[2]++;
                    iShape = 2;
                } else if(shape == ClusterAnalyzation.Shape.SCALENE){
                    sShape = "scalene";
                    iaShapeTypes[3]++;
                    iShape = 3;
                } else{
                    System.err.println("ERROR: It should be impossible to end up here. Notify the author please. ;-)");
                }

                // increment the correct bin
                iaBinningInfo[iaWhichBinForGeo[i]][iShape]++;

                String sTempOut = "\t " + i + "\t " + formInertias.format(daTempMoments[1]) + "\t\t " +
                        formInertias.format(daTempMoments[2]) + "\t\t " + sShape;
                llOutput.add(sTempOut);
            }

            // total analysis: ratios
            
            llOutput.add("\t spherical: \t" + formPercentages.format((double)iaShapeTypes[0]/ (double)iNumberOfStructs * 100) + " percent");
            llOutput.add("\t prolate: \t" + formPercentages.format((double)iaShapeTypes[1]/(double)iNumberOfStructs * 100) + " percent");
            llOutput.add("\t oblate: \t" + formPercentages.format((double)iaShapeTypes[2]/(double)iNumberOfStructs * 100) + " percent");
            llOutput.add("\t non-defined: \t" + formPercentages.format((double)iaShapeTypes[3]/(double)iNumberOfStructs * 100) + " percent");

            // create the binning output
            DecimalFormat formBins = new DecimalFormat("#.###");
            
            LinkedList<String> llBinOut = new LinkedList<>();
            
            llBinOut.add("#binoffset\t amounts: sph\t pro\t obl\t n-d\t ratios: sph\t pro\t obl\t n-d");
            for(int i = 0; i < iNoOfBins; i++){

                final StringBuffer sBuff = new StringBuffer();
                sBuff.append(formBins.format(daBins[i]));
                sBuff.append("\t\t");

                int iTotalInBin = 0;
                // first add the amounts
                for(int j = 0; j < 4; j++){
                    sBuff.append("\t");
                    sBuff.append(formBins.format(iaBinningInfo[i][j]));
                    iTotalInBin += iaBinningInfo[i][j];
                }

                sBuff.append("\t");

                // to avoid division by 0
                if(iTotalInBin == 0){
                    iTotalInBin = 1;
                }

                // then calculate the ratios
                for(int j = 0; j < 4; j++){
                    double dTempRatio = (double) iaBinningInfo[i][j] / (double)iTotalInBin * 100.0;
                    sBuff.append("\t");
                    sBuff.append(formPercentages.format(dTempRatio));
                }

                final String sLine = sBuff.toString();

                llBinOut.add(sLine);
            }

            // write the binning to a file
            try{
                OutputPrimitives.writeOut(sOutBinsShape, llBinOut, true);
            } catch(Exception e){
                System.err.println("ERROR: Failure to write binning data of shapes out. Continuing.");
                e.printStackTrace(System.err);
            }

            llOutput.add("gnuplot ready shapes binning output in " + sOutBinsShape);
        }

        /*
         * COM differences?
         */
        if(bCOMDiffs){

            final DecimalFormat formCOMs = new DecimalFormat("0.000");

            llOutput.add("ANALYZATION OF COMS");

            CartesianCoordinates cartesTemp;

            for(int i = 0; i < iNumberOfStructs; i++){
                llOutput.add("STRUCTURE NO. " + i);
                llOutput.add("\t COM-No\t x-value\t y-value\t z-value\t size");
                cartesTemp = pool.getIndividualAtPosition(i).getCartesians();

                final double[][] daCOMDiffs = ClusterAnalyzation.calculateDistanceToCOM(cartesTemp);

                for(int j = 0; j < daCOMDiffs.length; j++){
                    final double[] daTempDiff = daCOMDiffs[j];
                    final double dSize = Math.sqrt(
                            Math.pow(daTempDiff[0], 2) + Math.pow(daTempDiff[1], 2) +
                            Math.pow(daTempDiff[2], 2)
                            );
                    llOutput.add("\t " + j + "\t " + formCOMs.format(daTempDiff[0] * BOHRTOANG) + "\t\t "
                            + formCOMs.format(daTempDiff[1] * BOHRTOANG) + "\t\t "
                            + formCOMs.format(daTempDiff[2] * BOHRTOANG) + "\t\t "
                            + formCOMs.format(dSize * BOHRTOANG));
                }
            }
        }

        /*
         * dissociation?
         */
        if(bDissociations){
            CartesianCoordinates cartesTemp;

            llOutput.add("ANALYZING CLUSTERS FOR DISSOCIATIONS.");

            final int[][] iaBinningInfo = new int[iNoOfBins][2];

            final int[] iaDissocs = new int[2];

            for(int i = 0; i < iNumberOfStructs; i++){
                cartesTemp = pool.getIndividualAtPosition(i).getCartesians();

                final boolean bDissocTemp = ClusterAnalyzation.detectDissociations(cartesTemp, dBlowDissoc);

                String sDissoc = "";
                int iTempBin;
                if(bDissocTemp){
                    sDissoc = "dissociated";
                    iaDissocs[0]++;
                    iTempBin = 0;
                } else{
                    sDissoc = "non-dissociated";
                    iaDissocs[1]++;
                    iTempBin = 1;
                }

                // increment the correct bin
                iaBinningInfo[iaWhichBinForGeo[i]][iTempBin]++;

                String sTempOut = "\t" + i + "\t" + sDissoc;
                llOutput.add(sTempOut);
            }

            // total analysis: ratios

            llOutput.add("\t dissociated: \t" + formPercentages.format((double)iaDissocs[0]/(double)iNumberOfStructs * 100.0) + " percent");
            llOutput.add("\t non-dissociated: \t" + formPercentages.format((double)iaDissocs[1]/(double)iNumberOfStructs * 100.0) + " percent");

            // create the binning output
            final DecimalFormat formBins = new DecimalFormat("#.###");

            final LinkedList<String> llBinOut = new LinkedList<>();

            llBinOut.add("#binoffset\t amounts: diss\t non-diss\t ratios: diss\t non-diss");
            for(int i = 0; i < iNoOfBins; i++){

                final StringBuffer sBuff = new StringBuffer();
                sBuff.append(formBins.format(daBins[i]));
                sBuff.append("\t\t");

                int iTotalInBin = 0;
                // first add the amounts
                for(int j = 0; j < 2; j++){
                    sBuff.append("\t");
                    sBuff.append(formBins.format(iaBinningInfo[i][j]));
                    iTotalInBin += iaBinningInfo[i][j];
                }

                sBuff.append("\t");

                // to avoid division by 0
                if(iTotalInBin == 0){
                    iTotalInBin = 1;
                }

                // then calculate the ratios
                for(int j = 0; j < 2; j++){
                    double dTempRatio = (double) iaBinningInfo[i][j] /(double) iTotalInBin * 100.0;
                    sBuff.append("\t");
                    sBuff.append(formPercentages.format(dTempRatio));
                }

                llBinOut.add(sBuff.toString());
            }

            // write the binning to a file
            try{
                OutputPrimitives.writeOut(sOutBinsDiss, llBinOut, true);
            } catch(Exception e){
                System.err.println("ERROR: Failure to write binning data of dissociations out. Continuing.");
                e.printStackTrace(System.err);
            }

            llOutput.add("gnuplot ready dissociation binning output in " + sOutBinsDiss);
        }

        /*
         * calculate lj strains?
         */
        if(bLJStrains){
            CartesianCoordinates cartesTemp;

            llOutput.add("");
            llOutput.add("");
            llOutput.add("ANALYZING CLUSTERS FOR LJ STRAIN ENERGY.");
            llOutput.add("");
            llOutput.add("\tID\tNearest Energy \t Strain Energy\t Non Nearest Energy\t Total Energy");

            final DecimalFormat formEnergies = new DecimalFormat("0.000");

            final double dHartToKJ = org.ogolem.core.Constants.HARTREETOKJ;

            for(int i = 0; i < iNumberOfStructs; i++){
                cartesTemp = pool.getIndividualAtPosition(i).getCartesians();

                final double dNeighsE = ClusterAnalyzation.ljNearestNeighbors(cartesTemp);

                final double dStrainE = ClusterAnalyzation.ljStrainEnergy(cartesTemp);

                final double dNonNearestE = ClusterAnalyzation.ljEnergyNonNearest(cartesTemp);

                final double dTotalEnergy = dNeighsE + dStrainE + dNonNearestE;

                String sTempOut = "\t" + i + "\t" + formEnergies.format(dNeighsE * dHartToKJ)
                        + "\t\t\t" + formEnergies.format(dStrainE * dHartToKJ)
                        + "\t\t" + formEnergies.format(dNonNearestE * dHartToKJ)
                        + "\t\t" + formEnergies.format(dTotalEnergy * dHartToKJ);
                llOutput.add(sTempOut);
            }
        }

        // get niche distribution
        if(bGetNiches){
            
            if(pool.getNicheOfIndividualAtPos(0) != null){
            
                llOutput.add("Niche assignments of the pool: (ID-sorted)");
                llOutput.add("\tpool#\tIndivID\tfitness\t\t\tNicheName");
            
                List<String> listNicheNames = new ArrayList<>();
                List<List<Long>> listNicheMembers = new ArrayList<>();
                for(int i = 0; i < iNumberOfStructs; i++){
                    long myID = pool.getIndividualAtPosition(i).getID();
                    final Niche niche = pool.getNicheOfIndividualAtPos(i);
                    if(niche == null){
                        System.err.println("ERROR: No niche info for individual # " + myID + ". Exiting.");
                        System.exit(21);
                    }
                    final String myNicheName = niche.getID();
                    llOutput.add("\t" + i + "\t" + myID + "\t" + daFitnesses[i]+ "\t" + myNicheName);
                    if(listNicheNames.contains(myNicheName)){ 
                        listNicheMembers.get(listNicheNames.indexOf(myNicheName)).add(myID);
                    }else{
                        listNicheNames.add(myNicheName);
                        List<Long> myMembers = new ArrayList<>();
                        myMembers.add(myID);
                        listNicheMembers.add(myMembers);
                    }
                }
            
                llOutput.add("");
            
                llOutput.add("Niche assignments of the pool: (niche-sorted)");
                llOutput.add("NicheName: Population [memberIndivIDs]");
            
                List<String> listNicheNamesSorted = new ArrayList<>(listNicheNames);
                Collections.sort(listNicheNamesSorted);
                for(int i = 0; i < listNicheNames.size() && i < listNicheMembers.size(); i++){
                    final int j = listNicheNames.indexOf(listNicheNamesSorted.get(i));
                    llOutput.add(listNicheNames.get(j) + ": " 
                        + listNicheMembers.get(j).size() + " " + listNicheMembers.get(j).toString());
                }
            
                llOutput.add("");
            } else{
                System.err.println("ERROR: You asked for niche information but at least pool member(0) has none,");
                System.err.println("hence skipping niche analysis...");
            }
        }
        
        /*
         * print the results out
         */
        try{
            OutputPrimitives.writeOut(sOutputFile, llOutput, true);
        } catch(Exception e){
            System.err.println("ERROR: Failure to write results to file. Exiting.");
            e.printStackTrace(System.err);
            System.exit(69);
        }

        System.out.println(org.ogolem.helpers.Fortune.randomFortune());
    }
}
