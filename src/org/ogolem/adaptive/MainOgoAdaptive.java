/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015-2016, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Locale;
import org.ogolem.adaptive.genericfitness.GenericReferencePoint;
import org.ogolem.adaptive.genericfitness.ReferenceGeomData;
import org.ogolem.core.*;
import org.ogolem.generic.GenericSanityCheck;
import org.ogolem.generic.IndividualWriter;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.generichistory.GenericHistoryConfig;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.GenericPoolConfig;
import org.ogolem.generic.threading.GenericOGOLEMOptimization;
import org.ogolem.helpers.Tuple;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.properties.Energy;

/**
 * The main class of the adaptive package designed to do a paramterization.
 * @author Johannes Dieterich
 * @version 2016-07-20
 */
public class MainOgoAdaptive {

    private static final boolean SPACEEVAL = false;

    /**
     * The actual entry point to parametrization.
     * @param args first argument: the ogo file for the parametrization, second argument: the number of threads
     */
    public static void run(final String[] args) {
        
        if(args[0].equalsIgnoreCase("help")){
            System.out.println("This is the threading cluster global optimization functionality.");
            System.out.println("Required arguments:");
            System.out.println(" * the input file (.ogo format)");
            System.out.println(" * the number of threads to be used");
            System.exit(0);
        }
        
        final String configFile = args[0];
        if(!configFile.endsWith(".ogo")){
            System.err.println("ERROR: Please specify a correct input file. Exiting.");
            System.exit(1);
        }

        int noThreads;
        try {
            noThreads = Integer.parseInt(args[1]);
        } catch (Exception e) {
            System.err.println("WARNING: Couldn't read the number of threads that I am supposed to use, "
                    + "setting it to the number of available cores. " + e.toString());
            noThreads = Runtime.getRuntime().availableProcessors();
        }
        
        final Tuple<String,String> dirs = ManipulationPrimitives.outDirAndBaseName(configFile);
        final String outFolder = dirs.getObject1();
        final String baseName = dirs.getObject2();
        
        final String outFile = outFolder + File.separator + baseName + ".out";
        
        // parse config
        GlobalConfig globConf = null;
        try{
            globConf = org.ogolem.core.Input.ConfigureMe(configFile,false);
        } catch(Exception e){
            System.err.append("ERROR: Failure configure ogolem from input.");
            e.printStackTrace(System.err);
            System.exit(1);
        }
        
        globConf.setOutputFolder(outFolder);
        globConf.setOutputFile(outFile);
                
        // IDEALLY, we would only parse adapConf and set what we need ADDITIONALLY (restart?) in GlobalConfig
               
        /*
         * now we have the first branch: either read reference geometries in and/or
         * generate some yourself
         */
        final AdaptiveConf adapConf = globConf.getAdaptiveConf();

        if(adapConf == null){
            System.err.println("ERROR: The adaptive configuration is null, exiting ungracefully.");
            System.exit(1112);
        }
        adapConf.outputFolder = outFolder;

        // check whether all of our references should have the same spins and charges
        if(adapConf.allRefsSameChargesAndSpins){
            final List<GenericReferencePoint<Energy,ReferenceGeomData<Energy,CartesianCoordinates>>> referenceEnergies = adapConf.refPoints.<Energy,ReferenceGeomData<Energy,CartesianCoordinates>>retrieveReferencePointsForName("ENERGY");
            System.err.println("WARNING: allRefsSameChargesAndSpins: This currently only works in cases where we optimize all properties with the same set of reference geom data. Bug the author(s). :-)");
            final short[] spins = referenceEnergies.get(0).getReferenceInputData().c.getAllSpins().clone();
            final float[] charges = referenceEnergies.get(0).getReferenceInputData().c.getAllCharges().clone();
            for(final GenericReferencePoint<Energy,ReferenceGeomData<Energy, CartesianCoordinates>> point : referenceEnergies){
                point.getReferenceInputData().c.setAllSpins(spins);
                point.getReferenceInputData().c.setAllCharges(charges);
            }
        }

        /*
         * dirty little hack....
         */
        AdaptiveParameters pOne;
        if(args.length > 2 && args[2].startsWith("params=")){
            try{
                final String paramFile = args[2].substring(7);
                final String[] pD = Input.ReadFile(paramFile);
                pOne = new AdaptiveParameters(pD, -100);

                // make the fitness function
                org.ogolem.generic.GenericFitnessFunction<Double, AdaptiveParameters> fitFunc = null;
                try{
                    fitFunc = adapConf.getFitnessFunction();
                } catch(Exception e){
                    System.err.println("ERROR: No success in getting the paramter fitness function.");
                    e.printStackTrace(System.err);
                    System.exit(22);
                }
                final double fit = fitFunc.fitness(pOne, true).getFitness();
                System.out.println("Fitness of this parameter set is " + fit);
                System.exit(0);
            } catch(Exception e){
                e.printStackTrace(System.err);
                System.exit(23);
            }
        }


        /*
         * a difficult part: we need to figure out, how many parameters we need in total for all
         * reference structures and create a reference set
         */
        final AdaptiveGateway gate = new AdaptiveGateway(adapConf);

        final String sMethod = adapConf.whichMethod;
        
        // lets first check whether there is any stub file
        final String sStubFile = configFile.substring(0,configFile.indexOf(".")) + "-stub.aux";
        final boolean bStubFileExists = Input.doesFileExist(sStubFile);
        AdaptiveParameters params;
        if(!bStubFileExists){
            params = gate.createInitialParameterStub(adapConf.getAllReferenceGeoms(), sMethod);
        } else{
            try{
                final String[] data = Input.ReadFile(sStubFile);
                params = new AdaptiveParameters(data,-1);
            } catch(Exception e){
                System.err.println("ERROR: Custom stub file unreadable. Using default. " + e.toString());
                params = gate.createInitialParameterStub(adapConf.getAllReferenceGeoms(), sMethod);
            }
        }
        adapConf.paramExample = params.clone();

        // let's first check wether there is any borders file
        final String sBorderFile = configFile.substring(0,configFile.lastIndexOf(".")) + "-borders.aux";
        final boolean bBorderFileExists = Input.doesFileExist(sBorderFile);
        if (!bBorderFileExists){
            // use standard borders
            System.out.println("INFO: No borders auxiliary file exists. Using standard borders.");
            final double[][] daBorders = gate.minMaxBordersForParams(params);
            adapConf.lowerParameterBorder = daBorders[0];
            adapConf.upperParameterBorder = daBorders[1];
        } else {
            // populate the borders customly
            try {
                final double[][] daBorders = Input.readBordersIn(sBorderFile, params.getNumberOfParamters());
                adapConf.lowerParameterBorder = daBorders[0];
                adapConf.upperParameterBorder = daBorders[1];
                System.out.println("INFO: We are operating with custom parameter borders.");
            } catch (Exception e) {
                e.printStackTrace(System.err);
                System.out.println("INFO: We are operating without custom parameter borders.");
                final double[][] daBorders = gate.minMaxBordersForParams(params);
                adapConf.lowerParameterBorder = daBorders[0];
                adapConf.upperParameterBorder = daBorders[1];
            }
        }
        
        if (adapConf.printBorders) {
            System.out.println("INFO: Used borders for this run are:");
            for(int i = 0; i < adapConf.lowerParameterBorder.length; i++){
                System.out.println(" " + adapConf.lowerParameterBorder[i] + " \t " + adapConf.upperParameterBorder[i]);
            }
        }

        
        final double[] lower = adapConf.lowerParameterBorder;
        final double[] upper = adapConf.upperParameterBorder;
        final boolean checkBounds = adapConf.checkOutOfParamBounds;
        
        final GenericSanityCheck<Double,AdaptiveParameters> sanity = new ParameterSanityCheck(checkBounds,lower,upper);
        org.ogolem.generic.GenericFitnessFunction<Double,AdaptiveParameters> fitness = null;
        IndividualWriter<AdaptiveParameters> writer = null;
        try{
            fitness = adapConf.getFitnessFunction();
            writer = adapConf.getWriter(outFolder);
        } catch(Exception e){
            System.err.println("ERROR: Couldn't get the fitness function or the individual writer.");
            e.printStackTrace(System.err);
        }
        final double crossPoss = adapConf.crossPossibility;
        final double mutPoss = adapConf.mutationRatio;
        final boolean printBeforeFitness = adapConf.writeEveryParameterset;
        final int noOfTries = 1;
            
        // map the string to the global optimization algorithm
        final ParamGlobOptFactory globFactory = new ParamGlobOptFactory(sanity,fitness,writer,crossPoss,
            mutPoss,printBeforeFitness, noOfTries, lower, upper,adapConf.noOfGlobalSteps);

        try{
            adapConf.opter = globFactory.translateToGlobOpt(adapConf.globOptString);
        } catch(Exception e){
            System.err.println("ERROR: Couldn't translate the parameter glob opt string to a method.");
            e.printStackTrace(System.err);
        }

        /*
         * initialize the reference parameter locopt (always needed!),
         * very important that this happens now and not earlier!
         */
        org.ogolem.generic.GenericFitnessFunction<Double, AdaptiveParameters> fitFunc = null;
        try{
            fitFunc = adapConf.getFitnessFunction();
        } catch(Exception e){
            System.err.println("ERROR: No success in getting the paramter fitness function.");
            e.printStackTrace(System.err);
            System.exit(22);
        }
        
        assert(fitFunc != null);
        
        /*
         * we do one single energy evaluation to fill any caches that might potentially exist.
         * performance improvement and memory footprint reduction.
         */
        //XXX keyword
        if(true){
            final AdaptiveParameters adapt = new AdaptiveParameters(params);
            adapt.setRandomParams(adapConf.lowerParameterBorder, adapConf.upperParameterBorder);
            adapt.setID(-10);
            final long before = System.currentTimeMillis();
            final double fitne = fitFunc.fitness(adapt, true).getFitness();
            final long after = System.currentTimeMillis();
            System.out.println("INFO: Evaluated one fitness to fill the caches. Value (nonsense!) is: " + fitne + " took " + (after-before) + "ms.");
            final long before2 = System.currentTimeMillis();
            final double fitness2 = fitFunc.fitness(adapt, true).getFitness();
            final long after2 = System.currentTimeMillis();
            System.out.println("INFO: Evaluated one fitness after filling the caches. Value (nonsense!) is: " + fitness2 + " took " + (after2-before2) + "ms.");
            if(after-after2 < 1E-14){
                System.out.println("INFO: No caching induced error detected. Everything normal.");
            } else {
                System.err.println("ERROR: Caching induced error of " + (after-after2) + " detected. This is a problem!");
                System.err.println("ERROR: Aborting the run, disable caching!");
                System.exit(56);
            }
        }
        
        
        // check the input for sanity
        // TODO
        
        // serialize the globconf for restart purposes
	final String configBin = outFolder + "-config.bin";
        try{
            OutputPrimitives.writeObjToBinFile(configBin, globConf);
        } catch(Exception e){
            System.err.println("ERROR: Failure in binary writing the global configuration. " + e.toString());
            e.printStackTrace(System.err);
            System.exit(37);
        }


        final boolean checkMove = true; // XXX support perhaps?
        if(checkMove){
            try{
                ManipulationPrimitives.moveOldFolders(outFolder, 9, System.getProperty("user.dir"));
            } catch (Exception e) {
                System.err.println("ERROR: Failure in initial file/folder operations. " + e.toString());
                e.printStackTrace(System.err);
                System.exit(1);
            }
        }
        
        // create folder
        try{
            OutputPrimitives.createAFolder(outFolder);
        } catch(IOException e){
            System.err.println("ERROR: Couldn't create output folder!");
            e.printStackTrace(System.err);
            System.exit(2);
        }
                
        // assemble some helper objects
        final ParameterReader reader = new ParameterReader();
        GenericHistory<Double,AdaptiveParameters> history = null;
        GenericPool<Double,AdaptiveParameters> pool = null;
        try{
            assert(adapConf != null);
            final GenericHistoryConfig hisConf = adapConf.getGenericHistoryConfig();
            history = GenericHistory.getReference(hisConf);
            final GenericPoolConfig<Double,AdaptiveParameters> poolConf = adapConf.getGenericPoolConfig();
            final AdaptiveParameters example = adapConf.getExample();
            pool = GenericPool.getInstance(poolConf,example);
        } catch(Exception e){
            System.err.println("ERROR: Couldn't instantiate helper objects in the parameter optimizer.");
            e.printStackTrace(System.err);
            System.exit(4);
        }
        
        assert(pool != null);
        assert(history != null);
        
        // start our generic global optimization
        final GenericOGOLEMOptimization<Double,AdaptiveParameters> opter
                = new GenericOGOLEMOptimization<>(adapConf, pool, history, outFile, 
                        outFolder, writer, reader, adapConf.maxTasksToSubmit);
        
        // run it
        final AdaptiveParameters best = opter.globOpt(noThreads);
        
        // be (almost) done ;-)
        
        if(SPACEEVAL){
            // dump results
            final double[][] borders = new double[2][];
            borders[0] = adapConf.lowerParameterBorder;
            borders[1] = adapConf.upperParameterBorder;
            try{
                SearchspaceCapturer.dumpData(borders, adapConf.acceptableFitness);
            } catch(Exception e){
                System.err.println("ERROR: Couldn't plot gathered searchspace data.");
                e.printStackTrace(System.err);
            }
        }

        if(adapConf.structuralDataType() != AdaptiveConf.StructureDataType.Cartesian){
            // yeah, not happening yet
            System.err.println("WARNING: No final energy evaluation possible for non-cartesian structural data type.");
            System.exit(0);
        }
        
        // do a last set of energy evaluations
        System.out.println("INFO: Raw energies");
        System.out.println("INFO: Is for ID " + best.getID()
                + " with fitness " + best.getFitness());
        System.out.println("#No\t Ref energy\t Eq ref \t actual energy");
        final AdaptiveGateway adaptivable = new AdaptiveGateway(adapConf);
        final List<GenericReferencePoint<Energy,ReferenceGeomData<Energy,CartesianCoordinates>>> referenceEnergies = adapConf.refPoints.<Energy,ReferenceGeomData<Energy,CartesianCoordinates>>retrieveReferencePointsForName("ENERGY");
        final double first = referenceEnergies.get(0).getReferenceProperty().getValue();
        
        final List<ReferenceGeomData<Energy,CartesianCoordinates>> referenceData = adapConf.getReferenceGeomDataForEnergy();
        for (int i = 0; i < referenceEnergies.size(); i++) {
            final ReferenceGeomData<Energy,CartesianCoordinates> db = referenceData.get(i);
            final CartesianCoordinates cartes = db.c;
            final double energy = adaptivable.energyOfStructWithParams(cartes, best ,i, db.bonds);
            final double ref = referenceEnergies.get(i).getReferenceProperty().getValue();
            System.out.println(" " + i + "    " + String.format(Locale.US, "%20.9f", ref) + "    " + String.format(Locale.US, "%20.9f",(ref - first)) + "    " + String.format(Locale.US, "%20.9f",energy));
        }
    }
}
