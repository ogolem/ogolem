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
package org.ogolem.core;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import org.ogolem.adaptive.AdaptiveConf;
import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.generic.Configuration;
import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericGlobalOptimization;
import org.ogolem.generic.GenericInitializer;
import org.ogolem.generic.GenericLocOpt;
import org.ogolem.generic.IndividualReader;
import org.ogolem.generic.IndividualWriter;
import org.ogolem.generic.generichistory.GenericHistoryConfig;
import org.ogolem.generic.genericpool.*;
import org.ogolem.generic.threading.GenericGlobOptTask;
import org.ogolem.generic.threading.GenericInitTask;
import org.ogolem.generic.threading.TaskFactory;
import org.ogolem.md.MDConfig;
import org.ogolem.random.Lottery;

/**
 * All configuration data needed for the program to run is in this class. It is
 * reusable for the "constructor madness". Default values are provided.
 * @author Johannes Dieterich
 * @version 2020-04-25
 */
public class GlobalConfig implements Configuration<Molecule,Geometry> {

    // the ID
    private static final long serialVersionUID = (long) 20200425;

    // the pool size: PoolSize=
    int poolSize = 100;

    // whether the cell should grow or not: GrowCell=
    boolean growCell = false;
    
    /*
     * whether or not detailed stats are wanted at the expense of an additional
     * synchronization point (potential slowdown). only applicable for the SMP
     * case.
     * ClusterDetailedStats=
     */
    boolean enableDetailedStats = false;
    
    /**
     * To disable logging to either out or log file.
     * SilentMode=
     */
    boolean silentMode = false;

    /*
     * Chooses a generic parent selector.
     * GeometryChoice=
     */
    String whichGeomChoice = "fitnessrankbased:gausswidth=0.05";

    /**
     * whether we restart or not
     * restart=
     */
    boolean restart = false;
    
    /**
     * How many records shall be cached before the info is flushed into the output file.
     * GeneticRecordBufferSize=
     */
    int geneticRecordsToASCII = 1000;

    /**
     * Debug enables verbosity.
     * DebugLevel=
     */
    public static int DEBUGLEVEL = 0;
    
    /**
     * If the geometries shall be printed after the globopt step, before the locopt.
     * PrintGeomsBeforeLocOpt=
     */
    boolean printBeforeLocOpt = false;
    
    /**
     * Enable some niching.
     * DoGeometryNiching=
     */
    boolean doNiching = false;
    
    /**
     * How many adds till statistics are printed.
     * NichingAddsToStats=
     */
    int addsToNicheStats = 50;
    
    /**
     * Which Nicher should be used. The remaining String
     * WhichClusterNicher=
     * TODO doc
     */
    String nicherString = "";
    
    /**
     * How many individuals are allowed in one niche at maximum.
     * IndividualsPerNicheAtMax=
     */
    int noOfIndividualsPerNicheAtMax = 50;
    
    /**
     * the acceptable energy (only for benchmarking purposes!!!)
     * AcceptableFitness=
     */
    double acceptableFitness = Double.NEGATIVE_INFINITY;

    /*
     * all GeneticHistory related fields can also be made static since they are only used
     * on the main processor
     * GeneticRecordsToSerial=
     */
    int geneticRecordsToSerial = 1000;

    /*
     * how many geometries added till the pool gets serialized again.
     * GeometriesToSerial=
     */
    int geometriesToSerial = 99;

    /**
     * If the pool shall be serialized whenever a new best geometry is found.
     * SerializeAfterNewBest=
     */
    boolean serializePoolAfterBest = true;
    
    /**
     * If every geometry should be written out.
     * WriteEveryGeometry=
     */
    boolean writeEveryGeom = false;
    String geomDumpFolder = "allgeoms";

    /**
     * The statistics object, null by default. Needs to be filled in by the
     * main method before the pool is initialized.
     */
    GlobOptStatistics stats = null;

    /*
     * general configuration fields
     */
    // the output path
    String outputFolder = "default";

    // the output file
    String outputFile = "default" + System.getProperty("file.separator") + "default.out";

    /**
     * should the binary geometries be deleted after reading
     * DeleteOldGeoms=
     */
    boolean deleteOldGeoms = true;

    /**
     * Number of global optimization steps
     * NumberOfGlobIterations=
     */
    int noOfGlobalSteps = 1000;

    /**
     * offset for "naming" the resulting geometries
     * NamingOffset=
     */
    int namingOffsetGeoms = 0;

    /**
     * The possibility for a geometry crossover to take place in a given global step.
     * CrossoverPossibility=
     */
    double crossPossibility = 1.0;
    
    /**
     * The possibility for a molecular crossover to take place in a given global step. Applies only to fully flexible molecules.
     * MolCrossoverPossibility=
     */
    double molCrossPoss = 1.0;
    
    /**
     * The possibility for a geometry mutation to take in a given global step.
     * MutationPossibility=
     */
    double mutatePossibility = 0.05;
    
    /**
     * The possibility for a molecular mutation to take place in a given global step. Applies only to fully flexible molecules.
     * MolMutationPossibility=
     */
    double molMutatePoss = 0.05;
    
    /**
     * Options:
     * - energy diversity check
     * - energy and path based diversity check
     * - path based diversity check
     * - overlap based diversity check
     * - fitness-spiced diversity check
     * DiversityCheck=
     */
    DiversityChecker<Molecule,Geometry> diversityChecker;

    /**
     * Options are:
     * - simple pair wise check
     * - advanced pair wise check
     * - grid collision detection w/ exit on first clash
     * - grid collision detection w/o exit on first clash
     * default: simple pair wise check
     * CollisionDetection=
     */
    CollisionDetection.CDTYPE whichCollisionEngine = CollisionDetection.DEFAULTCD;

    /**
     * Options are:
     * - recursive DFS
     * - Warshall
     * default: dfs
     * DissociationDetection=
     */
    DissociationDetection.DDTYPE whichDissociationEngine = DissociationDetection.DEFAULTDD;

    /**
     * Options are:
     * - packing init
     * - packing init w/o considering size
     * - packing init fully randomized
     * - randomized init w/  DD
     * - randomized init w/o DD
     * default: packing init w/o size
     * InitialFillAlgo=
     */
    GeometryInit.INITSTYLE whichInitialFill = GeometryInit.DEFAULTINIT;

    /*
     * Where to find the seeding geometries. Also works as an identifier whether
     * to use seeding or not.
     * SeedingPath=
     */
    String seedingFolder = "N/A";

    @Override
    public String intermediatePoolFile() {
        return outputFolder + System.getProperty("file.separator") + "IntermediateClusterPool.bin";
    }

    static enum MOLMUTCHOICE{SOME,ONE};
    
     /*
      * which way to mutate molecular coordinates
      * some: mutate a couple of them
      * one: mutate just one of them
      * MolecularCoordinateMutation=
      */
    MOLMUTCHOICE whichMolecularMutation = MOLMUTCHOICE.SOME;
    
    /*
     * how many attempts the global optimization should make before giving up
     * GlobOptTries=
     */
    int howManyTries = 200;
    
    /**
     * Defines whether ogolem checks the sanity of the structure after each local optimization.
     * PostSanityCheck=
     */
    boolean doPostSanityCheck = true;
    
    /**
     * PostSanityCD=
     */
    boolean doPostCD = false;
 
    /**
     * PostSanityDD=
     */
    boolean doPostDD = false;
    
    /**
     * PreSanityCD=
     */
    boolean doPreCD = true;
 
    /**
     * PreSanityDD=
     */
    boolean doPreDD = true;

    /**
     * For setting randomized COMs: the cell size in which they are supposed to be.
     * CellSize=
     */
    double[] maxCellSize = {9, 9, 9};

    /**
     * What ratio of explicit degrees of freedom should be randomized in the
     * initialization per molecule.
     * RatioExplDoFInit=
     */
    float ratioExplDoFInit = 1.0f;

    /**
     * Whether in case of flexible molecules one attempts a collision detection
     * after each trial.
     * MolecularCDInInit=
     */
    boolean doMolecularCD = true;

    /**
     * The blow factor used to tell that there is a bond.
     * BlowBondDetect=
     */
    double blowFacBondDetect = 1.2;

    /**
     * the maximum number of iterations in the (internal) local optimization
     * MaxIterLocOpt=
     */
    int maxIterLocOpt = 2000;

    /**
     * the gradient threshold in the local optimization
     * ThreshLocOptGradient=
     */
    double threshLocOptGradient = 1E-8;

    /**
     * the coordinate threshold in the local optimization
     * ThreshLocOptCoord=
     */
    double threshLocOptCoord = 1E-8;

    /**
     * the maximum step size for the internal local optimization
     * IntLocOptMaxStep=
     */
    double maxStepSizeLocOpt = 1E0;

    /**
     * the blow factor for the dissociation detection
     * BlowFacDissoc=
     */
    double blowFacDissocDetect = 3.0;

    /**
     * the blow factor for the initial bond detection
     * BlowInitialBonds=
     */
    double blowFacInitialBondDetect = blowFacBondDetect;

    /**
     * the blow factor for the environment fitting detection:
     * TODO input keyword
     */
    double blowFacEnvClusterClashes = 0.8;

    /**
     * the maximum stretch factor for a bond
     * MaxBondStretch=
     */
    double maxBondStretch = 1.3;

    /**
     * The reference geometry configuration. Filled in using the GEOMETRY tags.
     */
    GeometryConfig geoConf = new GeometryConfig();

    /**
     * Is this a periodic optimization?
     * XXX keyword and usage
     */
    boolean isPeriodic = false;
    
    /**
     * What are the cell dimensions (as a parallelpiped) of the periodic cell?
     * XXX keyword and parsing
     */
    double[][] periodicCell = new double[3][3];

    /**
     * The number of tasks to submit in one chunk.
     * TasksToSubmit=
     */
    int tasksToSubmit = org.ogolem.generic.threading.GenericOGOLEMOptimization.DEFAULTSUBSTOWAIT;
    
    /**
     * for usage with adaptive methods, we also want/need an adaptive configuration in here
     * by default, it is null'd and just filled if specified in the input
     */
    org.ogolem.adaptive.AdaptiveConf adaptiveConf = null;
    
    /**
     * for usage with molecular dynamics, we also want/need a molecular dynamics
     * configuration. Null by default, filled in if in that mode.
     */
    org.ogolem.md.MDConfig mdConf = null;
    
    /**
     * for usage with adaptive methods, the parameters
     * AdaptiveParameters=
     */
    AdaptiveParameters parameters = null;

    Environment environment = null;

    AdditivePenaltyFunction penalty = null;

    /**
     * Gets automagically filled in during the input parsing.
     * LocOptAlgo=
     */
    GenericLocOpt<Molecule,Geometry> refNewton = null;

    /**
     * If filled in will be used for the initial local optimizations.
     * InitialLocOpt=
     */
    GenericLocOpt<Molecule,Geometry> initNewton = null;

    /**
     * What we optimize actually. Used to contruct a fitness function (potentially using the ref/initNewtons).
     * OptimizationTarget=
     */
    String fitnessFunctionConfig = "energy";
    
    /**
     * Gets automagically filled in during the input parsing.
     * GlobOptAlgo=
     */
    GenericGlobalOptimization<Molecule,Geometry> opter = null;

    HashMap<String,GenericBackend<Molecule,Geometry>> backendDefs = null;

    /*
     * Methods
     */

    /**
     * Creates a printable version of the configuration.
     * @return saContent The configuration as an array of Strings.
     */
    @Override
    public List<String> getFormattedConfig(){
        
        final List<String> configData = new ArrayList<>();

        configData.add("**************************");
        configData.add("");
        configData.add("CONFIGURATION");
        configData.add("");
        configData.add("**************************");
        
        configData.add("");
        configData.add("");

        /*
         * First the cartesians and bonds of all molecules.
         * Offset and endset are needed for the bond matrix printing.
         */
        int iOffSet = 0;
        int iEndSet = 0;

        for(int i = 0; i < geoConf.noOfParticles; i++){
            configData.add("CARTESIANS OF MOLECULE " + i + ":");
            configData.add("");
            String[] saTempCartes = geoConf.geomMCs.get(i).toCartesians().createPrintableCartesians();
            configData.addAll(Arrays.asList(saTempCartes));
            configData.add("");

            configData.add("BOND MATRIX OF MOLECULE " + i + ":");
            iEndSet += geoConf.geomMCs.get(i).noOfAtoms;
                for (int k = iOffSet; k < iEndSet; k++) {

                    configData.add("");

                    final StringBuffer sBuff = new StringBuffer();
                    for (int l = iOffSet; l < iEndSet; l++) {
                        sBuff.append("\t");
                        sBuff.append(geoConf.bonds.bondType(k, l));
                        assert(geoConf.bonds.bondType(k, l) == geoConf.bonds.bondType(l, k));
                    }

                    configData.add(sBuff.toString());
                }

            configData.add("");

            // so that we start next time from the right spot in the bond matrix
            iOffSet = iEndSet;
        }

        // now the degrees of freedom
        configData.add("------------------------------------------------------");
        configData.add("       EXPLICIT INTERNAL DEGREES OF FREEDOM           ");
        configData.add("------------------------------------------------------");
        for(int i = 0; i < geoConf.noOfParticles; i++){
            if(!geoConf.geomMCs.get(i).flexy){
                // molecule is not flexible
                configData.add("molecule " + i + " is not flexible");
            } else{
                // molecule is flexible
                configData.add("molecule " + i + " is flexible");
                configData.add("");
                configData.add("the zmatrix is:");
                String[] saTempZmat = geoConf.geomMCs.get(i).zmat.getPrintableZMatrix();
                configData.addAll(Arrays.asList(saTempZmat));
                configData.add("");

                // checking which coordinates of the molecule are supposed to be flexible
                boolean[][] baFlexCoords = geoConf.geomMCs.get(i).degreesOfFreedom;
                for(int k = 0; k < baFlexCoords.length; k++){
                    for(int l = 0; l< baFlexCoords[0].length; l++){
                        if(baFlexCoords[k][l]){
                            // the coordinate is flexible
                            configData.add("degree of freedom at coordinate " + k + "," + l);
                        } else{
                            // the coordinate is not flexible
                        }
                    }
                }
            }
        }
        
        // now the constraints
        configData.add("");
        configData.add("------------------------------------------------------");
        configData.add("                     CONSTRAINTS                      ");
        configData.add("------------------------------------------------------");
        configData.add("molecule \t constricted? \t atom \t coordinate");
        for(int i = 0; i < geoConf.noOfParticles; i++){
            if(geoConf.geomMCs.get(i).constricted){
                // some constraint on this molecule
                configData.add("\t " + i + " \t constricted");
                final boolean[][] constr = geoConf.geomMCs.get(i).constraints;
                for(int at = 0; at < geoConf.geomMCs.get(i).noOfAtoms; at++){
                    if(constr[0][at]){configData.add("\t\t\t\t" + at + " \t x-coord");}
                    if(constr[1][at]){configData.add("\t\t\t\t" + at + " \t y-coord");}
                    if(constr[2][at]){configData.add("\t\t\t\t" + at + " \t z-coord");}
                }
            } else{
                // no constraint
                configData.add("\t " + i + " \t unconstricted");
            }
        }
                    
        // charges
        configData.add("");
        configData.add("------------------------------------------------------");
        configData.add("                      CHARGES                         ");
        configData.add("------------------------------------------------------");

        float fTotalCharge = 0.0f;
        boolean bAnyCharge = false;
        ArrayList<float[]> alChargeArrays = new ArrayList<>();
        for(int i = 0; i < geoConf.noOfParticles; i++){
            float[] faCharges = geoConf.geomMCs.get(i).charges;
            alChargeArrays.add(faCharges);
            for(int j = 0; j < faCharges.length; j++){
                fTotalCharge += faCharges[j];
                if(faCharges[j] != 0){
                    bAnyCharge = true;
                }
            }
        }

        if(Math.round(fTotalCharge) == 0){
            // uncharged system
            configData.add("the total system is uncharged");
            configData.add("");
            if (bAnyCharge) {
                // charges add up to a total of zero
                configData.add("this stems from the partial charges summing up to zero");
                configData.add("------------------------------------------------------");
                configData.add("                   PARTIAL CHARGES                    ");
                configData.add("------------------------------------------------------");
                configData.add("molecule \t atom \t charge");
                for (int i = 0; i < alChargeArrays.size(); i++) {
                    final float[] faTempCharges = alChargeArrays.get(i);
                    for (int j = 0; j < faTempCharges.length; j++) {
                        if (faTempCharges[j] != 0.0) {
                            // found a charge
                            configData.add("\t " + i + " \t " +j + " \t " + faTempCharges[j]);
                        }
                    }
                }
                configData.add("------------------------------------------------------");
            }
        } else{
            // charged system
            configData.add("the total system is charged with a total of " + fTotalCharge);
            configData.add("------------------------------------------------------");
            configData.add("                   PARTIAL CHARGES                    ");
            configData.add("------------------------------------------------------");
            configData.add("molecule \t atom \t charge");
            for(int i = 0; i < alChargeArrays.size(); i++){
                final float[] faTempCharges = alChargeArrays.get(i);
                for(int j = 0; j < faTempCharges.length; j++){
                    if(faTempCharges[j] != 0.0){
                        // found a charge
                        configData.add("\t " + i + " \t " +j + " \t " + faTempCharges[j]);
                    }
                }
            }
            configData.add("------------------------------------------------------");
        }

        // general settings
        configData.add("");
        configData.add("------------------------------------------------------");
        configData.add("                  GENERAL SETTINGS                    ");
        configData.add("------------------------------------------------------");
        configData.add("");
        
        configData.add("random number generator is: " + Lottery.getInstance().getInformation());

        configData.add("pool size is " + poolSize);
        
        if(restart){
            configData.add("we are doing a restart");
        } else{
            configData.add("we are doing a fresh start");
        }

        configData.add("the number of global iterations is set to " + noOfGlobalSteps);
        configData.add("number of tries in globopt: " + howManyTries);
        
        configData.add("");
        configData.add("parent selector: " + whichGeomChoice);

        configData.add("");
        configData.add("algorithms used:");
        configData.add("*** global optimization: " + opter.getMyID());
        if(initNewton != null){
            configData.add("*** initial local optimization: " + initNewton.getMyID());
        }
        configData.add("*** local optimization: " + refNewton.getMyID());
        configData.add("*** diversity check: " + diversityChecker.getMyName());
        configData.add("*** molecular mutation: " + whichMolecularMutation.name());

        configData.add("");

        configData.add("threshhold settings for local optimization (not all might apply):");
        configData.add("\t in local optimization (coordinates): " + threshLocOptCoord);
        configData.add("\t in local optimization (gradient): " + threshLocOptGradient);
        configData.add("");

        configData.add("maximum number of iterations in local optimization: " + maxIterLocOpt);

        configData.add("");
        configData.add("informations related to dissociation and collision detection:");
        configData.add("\t blow factor in initial bond detection: " + blowFacInitialBondDetect);
        configData.add("\t blow factor in collision detection: " + blowFacBondDetect);
        configData.add("\t blow factor in dissociation detection: " + blowFacDissocDetect);
        configData.add("\t choice of collision detection engine: " + whichCollisionEngine);

        configData.add("");
        if(!Double.isInfinite(acceptableFitness)){
            configData.add("");
            configData.add("ACCEPTABLE FITNESS SET TO " + acceptableFitness);
            configData.add("BE CAREFUL!");
            configData.add("");
        }
        configData.add("**************************");
        configData.add("");
        configData.add("END OF CONFIGURATION");
        configData.add("");
        configData.add("**************************");

        return configData;
    }

    static DiversityChecker<Molecule,Geometry> mapStringToDiversityCheck(final String diverInput) throws Exception {
        
        DiversityChecker<Molecule,Geometry> diver;
        if(diverInput.startsWith("fitnessbased:")){
            final String s = diverInput.substring(13).trim();
            final double thresh = Double.parseDouble(s);
            diver = new GenericDiversityCheckers.FitnessDiversityChecker<>(thresh);
        } else if(diverInput.startsWith("percfitnessbased:")){
            final String s = diverInput.substring(17).trim();
            final double perc = Double.parseDouble(s)/100;
            diver = new GenericDiversityCheckers.PercentageFitnessDiversityChecker<>(perc);
        } else{
            throw new RuntimeException("Wrong input " + diverInput + " for diversity check. Please fix.");
        }
        
        return diver;
    }
    
    @Override
    public GenericPoolConfig<Molecule,Geometry> getGenericPoolConfig() throws Exception {
        
        final GenericPoolConfig<Molecule,Geometry> config = new GenericPoolConfig<>();
        
        config.setDoNiching(doNiching);
        config.setSerializeAfterNewBest(serializePoolAfterBest);
        config.setAcceptableFitness(acceptableFitness);
        config.setAddsToSerial(geometriesToSerial);
        config.setWriteEveryAdd(writeEveryGeom);
        config.setPoolSize(poolSize);
        config.setInterBinFile(intermediatePoolFile());
                
        if(silentMode){
            config.beSilent();
        }
            
        config.setDiversityChecker(diversityChecker);
        
        if(doNiching){
            config.setAddsToStats(this.addsToNicheStats);
            config.setNicher(new SimpleNicher<>(this.noOfIndividualsPerNicheAtMax));
        }
        
        ParentSelector<Molecule,Geometry> selec;
        try{
            selec = GenericParentSelectors.buildSelector(whichGeomChoice);
        } catch(Exception e){
            throw new RuntimeException("Error in creating parent selector.",e);
        }
        config.setSelector(selec);
        
        config.setWriter(new GeometryWriter(this.outputFolder));
        final Path outFolderPath = Paths.get(this.outputFolder);
        final String lastFolder = outFolderPath.getFileName().toString();
        config.setStats(new GenericStatistics(this.outputFolder + File.separator + lastFolder + ".log", geneticRecordsToSerial));

        return config;
    }
    
    private NicheComputer<Molecule,Geometry> mapStringToNicheComp(final String fullNicherString) throws Exception {
        
        
        /*
         * Syntax is as follows: NICHERNAME:NICHEROPT1;NICHEROPT2;NICHEROPT3
         */
        final int indCol = fullNicherString.indexOf(":");
        if(indCol <= 0){throw new RuntimeException("Syntax for nicher is nichername:nicheropt1;nicheropt2;nicheropt3");}
        final String whichNicher = fullNicherString.substring(0,indCol);
        final String nicherOptions = fullNicherString.substring(indCol+1);
        
        System.out.println("which nicher " + whichNicher);
        System.out.println("nicher options " + nicherOptions);
        
        NicheComputer<Molecule,Geometry> nicher;
        switch(whichNicher){
                case "hydrogenbond": nicher =  new HydrogenBondNicheComp(2.5*Constants.ANGTOBOHR, 3.5*Constants.ANGTOBOHR, 45.0); break;
                case "waterangle": 
                    final String angleOpt = nicherOptions.trim();
                    if(!angleOpt.isEmpty()){
                        if(angleOpt.startsWith("nichewidth=")){
                            final int nicheWidth = Integer.parseInt(angleOpt.substring(11).trim());
                            nicher = new WaterAngleNicheComp(nicheWidth);
                        } else{
                            throw new RuntimeException("WaterAngleNiche option " + angleOpt + " unknown!");
                        }
                    } else{
                        nicher = new WaterAngleNicheComp(); 
                    }
                    break;
                case "interiormolecule":
                    final SurfaceDetection.SURFDETECTTYPE type = SurfaceDetection.parseType(nicherOptions.trim());
                    final SurfaceDetectionEngine surfeng = new SurfaceDetection(type);
                    nicher = new InteriorMoleculeNicheComp(surfeng);
                    break;
                case "chained":
                    final String[] nicheStrings = nicherOptions.trim().split("\\|");
                    final List<NicheComputer<Molecule,Geometry>> nichers = new ArrayList<>();
                    for(final String nicheString : nicheStrings){
                        final NicheComputer<Molecule,Geometry> nic = mapStringToNicheComp(nicheString);
                        nichers.add(nic);
                    }
                    nicher = new ChainedNicheComp(nichers);
                    break;
                case "ljprojection":
                    double incrX = LJProjNicheComp.defXIncr;
                    double incrY = LJProjNicheComp.defYIncr;
                    double grid = LJProjNicheComp.defGridIncr;
                    int binStart = LJProjNicheComp.defBinStart;
                    int binWidth = LJProjNicheComp.defBinWidth;
                    
                    final String[] options = nicherOptions.trim().split("\\;");
                    for(final String opt : options){
                        if(opt.isEmpty()){
                            continue;
                        } else if(opt.trim().startsWith("incrx=")){
                            final String incr = opt.substring(6).trim();
                            incrX = Double.parseDouble(incr);
                        } else if(opt.trim().startsWith("incry=")){
                            final String incr = opt.substring(6).trim();
                            incrY = Double.parseDouble(incr);
                        } else if(opt.trim().startsWith("grid=")){
                            final String incr = opt.substring(5).trim();
                            grid = Double.parseDouble(incr);
                        } else if(opt.trim().startsWith("binStart=")){
                            final String incr = opt.substring(9).trim();
                            binStart = Integer.parseInt(incr);
                        } else if(opt.trim().startsWith("binWidth=")){
                            final String incr = opt.substring(9).trim();
                            binWidth = Integer.parseInt(incr);
                        } else {
                            throw new RuntimeException("Illegal option " + opt + " for ljprojection cluster nicher.");
                        }
                    }
                                        
                    final LJProjNicheComp ljproj = new LJProjNicheComp(incrX,incrY,grid,binStart,binWidth);
                    nicher = ljproj;
                    break;
                case "ljneighbor": 
                    int width = LJNeighborNicheComp.defWidth;
                    String mode = LJNeighborNicheComp.defMode;
                    
                    final String[] optionsLJN = nicherOptions.trim().split("\\;");
                    for(final String opt : optionsLJN){
                        if(opt.isEmpty()){
                            continue;
                        } else if(opt.trim().startsWith("mode=")){
                            mode = opt.substring(5).trim();
                        } else if(opt.trim().startsWith("witdh=")){
                            final String incr = opt.substring(6).trim();
                            width=Integer.parseInt(incr);
                        }
                    }
                    
                    final LJNeighborNicheComp ljneigh = new LJNeighborNicheComp(width,mode);
                    nicher = ljneigh; 
                    break;
                    
                default:
                    throw new RuntimeException("No such nicher " + whichNicher + ".");
            }
        
        return nicher;
        
    }

    @Override
    public int getNumberOfGlobalSteps() throws Exception {
        if(noOfGlobalSteps < 0){
            throw new RuntimeException("Number of global steps less than zero. Please specify zero or more steps.");
        }
        
        return noOfGlobalSteps;
    }

    @Override
    public GenericInitializer<Molecule, Geometry> getInitializer() throws Exception {
        
        final GeometryInitialization initer = new GeometryInit(whichInitialFill);
        final GenericFitnessFunction<Molecule,Geometry> fit = org.ogolem.core.FitnessFunctionFactory.build(this,(initNewton == null) ? refNewton : initNewton, fitnessFunctionConfig);
        
        final GeometryInitializationToGenericAdaptor adap = new GeometryInitializationToGenericAdaptor(initer,whichCollisionEngine,
                    whichDissociationEngine, maxCellSize.clone(),
                    blowFacDissocDetect, blowFacBondDetect, ratioExplDoFInit, doMolecularCD, fit);
        
        return adap;
    }

    @Override
    public <V extends GenericInitializer<Molecule, Geometry>> TaskFactory<Molecule, Geometry, V> getInitFactory() throws Exception {
        return new GenericInitTask<>();
    }

    @Override
    public GenericGlobalOptimization<Molecule, Geometry> getGlobalOptimization() throws Exception {
        
        assert(opter != null);
        return opter;
    }

    @Override
    public <V extends GenericGlobalOptimization<Molecule, Geometry>> TaskFactory<Molecule, Geometry, V> getGlobalFactory() throws Exception {
        return new GenericGlobOptTask<>();
    }

    @Override
    public GenericHistoryConfig getGenericHistoryConfig() throws Exception {
        
        final GenericHistoryConfig hisConf = new GenericHistoryConfig();
        hisConf.recordsToASCII = geneticRecordsToASCII;
        hisConf.recordsToSerial = geneticRecordsToSerial;
        hisConf.offset = (restart) ? namingOffsetGeoms + poolSize : poolSize;
        hisConf.binOut = outputFolder + File.separator + "genetic-history.bin";
        hisConf.asciiAppend = outputFile;
        hisConf.silentMode = silentMode;
        
        return hisConf;
    }

    @Override
    public String seedFolder() {
        return (seedingFolder.equalsIgnoreCase("N/A")) ? null : seedingFolder;
    }

    @Override
    public GenericFitnessFunction<Molecule, Geometry> getFitnessFunction() throws Exception {
        
        final GenericFitnessFunction<Molecule,Geometry> fit = (initNewton == null) ? 
                refNewton : initNewton;
        
        return fit;
    }

    @Override
    public Geometry getExample() throws Exception {
        return new Geometry(geoConf);
    }

    @Override
    public IndividualReader<Geometry> getReader() throws Exception {
        return new GeometryReader();
    }

    @Override
    public IndividualWriter<Geometry> getWriter(String folder) throws Exception {
        return new GeometryWriter(folder);
    }

    public void setOutputFolder(final String outputFolder) {
        assert(outputFolder != null);
        assert(!outputFolder.isEmpty());
        this.outputFolder = outputFolder;
    }

    public void setOutputFile(final String outputFile) {
        assert(outputFile != null);
        assert(!outputFile.isEmpty());
        this.outputFile = outputFile;
    }

    public String getOutputFolder() {
        return outputFolder;
    }

    public String getOutputFile() {
        return outputFile;
    }

    public boolean doesRestart() {
        return restart;
    }

    public Newton getRefNewton() {
        
        try{
            final Newton newton = LocOptFactory.translateLocalOpt(refNewton.clone());
            return newton;
        } catch(Exception e){
            e.printStackTrace(System.err);
            return null;
        }
    }

    public AdaptiveConf getAdaptiveConf() {
        return adaptiveConf;
    }

    public MDConfig getMdConf() {
        return mdConf;
    }
    
    public CollisionDetection.CDTYPE getWhichCollisionEngine() {
        return whichCollisionEngine;
    }

    public DissociationDetection.DDTYPE getWhichDissociationEngine() {
        return whichDissociationEngine;
    }
    
    public GeometryConfig geoConfCopy() {
        return geoConf.clone();
    }
    
    public double getBlowFacBondDetect() {
        return blowFacBondDetect;
    }

    public double getBlowFacDissocDetect() {
        return blowFacDissocDetect;
    }
    
    public int getTasksToSubmit(){
        return tasksToSubmit;
    }

    @Override
    public boolean wantsDetailedStats() throws Exception {
        return enableDetailedStats;
    }
    
    @Override
    public NicheComputer<Molecule, Geometry> getNicheComputer() throws Exception {
        if(doNiching){
            return mapStringToNicheComp(this.nicherString);
        }
        
        return null;
    }
}
