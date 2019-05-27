/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015-2017, J. M. Dieterich and B. Hartke
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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import org.ogolem.adaptive.genericfitness.BatchedPropertyCalculator;
import org.ogolem.adaptive.genericfitness.EnergyCalculator;
import org.ogolem.adaptive.genericfitness.FitnessTermConfig;
import org.ogolem.adaptive.genericfitness.GenericFitnessFunction;
import org.ogolem.adaptive.genericfitness.ReferenceGeomData;
import org.ogolem.adaptive.genericfitness.GenericFitnessTerm;
import org.ogolem.adaptive.genericfitness.GenericReferencePoint;
import org.ogolem.adaptive.genericfitness.PropertyCalculator;
import org.ogolem.adaptive.genericfitness.PseudoPropertyCalculator;
import org.ogolem.adaptive.genericfitness.RefBulkModulusData;
import org.ogolem.adaptive.genericfitness.RefCellVolumeData;
import org.ogolem.adaptive.genericfitness.ReferenceDeltaGaugeData;
import org.ogolem.adaptive.genericfitness.ReferenceDensityData;
import org.ogolem.adaptive.genericfitness.ReferenceEnergyOrderData;
import org.ogolem.adaptive.genericfitness.ReferenceForcesData;
import org.ogolem.adaptive.genericfitness.ReferenceInputData;
import org.ogolem.adaptive.genericfitness.ReferencePoint;
import org.ogolem.adaptive.genericfitness.SerialBatchedPropertyCalculator;
import org.ogolem.adaptive.genericfitness.SerialGenericFitnessTerm;
import org.ogolem.properties.Energy;
import org.ogolem.core.*;
import org.ogolem.generic.Configuration;
import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.GenericGlobalOptimization;
import org.ogolem.generic.GenericInitializer;
import org.ogolem.generic.IndividualReader;
import org.ogolem.generic.IndividualWriter;
import org.ogolem.generic.generichistory.GenericHistoryConfig;
import org.ogolem.generic.genericpool.*;
import org.ogolem.generic.threading.GenericGlobOptTask;
import org.ogolem.generic.threading.GenericInitTask;
import org.ogolem.generic.threading.TaskFactory;
import org.ogolem.helpers.Tuple;
import org.ogolem.io.InputPrimitives;
import org.ogolem.properties.BulkModulus;
import org.ogolem.properties.CellVolume;
import org.ogolem.properties.DeltaGauge;
import org.ogolem.properties.Density;
import org.ogolem.properties.EnergyOrder;
import org.ogolem.properties.Forces;
import org.ogolem.properties.Property;

/**
 * A configuration object for the adaptive package.
 * @author Johannes Dieterich
 * @version 2017-11-16
 */
public class AdaptiveConf implements Configuration<Double,AdaptiveParameters> {
    
    public static enum StructureDataType{Cartesian,Periodic};

    // the ID
    private static final long serialVersionUID = (long) 20170914;
    
    public AdaptiveConf(final String[] saAdaptiveContent, final GlobalConfig globconf) throws Exception{

        boolean bReferenceGeoms = false;

        final List<Integer> llReferenceTags = new LinkedList<>();
        final List<Integer> llReferenceCloseTags = new LinkedList<>();

        String sAdaptivable = null;
        
        final List<Tuple<Integer,Integer>> fitFuncDefTags = new LinkedList<>();

        for (int i = 0; i < saAdaptiveContent.length; i++) {
            final String lineAdapConf = saAdaptiveContent[i].trim();

            if(lineAdapConf.startsWith("#") || lineAdapConf.startsWith("//")){
                // well-known comment chars... ;-)
                continue;
            } if (lineAdapConf.equalsIgnoreCase("<REFERENCE>")) {
                llReferenceTags.add(i);
                bReferenceGeoms = true;
                while(i < saAdaptiveContent.length){
                    if (saAdaptiveContent[i].trim().equalsIgnoreCase("</REFERENCE>")) {
                        llReferenceCloseTags.add(i);
                        break;
                    }
                    i++;
                }
            } else if(lineAdapConf.startsWith("<PARAMFITNESSFUNCTION>")){
                final int start = i;
                while(i < saAdaptiveContent.length){
                    if (saAdaptiveContent[i].trim().equalsIgnoreCase("</PARAMFITNESSFUNCTION>")) {
                        final int end = i;
                        final Tuple<Integer,Integer> tags = new Tuple<>(start,end);
                        fitFuncDefTags.add(tags);
                        break;
                    }
                    i++;
                }
            } else if(lineAdapConf.startsWith("<FITNESSFUNCTIONTERMS>")){
                final int start = i;
                int end = -1;
                while(i < saAdaptiveContent.length){
                    if (saAdaptiveContent[i].trim().equalsIgnoreCase("</FITNESSFUNCTIONTERMS>")) {
                        end = i;
                        break;
                    }
                    i++;
                }
                this.allTermConfigs = new String[end-start-1];
                int c = 0;
                for(int x = start+1; x < end; x++){
                    allTermConfigs[c] = saAdaptiveContent[x];
                    c++;
                }
            } else if (lineAdapConf.startsWith("PopulationSize=")) {
                String sTemp2 = lineAdapConf.substring(15).trim();
                try {
                    int iPopSize = Integer.parseInt(sTemp2);
                    populationSize = iPopSize;
                } catch (Exception e) {
                    System.err.println("Wrong input for the parameter population size. Using default. "
                            + e.toString());
                }
            } else if(lineAdapConf.startsWith("ParameterDetailedStats=")){
                try{
                    enableDetailedStats = Boolean.parseBoolean(lineAdapConf.substring(23));
                } catch(Exception e){
                    System.err.println("Wrong input for ParameterDetailedStats: " + e.toString() + ". NOT ENABLED!");
                }
            } else if (lineAdapConf.startsWith("ParamSeedingFolder=")) {
                final String s = lineAdapConf.substring(19).trim();
                paramSeedFolder = s;
            } else if (lineAdapConf.startsWith("ParamGlobOptIter=")) {
                String sTemp2 = lineAdapConf.substring(17).trim();
                try {
                    int iGlobIter = Integer.parseInt(sTemp2);
                    noOfGlobalSteps = iGlobIter;
                } catch (Exception e) {
                    System.err.println("Wrong input for the parameter globopt iterations. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParamsToSerial=")) {
                final String sTemp2 = lineAdapConf.substring(15).trim();
                try {
                    final int iParamsToSerial = Integer.parseInt(sTemp2);
                    paramsToSerial = iParamsToSerial;
                } catch (Exception e) {
                    System.err.println("Wrong input for the parameter to serial option. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParamSerializeAfterNewBest=")) {
                final String s2 = lineAdapConf.substring(27).trim();
                try {
                    boolean b = Boolean.parseBoolean(s2);
                    serializeAfterNewBest = b;
                } catch (Exception e) {
                    System.err.println("Wrong input for ParamSerializeAfterNewBest=. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("DoNiching=")) {
                String sTemp2 = lineAdapConf.substring(10).trim();
                try {
                    boolean b = Boolean.parseBoolean(sTemp2);
                    doNiching = b;
                } catch (Exception e) {
                    System.err.println("Wrong input for niching. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("WhichNicher=")) {
                final String[] sa = lineAdapConf.substring(12).trim().split("\\;");
                try {
                    whichNicher = Integer.parseInt(sa[0]);
                    nicherString = (sa.length > 1) ? sa[1] : "";
                } catch (Exception e) {
                    System.err.println("Wrong input for nicher. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("WhichDiversityCheck=")) {
                final String sTemp2 = lineAdapConf.substring(20).trim();
                try {
                    final int iTemp = Integer.parseInt(sTemp2);
                    whichDiversityCheck = iTemp;
                } catch (Exception e) {
                    System.err.println("Wrong input for diversity check. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("NichesPerDim=")) {
                final String sTemp2 = lineAdapConf.substring(13).trim();
                try {
                    final int iTemp = Integer.parseInt(sTemp2);
                    noOfNichesPerDim = iTemp;
                } catch (Exception e) {
                    System.err.println("Wrong input for number of niches per dimension. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("MaxIndividualsPerNiche=")) {
                final String sTemp2 = lineAdapConf.substring(23).trim();
                try {
                    final int iTemp = Integer.parseInt(sTemp2);
                    noOfIndividualsPerNicheMax = iTemp;
                } catch (Exception e) {
                    System.err.println("Wrong input for number of individuals per niche. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParamBorderPrint=")) {
                final String sTemp2 = lineAdapConf.substring(17).trim();
                try {
                    final boolean b = Boolean.parseBoolean(sTemp2);
                    printBorders = b;
                } catch (Exception e) {
                    System.err.println("Wrong input for ParamBorderPrint=. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParamThreshDiv=")) {
                String sTemp2 = lineAdapConf.substring(15).trim();
                try {
                    final double d = Double.parseDouble(sTemp2);
                    paramDiversityThresh = d;
                } catch (Exception e) {
                    System.err.println("Wrong input for the parameter diversity threshhold. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParamEnergyDiv=")) {
                String sTemp2 = lineAdapConf.substring(15).trim();
                try {
                    double dEnergDiv = Double.parseDouble(sTemp2);
                    energyDiversity = dEnergDiv;
                } catch (Exception e) {
                    System.err.println("Wrong input for the parameter energy diversity. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParamsToBorders=")) {
                final String sTemp2 = lineAdapConf.substring(16).trim();
                try {
                    final boolean b = Boolean.parseBoolean(sTemp2);
                    normParams = b;
                } catch (Exception e) {
                    System.err.println("Wrong input for ParamsToBorders. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParamThreshLocOptGradient=")) {
                final String sTemp2 = lineAdapConf.substring(26).trim();
                try {
                    final double d = Double.parseDouble(sTemp2);
                    threshLocOptGradient = d;
                } catch (Exception e) {
                    System.err.println("Wrong input for ParamThreshLocOptGradient. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParamThreshLocOptParam=")) {
                final String sTemp2 = lineAdapConf.substring(23).trim();
                try {
                    final double d = Double.parseDouble(sTemp2);
                    threshLocOptParam = d;
                } catch (Exception e) {
                    System.err.println("Wrong input for ParamThreshLocOptParam. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParamLinesearchThreshParam=")) {
                final String sTemp2 = lineAdapConf.substring(27).trim();
                try {
                    final double d = Double.parseDouble(sTemp2);
                    linesearchThreshParams = d;
                } catch (Exception e) {
                    System.err.println("Wrong input for ParamLinesearchThreshParam. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParamLocOptMaxStep=")) {
                final String sTemp2 = lineAdapConf.substring(19).trim();
                try {
                    final double d = Double.parseDouble(sTemp2);
                    maxStepSizeLinesearch = d;
                } catch (Exception e) {
                    System.err.println("Wrong input for ParamLocOptMaxStep. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParamLinesearchEnergyDecr=")) {
                final String sTemp2 = lineAdapConf.substring(26).trim();
                try {
                    final double d = Double.parseDouble(sTemp2);
                    linesearchEnergyDec = d;
                } catch (Exception e) {
                    System.err.println("Wrong input for ParamLinesearchEnergyDecr. Using default. "
                            + e.toString());
                }
            } else if (lineAdapConf.startsWith("AllRefsSameChargesAndSpins=")) {
                String sTemp2 = lineAdapConf.substring(27).trim();
                try {
                    boolean b = Boolean.parseBoolean(sTemp2);
                    allRefsSameChargesAndSpins = b;
                } catch (Exception e) {
                    System.err.println("Wrong input in AllRefsSameChargesAndSpins= . Using default. " + e.toString());
                }
            } else if (lineAdapConf.startsWith("AnyParamHistory=")) {
                final String sTemp2 = lineAdapConf.substring(16);
                try {
                    boolean bHistory = Boolean.parseBoolean(sTemp2);
                    anyHistory = bHistory;
                } catch (Exception e) {
                    System.err.println("Wrong input in AnyParamHistory= . Using default. " + e.toString());
                }
            } else if (lineAdapConf.startsWith("MaxTasksToSubmit=")) {
                final String sTemp2 = lineAdapConf.substring(17).trim();
                try {
                    final int iMaxTasks = Integer.parseInt(sTemp2);
                    maxTasksToSubmit = iMaxTasks;
                } catch (Exception e) {
                    System.err.println("Wrong input in MaxTasksToSubmit= . Using default. " + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParameterLocOpt=")) {
                final String sTemp2 = lineAdapConf.substring(16).trim();
                whichParamLocOpt = sTemp2;
            } else if (lineAdapConf.startsWith("AdaptivableChoice=")) {

                sAdaptivable = lineAdapConf.substring(18).trim();

            } else if (lineAdapConf.startsWith("DoubleMorseDMax=")) {
                String sWhichDMax = lineAdapConf.substring(16).trim();
                try {
                    double dMax = Double.parseDouble(sWhichDMax);
                    maxDDoubleMorse = dMax;
                } catch (Exception e) {
                    System.err.println("WARNING: Couldn't cast the double choice"
                            + " of DoubleMorseDMax=. Using default. " + e.toString());
                }
            } else if (lineAdapConf.startsWith("AcceptableFitness=")) {
                final String sAccFit = lineAdapConf.substring(18).trim();
                try {
                    double dAccFit = Double.parseDouble(sAccFit);
                    acceptableFitness = dAccFit;
                } catch (Exception e) {
                    System.err.println("WARNING: Couldn't cast the double choice"
                            + " of DoubleMorseDMax=. Using default. " + e.toString());
                }
            } else if (lineAdapConf.startsWith("ParameterGlobOpt=")) {
                globOptString = lineAdapConf.substring(17);
            } else if (lineAdapConf.startsWith("MaxNonConvSteps=")) {
                String sTemp2 = lineAdapConf.substring(16).trim();

                try {
                    int iAmount = Integer.parseInt(sTemp2);
                    maxNonConvSteps = iAmount;
                } catch (Exception e) {
                    System.err.println("WARNING: Couldn't cast the integer choice"
                            + " of MaxNonConvSteps=. Using default. " + e.toString());
                }
            } else if (lineAdapConf.startsWith("UsePercentageScaledDiffs=")) {
                String sTemp2 = lineAdapConf.substring(25).trim();

                try {
                    boolean bUsePercentages = Boolean.parseBoolean(sTemp2);
                    bUsePercentageDiffs = bUsePercentages;
                } catch (Exception e) {
                    System.err.println("WARNING: Couldn't cast the boolean choice"
                            + " of UsePercentageScaledDiffs=. Using default. " + e.toString());
                }
            } else if (lineAdapConf.startsWith("PercentageForDiffs=")) {
                String sTemp2 = lineAdapConf.substring(19).trim();

                try {
                    double dPerc = Double.parseDouble(sTemp2);
                    percentageForDiff = dPerc;
                } catch (Exception e) {
                    System.err.println("WARNING: Couldn't cast the double choice"
                            + " of PercentageForDiffs=. Using default. " + e.toString());
                }
            } else if (lineAdapConf.equalsIgnoreCase("PrintFitnessContributions")){
                printContributions = true;
            } else {
                throw new RuntimeException("Unknown line " + lineAdapConf + " exiting.");
            }
        }

        if (sAdaptivable == null) {
            throw new InitIOException("Adaptivable not set.");
        }
        try {
            org.ogolem.helpers.Tuple<org.ogolem.adaptive.Adaptivable, String> tupel
                    = mapStringToAdaptivable(sAdaptivable, globconf);
            whichMethod = tupel.getObject2();
            refAdaptivable = tupel.getObject1();
        } catch (Exception e) {
            throw new CastException("Error in adaptivable parsing.", e);
        }

        if (bReferenceGeoms) {

            if (llReferenceTags.size() != llReferenceCloseTags.size()) {
                System.err.println("ERROR: Detected unequal amounts of closing and opening <REFERENCE> tags!");
                throw new CastException("ERROR: Detected unequal amounts of closing and opening <REFERENCE> tags!");
            }

            for (int tags = 0; tags < llReferenceTags.size(); tags++) {

                final int refStart = llReferenceTags.get(tags);
                final int refEnd = llReferenceCloseTags.get(tags);

                // copy the content to a new string[]
                final String[] refTagContent = new String[refEnd - refStart - 1];

                for (int i = refStart + 1; i < refEnd; i++) {
                    refTagContent[i - refStart - 1] = saAdaptiveContent[i].trim();
                }

                // might throw an exception if parsing fails
                try {
                    addReferenceInformationsCartesians(refTagContent, tags);
                } catch (Exception e) {
                    throw new InitIOException("Error in parsing references.", e);
                }
            }
            // put it into the adaptive configuration
        }

        /*
         * Using this has the advantage that the modelling of references
         * with smaller relative energies is more important than bigger
         * ones. This is rather important for an exact globopt modelling.
         */
        if (bUsePercentageDiffs) {
            // use percentage scaled maximal allowed differences
            for (int i = 0; i < referenceEnergies.size(); i++) {
                final double dRefEnergy = referenceEnergies.get(i).getReferenceProperty().getValue();
                
                final double maxAllowedDiff = (Math.abs(dRefEnergy) <= 1E-5) ? 1E-6 : Math.abs(percentageForDiff * 0.01 * dRefEnergy); // first: cutoff
                final GenericReferencePoint<Energy,ReferenceGeomData<Energy,?>> point = referenceEnergies.get(i);
                final GenericReferencePoint<Energy,ReferenceGeomData<Energy,?>> newPoint = new ReferencePoint<Energy,ReferenceGeomData<Energy,?>>(point.getReferenceProperty(),
                        point.getReferenceInputData(),point.getReferenceID(),point.getRefWeight(),maxAllowedDiff);
                // replace old with new
                referenceEnergies.set(i, newPoint);
            }
        }
        
        // construct all the fitness functions, since these are embedded into tags, it makes adding more involved stuff easier, I hope
        if(!fitFuncDefTags.isEmpty()){
        for(final Tuple<Integer,Integer> t : fitFuncDefTags){
            
            String tag = null;
            GenericFitnessFunction fitFunc = null;
            String fitFuncStyle = "full";
            
            for(int x = t.getObject1()+1; x < t.getObject2(); x++){
                final String line = saAdaptiveContent[x].trim();
                if(line.startsWith("#") || line.startsWith("//")){
                    continue;
                } else if (line.startsWith("FitFunctionTag=")){
                    tag = line.substring(15).trim();
                } else if(line.startsWith("FitnessFunction=")){
                    
                    final String def = line.substring(16).trim();
                    switch(def){
                        case "default":
                            fitFunc = getCartesianFitnessFunc();
                            break;
                        default:
                            throw new RuntimeException("Unknown fitness function " + def);
                    }
                } else if(line.startsWith("FitnessFunctionStyle=")){
                    fitFuncStyle = line.substring(21).trim();
                } else {
                    throw new RuntimeException("Unknown line: " + line + " in fitness function definition block.");
                }
            }
            
            if(tag == null){throw new RuntimeException("No tag specified in fitness function definition block.");}
            if(fitFunc == null){throw new RuntimeException("No fitness function specified in fitness function definition block.");}
            
            GenericBackend<Double,AdaptiveParameters> func = null;
            if(fitFuncStyle.equalsIgnoreCase("full")){
                func = new FitFuncToBackend(fitFunc, normParams, lowerParameterBorder, upperParameterBorder);
            } else if(fitFuncStyle.startsWith("ranged:")){
                
                final String rangeDef = fitFuncStyle.substring(7).trim();
                final String[] rangesDef = rangeDef.split("\\/");
                final List<RangedFitFuncToBackend.Range> ranges = new ArrayList<>();
                for(final String rad : rangesDef){
                    if(rad.contains("-")){
                        final String[] def = rad.split("\\-");
                        if(def.length != 2){
                            throw new RuntimeException("Range definition wrong. Must be X-Y (both inclusive).");
                        }
                        final int start = Integer.parseInt(def[0].trim());
                        final int end = Integer.parseInt(def[1].trim());
                        final RangedFitFuncToBackend.Range range = new RangedFitFuncToBackend.Range(start, end);
                        ranges.add(range);
                    } else {
                        final int param = Integer.parseInt(rad.trim());
                        final RangedFitFuncToBackend.Range range = new RangedFitFuncToBackend.Range(param, param);
                        ranges.add(range);
                    }
                }
                
                func = new RangedFitFuncToBackend(fitFunc, normParams, ranges, lowerParameterBorder, upperParameterBorder);
            } else {
                throw new RuntimeException("Unknown fitness function style: " + fitFuncStyle);
            }
            
            backendDict.put(tag, func);
        }
        }
    }

    /**
     * the poolsize of parameters
     * PopulationSize=
     */
    int populationSize = 100;

    /**
     * the number of global optimization iterations in the parameter fit
     * ParamGlobOptIter=
     */
    int noOfGlobalSteps = 9900;

    /**
     * the energy diversity supposed to be maintained in the pool
     * ParamEnergyDiv=
     */
    double energyDiversity = 1E-8;

    /**
     * the threshhold for the parameter diversity (used in case of the corresponding check)
     * ParamThreshDiv=
     */
    double paramDiversityThresh = 1E-8;

    /**
     * Which diversity check to use.
     * WhichDiversityCheck=
     */
    int whichDiversityCheck = 0;

    /**
     * How many parameter sets get added to the pool before serializing it.
     * ParamsToSerial=
     */
    int paramsToSerial = 10000;

    /**
     * Defines what an acceptable fitness is
     * AcceptableFitness=
     */
    double acceptableFitness = 0.0;

    /**
     * Whether one wants to do niching or not.
     * DoNiching=
     */
    boolean doNiching = false;

    /**
     * Which nicher to take.
     * WhichNicher=
     */
    int whichNicher = 0;
    
    /**
     * Which Nicher should be used. The remaining String
     * WhichNicher=
     * TODO doc
     */
    String nicherString = "";

    /**
     * How many niches per dimension.
     * NichesPerDim=
     */
    int noOfNichesPerDim = 10;

    /**
     * How many individuals per niche at the maximum.
     * MaxIndividualsPerNiche=
     */
    int noOfIndividualsPerNicheMax = 50;
    
    /**
     * Where intermediate parameters should be dumped to.
     * TODO keyword and doc
     */
    String paramsDumpFolder = "allparams";
    
    /**
     * Print the borders at the intial stage.
     * ParamBorderPrint=
     * TODO doc
     */
    boolean printBorders = true;
    
    /*
     * whether or not detailed stats are wanted at the expense of an additional
     * synchronization point (potential slowdown). only applicable for the SMP
     * case.
     * ParameterDetailedStats=
     */
    boolean enableDetailedStats = false;
    
    /**
     * TODO keyword and doc
     */
    int addsToNicheStats = 100000;
    
    /**
     * TODO keyword and doc
     */
    String whichParentsChoice = "fitnessrankbased:gausswidth=0.05";
    
    /**
     * TODO doc
     * ParamSerializeAfterNewBest=
     */
    boolean serializeAfterNewBest = true;

    /**
     * TODO keyword and doc
     */
    boolean writeEveryParameterset = false;
    
    /**
     * highly double morse specific: the maximum d factor
     * DoubleMorseDMax=
     */
    double maxDDoubleMorse = 0.3381355712;

    /**
     * for the output. is supposed to be automatically set by the configurating method.
     */
    String whichMethod = "UNSET";


    /**
     * Whether the parameters should be put back into their borders in the local
     * optimizations.
     * ParamsToBorders=
     */
    boolean normParams = false;

    /**
     * Defines if there should be any genetic history. Since this is a potential
     * cause of a high memory need, it is disabled by default.
     * AnyParamHistory=
     */
    boolean anyHistory = false;
    
    /**
     * TODO keyword and doc
     */
    double mutationRatio = 0.05;
    
    /**
     * TODO keyword and doc
     */
    double crossPossibility = 1.0;

    /**
     * for the energy property
     */
    List<GenericReferencePoint<Energy,ReferenceGeomData<Energy,?>>> referenceEnergies = new ArrayList<>();
    
    /**
     * for the bulk modulus property
     */
    List<GenericReferencePoint<BulkModulus,RefBulkModulusData<?>>> referenceBulkMod = new ArrayList<>();
    
    /**
     * for the forces property
     */
    List<GenericReferencePoint<Forces,ReferenceForcesData<?>>> referenceForces = new ArrayList<>();
    
    /**
     * for the cell volume property
     */
    List<GenericReferencePoint<CellVolume,RefCellVolumeData<?>>> referenceCellVolume = new ArrayList<>();
    
    /**
     * for the energy order property
     */
    List<GenericReferencePoint<EnergyOrder,ReferenceEnergyOrderData<?>>> referenceEnergyOrder = new ArrayList<>();
    
    /**
     * for the delta gauge property
     */
    List<GenericReferencePoint<DeltaGauge,ReferenceDeltaGaugeData<?>>> referenceDeltaGauge = new ArrayList<>();
    
    /**
     * for the density property
     */
    List<GenericReferencePoint<Density,ReferenceDensityData<?>>> referenceDensity = new ArrayList<>();
    
    //XXX PUT PROPERTIES HERE

    /**
     * Whether all reference geometries shall have the same charges and spins,
     * defined by the first reference.
     * AllRefsSameChargesAndSpins=
     */
    boolean allRefsSameChargesAndSpins = false;

    /**
     * If we want to use percentage scaled differences for the maximum allowed
     * differences.
     * Attention: If this is set, the reference data should be 1-body free!
     * UsePercentageScaledDiffs=
     */
    boolean bUsePercentageDiffs = false;

    /**
     * The percentage for the scaled maximum differences.
     * PercentageForDiffs=
     */
    double percentageForDiff = 1;


    /*
     * Boundaries for the parameters. Supposed to be set by the adaptive locopt/
     * backend.
     */
    double[] lowerParameterBorder = null;

    double[] upperParameterBorder = null;

    String[] allTermConfigs = new String[0];
    
    /*
     *
     * 
     * configuration options for the internal local optimization engines
     *
     *
     */
    
    /**
     * maximum number of iterations in the local optimization before it
     * exits the locopt gracefully
     * ParamLocOptIter=
     */
    int maxIterLocOpt = 5000;

    /**
     * threshold for the parameter convergence
     * ParamThreshLocOptParam=
     */
    double threshLocOptParam = 1E-8;

    /**
     * threshold for the gradient convergence
     * ParamThreshLocOptGradient=
     */
    double threshLocOptGradient = 1E-8;

    /**
     * threshold for parameter convergence in linesearch
     * ParamLinesearchThreshParam=
     */
    double linesearchThreshParams=1E-8;

    /**
     * maximum stepsize in linesearch
     * ParamLocOptMaxStep=
     */
    double maxStepSizeLinesearch=1E-5;

    /**
     * threshold for energy decrease in linesearch
     * ParamLinesearchEnergyDecr=
     */
    double linesearchEnergyDec=1E-10;

    /**
     * how many non-converged steps in local optimization are allowed before
     * backing out
     * MaxNonConvSteps=
     */
    int maxNonConvSteps = 100;

    /**
     * how many tasks we submit to the threadpool before waiting for their finishing
     * MaxTasksToSubmit=
     */
    int maxTasksToSubmit = org.ogolem.generic.threading.GenericOGOLEMOptimization.DEFAULTSUBSTOWAIT;

    /**
     * A reference adaptivable which is automatically set during the
     * input reading.
     * AdaptivableChoice=
     */
    Adaptivable refAdaptivable = null;

    /**
     * Which local optimization to be used for the parameters.
     * Defaults to lbfgs.
     */
    String whichParamLocOpt = "lbfgs:backend=kakapo";

    /**
     * A reference object for the parameter locopt. Needs to be de-null'd.
     */
    org.ogolem.generic.GenericFitnessFunction<Double, AdaptiveParameters> refNewton = null;
    
    /**
     * The seeding folder. ParamSeedingFolder=
     */
    String paramSeedFolder = null;
    
    /**
     * If we should check the parameters for out-of-bounds cases and set those to a high fitness.
     * TODO keyword, doc
     */
    boolean checkOutOfParamBounds = false;
    
    String outputFolder = null;
        
    AdaptiveParameters paramExample = null;
    
    String globOptString = "NO GLOBOPT SPECIFIED";
    GenericGlobalOptimization<Double, AdaptiveParameters> opter = null;
    
    /**
     * Debug option to print fitness contributions.
     * PrintFitnessContributions
     * TODO doc
     */
    boolean printContributions = false;
    
    private HashMap<String,GenericBackend<Double,AdaptiveParameters>> backendDict = null;
    
    public static org.ogolem.helpers.Tuple<Adaptivable,String> mapStringToAdaptivable(final String s, final GlobalConfig globConf)
        throws Exception {

        Adaptivable adapt = null;
        String sWhichMethod;

        if (s.startsWith("orca:")) {
            String sTemp = s.substring(5).trim();
            sWhichMethod = sTemp;
            if (sWhichMethod.equalsIgnoreCase("AM1")) {
                adapt = new AdaptiveOrcaCaller(globConf.getBlowFacBondDetect(), false);
            } else {
                System.err.println("WARNING: No valid adaptivable choice found in selection orca."
                        + "Using orca:am1 instead.");
                sWhichMethod = "AM1";
                adapt = new AdaptiveOrcaCaller(globConf.getBlowFacBondDetect(), false);
            }
        } else if (s.startsWith("mopac:")) {
            String sTemp = s.substring(6).trim();
            sWhichMethod = sTemp;
            if (sWhichMethod.equalsIgnoreCase("AM1")) {
                adapt = new AdaptiveMopacCaller(true,false, globConf.getBlowFacBondDetect(), true);
            } else if (sWhichMethod.equalsIgnoreCase("azobenzene")) {
                adapt = new AdaptiveMopacCaller(true,true, globConf.getBlowFacBondDetect(), true);
            } else {
                System.err.println("WARNING: No valid adaptivable choice found in selection mopac."
                        + "Using Mopac:am1 instead.");
                sWhichMethod = "AM1";
                adapt = new AdaptiveMopacCaller(true,false, globConf.getBlowFacBondDetect(), true);
            }
        } else if (s.equalsIgnoreCase("adaptiveUFF")) {
            sWhichMethod = "AdaptiveUFF";
            adapt = new AdaptiveUFF(true);
        } else if (s.startsWith("doublemorse")) {
            if (s.endsWith("symmetric")) {
                adapt = new DoubleMorse(true,globConf.getAdaptiveConf().maxDDoubleMorse);
                sWhichMethod = "doublemorse:symmetric";
            } else {
                adapt = new DoubleMorse(false,globConf.getAdaptiveConf().maxDDoubleMorse);
                sWhichMethod = "doublemorse";
            }
        } else if (s.startsWith("adaptiveLJFF")) {
            if (s.endsWith("easy;6,12,6")) {
                adapt = new AdaptiveLJFF(true,true,6,12,6,null,false,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:easy;6,12,6";
            } else if (s.endsWith("easy;2,12,2")) {
                adapt = new AdaptiveLJFF(true,true,2,12,2,null,false,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:easy;2,12,2";
            } else if (s.endsWith("easy;6,16,2")) {
                adapt = new AdaptiveLJFF(true,true,6,16,2,null,false,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:easy;6,16,2";
            } else if (s.endsWith("nomix;6,12,6")) {
                adapt = new AdaptiveLJFF(true,false,6,12,6,null,false,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:nomix;6,12,6";
            } else if (s.endsWith("nomix;2,12,2")) {
                adapt = new AdaptiveLJFF(true,false,2,12,2,null,false,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:nomix;2,12,2";
            } else if (s.endsWith("nomix;6,16,2")) {
                adapt = new AdaptiveLJFF(true,false,6,16,2,null,false,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:nomix;6,16,2";
            } else if (s.endsWith("easy;6,12,6;3body")) {
                adapt = new AdaptiveLJFF(true,true,6,12,6,null,true,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:easy;6,12,6;3body";
            } else if (s.endsWith("easy;2,12,2;3body")) {
                adapt = new AdaptiveLJFF(true,true,2,12,2,null,true,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:easy;2,12,2;3body";
            } else if (s.endsWith("easy;6,16,2;3body")) {
                adapt = new AdaptiveLJFF(true,true,6,16,2,null,true,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:easy;6,16,2;3body";
            } else if (s.endsWith("nomix;6,12,6;3body")) {
                adapt = new AdaptiveLJFF(true,false,6,12,6,null,true,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:nomix;6,12,6;3body";
            } else if (s.endsWith("nomix;2,12,2;3body")) {
                adapt = new AdaptiveLJFF(true,false,2,12,2,null,true,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:nomix;2,12,2;3body";
            } else if (s.endsWith("nomix;6,16,2;3body")) {
                adapt = new AdaptiveLJFF(true,false,6,16,2,null,true,1.2,20.0,false,false);
                sWhichMethod = "adaptiveLJFF:nomix;6,16,2;3body";
            } else {
                boolean mix;
                int start;
                int end;
                int incr;
                boolean threeB;
                double blowClose;
                double blowDist;
                boolean cache;
                boolean disc2body;
                try{
                    final String s2 = s.substring(12).trim();
                    final String[] sa = s2.split("\\;");
                    mix = Boolean.parseBoolean(sa[0].trim());
                    final String[] sa2 = sa[1].split("\\,");
                    start = Integer.parseInt(sa2[0].trim());
                    end = Integer.parseInt(sa2[1].trim());
                    incr = Integer.parseInt(sa2[2].trim());
                    threeB = Boolean.parseBoolean(sa[2].trim());
                    blowClose = Double.parseDouble(sa[3].trim());
                    blowDist = Double.parseDouble(sa[4].trim());
                    cache = Boolean.parseBoolean(sa[5].trim());
                    disc2body = Boolean.parseBoolean(sa[6].trim());
                    adapt = new AdaptiveLJFF(true,mix,start,end,incr,null,threeB,blowClose,blowDist,cache,disc2body);
                } catch(Exception e){
                    throw new Exception("ERROR: Couldn't parse adaptiveLJ Input.", e);
                }
                sWhichMethod = s;
            }
        } else if(s.startsWith("adaptivegupta")){
            final String s2 = s.substring(14).trim();
            final String[] sa = s2.split("\\,");
            try{
                final double distBlow = Double.parseDouble(sa[0]);
                final boolean cache = Boolean.parseBoolean(sa[1]);
                int expF = 0;
                if(sa.length == 3) {expF = Integer.parseInt(sa[2]);}
                adapt = new AdaptiveGUPTA(true, null, distBlow, cache, expF);
            } catch(Exception e){
                System.err.println("ERROR: Couldn't parse required information for adaptivegupta. " + e.toString());
                return null;
            }
            sWhichMethod = s;
        } else if(s.startsWith("adaptivemorse:")){
            adapt = new AdaptiveMorse(true,s,null);
            sWhichMethod = "adaptivemorse";
        } else if (s.startsWith("amberff")){
            // XXX this parsing is crap!
            boolean useCache = true;
            int whichAcos = 0;
            int whichCos = 0;
            int whichSin = 0;
            boolean totalShift = false;
            if(s.contains("cos=0")) {whichCos=0;}
            if(s.contains("cos=1")) {whichCos=1;}
            if(s.contains("cos=2")) {whichCos=2;}
            if(s.contains("cos=3")) {whichCos=3;}
            if(s.contains("cos=4")) {whichCos=4;}
            if(s.contains("cos=5")) {whichCos=5;}
            if(s.contains("cos=10")) {whichCos=10;}
            if(s.contains("cos=11")) {whichCos=11;}
            if(s.contains("cos=12")) {whichCos=12;}
            if(s.contains("cos=13")) {whichCos=13;}
            if(s.contains("cos=14")) {whichCos=14;}
            if(s.contains("cos=15")) {whichCos=15;}
            if(s.contains("acos=0")) {whichAcos=0;}
            if(s.contains("acos=1")) {whichAcos=1;}
            if(s.contains("acos=2")) {whichAcos=2;}
            if(s.contains("acos=3")) {whichAcos=3;}
            if(s.contains("acos=4")) {whichAcos=4;}
            if(s.contains("acos=5")) {whichAcos=5;}
            if(s.contains("acos=10")) {whichAcos=10;}
            if(s.contains("acos=11")) {whichAcos=11;}
            if(s.contains("acos=12")) {whichAcos=12;}
            if(s.contains("acos=13")) {whichAcos=13;}
            if(s.contains("acos=14")) {whichAcos=14;}
            if(s.contains("acos=15")) {whichAcos=15;}
            if(s.contains("fast")){
                whichCos = 4;
                whichAcos = 4;
                whichSin = 4;
            }
            if(s.contains("totalshift")) {totalShift = true;}
            if(s.endsWith("nocache")) {useCache = false;}
            try{
                adapt = new AdaptiveAmberFF(true, null, useCache, whichAcos, whichCos, whichSin, totalShift);
                sWhichMethod = s;
            } catch(Exception e){
                System.err.println("ERROR: Couldn't setup adaptive AMBER FF. Stack coming...");
                e.printStackTrace(System.err);
                return null;
            }
        } else if (s.startsWith("skalevala:")){
            boolean cache = true;
            final String s2 = s.substring(10).trim();
            final String[] sa = s2.split("\\,");
            final int whichMethod = Integer.parseInt(sa[0]);
            if(sa.length > 1) {cache = Boolean.parseBoolean(sa[1]);}
            adapt = new AdaptiveSkalevalaCaller(whichMethod, true, null, cache);
            sWhichMethod = s;
        } else if (s.startsWith("schwefelbench:")) {
            // try to parse the dimensionality
            String s2 = s.substring(14).trim();

            if(!s2.endsWith("D") || !s2.endsWith("D")){
                System.err.println("ERROR: Dimensionality choice doesn't end with D/d. Using 15 D now.");
                adapt = new BenchSchwefels(15);
                sWhichMethod = "schwefelbench:15D";
            } else {

                // chop the D/d off
                s2 = s2.substring(0, s2.length()-1).trim();
                try{
                    int dims = Integer.parseInt(s2);
                    adapt = new BenchSchwefels(dims);
                    sWhichMethod = "schwefelbench:"+dims+"D";
                } catch(Exception e){
                    System.err.println("WARNING: No explicit choice found for the Rastrigin benchmark "
                        + "using 15D now. " + e.toString());
                    adapt = new BenchSchwefels(15);
                    sWhichMethod = "Schwefelbench:15D";
                }
            }
        } else if (s.startsWith("berndbench")) {
            if (s.trim().endsWith(":10D")) {
                adapt = new BenchBerndGauss10D();
                sWhichMethod = "berndbench:10D";
            } else {
                System.err.println("WARNING: No explicit choice found for the Bernd benchmark "
                        + "using 10D now. " + s);
                adapt = new BenchBerndGauss10D();
                sWhichMethod = "berndbench:10D";
            }
        } else if (s.startsWith("schafferf6bench")) {
            if (s.trim().endsWith(":2D")) {
                adapt = new BenchSchafferF6();
                sWhichMethod = "schafferf6bench:2D";
            } else {
                System.err.println("WARNING: No explicit choice found for the Schaffer F6 benchmark "
                        + "using 2D now. " + s);
                adapt = new BenchSchafferF6();
                sWhichMethod = "schafferf6bench:2D";
            }
        } else if (s.startsWith("rastriginbench:")) {
            // try to parse the dimensionality
            String s2 = s.substring(15).trim();

            if(!s2.endsWith("D") || !s2.endsWith("D")){
                System.err.println("ERROR: Dimensionality choice doesn't end with D/d. Using 15 D now.");
                adapt = new BenchRastrigin(15);
                sWhichMethod = "rastriginbench:15D";
            } else {

                // chop the D/d off
                s2 = s2.substring(0, s2.length()-1).trim();
                try{
                    int dims = Integer.parseInt(s2);
                    adapt = new BenchRastrigin(dims);
                    sWhichMethod = "rastriginbench:"+dims+"D";
                } catch(Exception e){
                    System.err.println("WARNING: No explicit choice found for the Rastrigin benchmark "
                        + "using 15D now. " + e.toString());
                    adapt = new BenchRastrigin(15);
                    sWhichMethod = "rastriginbench:15D";
                }
            }
        } else if (s.startsWith("ackleybench:")) {
            // try to parse the dimensionality
            String s2 = s.substring(12).trim();

            if(!s2.endsWith("D") || !s2.endsWith("D")){
                System.err.println("ERROR: Dimensionality choice doesn't end with D/d. Using 15 D now.");
                adapt = new BenchAckley(15, true);
                sWhichMethod = "ackleybench:15D";
            } else {

                // chop the D/d off
                s2 = s2.substring(0, s2.length()-1).trim();
                try{
                    int dims = Integer.parseInt(s2);
                    adapt = new BenchAckley(dims, true);
                    sWhichMethod = "ackleybench:"+dims+"D";
                } catch(Exception e){
                    System.err.println("WARNING: No explicit choice found for the Ackley benchmark "
                        + "using 15D now. " + e.toString());
                    adapt = new BenchAckley(15, true);
                    sWhichMethod = "ackleybench:15D";
                }
            }
        }  else if (s.startsWith("ackleybenchoptgrad:")) {
            // try to parse the dimensionality
            String s2 = s.substring(12).trim();

            if(!s2.endsWith("D") || !s2.endsWith("D")){
                System.err.println("ERROR: Dimensionality choice doesn't end with D/d. Using 15 D now.");
                adapt = new BenchAckley(15, false);
                sWhichMethod = "ackleybench:15D";
            } else {

                // chop the D/d off
                s2 = s2.substring(0, s2.length()-1).trim();
                try{
                    int dims = Integer.parseInt(s2);
                    adapt = new BenchAckley(dims, false);
                    sWhichMethod = "ackleybench:"+dims+"D";
                } catch(Exception e){
                    System.err.println("WARNING: No explicit choice found for the Ackley benchmark "
                        + "using 15D now. " + e.toString());
                    adapt = new BenchAckley(15, false);
                    sWhichMethod = "ackleybench:15D";
                }
            }
        } else if (s.startsWith("lunacekbench:")) {

            // try to parse the dimensionality
            String s2 = s.substring(13).trim();

            if(!s2.endsWith("D") || !s2.endsWith("D")){
                System.err.println("ERROR: Dimensionality choice doesn't end with D/d. Using 15 D now.");
                adapt = new BenchLunacek(15);
                sWhichMethod = "lunacekbench:15D";
            } else {

                // chop the D/d off
                s2 = s2.substring(0, s2.length()-1).trim();
                try{
                    int dims = Integer.parseInt(s2);
                    adapt = new BenchLunacek(dims);
                    sWhichMethod = "lunacekbench:"+dims+"D";
                } catch(Exception e){
                    System.err.println("WARNING: No explicit choice found for the Lunacek benchmark "
                        + "using 15D now. " + e.toString());
                    adapt = new BenchLunacek(15);
                    sWhichMethod = "lunacekbench:15D";
                }
            }
        } else if (s.startsWith("griewangkrosenbrockbench:")) {
            
            // try to parse the dimensionality
            String s2 = s.substring(25).trim();

            if(!s2.endsWith("D") || !s2.endsWith("D")){
                System.err.println("ERROR: Dimensionality choice doesn't end with D/d. Using 15 D now.");
                adapt = new BenchGriewangkRosenbrock(15);
                sWhichMethod = "griewangkrosenbrockbench:15D";
            } else {

                // chop the D/d off
                s2 = s2.substring(0, s2.length()-1).trim();
                try{
                    int dims = Integer.parseInt(s2);
                    adapt = new BenchGriewangkRosenbrock(dims);
                    sWhichMethod = "griewangkrosenbrockbench:"+dims+"D";
                } catch(Exception e){
                    System.err.println("WARNING: No explicit choice found for the Griewangk-Rosenbrock benchmark "
                        + "using 15D now. " + e.toString());
                    adapt = new BenchGriewangkRosenbrock(15);
                    sWhichMethod = "griewangkrosenbrockbench:15D";
                }
            }
        } else if (s.startsWith("schafferf7bench:")) {
            // try to parse the dimensionality
            String s2 = s.substring(16).trim();

            if(!s2.endsWith("D") || !s2.endsWith("D")){
                System.err.println("ERROR: Dimensionality choice doesn't end with D/d. Using 15 D now.");
                adapt = new BenchSchaffer(15);
                sWhichMethod = "schafferf7bench:15D";
            } else {

                // chop the D/d off
                s2 = s2.substring(0, s2.length()-1).trim();
                try{
                    int dims = Integer.parseInt(s2);
                    adapt = new BenchSchaffer(dims);
                    sWhichMethod = "schafferf7bench:"+dims+"D";
                } catch(Exception e){
                    System.err.println("WARNING: No explicit choice found for the Schaffer F7 benchmark "
                        + "using 15D now. " + e.toString());
                    adapt = new BenchSchaffer(15);
                    sWhichMethod = "schafferf7bench:15D";
                }
            }
        } else if(s.startsWith("polaris:")){
            final String s2 = s.substring(8).trim();
            if(s2.equalsIgnoreCase("totalenergy")){
                adapt = new AdaptivePolaris(true);
                sWhichMethod = "polaris:totalenergy";
            } else if(s2.equalsIgnoreCase("waterinteraction")){
                adapt = new AdaptivePolaris(false);
                sWhichMethod = "polaris:waterinteraction";
            } else{
                System.err.println("WARNING: Unknown choice " + s2 + " for polaris. Fitting total energy.");
                adapt = new AdaptivePolaris(true);
                sWhichMethod = "polaris:totalenergy";
            }
        } else if(s.startsWith("spillage:")){
            final String s2 = s.substring(9);
            final String[] sa = s2.trim().split("\\;");
            System.out.println("DEBUG: " + s2 + "  " + sa[0] + "  " + sa[1] + "  " + sa[2]);
            final String[] gs = sa[1].trim().split("\\,");
            final String[] cs = sa[2].trim().split("\\,");
            final int noChannels = Integer.parseInt(sa[0].trim());
            final int[] nog = new int[noChannels];
            final int[] noc = new int[noChannels];
            for(int i = 0; i < noChannels; i++){
                nog[i] = Integer.parseInt(gs[i].trim());
                noc[i] = Integer.parseInt(cs[i].trim());
            }
            
            adapt = new AdaptiveCustomSpillage(noChannels,nog,noc);
            sWhichMethod = "spillage";
        } else if(s.startsWith("forcematchingprofess:")){

            final String s2 = s.substring(21);
            final String[] sa = s2.trim().split("\\;");
            
            boolean doOverallGauss = true;
            boolean doMirrorGauss = true;
            
            int noGaussians = 3;
            double potCutoff = -1;
            double overallBeta = -1;
            String origPPFile = "non-fm-pp.pot";
            final List<Short> nonOptAtoms = new ArrayList<>();
            for(final String token : sa){
                if(token.startsWith("nogauss=")){
                    final String x = token.substring(8).trim();
                    noGaussians = Integer.parseInt(x);
                } else if(token.startsWith("potcutoff=")){
                    final String x = token.substring(10).trim();
                    potCutoff = Double.parseDouble(x);
                } else if(token.startsWith("origppfile=")){
                    origPPFile = token.substring(11).trim();
                } else if(token.startsWith("bigbeta=")){
                    overallBeta = Double.parseDouble(token.substring(8).trim());
                } else if(token.startsWith("noppoptfor=")){
                    final String[] atoms = token.substring(11).trim().split("\\,");
                    for(final String atom : atoms){
                        final short no = AtomicProperties.giveAtomicNumber(atom.trim());
                        nonOptAtoms.add(no);
                    }
                } else{
                    throw new RuntimeException("Illegal option " + token + " for forcematchingprofess.");
                }
            }
            
            if(potCutoff <= 0){throw new RuntimeException("Potential cutoff must be specified for force matching optimization.");}
            if(overallBeta <= 0.0){throw new RuntimeException("Overall beta must be larger than 0.0.");}
            if(noGaussians <= 0){throw new RuntimeException("Must optimize at least one Gaussian.");}
            
            adapt = new FMProfessPseudoPotential(noGaussians, potCutoff, origPPFile, doMirrorGauss, doOverallGauss, overallBeta, nonOptAtoms);
            sWhichMethod = "forcematchingprofess";
        } else {
            throw new RuntimeException("ERROR: No valid adaptivable choice found.");
        }

        Tuple<Adaptivable,String> tupel = new Tuple<>(adapt,sWhichMethod);

        return tupel;
    }
    
    @Override
    public GenericPoolConfig<Double,AdaptiveParameters> getGenericPoolConfig(){
        
        final GenericPoolConfig<Double,AdaptiveParameters> config = new GenericPoolConfig<>();
        
        config.setDoNiching(doNiching);
        config.setSerializeAfterNewBest(this.serializeAfterNewBest);
        config.setAcceptableFitness(acceptableFitness);
        config.setAddsToSerial(paramsToSerial);
        config.setWriteEveryAdd(this.writeEveryParameterset);
        config.setPoolSize(populationSize);
        config.setInterBinFile(intermediatePoolFile());
        
        DiversityChecker<Double,AdaptiveParameters> diver;
        switch(whichDiversityCheck){
            case 0: diver = new GenericDiversityCheckers.FitnessDiversityChecker<>(energyDiversity); break;
            case 2: diver = new ParameterDiversityCheckers.ParamsDiversityChecker(paramDiversityThresh); break;
                //TODO more!!!
            default: diver = new GenericDiversityCheckers.FitnessDiversityChecker<>(energyDiversity); break;
        }
        config.setDiversityChecker(diver);
        
        if(doNiching){
            config.setAddsToStats(this.addsToNicheStats);
            config.setNicher(new SimpleNicher<Double,AdaptiveParameters>(noOfIndividualsPerNicheMax));
        }
        
        ParentSelector<Double,AdaptiveParameters> selec;
        try{
            selec = GenericParentSelectors.buildSelector(whichParentsChoice);
        } catch(Exception e){
            throw new RuntimeException("Error in creating parent selector.",e);
        }
        config.setSelector(selec);
        
        config.setWriter(new ParameterWriter(paramsDumpFolder));
        //XXX config.setStats(new GenericStatistics(this.sOutputPath + System.getProperty("file.separator") + this.sOutputPath + ".log"));
        config.setStats(new GenericStatistics(outputFolder + File.separator + "parameterglobopt.log",10000l)); //XXX hard coded...
        
        return config;
    }
    
    private void addReferenceInformationsCartesians(final String[] data, final int id) throws Exception{
        
        CartesianCoordinates currentCartes = null;        
        BondInfo bonds = null;
        double dRefWeight = 1.0;
        double dRefMaxDiff = Double.MAX_VALUE;
        String pathToBonds = null;
        Energy energy = null;
        Forces forces = null;
        BulkModulus modulus = null;
        CellVolume cellVolume = null;
        EnergyOrder energyOrder = null;
        DeltaGauge deltaGauge = null;
        Density density = null;
        RefBulkModulusData<CartesianCoordinates> bulkModData = null;
        ReferenceForcesData<CartesianCoordinates> forcesData = null;
        RefCellVolumeData<CartesianCoordinates> cellVolumeData = null;
        ReferenceEnergyOrderData<CartesianCoordinates> energyOrderData = null;
        ReferenceDeltaGaugeData<CartesianCoordinates> deltaGaugeData = null;
        ReferenceDensityData<CartesianCoordinates> densityData = null;
        String tag = "N/A";

        reference:
        for (int i = 0; i < data.length; i++) {
            if (data[i].equalsIgnoreCase("<PATH>")) {
                if (!data[i + 2].equalsIgnoreCase("</PATH>")) {
                    throw new Exception("ERROR: No closing </PATH> tag found.");
                } else {
                    final String sPath = data[i + 1].trim();
                            System.out.println("DEBUG: Tag for geom " + id + "\t" + sPath);
                    if (sPath.endsWith(".xyz")) {
                        currentCartes = org.ogolem.core.Input.readCartesFromFile(sPath);
                    } else if (sPath.endsWith(".zmat")) {
                        ZMatrix zmat = org.ogolem.core.Input.readZmatFromFile(sPath);
                        // translate it
                        currentCartes = CoordTranslation.zMatToCartesians(zmat);
                    } else if (sPath.endsWith(".pxyz")) {
                        System.out.println("INFO: You've specified a .pxyz file. "
                                + "Be aware that this is only OK for double morse fitting.");
                        currentCartes = org.ogolem.core.Input.readCartesFromFile(sPath);
                    } else {
                        System.err.println("WARNING: Unknown extension for the reference path. Assuming it to be xyz, continuing on own risk.");
                        currentCartes = org.ogolem.core.Input.readCartesFromFile(sPath);
                    }

                    // we null spins and charges
                    currentCartes.setAllSpins(new short[currentCartes.getNoOfAtoms()]);
                    currentCartes.setAllCharges(new float[currentCartes.getNoOfAtoms()]);
                    
                    i += 2;
                    continue reference;
                }
            } else if (data[i].equalsIgnoreCase("<ENERGY>")) {
                if (!data[i + 2].equalsIgnoreCase("</ENERGY>")) {
                    throw new RuntimeException("ERROR: No closing </ENERGY> tag found.");
                } else {
                    try {
                        final double e = Double.parseDouble(data[i + 1].trim());
                        energy = new Energy(e);
                    } catch (Exception e) {
                        throw new RuntimeException("Couldn't cast reference energy. ",e);
                    }
                    i += 2;
                    continue reference;
                }
            } else if (data[i].equalsIgnoreCase("<BULKMODULUS>")){
                if (!data[i + 2].equalsIgnoreCase("</BULKMODULUS>")) {
                    throw new RuntimeException("ERROR: No closing </BULKMODULUS> tag found.");
                } else {
                    final String[] line = data[i + 1].trim().split("\\:");
                    final String crystStruct = line[0].trim();
                    final int indexBra = crystStruct.indexOf("(");
                    final String atom = crystStruct.substring(0, indexBra).trim();
                    final short atomNo = AtomicProperties.giveAtomicNumber(atom);
                    final String sym = crystStruct.substring(indexBra+1,crystStruct.length()-1).trim();
                    ReferenceGeomData<BulkModulus,CartesianCoordinates> geom = new ReferenceGeomData<>(currentCartes, bonds, id);
                    bulkModData = new RefBulkModulusData<>(id,sym,atomNo, geom);
                    try {
                        final double bulkMod = Double.parseDouble(line[1].trim());
                        modulus = new BulkModulus(bulkMod);
                    } catch (Exception e) {
                        throw new RuntimeException("Couldn't cast reference bulk modulus.",e);
                    }
                    i += 2;
                    continue reference;
                }
            } else if (data[i].equalsIgnoreCase("<CELLVOLUME>")){
                if (!data[i + 2].equalsIgnoreCase("</CELLVOLUME>")) {
                    throw new RuntimeException("ERROR: No closing </CELLVOLUME> tag found.");
                } else {
                    final String[] line = data[i + 1].trim().split("\\:");
                    final String crystStruct = line[0].trim();
                    final int indexBra = crystStruct.indexOf("(");
                    final String atom = crystStruct.substring(0, indexBra).trim();
                    final short atomNo = AtomicProperties.giveAtomicNumber(atom);
                    final String sym = crystStruct.substring(indexBra+1,crystStruct.length()-1).trim();
                    ReferenceGeomData<CellVolume,CartesianCoordinates> geom = new ReferenceGeomData<>(currentCartes, bonds, id);;
                    cellVolumeData = new RefCellVolumeData<>(id,sym,atomNo,geom);
                    try {
                        final double cellVol = Double.parseDouble(line[1].trim());
                        cellVolume = new CellVolume(cellVol);
                    } catch (Exception e) {
                        throw new RuntimeException("Couldn't cast reference cell volume.",e);
                    }
                    i += 2;
                    continue reference;
                }
            } else if (data[i].equalsIgnoreCase("<FORCES>")){
                if (!data[i + 2].equalsIgnoreCase("</FORCES>")) {
                    throw new RuntimeException("ERROR: No closing </FORCES> tag found.");
                } else {

                    final String forceLine = data[i + 1].trim();
                    if(forceLine.startsWith("ForcesFile=")){
                        final String file = forceLine.substring(11).trim();
                        final String[] fileData = InputPrimitives.readFileIn(file);
                        final int noAtoms = Integer.parseInt(fileData[0].trim());
                        final double[][] forceVals = new double[3][noAtoms];
                        for(int at = 0; at < noAtoms; at++){
                            final String[] line = fileData[at+2].trim().split("\\s+");
                            forceVals[0][at] = Double.parseDouble(line[1]);
                            forceVals[1][at] = Double.parseDouble(line[2]);
                            forceVals[2][at] = Double.parseDouble(line[3]);
                        }
                        
                        forces = new Forces(forceVals);
                    } else {
                        throw new RuntimeException("<FORCES> tag must contain ForcesFile= followed by the path to the file containing the reference forces (in pseudo-xyz format).");
                    }
                    
                    i += 2;
                    continue reference;
                }
            } else if (data[i].equalsIgnoreCase("<DENSITY>")){
                if (!data[i + 3].equalsIgnoreCase("</DENSITY>")) {
                    throw new RuntimeException("ERROR: No closing </DENSITY> tag found.");
                } else {

                    final String densityLine = data[i + 1].trim();
                    final String densityLine2 = data[i + 2].trim();
                    if(densityLine.startsWith("DensityFile=")){
                        final String file = densityLine.substring(12).trim();
                        final String[] fileData = InputPrimitives.readFileIn(file);
                        
                        // first line: grid dimension
                        try{
                            final String[] sa = fileData[0].trim().split("\\s+");
                            final int xDim = Integer.parseInt(sa[0]);
                            final int yDim = Integer.parseInt(sa[1]);
                            final int zDim = Integer.parseInt(sa[2]);
                            
                            final double[][][] densityValues = new double[xDim][yDim][zDim];
                        
                            int count = 1;
                            for(int x = 0; x < xDim; x++){
                                for(int y = 0; y < yDim; y++){
                                    for(int z = 0; z < zDim; z++){
                                        final String val = fileData[count].trim();
                                        densityValues[x][y][z] = Double.parseDouble(val);
                                        count++;
                                    }
                                }
                            }
                            
                            if(false){
                                System.out.println("DEBUG: Read in density with dimensions " + xDim + " " + yDim + " " + zDim);
                            }
                            
                            density = new Density(densityValues);
                        } catch (Exception e){
                            throw new RuntimeException("<DENSITY> tag must contain DensityFile= followed by the path to the file containing the reference density (first three grid dimensions followed by one density element per line).", e);
                        }
                    } else {
                        throw new RuntimeException("<DENSITY> tag must contain DensityFile= in the first line followed by the path to the file containing the reference density (first three grid dimensions followed by one density element per line).");
                    }
                    
                    if(densityLine2.startsWith("ReferenceTag=")){
                        tag = densityLine2.substring(13).trim();
                    } else {
                        throw new RuntimeException("<DENSITY> tag must contain ReferenceTag= in the second line followed by a unique identifier for this density.");
                    }
                    
                    i += 2;
                    continue reference;
                }
            } else if(data[i].equalsIgnoreCase("<ENERGYORDER>")){
                
                final List<Tuple<String,ReferenceGeomData<EnergyOrder,CartesianCoordinates>>> ll = new ArrayList<>();
                final List<Tuple<String,Double>> eos = new ArrayList<>();
                i++;
                
                boolean considerRelEnergies = false;
                if(data[i].equalsIgnoreCase("ConsiderRelativeEnergies")){
                    considerRelEnergies = true;
                    i++;
                }
                
                while(!data[i].trim().equalsIgnoreCase("</ENERGYORDER>")){
                    
                    final String[] sa = data[i].trim().split("\\;");
                    // format is: key, file, energy
                    final String key = sa[0].trim();
                    final String file = sa[1].trim();
                    final double e = Double.parseDouble(sa[2].trim());
                    
                    CartesianCoordinates cartes = null;                    
                    if (file.endsWith(".xyz")) {
                        cartes = org.ogolem.core.Input.readCartesFromFile(file);
                    } else if (file.endsWith(".zmat")) {
                        ZMatrix zmat = org.ogolem.core.Input.readZmatFromFile(file);
                        // translate it
                        cartes = CoordTranslation.zMatToCartesians(zmat);
                    } else if (file.endsWith(".pxyz")) {
                        System.out.println("INFO: You've specified a .pxyz file. "
                                + "Be aware that this is only OK for double morse fitting.");
                        cartes = org.ogolem.core.Input.readCartesFromFile(file);
                    } else {
                        System.err.println("WARNING: Unknown extension for the reference path. Assuming it to be xyz, continuing on own risk.");
                        cartes = org.ogolem.core.Input.readCartesFromFile(file);
                    }

                    // we null spins and charges
                    cartes.setAllSpins(new short[cartes.getNoOfAtoms()]);
                    cartes.setAllCharges(new float[cartes.getNoOfAtoms()]);
                    
                    System.out.println("INFO: We currently auto-create the bond info for energy order item " + key + ". If this is a problem for you, contact the author(s).");
                    final BondInfo bondInfo = CoordTranslation.checkForBonds(cartes, 1.3);
                    
                    final ReferenceGeomData<EnergyOrder,CartesianCoordinates> geom = new ReferenceGeomData<>(cartes,bondInfo,id);
                    final Tuple<String,ReferenceGeomData<EnergyOrder,CartesianCoordinates>> tup = new Tuple<>(key,geom);
                    ll.add(tup);
                    
                    final Tuple<String,Double> eTup = new Tuple<>(key,e);
                    eos.add(eTup);
                    
                    i++;
                }
                
                if(ll.isEmpty() || ll.size() < 2){
                    throw new RuntimeException("<ENERGYORDER> should contain two or more reference points. Otherwise it is not an order, is it?");
                }
                
                energyOrderData = new ReferenceEnergyOrderData<>(id,ll);
                energyOrder = new EnergyOrder(eos, considerRelEnergies);
            } else if (data[i].equalsIgnoreCase("<DELTAGAUGE>")){
                if (!data[i + 2].equalsIgnoreCase("</DELTAGAUGE>")) {
                    throw new RuntimeException("ERROR: No closing </DELTAGAUGE> tag found.");
                } else {
                    final String[] line = data[i + 1].trim().split("\\:");
                    final String crystStruct = line[0].trim();
                    final int indexBra = crystStruct.indexOf("(");
                    final String atom = crystStruct.substring(0, indexBra).trim();
                    final short atomNo = AtomicProperties.giveAtomicNumber(atom);
                    final String sym = crystStruct.substring(indexBra+1,crystStruct.length()-1).trim();
                    ReferenceGeomData<DeltaGauge,CartesianCoordinates> geom = null;
                    deltaGaugeData = new ReferenceDeltaGaugeData<>(id,sym,atomNo, geom);
                    try {
                        final double dg = Double.parseDouble(line[1].trim());
                        deltaGauge = new DeltaGauge(dg);
                    } catch (Exception e) {
                        throw new RuntimeException("Couldn't cast reference delta gauge.",e);
                    }
                    i += 2;
                    continue reference;
                }
            } else if (data[i].equalsIgnoreCase("<REFCHARGES>")) {
                // try to find out where the charges end
                int endChargeBlock = -1;
                for (int j = 0; j < data.length; j++) {
                    if (data[j].equalsIgnoreCase("</REFCHARGES>")) {
                        endChargeBlock = j;
                        break;
                    }
                }
                if (endChargeBlock == -1) {
                    throw new RuntimeException("ERROR: No closing </REFCHARGES> tag found.");
                }

                // get the content out
                final String[] chargeCont = new String[endChargeBlock - i - 1];
                for (int cont = i + 1; cont < endChargeBlock; cont++) {
                    chargeCont[cont - i - 1] = data[cont];
                }

                // get the charges
                final float[][] tmpCharges = org.ogolem.core.Input.readChargesIn(chargeCont, 1, currentCartes.getAllAtomsPerMol(), new int[]{0});
                currentCartes.setAllCharges(tmpCharges[0]);

                // continue the loop
                i = endChargeBlock;
                continue reference;
            } else if (data[i].equalsIgnoreCase("<REFSPINS>")) {
                // try to find out where the spins end
                int endSpinBlock = -1;
                for (int j = 0; j < data.length; j++) {
                    if (data[j].equalsIgnoreCase("</REFSPINS>")) {
                        endSpinBlock = j;
                        break;
                    }
                }
                if (endSpinBlock == -1) {
                    throw new RuntimeException("ERROR: No closing </REFSPINS> tag found.");
                }

                // get the content out
                final String[] spinCont = new String[endSpinBlock - i - 1];
                for (int cont = i + 1; cont < endSpinBlock; cont++) {
                    spinCont[cont - i - 1] = data[cont];
                }

                // get the spins
                final short[][] tmpSpins = org.ogolem.core.Input.readSpinsIn(spinCont, 1, currentCartes.getAllAtomsPerMol(), new int[]{0});
                currentCartes.setAllSpins(tmpSpins[0]);

                // continue the loop
                i = endSpinBlock;
                continue reference;
            } else if (data[i].equalsIgnoreCase("<REFWEIGHT>")) {
                if (!data[i + 2].equalsIgnoreCase("</REFWEIGHT>")) {
                    throw new RuntimeException("ERROR: No closing </REFWEIGHT> tag found.");
                } else {
                    try {
                        dRefWeight = Double.parseDouble(data[i + 1].trim());
                    } catch (Exception e) {
                        throw new RuntimeException("Couldn't cast reference weight.",e);
                    }
                    i = i + 2;
                    continue reference;
                }
            } else if (data[i].equalsIgnoreCase("<REFMAXDIFF>")) {
                if (!data[i + 2].equalsIgnoreCase("</REFMAXDIFF>")) {
                    throw new RuntimeException("ERROR: No closing </REFMAXDIFF> tag found.");
                } else {
                    try {
                        dRefMaxDiff = Double.parseDouble(data[i + 1].trim());
                    } catch (Exception e) {
                        throw new RuntimeException("Couldn't cast maximum allowed difference.",e);
                    }
                    i += 2;
                    continue reference;
                }
            } else if (data[i].equalsIgnoreCase("<REFBONDS>")) {
                if (!data[i + 2].equalsIgnoreCase("</REFBONDS>")) {
                    throw new RuntimeException("ERROR: No closing </REFBONDS> tag found.");
                } else {
                    try {
                        pathToBonds = data[i + 1].trim();
                    } catch (Exception e) {
                        throw new RuntimeException("Couldn't cast maximum allowed difference.",e);
                    }
                    i += 2;
                    continue reference;
                }
            }
        }
        
        if(currentCartes != null){
            // now we check whether or not spins and charges were set
            if (currentCartes.getAllCharges() == null) {
                // set zero'd charges in
                currentCartes.setAllCharges(new float[currentCartes.getNoOfAtoms()]);
            }

            if (currentCartes.getAllSpins() == null) {
                // set zero'd spins in
                currentCartes.setAllSpins(new short[currentCartes.getNoOfAtoms()]);
            }

            // deal with the bonds
            if (pathToBonds == null) {
                bonds = CoordTranslation.checkForBonds(currentCartes, 1.3); // hardcoded!
            } else {
                try {
                    if (pathToBonds.endsWith(".list")) {
                        bonds = org.ogolem.core.Input.readBondsListIn(pathToBonds, currentCartes.getNoOfAtoms());
                    } else {
                        bonds = org.ogolem.core.Input.readBondsIn(pathToBonds);
                    }
                } catch (Exception e) {
                    throw new RuntimeException("ERROR: Couldn't read bonds in. ",e);
                }
            }
        }
        
        // compile reference data and see what we got
        
        /* for the energy property */
        if(energy != null){
            if(currentCartes == null || bonds == null){
                throw new RuntimeException("Energy reference specified but no reference input data for either cartesians or bonds!");
            }
            final ReferenceGeomData<Energy,CartesianCoordinates> geomD = new ReferenceGeomData<>(currentCartes, bonds, id);
            final GenericReferencePoint<Energy,ReferenceGeomData<Energy,?>> point = new ReferencePoint<Energy,ReferenceGeomData<Energy,?>>(energy,geomD,id,dRefWeight,dRefMaxDiff);
            referenceEnergies.add(point);
        }
        
        /* for the forces property */
        if(forces != null){
            forcesData = new ReferenceForcesData(id, new ReferenceGeomData<>(currentCartes, bonds, id));
            final GenericReferencePoint<Forces,ReferenceForcesData<?>> point = new ReferencePoint<Forces,ReferenceForcesData<?>>(forces,forcesData,id,dRefWeight,dRefMaxDiff);
            referenceForces.add(point);
        }
        
        /* for the density property */
        if(density != null){
            densityData = new ReferenceDensityData(id, tag, new ReferenceGeomData<>(currentCartes, bonds, id));
            final GenericReferencePoint<Density,ReferenceDensityData<?>> point = new ReferencePoint<Density,ReferenceDensityData<?>>(density,densityData,id,dRefWeight,dRefMaxDiff);
            referenceDensity.add(point);
        }
        
        /* for the bulk modulus property */
        if(modulus != null){
            if(bulkModData == null){
                throw new RuntimeException("Bulk modulus reference specified but no reference input data!");
            }
            final GenericReferencePoint<BulkModulus,RefBulkModulusData<?>> point = new ReferencePoint<BulkModulus,RefBulkModulusData<?>>(modulus,bulkModData,id,dRefWeight,dRefMaxDiff);
            referenceBulkMod.add(point);
        }
        
        /* for the equillibrium cell volume */
        if(cellVolume != null){
            if(cellVolumeData == null){
                throw new RuntimeException("Cell volume reference specified but no reference input data!");
            }
            final GenericReferencePoint<CellVolume,RefCellVolumeData<?>> point = new ReferencePoint<CellVolume,RefCellVolumeData<?>>(cellVolume,cellVolumeData,id,dRefWeight,dRefMaxDiff);
            referenceCellVolume.add(point);
        }
        
        /* for the energy order */
        if(energyOrder != null){
            if(energyOrderData == null){
                throw new RuntimeException("Energy order reference specified but no reference input data!");
            }
            final GenericReferencePoint<EnergyOrder,ReferenceEnergyOrderData<?>> point = new ReferencePoint<EnergyOrder,ReferenceEnergyOrderData<?>>(energyOrder,energyOrderData,id,dRefWeight,dRefMaxDiff);
            referenceEnergyOrder.add(point);
        }
        
        /* for the delta gauge property */
        if(deltaGauge != null){
            if(deltaGaugeData == null){
                throw new RuntimeException("Delta gauge reference specified but no reference input data!");
            }
            deltaGaugeData.setGeomData(new ReferenceGeomData<DeltaGauge,CartesianCoordinates>(currentCartes, bonds, id));
            final GenericReferencePoint<DeltaGauge,ReferenceDeltaGaugeData<?>> point = new ReferencePoint<DeltaGauge,ReferenceDeltaGaugeData<?>>(deltaGauge,deltaGaugeData,id,dRefWeight,dRefMaxDiff);
            referenceDeltaGauge.add(point);
        }
    }
    
    final GenericFitnessFunction getCartesianFitnessFunc(){
        
        final List<GenericFitnessTerm<?>> terms = new ArrayList<>();
        final List<Double> weights = new ArrayList<>();
        
        final List<SerialBatchedPropertyCalculator.PropertyBatch> batches = createBatches();
        final BatchedPropertyCalculator batcher = new SerialBatchedPropertyCalculator(refAdaptivable.clone(),batches);
        
        if(!referenceEnergies.isEmpty()){
            // XXX: ideally deprecate eventually!
            final EnergyCalculator enerCalc = new EnergyCalculator(refAdaptivable.clone());
            final PseudoPropertyCalculator<Energy,ReferenceGeomData<Energy,CartesianCoordinates>> wrapped = new PseudoPropertyCalculator<>(enerCalc,batcher,new Energy(-42));
                        
            FitnessTermConfig<Energy> conf;
            try{
                final String[] block = FitnessTermConfig.findBlockFor(allTermConfigs, "ENERGYTERM");
                conf = new FitnessTermConfig<>(block);
            } catch(Exception e){
                throw new RuntimeException("Couldn't create config for energy term.",e);
            }
            
            List<GenericReferencePoint<Energy,ReferenceGeomData<Energy,CartesianCoordinates>>> referenceEs = new ArrayList<>();
            for(final GenericReferencePoint<Energy,ReferenceGeomData<Energy,?>> p : referenceEnergies){
                // let's just cast this around...
                // XXX god is this ugly
                referenceEs.add((GenericReferencePoint<Energy,ReferenceGeomData<Energy,CartesianCoordinates>>) (Object)p);
            }
            
            final GenericFitnessTerm<Energy> energyTerm = new SerialGenericFitnessTerm<>(referenceEs,
                wrapped, conf, printContributions);
        
            terms.add(energyTerm);
            weights.add(conf.getTermWeight());
        }
        
        if(!referenceForces.isEmpty()){
            final Forces forces = new Forces(new double[0][0]);
            final ReferenceForcesData<CartesianCoordinates> dummy = new ReferenceForcesData<>(-1,new ReferenceGeomData<Forces,CartesianCoordinates>(null,null,-1));
            final PropertyCalculator<Forces,ReferenceForcesData<CartesianCoordinates>> propCalc = refAdaptivable.clone().<Forces,ReferenceForcesData<CartesianCoordinates>>getCalculatorForProperty(forces, dummy);
            if(propCalc == null){
                throw new RuntimeException("No property calculator available for forces and reference forces data, although term was requested!");
            }
            final PseudoPropertyCalculator<Forces,ReferenceForcesData<CartesianCoordinates>> wrapped = new PseudoPropertyCalculator<Forces,ReferenceForcesData<CartesianCoordinates>>(propCalc,batcher,forces);
            
            List<GenericReferencePoint<Forces,ReferenceForcesData<CartesianCoordinates>>> referenceFs = new ArrayList<>();
            for(final GenericReferencePoint<Forces,ReferenceForcesData<?>> p : referenceForces){
                // let's just cast this around...
                // XXX god is this ugly
                referenceFs.add((GenericReferencePoint<Forces,ReferenceForcesData<CartesianCoordinates>>) (Object)p);
            }
            
            FitnessTermConfig<Forces> conf;
            try{
                final String[] block = FitnessTermConfig.findBlockFor(allTermConfigs, "FORCESTERM");
                conf = new FitnessTermConfig<>(block);
            } catch(Exception e){
                throw new RuntimeException("Couldn't create config for forces term.",e);
            }
            
            final GenericFitnessTerm<Forces> forcesTerm = new SerialGenericFitnessTerm<>(referenceFs,
                wrapped, conf, printContributions);
            
            terms.add(forcesTerm);
            weights.add(conf.getTermWeight());
        }
        
        if(!referenceBulkMod.isEmpty()){
            final BulkModulus modulus = new BulkModulus(42.0);
            final RefBulkModulusData<CartesianCoordinates> dummy = new RefBulkModulusData<>(-1,"somesymmetry", (short)0, null);
            final PropertyCalculator<BulkModulus,RefBulkModulusData<CartesianCoordinates>> propCalc = refAdaptivable.clone().<BulkModulus,RefBulkModulusData<CartesianCoordinates>>getCalculatorForProperty(modulus,dummy);
            if(propCalc == null){
                throw new RuntimeException("No property calculator available for bulk modulus and reference bulk modulus data, although term was requested!");
            }
            final PseudoPropertyCalculator<BulkModulus,RefBulkModulusData<CartesianCoordinates>> wrapped = new PseudoPropertyCalculator<>(propCalc,batcher,modulus);
            
            List<GenericReferencePoint<BulkModulus,RefBulkModulusData<CartesianCoordinates>>> referenceBs = new ArrayList<>();
            for(final GenericReferencePoint<BulkModulus,RefBulkModulusData<?>> p : referenceBulkMod){
                // let's just cast this around...
                // XXX god is this ugly
                referenceBs.add((GenericReferencePoint<BulkModulus,RefBulkModulusData<CartesianCoordinates>>) (Object)p);
            }
            
            FitnessTermConfig<BulkModulus> conf;
            try{
                final String[] block = FitnessTermConfig.findBlockFor(allTermConfigs, "BULKMODULUSTERM");
                conf = new FitnessTermConfig<>(block);
            } catch(Exception e){
                throw new RuntimeException("Couldn't create config for bulk modulus term.",e);
            }
            
            final GenericFitnessTerm<BulkModulus> modulusTerm = new SerialGenericFitnessTerm<>(referenceBs,
                wrapped, conf, printContributions);
            
            terms.add(modulusTerm);
            weights.add(conf.getTermWeight());
        }
        
        if(!referenceCellVolume.isEmpty()){
            final CellVolume volume = new CellVolume(42.0);
            final RefCellVolumeData<CartesianCoordinates> dummy = new RefCellVolumeData<>(-1,"somesymmetry", (short)0, null);
            final PropertyCalculator<CellVolume,RefCellVolumeData<CartesianCoordinates>> propCalc = refAdaptivable.clone().<CellVolume,RefCellVolumeData<CartesianCoordinates>>getCalculatorForProperty(volume,dummy);
            if(propCalc == null){
                throw new RuntimeException("No property calculator available for cell volume and reference cell volume data, although term was requested!");
            }
            final PseudoPropertyCalculator<CellVolume,RefCellVolumeData<CartesianCoordinates>> wrapped = new PseudoPropertyCalculator<>(propCalc,batcher,volume);
            
            List<GenericReferencePoint<CellVolume,RefCellVolumeData<CartesianCoordinates>>> referenceVs = new ArrayList<>();
            for(final GenericReferencePoint<CellVolume,RefCellVolumeData<?>> p : referenceCellVolume){
                // let's just cast this around...
                // XXX god is this ugly
                referenceVs.add((GenericReferencePoint<CellVolume,RefCellVolumeData<CartesianCoordinates>>) (Object)p);
            }
            
            FitnessTermConfig<CellVolume> conf;
            try{
                final String[] block = FitnessTermConfig.findBlockFor(allTermConfigs, "CELLVOLUMETERM");
                conf = new FitnessTermConfig<>(block);
            } catch(Exception e){
                throw new RuntimeException("Couldn't create config for cell volume term.",e);
            }
            
            final GenericFitnessTerm<CellVolume> cellVolumeTerm = new SerialGenericFitnessTerm<>(referenceVs,
                wrapped, conf, printContributions);
            
            terms.add(cellVolumeTerm);
            weights.add(conf.getTermWeight());
        }
        
        if(!referenceEnergyOrder.isEmpty()){
            final Tuple<String,Double> t1 = new Tuple<>("bcc",21.0);
            final Tuple<String,Double> t2 = new Tuple<>("fcc",42.0);
            final Tuple<String,Double> t3 = new Tuple<>("hcp",84.0);
            final List<Tuple<String,Double>> ll = new ArrayList<>();
            ll.add(t1);
            ll.add(t2);
            ll.add(t3);            
            
            final boolean consRelOrder = referenceEnergyOrder.get(0).getReferenceProperty().shouldConsiderRelEnergies();
            final EnergyOrder order = new EnergyOrder(ll, consRelOrder);
            
            final List<Tuple<String,ReferenceGeomData<EnergyOrder,CartesianCoordinates>>> refLi = new ArrayList<>();
            final Tuple<String,ReferenceGeomData<EnergyOrder,CartesianCoordinates>> tt = new Tuple<>("bcc", null);
            refLi.add(tt);
            
            final ReferenceEnergyOrderData<CartesianCoordinates> dummy = new ReferenceEnergyOrderData<>(-1,refLi);
            final PropertyCalculator<EnergyOrder,ReferenceEnergyOrderData<CartesianCoordinates>> propCalc = refAdaptivable.clone().<EnergyOrder,ReferenceEnergyOrderData<CartesianCoordinates>>getCalculatorForProperty(order,dummy);
            if(propCalc == null){
                throw new RuntimeException("No property calculator available for energy order and reference energy order data, although term was requested!");
            }
            final PseudoPropertyCalculator<EnergyOrder,ReferenceEnergyOrderData<CartesianCoordinates>> wrapped = new PseudoPropertyCalculator<>(propCalc,batcher,order);
            
            List<GenericReferencePoint<EnergyOrder,ReferenceEnergyOrderData<CartesianCoordinates>>> referenceEOs = new ArrayList<>();
            for(final GenericReferencePoint<EnergyOrder,ReferenceEnergyOrderData<?>> p : referenceEnergyOrder){
                // let's just cast this around...
                // XXX god is this ugly
                referenceEOs.add((GenericReferencePoint<EnergyOrder,ReferenceEnergyOrderData<CartesianCoordinates>>) (Object)p);
            }
            
            FitnessTermConfig<EnergyOrder> conf;
            try{
                final String[] block = FitnessTermConfig.findBlockFor(allTermConfigs, "ENERGYORDERTERM");
                conf = new FitnessTermConfig<>(block);
            } catch(Exception e){
                throw new RuntimeException("Couldn't create config for energy order term.",e);
            }
            
            final GenericFitnessTerm<EnergyOrder> energyOrderTerm = new SerialGenericFitnessTerm<>(referenceEOs,
                wrapped, conf, printContributions);
            
            terms.add(energyOrderTerm);
            weights.add(conf.getTermWeight());
        }
        
        if(!referenceDeltaGauge.isEmpty()){
            final DeltaGauge deltaGauge = new DeltaGauge(42.0);
            final ReferenceDeltaGaugeData<CartesianCoordinates> dummy = new ReferenceDeltaGaugeData<>(-1,"somesymmetry", (short)0, null);
            final PropertyCalculator<DeltaGauge,ReferenceDeltaGaugeData<CartesianCoordinates>> propCalc = refAdaptivable.clone().<DeltaGauge,ReferenceDeltaGaugeData<CartesianCoordinates>>getCalculatorForProperty(deltaGauge,dummy);
            if(propCalc == null){
                throw new RuntimeException("No property calculator available for delta gauge and reference delta gauge data, although term was requested!");
            }
            final PseudoPropertyCalculator<DeltaGauge,ReferenceDeltaGaugeData<CartesianCoordinates>> wrapped = new PseudoPropertyCalculator<>(propCalc,batcher,deltaGauge);
            
            List<GenericReferencePoint<DeltaGauge,ReferenceDeltaGaugeData<CartesianCoordinates>>> referenceDGs = new ArrayList<>();
            for(final GenericReferencePoint<DeltaGauge,ReferenceDeltaGaugeData<?>> p : referenceDeltaGauge){
                // let's just cast this around...
                // XXX god is this ugly
                referenceDGs.add((GenericReferencePoint<DeltaGauge,ReferenceDeltaGaugeData<CartesianCoordinates>>) (Object)p);
            }
            
            FitnessTermConfig<DeltaGauge> conf;
            try{
                final String[] block = FitnessTermConfig.findBlockFor(allTermConfigs, "DELTAGAUGETERM");
                conf = new FitnessTermConfig<>(block);
            } catch(Exception e){
                throw new RuntimeException("Couldn't create config for delta gauge term.",e);
            }
            
            final GenericFitnessTerm<DeltaGauge> modulusTerm = new SerialGenericFitnessTerm<>(referenceDGs,
                wrapped, conf, printContributions);
            
            terms.add(modulusTerm);
            weights.add(conf.getTermWeight());
        }
        
        if(!referenceDensity.isEmpty()){
            
            final GenericReferencePoint<Density,ReferenceDensityData<?>> zeroth = referenceDensity.get(0);
            
            final Density density = new Density(null);
            final ReferenceDensityData<CartesianCoordinates> dummy = new ReferenceDensityData<>(-1,"dummy",new ReferenceGeomData<Density,CartesianCoordinates>(null,null,-1));
            final PropertyCalculator<Density,ReferenceDensityData<CartesianCoordinates>> propCalc = refAdaptivable.clone().<Density,ReferenceDensityData<CartesianCoordinates>>getCalculatorForProperty(density, dummy);
            if(propCalc == null){
                throw new RuntimeException("No property calculator available for forces and reference forces data, although term was requested!");
            }
            final PseudoPropertyCalculator<Density,ReferenceDensityData<CartesianCoordinates>> wrapped = new PseudoPropertyCalculator<Density,ReferenceDensityData<CartesianCoordinates>>(propCalc,batcher,density);
            
            List<GenericReferencePoint<Density,ReferenceDensityData<CartesianCoordinates>>> referenceDs = new ArrayList<>();
            for(final GenericReferencePoint<Density,ReferenceDensityData<?>> p : referenceDensity){
                // let's just cast this around...
                // XXX god is this ugly
                referenceDs.add((GenericReferencePoint<Density,ReferenceDensityData<CartesianCoordinates>>) (Object)p);
            }
            
            FitnessTermConfig<Density> conf;
            try{
                final String[] block = FitnessTermConfig.findBlockFor(allTermConfigs, "DENSITYTERM");
                conf = new FitnessTermConfig<>(block);
            } catch(Exception e){
                throw new RuntimeException("Couldn't create config for density term.",e);
            }
            
            final GenericFitnessTerm<Density> densityTerm = new SerialGenericFitnessTerm<>(referenceDs,
                wrapped, conf, printContributions);
            
            terms.add(densityTerm);
            weights.add(conf.getTermWeight());
        }
        
        if(terms.isEmpty()){
            throw new RuntimeException("No reference terms to optimize for!");
        }
                
        final GenericFitnessFunction genFunc = new GenericFitnessFunction(batcher, terms, weights, printContributions);
                
        return genFunc;
    }
    
    private List<SerialBatchedPropertyCalculator.PropertyBatch> createBatches(){
        
        int maxRefID = -1;
        if(!referenceEnergies.isEmpty()){
            for(final GenericReferencePoint<Energy,ReferenceGeomData<Energy,?>> point : referenceEnergies){
                maxRefID = Math.max(maxRefID, point.getReferenceID());
            }
        }
        
        if(!referenceForces.isEmpty()){
            for(final GenericReferencePoint<Forces,ReferenceForcesData<?>> point : referenceForces){
                maxRefID = Math.max(maxRefID, point.getReferenceID());
            }
        }
        
        if(!referenceBulkMod.isEmpty()){
            for(final GenericReferencePoint<BulkModulus,RefBulkModulusData<?>> point : referenceBulkMod){
                maxRefID = Math.max(maxRefID, point.getReferenceID());
            }
        }
        
        if(!referenceCellVolume.isEmpty()){
            for(final GenericReferencePoint<CellVolume,RefCellVolumeData<?>> point : referenceCellVolume){
                maxRefID = Math.max(maxRefID, point.getReferenceID());
            }
        }
        
        if(!referenceEnergyOrder.isEmpty()){
            for(final GenericReferencePoint<EnergyOrder,ReferenceEnergyOrderData<?>> point : referenceEnergyOrder){
                maxRefID = Math.max(maxRefID, point.getReferenceID());
            }
        }
        
        if(!referenceDeltaGauge.isEmpty()){
            for(final GenericReferencePoint<DeltaGauge,ReferenceDeltaGaugeData<?>> point : referenceDeltaGauge){
                maxRefID = Math.max(maxRefID, point.getReferenceID());
            }
        }
        
        if(!referenceDensity.isEmpty()){
            for(final GenericReferencePoint<Density,ReferenceDensityData<?>> point : referenceDensity){
                maxRefID = Math.max(maxRefID, point.getReferenceID());
            }
        }
        
        // to avoid looping over the references too often
        int indEn = 0;
        int indFo = 0;
        int indBu = 0;
        int indCe = 0;
        int indDG = 0;
        int indDe = 0;
        
        final List<SerialBatchedPropertyCalculator.PropertyBatch> batches = new ArrayList<>();
        for(int id = 0; id <= maxRefID; id++){
            
            final List<GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>>> points = new ArrayList<>();
            
            // this works, because we know that the references are sorted in ascending order
            if(!referenceEnergies.isEmpty()){
                for(int i = indEn; i < referenceEnergies.size(); i++){
                    if(referenceEnergies.get(i).getReferenceID() == id){
                        // found
                        points.add(referenceEnergies.get(i));
                        indEn = i+1;
                        break;
                    }
                }
            }
            if(!referenceForces.isEmpty()){
                for(int i = indFo; i < referenceForces.size(); i++){
                    if(referenceForces.get(i).getReferenceID() == id){
                        // found
                        points.add(referenceForces.get(i));
                        indFo = i+1;
                        break;
                    }
                }
            }
            if(!referenceBulkMod.isEmpty()){
                for(int i = indBu; i < referenceBulkMod.size(); i++){
                    if(referenceBulkMod.get(i).getReferenceID() == id){
                        // found
                        points.add(referenceBulkMod.get(i));
                        indBu = i+1;
                        break;
                    }
                }
            }
            if(!referenceCellVolume.isEmpty()){
                for(int i = indCe; i < referenceCellVolume.size(); i++){
                    if(referenceCellVolume.get(i).getReferenceID() == id){
                        // found
                        points.add(referenceCellVolume.get(i));
                        indCe = i+1;
                        break;
                    }
                }
            }
            if(!referenceEnergyOrder.isEmpty()){
                for(int i = indCe; i < referenceEnergyOrder.size(); i++){
                    if(referenceEnergyOrder.get(i).getReferenceID() == id){
                        // found
                        points.add(referenceEnergyOrder.get(i));
                        indCe = i+1;
                        break;
                    }
                }
            }
            if(!referenceDeltaGauge.isEmpty()){
                for(int i = indDG; i < referenceDeltaGauge.size(); i++){
                    if(referenceDeltaGauge.get(i).getReferenceID() == id){
                        // found
                        points.add(referenceDeltaGauge.get(i));
                        indDG = i+1;
                        break;
                    }
                }
            }
            if(!referenceDensity.isEmpty()){
                for(int i = indDe; i < referenceDensity.size(); i++){
                    if(referenceDensity.get(i).getReferenceID() == id){
                        // found
                        points.add(referenceDensity.get(i));
                        indDe = i+1;
                        break;
                    }
                }
            }
            
            final SerialBatchedPropertyCalculator.PropertyBatch batch = new SerialBatchedPropertyCalculator.PropertyBatch(points,id);
            batches.add(batch);
        }
        
        return batches;
    }
    
    ArrayList<CartesianCoordinates> getAllReferenceGeoms(){
        
        final ArrayList<CartesianCoordinates> all = new ArrayList<>(referenceEnergies.size());
        
        // energy ones
        for(int i = 0; i < referenceEnergies.size(); i++){
            all.add(referenceEnergies.get(i).getReferenceInputData().c.getCartesianCoordinates());
        }
        
        // forces ones
        for(int i = 0; i < referenceForces.size(); i++){
            all.add(referenceForces.get(i).getReferenceInputData().getGeomData().c.getCartesianCoordinates());
        }
        
        // bulk modulus ones
        for(int i = 0; i < referenceBulkMod.size(); i++){
            all.add(referenceBulkMod.get(i).getReferenceInputData().getGeomData().c.getCartesianCoordinates());
        }
        
        // cell volume ones
        for(int i = 0; i < referenceCellVolume.size(); i++){
            all.add(referenceCellVolume.get(i).getReferenceInputData().getGeomData().c.getCartesianCoordinates());
        }
        
        // energy order ones
        for(int i = 0; i < referenceEnergyOrder.size(); i++){
            
            // *cough* *cough* XXX fugly
            final List<Tuple<String,ReferenceGeomData<EnergyOrder,? extends StructuralData>>> data = (List<Tuple<String,ReferenceGeomData<EnergyOrder,? extends StructuralData>>>) (Object) referenceEnergyOrder.get(i).getReferenceInputData().getAllGeomData();
            for(final Tuple<String,ReferenceGeomData<EnergyOrder,? extends StructuralData>> p : data){
                all.add(p.getObject2().c.getCartesianCoordinates());
            }
        }
        
        // delta gauge ones
        for(int i = 0; i < referenceDeltaGauge.size(); i++){
            all.add(referenceDeltaGauge.get(i).getReferenceInputData().getGeomData().c.getCartesianCoordinates());
        }
        
        // density ones
        for(int i = 0; i < referenceForces.size(); i++){
            all.add(referenceDensity.get(i).getReferenceInputData().getGeomData().c.getCartesianCoordinates());
        }
        
        return all;
    }
    
/*    <V extends StructuralData>
    ArrayList<? extends StructuralData> getAllReferenceStructuralData(){
        
        final ArrayList<V> all = new ArrayList<>(referenceEnergies.size());
        for(final GenericReferencePoint<Energy,ReferenceGeomData<Energy,?>> point : referenceEnergies){
            // we know we can cast (assuming a sane client). if not, it'll go bang.
            // XXX god is this ugly
            final GenericReferencePoint<Energy,ReferenceGeomData<Energy,V>> px = (GenericReferencePoint<Energy,ReferenceGeomData<Energy,V>>) (Object) point;
            all.add(px.getReferenceInputData().c);
        }
        
        return all;
    }*/

    @Override
    public List<String> getFormattedConfig() {
        
        final List<String> output = new LinkedList<>();
        output.add("ADAPTIVE PARAMETER CONFIG");
        output.add("output not supported yet, sorry!");
        //TODO output
        
        return output;
    }

    @Override
    public String seedFolder() {
        return paramSeedFolder;
    }

    @Override
    public org.ogolem.generic.GenericFitnessFunction<Double, AdaptiveParameters> getFitnessFunction() throws Exception {
        
        if(this.refNewton != null){
            // great, we have been initialized, return that! :-)
            return refNewton;
        }
                
        // slightly more involved: wrap a locopt around it
        final ParamLocOptFactory fac = new ParamLocOptFactory(threshLocOptParam,maxIterLocOpt,
                normParams,lowerParameterBorder, upperParameterBorder, backendDict, this);
        
        final org.ogolem.generic.GenericFitnessFunction<Double, AdaptiveParameters> locopt = fac.buildLocalOpt(whichParamLocOpt);
        
        this.refNewton = locopt;
        
        return refNewton;
    }

    @Override
    public GenericHistoryConfig getGenericHistoryConfig() throws Exception {
        
        final GenericHistoryConfig histconf = new GenericHistoryConfig();
        histconf.offset = populationSize;
        histconf.recordsToASCII = 10000; //XXX should be tunable
        histconf.recordsToSerial = 10000; //XXX should be tunable
        
        return histconf;
    }

    @Override
    public int getNumberOfGlobalSteps() throws Exception {
        return noOfGlobalSteps;
    }

    @Override
    public GenericInitializer<Double, AdaptiveParameters> getInitializer() throws Exception {
        
        final org.ogolem.generic.GenericFitnessFunction<Double,AdaptiveParameters> fit = getFitnessFunction();
        final ParameterInit init = new ParameterInit(fit,lowerParameterBorder,upperParameterBorder);
        
        return init;
    }

    @Override
    public <V extends GenericInitializer<Double, AdaptiveParameters>> TaskFactory<Double, AdaptiveParameters, V> getInitFactory() throws Exception {
        return new GenericInitTask<>();
    }

    @Override
    public GenericGlobalOptimization<Double, AdaptiveParameters> getGlobalOptimization() throws Exception {
        assert(opter != null);
        return opter;
    }

    @Override
    public <V extends GenericGlobalOptimization<Double, AdaptiveParameters>> TaskFactory<Double, AdaptiveParameters, V> getGlobalFactory() throws Exception {
        return new GenericGlobOptTask<>();
    }

    @Override
    public AdaptiveParameters getExample() throws Exception {
        assert(paramExample != null);
        return new AdaptiveParameters(paramExample);
    }

    @Override
    public IndividualReader<AdaptiveParameters> getReader() throws Exception {
        return new ParameterReader();
    }

    @Override
    public IndividualWriter<AdaptiveParameters> getWriter(String folder) throws Exception {
        return new ParameterWriter(paramsDumpFolder);
    }
    
    public
    <V extends StructuralData>
    List<ReferenceGeomData<Energy,V>> getReferenceGeomDataForEnergy(){
        
        final List<ReferenceGeomData<Energy,V>> allReferenceData = new ArrayList<>();
        for(final GenericReferencePoint<Energy,ReferenceGeomData<Energy,?>> point : referenceEnergies){
            final ReferenceGeomData<Energy,?> geom = point.getReferenceInputData();
            allReferenceData.add((ReferenceGeomData<Energy,V>) (Object)geom);
        }
        
        return allReferenceData;
    }
    
    public List<Energy> getReferenceEnergies(){
        
        final List<Energy> allEnergies = new ArrayList<>();
        for(final GenericReferencePoint<Energy,ReferenceGeomData<Energy,?>> point : referenceEnergies){
            final Energy e = point.getReferenceProperty();
            allEnergies.add(e);
        }
        
        return allEnergies;
    }

    @Override
    public boolean wantsDetailedStats() throws Exception {
        return enableDetailedStats;
    }

    @Override
    public NicheComputer<Double,AdaptiveParameters> getNicheComputer() throws Exception {
        if(doNiching){
            switch(whichNicher){
                case 0: return new StaticGridNiches(noOfNichesPerDim, this.lowerParameterBorder, this.upperParameterBorder, false);
                case 1: return new StaticGridNiches(noOfNichesPerDim, this.lowerParameterBorder, this.upperParameterBorder, true);
                default: throw new RuntimeException("Illegal nicher choice " + whichNicher);
            }
        }
        
        return null;
    }

    @Override
    public String intermediatePoolFile() {
        return outputFolder + System.getProperty("file.separator") + "IntermediateParameterPool.bin";
    }

    private static StructureDataType type = StructureDataType.Cartesian;
    
    public StructureDataType structuralDataType(){
        return type;
    }
}
