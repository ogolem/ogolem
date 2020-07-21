/**
Copyright (c) 2014, J. M. Dieterich
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

import java.util.ArrayList;
import java.util.List;
import static org.ogolem.core.Constants.ANGTOBOHR;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericGlobalOptimization;
import org.ogolem.generic.GenericGlobalOptimizationFactory;
import org.ogolem.generic.GenericLocOpt;
import org.ogolem.generic.GenericMultipleGlobOpt;
import org.ogolem.generic.GenericMutation;
import org.ogolem.generic.GenericSanityCheck;
import org.ogolem.generic.IndividualWriter;
import org.ogolem.heat.LocalHeatPulses;
import org.ogolem.helpers.Tuple;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Factory to build the global optimization from input.
 * @author Johannes Dieterich
 * @version 2020-07-19
 */
public class GlobOptAlgoFactory extends GenericGlobalOptimizationFactory<Molecule,Geometry>{
    
    private static final long serialVersionUID = (long) 20200225;
    private static final Logger LOG = LoggerFactory.getLogger(GlobOptAlgoFactory.class);
    private final GenericGlobalOptimizationFactory<Double,Molecule> molecularFactory;
    private final GenericFitnessFunction<Molecule,Geometry> fitness;
    private final GlobalConfig globConf;
    private final GenericSanityCheck<Molecule,Geometry> sanity;
    private final IndividualWriter<Geometry> writer;
    private final int noTries;
    private final double crossPoss;
    private final double mutPoss;
    private final boolean printBeforeFitness;
    private final double molXOverProb;
    private final double molMutProb;
    
    
    public GlobOptAlgoFactory(final GenericSanityCheck<Molecule,Geometry> sanity,
            final GenericFitnessFunction<Molecule,Geometry> fitness,
            final IndividualWriter<Geometry> writer, final double crossPoss, final double mutPoss,
            final boolean printBeforeFitness, final int noOfTries, final double molXOverProb, final double molMutProb,
            final GlobalConfig globConf) {
        
        super(globConf.noOfGlobalSteps);
        this.molecularFactory = new GlobOptAlgoFactory.MolGlobOptAlgoFactory(globConf.noOfGlobalSteps);
        this.crossPoss = crossPoss;
        this.fitness = fitness;
        this.molMutProb = molMutProb;
        this.noTries = noOfTries;
        this.sanity = sanity;
        this.printBeforeFitness = printBeforeFitness;
        this.globConf = globConf;
        this.writer = writer;
        this.mutPoss = mutPoss;
        this.molXOverProb = molXOverProb;
    }
    
    @Override
    public GenericGlobalOptimization<Molecule,Geometry> translateToGlobOpt(final String globOptString) throws Exception {
        
        LOG.debug("Trying to parse: " + globOptString + " as global optimization for cluster.");
        if(globOptString.startsWith("multiple{")){
            // parse mutliple global optimizations
            try{
                final GenericGlobalOptimization<Molecule,Geometry> multiple = parseMultiple(globOptString);
                return multiple;
            } catch(Exception e){
                throw new RuntimeException("Failure to parse multiple input.",e);
            }
        } else {
        
            // syntax: cluster{xover() mutation()} molecules{xover() mutation()}
            final String w = globOptString.trim();
            final String geomString = w.substring(0, w.indexOf("}"));
            final int start = w.indexOf("}");
            final int end = w.lastIndexOf("}");
            String molString = "";
            if (start < end) {
                molString = w.substring(start + 1, end);
            }

            if (molString.trim().isEmpty()) {
                // we do not always need it
                molString = "molecules{xover(germany:) mutation(germany:)}";
            }

            if (!geomString.startsWith("cluster{") || !molString.startsWith("molecules{")) {
                throw new RuntimeException("Syntax: cluster{xover()mutation()} molecules{xover() mutation()} not met!"
                        + " Cluster: " + geomString + " Molecules: " + molString);
            }

            final String wGeom = geomString.substring(8).trim();
            final String wMol = molString.substring(10).trim();

            if(!wGeom.startsWith("xover(")){
                throw new RuntimeException("Syntax: cluster{xover()mutation()} molecules{xover() mutation()} not met!");
            }
            final String gXStr = wGeom.substring(6,wGeom.indexOf(")"));
            final String gMStr = wGeom.substring(wGeom.lastIndexOf("(")+1,wGeom.lastIndexOf(")"));
        
            if(!wMol.startsWith("xover(")){
                throw new RuntimeException("Syntax: cluster{xover()mutation()} molecules{xover() mutation()} not met!");
            }
            final String mXStr = wMol.substring(6,wMol.indexOf(")"));
            final String mMStr = wMol.substring(wMol.lastIndexOf("(")+1,wMol.lastIndexOf(")"));
        
            try{
                final GenericCrossover<Molecule,Geometry> gX = getXOver(gXStr);
                final GenericCrossover<Double,Molecule> mX = molecularFactory.getXOver(mXStr);
                
                final GenericMutation<Molecule,Geometry> gM = getMutation(gMStr);
                final GenericMutation<Double,Molecule> mM = molecularFactory.getMutation(mMStr);
            
                final GenericGeometryDarwin globOpt = new GenericGeometryDarwin(gX,gM,sanity,fitness,
                    writer, crossPoss, mutPoss, printBeforeFitness, noTries, molXOverProb, molMutProb,
                    mX, mM);
            
                return globOpt;
                
            } catch(Exception e){    
                throw new RuntimeException("No global optimization engine found for input: " + globOptString + "."
                    + " Tried both specialized and generic ones, this is FATAL!",e);
            }
        }
    }
    
    private GenericGlobalOptimization<Molecule,Geometry> parseMultiple(final String globOptString) throws Exception {
        
        // syntax: multiple{XX%[GLOBOPTDEF1]YY%[GLOBOPTDEF2]}
        String w = globOptString.trim().substring(9, globOptString.lastIndexOf("}")).trim();
        System.out.println("Working with multiple string : " + w);
        
        boolean cont = true;
        final List<GenericGlobalOptimization<Molecule,Geometry>> allOpts = new ArrayList<>();
        final List<Double> allProbs = new ArrayList<>();
        while(cont){
            
            // get this work string
            final int percIndex = w.indexOf("%");
            final String percString = w.substring(0,percIndex).trim();
            final double perc = Double.parseDouble(percString);
            allProbs.add(perc/100);
            
            LOG.debug("Parsed " + perc + " %.");
            
            final int bracIndex = w.indexOf("]");
            
            final String thisGlobOpt = w.substring(percIndex+2,bracIndex).trim();
            LOG.debug("Parsing: " + thisGlobOpt + " as global optimization string.");
            final GenericGlobalOptimization<Molecule,Geometry> thisOpt = this.translateToGlobOpt(thisGlobOpt);
            allOpts.add(thisOpt);
            
            // store the rest
            w = w.substring(bracIndex+1).trim();
            LOG.debug("Leftover global optimization definition is " + w);
            
            // check if there is more
            if(w.isEmpty() || !w.contains("]")){cont = false;}
        }
        
        return new GenericMultipleGlobOpt<>(allOpts,allProbs);
    }
    
    @Override
    protected GenericCrossover<Molecule,Geometry> specializedXOver(final String xOverString) throws Exception {

        LOG.debug("working on crossover string " + xOverString);
        
        if(xOverString.startsWith("germany:")){
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(8));
            double gaussWidth = 0.3;
            for (final String token : tokens) {
                if (token.startsWith("gausswidth=")) {
                    gaussWidth = doubleToken("gausswidth=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (germany).");
                }
            }

            return new GermanyGeometryXOver(gaussWidth);
        } else if(xOverString.startsWith("foehr:")){
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(6));
            
            int howManyCuts = 1;
            for (final String token : tokens) {
                if (token.startsWith("noofcuts=")) {
                    howManyCuts = integerToken("noofcuts=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (foehr).");
                }
            }
            
            return new FoehrGeometryXOver(howManyCuts);
        } else if(xOverString.startsWith("lapland:")){
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(8));
            GlobOptAtomics.CUTTINGMODE cutStyle = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
            
            for (final String token : tokens) {
                if (token.startsWith("cutstyle=")) {
                    final int cut = integerToken("cutstyle=",token);
                    switch(cut){
                        case 0:
                            cutStyle = GlobOptAtomics.CUTTINGMODE.ZEROZ;
                            break;
                        case 1:
                            cutStyle = GlobOptAtomics.CUTTINGMODE.GAUSSDISTR;
                            break;
                        case 2:
                            cutStyle = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
                            break;
                        default:
                            throw new RuntimeException("Cut style must be 0 to 2. Is: " + cut);
                    }
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (lapland).");
                }
            }
            
            return new LaplandGeometryXOver(cutStyle);            
        } else if(xOverString.startsWith("aland:")){
            
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(6));
            
            final MergingPhenoXOver.MergingPhenoConfig mergeConfig = new MergingPhenoXOver.MergingPhenoConfig(true);
            CollisionDetection.CDTYPE whichCollDetect = globConf.getWhichCollisionEngine();
            GlobOptAtomics.CUTTINGMODE cutStyle = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
            double blowColl = globConf.blowFacBondDetect;
            double blowDiss = globConf.blowFacDissocDetect;
            boolean doMerging = false;
            boolean doCDDD = true;
            boolean doRandomTrans = false;

            for (final String token : tokens) {
                if (token.startsWith("colldetect=")) {
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.startsWith("blowcoll=")){
                    blowColl = doubleToken("blowcoll=",token);
                } else if(token.startsWith("blowdiss=")){
                    blowDiss = doubleToken("blowdiss=",token);
                } else if(token.startsWith("inittrust=")){
                    mergeConfig.initialTrustRegionRadius = doubleToken("inittrust=",token);
                } else if(token.startsWith("stoptrust=")){
                    mergeConfig.stoppingTrustRegionRadius = doubleToken("stoptrust=",token);
                } else if(token.startsWith("mergeiter=")){
                    mergeConfig.noIterations = integerToken("mergeiter=",token);
                } else if(token.startsWith("cutstyle=")){
                    final int cut = integerToken("cutstyle=",token);
                    switch(cut){
                        case 0:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.ZEROZ;
                            break;
                        case 1:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.GAUSSDISTR;
                            break;
                        case 2:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
                            break;
                        default:
                            throw new RuntimeException("Cut style must be 0 to 2. Is: " + cut);
                    }
                    cutStyle = mergeConfig.planeMode;
                } else if(token.startsWith("nointerpoints=")){
                    mergeConfig.numberOfInterpolationPoints = integerToken("nointerpoints=",token);
                } else if(token.startsWith("optbounds=")){
                    final String[] sa = tokenizeFourthLevel(token);
                    mergeConfig.optBounds[0][0] = Double.parseDouble(sa[0])*ANGTOBOHR;
                    mergeConfig.optBounds[1][0] = Double.parseDouble(sa[1])*ANGTOBOHR;
                    mergeConfig.optBounds[0][1] = Double.parseDouble(sa[2]);
                    mergeConfig.optBounds[1][1] = Double.parseDouble(sa[3]);
                } else if(token.equalsIgnoreCase("mergingpheno")){
                    doMerging = true;
                }  else if(token.startsWith("docddd=")){
                    doCDDD = booleanToken("docddd=",token);
                } else if(token.startsWith("dorandomtrans=")){
                    doRandomTrans = booleanToken("dorandomtrans=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (aland).");
                }
            }
            
            final CollisionDetectionEngine collDetect = new CollisionDetection(whichCollDetect);
            return new AlandGeometryXOver(doMerging,globConf.getRefNewton().getBackend(),mergeConfig,blowColl,blowDiss,collDetect,cutStyle,doCDDD,doRandomTrans);
        } else if(xOverString.startsWith("simplealand:")){
            
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(12));
            
            final MergingPhenoXOver.MergingPhenoConfig mergeConfig = new MergingPhenoXOver.MergingPhenoConfig(false);
            CollisionDetection.CDTYPE whichCollDetect = globConf.getWhichCollisionEngine();
            GlobOptAtomics.CUTTINGMODE cutStyle = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
            double blowColl = globConf.blowFacBondDetect;
            double blowDiss = globConf.blowFacDissocDetect;
            boolean doMerging = false;
            boolean doCDDD = false;

            for (final String token : tokens) {
                if (token.startsWith("colldetect=")) {
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.startsWith("blowcoll=")){
                    blowColl = doubleToken("blowcoll=",token);
                } else if(token.startsWith("blowdiss=")){
                    blowDiss = doubleToken("blowdiss=",token);
                } else if(token.startsWith("inittrust=")){
                    mergeConfig.initialTrustRegionRadius = doubleToken("inittrust=",token);
                } else if(token.startsWith("stoptrust=")){
                    mergeConfig.stoppingTrustRegionRadius = doubleToken("stoptrust=",token);
                } else if(token.startsWith("mergeiter=")){
                    mergeConfig.noIterations = integerToken("mergeiter=",token);
                } else if(token.startsWith("cutstyle=")){
                    final int cut = integerToken("cutstyle=",token);
                    switch(cut){
                        case 0:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.ZEROZ;
                            break;
                        case 1:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.GAUSSDISTR;
                            break;
                        case 2:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
                            break;
                        default:
                            throw new RuntimeException("Cut style must be 0 to 2. Is: " + cut);
                    }
                    cutStyle = mergeConfig.planeMode;
                } else if(token.startsWith("nointerpoints=")){
                    mergeConfig.numberOfInterpolationPoints = integerToken("nointerpoints=",token);
                } else if(token.startsWith("optbounds=")){
                    final String[] sa = tokenizeFourthLevel(token);
                    mergeConfig.optBounds[0][0] = Double.parseDouble(sa[0])*ANGTOBOHR;
                    mergeConfig.optBounds[1][0] = Double.parseDouble(sa[1])*ANGTOBOHR;
                    mergeConfig.optBounds[0][1] = Double.parseDouble(sa[2]);
                    mergeConfig.optBounds[1][1] = Double.parseDouble(sa[3]);
                } else if(token.equalsIgnoreCase("mergingpheno")){
                    doMerging = true;
                } else if(token.startsWith("docddd=")){
                    doCDDD = booleanToken("docddd=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (aland).");
                }
            }
            
            final CollisionDetectionEngine collDetect = new CollisionDetection(whichCollDetect);
            return new AlandGeometryXOver(doMerging,globConf.getRefNewton().getBackend(),mergeConfig,blowColl,blowDiss,collDetect,cutStyle,doCDDD,false);
        } else if(xOverString.startsWith("sweden:")){
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(7));
            GlobOptAtomics.CUTTINGMODE cutStyle = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
            
            for (final String token : tokens) {
                if (token.startsWith("cutstyle=")) {
                    final int cut = integerToken("cutstyle=",token);
                    switch(cut){
                        case 0:
                            cutStyle = GlobOptAtomics.CUTTINGMODE.ZEROZ;
                            break;
                        case 1:
                            cutStyle = GlobOptAtomics.CUTTINGMODE.GAUSSDISTR;
                            break;
                        case 2:
                            cutStyle = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
                            break;
                        default:
                            throw new RuntimeException("Cut style must be 0 to 2. Is: " + cut);
                    }
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (sweden).");
                }
            }
            
            return new SwedenGeometryXOver(cutStyle);
        } else if(xOverString.startsWith("norrbotten:")){
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(11));
            GlobOptAtomics.CUTTINGMODE cutStyle = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
            NorrbottenGeometryXOver.ROTATIONPLANE plane = NorrbottenGeometryXOver.ROTATIONPLANE.XY;

            for (final String token : tokens) {
                if (token.startsWith("cutstyle=")) {
                    final int cut = integerToken("cutstyle=",token);
                    switch(cut){
                        case 0:
                            cutStyle = GlobOptAtomics.CUTTINGMODE.ZEROZ;
                            break;
                        case 1:
                            cutStyle = GlobOptAtomics.CUTTINGMODE.GAUSSDISTR;
                            break;
                        case 2:
                            cutStyle = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
                            break;
                        default:
                            throw new RuntimeException("Cut style must be 0 to 2. Is: " + cut);
                    }
                } else if(token.startsWith("inplane=")){
                    final String pl = token.substring("inplane=".length()).trim();
                    if(pl.equalsIgnoreCase("xy")){
                        plane = NorrbottenGeometryXOver.ROTATIONPLANE.XY;
                    } else if(pl.equalsIgnoreCase("xz")){
                        plane = NorrbottenGeometryXOver.ROTATIONPLANE.XZ;
                    } else if(pl.equalsIgnoreCase("yz")){
                        plane = NorrbottenGeometryXOver.ROTATIONPLANE.YZ;
                    } else {
                        throw new RuntimeException("inplane= choice must be xy or xz or yz. Is: " + plane);
                    }
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (norrbotten).");
                }
            }

            return new NorrbottenGeometryXOver(cutStyle, plane);
        } else if(xOverString.startsWith("iceland:")){
            
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(8));
            
            final MergingPhenoXOver.MergingPhenoConfig mergeConfig = new MergingPhenoXOver.MergingPhenoConfig(true);
            CollisionDetection.CDTYPE whichCollDetect = globConf.getWhichCollisionEngine();
            double blowColl = globConf.blowFacBondDetect;
            double blowDiss = globConf.blowFacDissocDetect;
            boolean doCDDD = true;
            boolean doRandomTrans = false;

            for (final String token : tokens) {
                if (token.startsWith("colldetect=")) {
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.startsWith("blowcoll=")){
                    blowColl = doubleToken("blowcoll=",token);
                } else if(token.startsWith("blowdiss=")){
                    blowDiss = doubleToken("blowdiss=",token);
                } else if(token.startsWith("inittrust=")){
                    mergeConfig.initialTrustRegionRadius = doubleToken("inittrust=",token);
                } else if(token.startsWith("stoptrust=")){
                    mergeConfig.stoppingTrustRegionRadius = doubleToken("stoptrust=",token);
                } else if(token.startsWith("mergeiter=")){
                    mergeConfig.noIterations = integerToken("mergeiter=",token);
                } else if(token.startsWith("cutstyle=")){
                    final int cut = integerToken("cutstyle=",token);
                    switch(cut){
                        case 0:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.ZEROZ;
                            break;
                        case 1:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.GAUSSDISTR;
                            break;
                        case 2:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
                            break;
                        default:
                            throw new RuntimeException("Cut style must be 0 to 2. Is: " + cut);
                    }
                } else if(token.startsWith("nointerpoints=")){
                    mergeConfig.numberOfInterpolationPoints = integerToken("nointerpoints=",token);
                } else if(token.startsWith("optbounds=")){
                    final String[] sa = tokenizeFourthLevel(token);
                    mergeConfig.optBounds[0][0] = Double.parseDouble(sa[0])*ANGTOBOHR;
                    mergeConfig.optBounds[1][0] = Double.parseDouble(sa[1])*ANGTOBOHR;
                    mergeConfig.optBounds[0][1] = Double.parseDouble(sa[2]);
                    mergeConfig.optBounds[1][1] = Double.parseDouble(sa[3]);
                } else if(token.startsWith("docddd=")){
                    doCDDD = booleanToken("docddd=",token);
                } else if(token.startsWith("dorandomtrans=")){
                    doRandomTrans = booleanToken("dorandomtrans=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (sweden).");
                }
            }
            
            final CollisionDetectionEngine collDetect = new CollisionDetection(whichCollDetect);
            return new IcelandGeometryXOver(globConf.getRefNewton().getBackend(),blowColl,blowDiss,collDetect,mergeConfig,doCDDD,doRandomTrans);
            
        } else if(xOverString.startsWith("simpleiceland:")){
            
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(14));
            
            final MergingPhenoXOver.MergingPhenoConfig mergeConfig = new MergingPhenoXOver.MergingPhenoConfig(false);
            CollisionDetection.CDTYPE whichCollDetect = globConf.getWhichCollisionEngine();
            double blowColl = globConf.blowFacBondDetect;
            double blowDiss = globConf.blowFacDissocDetect;
            boolean doCDDD = true;

            for (final String token : tokens) {
                if (token.startsWith("colldetect=")) {
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.startsWith("blowcoll=")){
                    blowColl = doubleToken("blowcoll=",token);
                } else if(token.startsWith("blowdiss=")){
                    blowDiss = doubleToken("blowdiss=",token);
                } else if(token.startsWith("inittrust=")){
                    mergeConfig.initialTrustRegionRadius = doubleToken("inittrust=",token);
                } else if(token.startsWith("stoptrust=")){
                    mergeConfig.stoppingTrustRegionRadius = doubleToken("stoptrust=",token);
                } else if(token.startsWith("mergeiter=")){
                    mergeConfig.noIterations = integerToken("mergeiter=",token);
                } else if(token.startsWith("cutstyle=")){
                    final int cut = integerToken("cutstyle=",token);
                    switch(cut){
                        case 0:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.ZEROZ;
                            break;
                        case 1:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.GAUSSDISTR;
                            break;
                        case 2:
                            mergeConfig.planeMode = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
                            break;
                        default:
                            throw new RuntimeException("Cut style must be 0 to 2. Is: " + cut);
                    }
                } else if(token.startsWith("nointerpoints=")){
                    mergeConfig.numberOfInterpolationPoints = integerToken("nointerpoints=",token);
                } else if(token.startsWith("optbounds=")){
                    final String[] sa = tokenizeFourthLevel(token);
                    mergeConfig.optBounds[0][0] = Double.parseDouble(sa[0])*ANGTOBOHR;
                    mergeConfig.optBounds[1][0] = Double.parseDouble(sa[1])*ANGTOBOHR;
                    mergeConfig.optBounds[0][1] = Double.parseDouble(sa[2]);
                    mergeConfig.optBounds[1][1] = Double.parseDouble(sa[3]);
                } else if(token.startsWith("docddd=")){
                    doCDDD = booleanToken("docddd=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (sweden).");
                }
            }
            
            final CollisionDetectionEngine collDetect = new CollisionDetection(whichCollDetect);
            return new IcelandGeometryXOver(globConf.getRefNewton().getBackend(),blowColl,blowDiss,collDetect,mergeConfig,doCDDD,false);
            
        } else if(xOverString.startsWith("snaefellsjoekull:")) {
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(17));
            int cutStyle = 0;
            double gaussWidth = 1.0;
            boolean adjustRadius = false;
            boolean doRigidOpt = false;
            boolean checkForCollAndInflate = false;
            SnaefellsjoekullGeometryXOver.OptConfig config = new SnaefellsjoekullGeometryXOver.OptConfig();
            boolean doCD = false;
            boolean doDD = false;
            double blowFacCD = globConf.blowFacBondDetect;
            double blowFacDD = globConf.blowFacDissocDetect;
            CollisionDetection.CDTYPE whichCollDetect = globConf.whichCollisionEngine;
            boolean intermediateSanityChecks = true;
            
            for (final String token : tokens) {
                if (token.startsWith("cutstyle=")) {
                    cutStyle = integerToken("cutstyle=",token);
                    if(cutStyle < 0 || cutStyle > 2){throw new RuntimeException("cutsyle= specification out of allowed bounds (0-2).");}
                } else if(token.startsWith("gausswidth=")){
                    gaussWidth = doubleToken("gausswidth=",token);
                } else if(token.startsWith("adjustradius=")){
                    adjustRadius = booleanToken("adjustradius=",token);
                } else if(token.startsWith("dorigidopt=")){
                    doRigidOpt = booleanToken("dorigidopt=",token);
                } else if(token.startsWith("doinflate=")){
                    checkForCollAndInflate = booleanToken("doinflate=",token);
                } else if(token.startsWith("maxinflate=")){
                    config.maxInflate = doubleToken("maxinflate=",token);
                } else if(token.startsWith("inflatestep=")){
                    config.incrInflate = doubleToken("inflatestep=",token);
                } else if(token.startsWith("inittrust=")){
                    config.initialTrustRegionRadius = doubleToken("inittrust=",token);
                } else if(token.startsWith("stoptrust=")){
                    config.stoppingTrustRegionRadius = doubleToken("stoptrust=",token);
                } else if(token.startsWith("rigiditer=")){
                    config.noIterations = integerToken("rigiditer=",token);
                } else if(token.startsWith("nointerpoints=")){
                    config.numberOfInterpolationPoints = integerToken("nointerpoints=",token);
                } else if (token.startsWith("colldetect=")) {
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.startsWith("dointermediatesanity=")){
                    intermediateSanityChecks = booleanToken("dointermediatesanity=",token);
                } else if(token.startsWith("docd=")){
                    doCD = booleanToken("docd=",token);
                } else if(token.startsWith("dodd=")){
                    doDD = booleanToken("dodd=",token);
                } else if(token.startsWith("blowcoll=")){
                    blowFacCD = doubleToken("blowcoll=",token);
                } else if(token.startsWith("blowdiss=")){
                    blowFacDD = doubleToken("blowdiss=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (snaefellsjoekull).");
                }
            }
            
            SnaefellsjoekullGeometryXOver.CUTTYPE cut = SnaefellsjoekullGeometryXOver.CUTTYPE.FULLRANDOM; // == 0
            if(cutStyle == 1){
                cut = SnaefellsjoekullGeometryXOver.CUTTYPE.GAUSSDISTRIBUTED;
            } else if(cutStyle == 2){
                cut = SnaefellsjoekullGeometryXOver.CUTTYPE.INVERTEDGAUSSDISTRIBUTED;
            }
            
            final CollisionDetectionEngine cd = new CollisionDetection(whichCollDetect);
                    
            return new SnaefellsjoekullGeometryXOver(adjustRadius,cut,gaussWidth,doRigidOpt,
                    checkForCollAndInflate, doCD, doDD, blowFacCD, blowFacDD, 
                    cd, globConf.getRefNewton().getBackend(), config, intermediateSanityChecks);
        } else if(xOverString.startsWith("vinland:")){

            final double blowFac = globConf.blowFacBondDetect;

            final String subXOver = xOverString.substring(8).trim();
            final String[] major = tokenizeSecondLevel(subXOver);

            final GenericCrossover<Molecule,Geometry> xover = getXOver(major[1].trim());

            final String[] tokens = tokenizeThirdLevel(major[0].trim());

            VinlandGeometryXOver.MOVEMODE mode = VinlandGeometryXOver.MOVEMODE.FULLRANDOM;
            double movePerStep = VinlandGeometryXOver.DEFAULTINCRBOHR;
            int maxMoveTries = VinlandGeometryXOver.DEFAULTMAXMOVETRIES;
            int maxTries = VinlandGeometryXOver.DEFAULTMAXTRIES;

            for(final String token : tokens){
                if(token.startsWith("mode=")){
                    final String sub = token.substring("mode=".length()).trim();
                    if(sub.equalsIgnoreCase("fullrandom")){
                        mode = VinlandGeometryXOver.MOVEMODE.FULLRANDOM;
                    } else if(sub.equalsIgnoreCase("surface")){
                        mode = VinlandGeometryXOver.MOVEMODE.SURFACESTYLE;
                    } else {
                        throw new RuntimeException("Illegal choice " + sub + " for move mode in Vinland X-Over.");
                    }
                } else if(token.startsWith("moveperstep=")){
                    movePerStep = doubleToken("moveperstep=", token);
                } else if(token.startsWith("maxmoves=")){
                    maxMoveTries = integerToken("maxmoves=", token);
                } else if(token.startsWith("maxtries=")){
                    maxTries = integerToken("maxtries=", token);
                } else {
                    throw new RuntimeException("Unknown token: " + token + " in Vinland X-Over.");
                }
            }

            return new VinlandGeometryXOver(xover, blowFac, mode, movePerStep, maxMoveTries, maxTries);
        } else if(xOverString.startsWith("portugal:")){
            
            final String[] tokens = tokenizeThirdLevel(xOverString.substring(9));
            int noCuts = 1;
            for (final String token : tokens) {
                if (token.startsWith("nocuts=")) {
                    noCuts = integerToken("nocuts=",token);
                    if(noCuts < 0){throw new RuntimeException("Number of cuts must be positive!");}
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized xover (portugal).");
                }
            }

            return new PortugalGeometryXOver(noCuts);
        }
        
        // no specialized crossover found 
        return null;
    }
    
    @Override
    public GenericMutation<Molecule,Geometry> specializedMutation(final String mutString) throws Exception {
        
        LOG.debug("Working on mutation string " + mutString);
        
        if(mutString.startsWith("germany:")){
            
            final String[] tokens = tokenizeThirdLevel(mutString.substring(8));
                
            int whichMolMut = 0;
            double lowCOM = GermanyGeometryMutation.DEFAULTLOWCOM;
            double highCOM = GermanyGeometryMutation.DEFAULTHIGHCOM;
            double widthCOM = GermanyGeometryMutation.DEFAULTWIDTHCOM;
            double lowEuler = GermanyGeometryMutation.DEFAULTLOWEULER;
            double highEuler = GermanyGeometryMutation.DEFAULTHIGHEULER;
            double widthEuler = GermanyGeometryMutation.DEFAULTWIDTHEULER;
            
            for (final String token : tokens) {
                if (token.startsWith("lowcom=")) {
                    lowCOM = doubleToken("lowcom=",token);
                } else if (token.startsWith("highcom=")) {
                    highCOM = doubleToken("highcom=",token);
                } else if (token.startsWith("widthcom=")) {
                    widthCOM = doubleToken("widthcom=",token);
                } else if (token.startsWith("loweuler=")) {
                    lowEuler = doubleToken("loweuler=",token);
                } else if (token.startsWith("higheuler=")) {
                    highEuler = doubleToken("higheuler=",token);
                } else if (token.startsWith("widtheuler=")) {
                    widthEuler = doubleToken("widtheuler=",token);
                } else if (token.startsWith("molmut=")){
                    whichMolMut = integerToken("molmut=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (germany).");
                }
            }

            return new GermanyGeometryMutation(whichMolMut, lowCOM, highCOM, widthCOM, lowEuler, highEuler, widthEuler);
            
        } else if(mutString.startsWith("foehr:")){
            
            FoehrGeometryMutation.HOWMANY howmany = FoehrGeometryMutation.HOWMANY.MULTIPLEGAUSS;
            FoehrGeometryMutation.MODUS modus = FoehrGeometryMutation.MODUS.INCREMENT;
            
            final String[] tokens = tokenizeThirdLevel(mutString.substring(6));
            for (final String token : tokens) {
                if (token.startsWith("modus=")) {
                    final String s = token.substring(6).trim();
                    if(s.equalsIgnoreCase("increment")){
                        modus = FoehrGeometryMutation.MODUS.INCREMENT;
                    } else if(s.equalsIgnoreCase("fullyrandom")){
                        modus = FoehrGeometryMutation.MODUS.FULLYRANDOM;
                    } else {
                        throw new RuntimeException("Unknown modus " + s + " in specialized mutation (foehr).");
                    }
                } else if(token.startsWith("howmany=")){
                    final String s = token.substring(8).trim();
                    if(s.equalsIgnoreCase("one")){
                        howmany = FoehrGeometryMutation.HOWMANY.SINGLE;
                    } else if(s.equalsIgnoreCase("multiple")){
                        howmany = FoehrGeometryMutation.HOWMANY.MULTIPLE;
                    } else if(s.equalsIgnoreCase("multiplegauss")){
                        howmany = FoehrGeometryMutation.HOWMANY.MULTIPLEGAUSS;
                    } else if(s.equalsIgnoreCase("all")){
                        howmany = FoehrGeometryMutation.HOWMANY.ALL;
                    } else {
                        throw new RuntimeException("Unknown howmany " + s + " in specialized mutation (foehr).");
                    }
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (foehr).");
                }
            }
            
            return new FoehrGeometryMutation(modus,howmany);
            
        } else if(mutString.startsWith("montecarlo:")){
            
            final String[] tokens = tokenizeThirdLevel(mutString.substring(11));
            MonteCarloMutation.MOVEMODE mode = MonteCarloMutation.MOVEMODE.ONE;
            double maxMove = 0.2;
            int maxGauss = globConf.geoConfCopy().noOfParticles/2;
            double gaussWidth = 1.0;

            for (final String token : tokens) {
                if (token.startsWith("mode=")) {
                    final String s = token.substring(5).trim();
                    if (s.equalsIgnoreCase("all")) {
                        mode = MonteCarloMutation.MOVEMODE.ALL;
                    } else if (s.equalsIgnoreCase("one")) {
                        mode = MonteCarloMutation.MOVEMODE.ONE;
                    } else if (s.equalsIgnoreCase("some")) {
                        mode = MonteCarloMutation.MOVEMODE.SOME;
                    } else if (s.equalsIgnoreCase("gaussian")) {
                        mode = MonteCarloMutation.MOVEMODE.GAUSSIAN;
                    } else if (s.equalsIgnoreCase("molall")) {
                        mode = MonteCarloMutation.MOVEMODE.ALLMOL;
                    } else if (s.equalsIgnoreCase("molone")) {
                        mode = MonteCarloMutation.MOVEMODE.ONEMOL;
                    } else if (s.equalsIgnoreCase("molsome")) {
                        mode = MonteCarloMutation.MOVEMODE.SOMEMOL;
                    } else if (s.equalsIgnoreCase("molgaussian")) {
                        mode = MonteCarloMutation.MOVEMODE.GAUSSIANMOL;
                    } else {
                        throw new RuntimeException("Unknown move mode " + mode + " in MC mutation for geometries.");
                    }
                } else if (token.startsWith("maxmove=")) {
                    maxMove = doubleToken("maxmove=", token);
                } else if (token.startsWith("gaussmax=")){
                    maxGauss = integerToken("gaussmax=", token);
                } else if (token.startsWith("gausswidth=")){
                    gaussWidth = doubleToken("gausswidth=", token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (MonteCarlo).");
                }
            }

            return new MonteCarloGeometryMutation(mode, maxMove, maxGauss, gaussWidth);
                
        } else if(mutString.startsWith("extcoordmc:")){
            
            final String[] tokens = tokenizeThirdLevel(mutString.substring(11));
            MonteCarloExtOnlyGeometryMutation.MOVEMODE mode = MonteCarloExtOnlyGeometryMutation.MOVEMODE.ONE;
            double maxMoveCOM = 0.2; // 0.2 bohr max move
            double maxMoveEuler = 0.2; // 0.2 of all Eulers max
            int maxGauss = globConf.geoConfCopy().noOfParticles/2;
            double gaussWidth = 1.0;

            for (final String token : tokens) {
                if (token.startsWith("mode=")) {
                    final String s = token.substring(5).trim();
                    if (s.equalsIgnoreCase("all")) {
                        mode = MonteCarloExtOnlyGeometryMutation.MOVEMODE.ALL;
                    } else if (s.equalsIgnoreCase("one")) {
                        mode = MonteCarloExtOnlyGeometryMutation.MOVEMODE.ONE;
                    } else if (s.equalsIgnoreCase("some")) {
                        mode = MonteCarloExtOnlyGeometryMutation.MOVEMODE.SOME;
                    } else if (s.equalsIgnoreCase("gaussian")){
                        mode = MonteCarloExtOnlyGeometryMutation.MOVEMODE.GAUSSIAN;
                    } else {
                        throw new RuntimeException("Unknown move mode " + mode + " in MC extonly mutation for geometries.");
                    }
                } else if (token.startsWith("maxmovecom=")) {
                    maxMoveCOM = doubleToken("maxmovecom=", token);
                } else if (token.startsWith("maxmoveeuler=")) {
                    maxMoveEuler = doubleToken("maxmoveeuler=", token);
                } else if (token.startsWith("gaussmax=")){
                    maxGauss = integerToken("gaussmax=", token);
                } else if (token.startsWith("gausswidth=")){
                    gaussWidth = doubleToken("gausswidth=", token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (ExtCoordMonteCarlo).");
                }
            }

            return new MonteCarloExtOnlyGeometryMutation(mode, maxMoveCOM, maxMoveEuler, maxGauss, gaussWidth);
                
        } else if(mutString.startsWith("2d-extcoordmc:")){

            final String[] tokens = tokenizeThirdLevel(mutString.substring(14));
            MonteCarlo2DExtOnlyGeometryMutation.MOVEMODE mode = MonteCarlo2DExtOnlyGeometryMutation.MOVEMODE.ONE;
            double maxMoveCOM = 0.2; // 0.2 bohr max move
            double maxMoveEuler = 0.2; // 0.2 of all Eulers max
            int maxGauss = globConf.geoConfCopy().noOfParticles/2;
            double gaussWidth = 1.0;

            for (final String token : tokens) {
                if (token.startsWith("mode=")) {
                    final String s = token.substring(5).trim();
                    if (s.equalsIgnoreCase("all")) {
                        mode = MonteCarlo2DExtOnlyGeometryMutation.MOVEMODE.ALL;
                    } else if (s.equalsIgnoreCase("one")) {
                        mode = MonteCarlo2DExtOnlyGeometryMutation.MOVEMODE.ONE;
                    } else if (s.equalsIgnoreCase("some")) {
                        mode = MonteCarlo2DExtOnlyGeometryMutation.MOVEMODE.SOME;
                    } else if (s.equalsIgnoreCase("gaussian")){
                        mode = MonteCarlo2DExtOnlyGeometryMutation.MOVEMODE.GAUSSIAN;
                    } else {
                        throw new RuntimeException("Unknown move mode " + mode + " in MC 2d-extonly mutation for geometries.");
                    }
                } else if (token.startsWith("maxmovecom=")) {
                    maxMoveCOM = doubleToken("maxmovecom=", token);
                } else if (token.startsWith("maxmoveeuler=")) {
                    maxMoveEuler = doubleToken("maxmoveeuler=", token);
                } else if (token.startsWith("gaussmax=")){
                    maxGauss = integerToken("gaussmax=", token);
                } else if (token.startsWith("gausswidth=")){
                    gaussWidth = doubleToken("gausswidth=", token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (2D-ExtCoordMonteCarlo).");
                }
            }

            return new MonteCarlo2DExtOnlyGeometryMutation(mode, maxMoveCOM, maxMoveEuler, maxGauss, gaussWidth);

        } else if(mutString.startsWith("finland:")){
            
            final String[] tokens = tokenizeThirdLevel(mutString.substring(8));
            
            double blowDiss = globConf.blowFacDissocDetect;
            double blowColl = globConf.blowFacBondDetect;
            GlobOptAtomics.CUTTINGMODE cutType = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
            for (final String token : tokens) {
                if(token.startsWith("cutstyle=")){
                    final int cut = integerToken("cutstyle=",token);
                    switch(cut){
                        case 0:
                            cutType = GlobOptAtomics.CUTTINGMODE.ZEROZ;
                            break;
                        case 1:
                            cutType = GlobOptAtomics.CUTTINGMODE.GAUSSDISTR;
                            break;
                        case 2:
                            cutType = GlobOptAtomics.CUTTINGMODE.NORMDISTR;
                            break;
                        default:
                            throw new RuntimeException("Cut style must be 0 to 2. Is: " + cut);
                    }
                } else if(token.startsWith("blowdiss=")){
                    blowDiss = doubleToken("blowdiss=",token);
                } else if(token.startsWith("blowcoll=")){
                    blowColl = doubleToken("blowcoll=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (finland).");
                }
            }
            
            return new FinlandGeometryMutation(cutType,blowDiss,blowColl);
        } else if(mutString.startsWith("norway:")){
            
            final String[] tokens = tokenizeThirdLevel(mutString.substring(7));
            
            double blowDiss = globConf.blowFacDissocDetect;
            double blowColl = globConf.blowFacBondDetect;
            CollisionDetection.CDTYPE whichCollDetect = globConf.getWhichCollisionEngine();
            DissociationDetection.DDTYPE whichDissDetect = globConf.whichDissociationEngine;
            NorwayGeometryMutation.MUTMODE mode = NorwayGeometryMutation.MUTMODE.ASCENDING; // ascending by default
            for (final String token : tokens) {
                if(token.startsWith("colldetect=")){
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.startsWith("dissdetect=")){
                    final String sub = stringToken("dissdetect=",token);
                    whichDissDetect = DissociationDetection.parseType(sub);
                } else if(token.startsWith("blowdiss=")){
                    blowDiss = doubleToken("blowdiss=",token);
                } else if(token.startsWith("blowcoll=")){
                    blowColl = doubleToken("blowcoll=",token);
                } else if(token.startsWith("mode=")){
                    final String s = token.substring(5).trim();
                    if(s.equalsIgnoreCase("ascending")){
                        mode = NorwayGeometryMutation.MUTMODE.ASCENDING;
                    } else if(s.equalsIgnoreCase("random")){
                        mode = NorwayGeometryMutation.MUTMODE.RANDOM;
                    } else if(s.equalsIgnoreCase("size")){
                        mode = NorwayGeometryMutation.MUTMODE.BYSIZE;
                    } else{
                        throw new RuntimeException("Illegal mode " + mode + " for setting up norway.");
                    }
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (norway).");
                }
            }
            
            return new NorwayGeometryMutation(whichCollDetect,blowColl,blowDiss,whichDissDetect,mode);
        } else if(mutString.startsWith("norway2D:")){
            
            final String[] tokens = tokenizeThirdLevel(mutString.substring(7));
            
            double blowDiss = globConf.blowFacDissocDetect;
            double blowColl = globConf.blowFacBondDetect;
            CollisionDetection.CDTYPE whichCollDetect = globConf.getWhichCollisionEngine();
            DissociationDetection.DDTYPE whichDissDetect = globConf.whichDissociationEngine;
            Norway2DGeometryMutation.MUTMODE mode = Norway2DGeometryMutation.MUTMODE.ASCENDING; // ascending by default
            for (final String token : tokens) {
                if(token.startsWith("colldetect=")){
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.startsWith("dissdetect=")){
                    final String sub = stringToken("dissdetect=",token);
                    whichDissDetect = DissociationDetection.parseType(sub);
                } else if(token.startsWith("blowdiss=")){
                    blowDiss = doubleToken("blowdiss=",token);
                } else if(token.startsWith("blowcoll=")){
                    blowColl = doubleToken("blowcoll=",token);
                } else if(token.startsWith("mode=")){
                    final String s = token.substring(5).trim();
                    if(s.equalsIgnoreCase("ascending")){
                        mode = Norway2DGeometryMutation.MUTMODE.ASCENDING;
                    } else if(s.equalsIgnoreCase("random")){
                        mode = Norway2DGeometryMutation.MUTMODE.RANDOM;
                    } else if(s.equalsIgnoreCase("size")){
                        mode = Norway2DGeometryMutation.MUTMODE.BYSIZE;
                    } else{
                        throw new RuntimeException("Illegal mode " + mode + " for setting up 2D-norway.");
                    }
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (2D-norway).");
                }
            }
            
            return new Norway2DGeometryMutation(whichCollDetect,blowColl,blowDiss,whichDissDetect,mode);
        } else if(mutString.startsWith("xchangemut:")){
            
            final String[] tokens = tokenizeThirdLevel(mutString.substring(11));
            int mode = XChangeGeometryMutation.SINGLEEXCHANGEMODE;
            double gaussWidth = XChangeGeometryMutation.DEFAULTGAUSSWIDTH;
            
            for (final String token : tokens) {
                if(token.startsWith("mode=")){
                    final String s = token.substring(5).trim();
                    if(s.equalsIgnoreCase("single")){
                        mode = XChangeGeometryMutation.SINGLEEXCHANGEMODE;
                    } else if(s.equalsIgnoreCase("multiple")){
                        mode = XChangeGeometryMutation.MULTIPLEEXCHANGEMODE;
                    } else {
                        throw new RuntimeException("Unknown mode " + mode + " in xchange mutation.");
                    }
                } else if(token.startsWith("gausswidth=")){
                    gaussWidth = doubleToken("gausswidth=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (xchange).");
                }
            }
            
            return new XChangeGeometryMutation(mode,gaussWidth);
        } else if(mutString.startsWith("graphbasedmut:")){
            
            final String[] tokens = tokenizeThirdLevel(mutString.substring(14));
            
            double blowColl = globConf.blowFacBondDetect;
            double blowBonds = 3.05; //NOT the same as the one above
            CollisionDetection.CDTYPE whichCollDetect = globConf.whichCollisionEngine;
            double gridHalfLength = GraphBasedDirMut.DEFAULTGRIDSIZE;
            double gridIncr = GraphBasedDirMut.DEFAULTGRIDINCR;
            double eulerIncr = GraphBasedDirMut.DEFAULTEULERINCR;
            boolean useGridApproach = true;
            boolean useOptApproach = false;
            boolean fullyRelaxed = false;
            boolean doLocOpt = false;
            boolean doCDFirst = false;
            
            for (final String token : tokens) {
                if(token.startsWith("blowbonds=")){
                    blowBonds = doubleToken("blowBonds=",token);
                } else if(token.startsWith("blowcoll=")){
                    blowColl = doubleToken("blowcoll=",token);
                } else if(token.startsWith("colldetect=")){
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.equalsIgnoreCase("nogrid")){
                    useGridApproach = false;
                } else if(token.equalsIgnoreCase("doopt")){
                    useGridApproach = false;
                    useOptApproach = true;
                } else if(token.equalsIgnoreCase("docdfirst")){
                    doCDFirst = true;
                } else if(token.equalsIgnoreCase("fullyrelaxed")){
                    fullyRelaxed = true;
                } else if(token.startsWith("gridhalf=")){
                    gridHalfLength = doubleToken("gridhalf=",token);
                } else if(token.startsWith("gridincr=")){
                    gridIncr = doubleToken("gridincr=",token);
                } else if(token.startsWith("eulerincr=")){
                    eulerIncr = doubleToken("eulerincr=",token);
                } else if(token.equalsIgnoreCase("dolocopt")){
                    doLocOpt = true;
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (graphdirmut).");
                }
            }
            
            final CollisionDetectionEngine collDetect = new CollisionDetection(whichCollDetect);
            return new GraphBasedDirMut(blowBonds, collDetect, blowColl, globConf.getRefNewton(),
                gridHalfLength, gridIncr,eulerIncr,useGridApproach, useOptApproach,
                globConf.blowFacBondDetect, globConf.blowFacDissocDetect, doLocOpt,
                fullyRelaxed, doCDFirst);
        } else if(mutString.startsWith("advancedgdm:")){
            
            final Tuple<String,String[]> allOpts = tokenizeThirdLevelWithBraces(mutString.substring(12));
            final String[] tokens = allOpts.getObject2();
            
            final AdvancedGraphBasedDirMut.GDMConfiguration config = new AdvancedGraphBasedDirMut.GDMConfiguration();
            DirMutPointProviders.PointProvider provider = null;
            DirMutOptStrategies.PointOptStrategy strategy = null;
            
            config.blowBonds = 4.5;
            config.blowCD = globConf.blowFacBondDetect;
            
            CollisionDetection.CDTYPE whichCollDetect = globConf.whichCollisionEngine;
            boolean doLocOpt = false;
            boolean doCDFirst = false;
            boolean stratFullyRelaxed = false;
            boolean envAware = true;
            GenericLocOpt<Molecule,Geometry> locopt = globConf.refNewton;
            
            String strategyString = null;
            
            for (final String token : tokens) {
                if(token.startsWith("blowbonds=")){
                    config.blowBonds = doubleToken("blowBonds=",token);
                } else if(token.startsWith("blowcoll=")){
                    config.blowCD = doubleToken("blowcoll=",token);
                } else if(token.startsWith("colldetect=")){
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.startsWith("nomolstomove=")){
                    config.noOfMolsToMove = integerToken("nomolstomove=",token);
                } else if(token.equalsIgnoreCase("docdfirst")){
                    doCDFirst = true;
                } else if(token.equalsIgnoreCase("fullyrelaxed")){
                    config.fullOptEachMove = true;
                } else if(token.equalsIgnoreCase("movedmolsaremoveable")){
                    config.markMovedMolsUnmoveable = false;
                } else if(token.equalsIgnoreCase("dontcdeverymove")){
                    config.doCDEveryTrial = false;
                } else if(token.equalsIgnoreCase("dolocopt")){
                    doLocOpt = true;
                } else if(token.startsWith("gdmprovider=")){
                    final String sub = token.substring(12).trim();
                    provider = DirMutPointProviders.parseProvider(sub);
                } else if(token.startsWith("gdmstrategy=")){
                    final String sub = token.substring(12).trim();
                    strategyString = sub;
                } else if(token.equalsIgnoreCase("strategyfullyrelaxed=")){
                    stratFullyRelaxed = true;
                } else if(token.equalsIgnoreCase("notenvaware")){
                    envAware = false;
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (graphdirmut).");
                }
            }
            
            final String intLocString = allOpts.getObject1();
            if(intLocString != null){
                if(intLocString.startsWith("internallocopt=[")){
                    final String locoptString = intLocString.substring(16,intLocString.length()-1);
                    final LocOptFactory factory = new LocOptFactory(globConf,globConf.parameters,globConf.backendDefs); // XXX parameters may not be suitable!
                    locopt = factory.buildLocalOpt(locoptString);
                } else {
                    throw new RuntimeException("Unknown token " + intLocString + " in specialized mutation (graphdirmut) for local optimization.");
                }
            }
            
            if(provider == null || strategyString == null){
                throw new RuntimeException("Both a provider and strategy must be specified for AdvancedGDM.");
            }
            
            final CollisionDetectionEngine collDetect = new CollisionDetection(whichCollDetect);
            
            strategy = DirMutOptStrategies.parseStrategy(strategyString, collDetect, config.blowCD,
                    stratFullyRelaxed);
            
            return new AdvancedGraphBasedDirMut(config, collDetect,
                locopt, globConf.blowFacBondDetect, globConf.blowFacDissocDetect,
                doLocOpt, doCDFirst, provider,
                strategy, envAware);
        } else if(mutString.startsWith("pluggabledm:")){
            
            final Tuple<String,String[]> allOpts = tokenizeThirdLevelWithBraces(mutString.substring(12));
            final String[] tokens = allOpts.getObject2();
            
            final AdvancedGraphBasedDirMut.GDMConfiguration config = new AdvancedGraphBasedDirMut.GDMConfiguration();
            DirMutPointProviders.PointProvider provider = null;
            DirMutOptStrategies.PointOptStrategy strategy = null;
            
            config.blowBonds = 4.5;
            config.blowCD = globConf.blowFacBondDetect;
            
            CollisionDetection.CDTYPE whichCollDetect = globConf.whichCollisionEngine;
            boolean doLocOpt = false;
            boolean doCDFirst = false;
            boolean stratFullyRelaxed = false;
            GenericLocOpt<Molecule,Geometry> locopt = globConf.refNewton;
            
            String strategyString = null;
            
            for (final String token : tokens) {
                if(token.startsWith("blowbonds=")){
                    config.blowBonds = doubleToken("blowBonds=",token);
                } else if(token.startsWith("blowcoll=")){
                    config.blowCD = doubleToken("blowcoll=",token);
                } else if(token.startsWith("colldetect=")){
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.startsWith("nomolstomove=")){
                    config.noOfMolsToMove = integerToken("nomolstomove=",token);
                } else if(token.equalsIgnoreCase("docdfirst")){
                    doCDFirst = true;
                } else if(token.equalsIgnoreCase("fullyrelaxed")){
                    config.fullOptEachMove = true;
                } else if(token.equalsIgnoreCase("movedmolsaremoveable")){
                    config.markMovedMolsUnmoveable = false;
                } else if(token.equalsIgnoreCase("dontcdeverymove")){
                    config.doCDEveryTrial = false;
                } else if(token.equalsIgnoreCase("dolocopt")){
                    doLocOpt = true;
                } else if(token.startsWith("gdmprovider=")){
                    final String sub = token.substring(12).trim();
                    provider = DirMutPointProviders.parseProvider(sub);
                } else if(token.startsWith("gdmstrategy=")){
                    final String sub = token.substring(12).trim();
                    strategyString = sub;
                } else if(token.equalsIgnoreCase("strategyfullyrelaxed=")){
                    stratFullyRelaxed = true;
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (graphdirmut).");
                }
            }
            
            final String intLocString = allOpts.getObject1();
            if(intLocString != null){
                if(intLocString.startsWith("internallocopt=[")){
                    final String locoptString = intLocString.substring(16,intLocString.length()-1);
                    final LocOptFactory factory = new LocOptFactory(globConf,globConf.parameters,globConf.backendDefs); // XXX parameters may not be suitable!
                    locopt = factory.buildLocalOpt(locoptString);
                } else {
                    throw new RuntimeException("Unknown token " + intLocString + " in specialized mutation (graphdirmut) for local optimization.");
                }
            }
            
            if(provider == null || strategyString == null){
                throw new RuntimeException("Both a provider and strategy must be specified for PluggableDM.");
            }
            
            final CollisionDetectionEngine collDetect = new CollisionDetection(whichCollDetect);
            
            strategy = DirMutOptStrategies.parseStrategy(strategyString, collDetect, config.blowCD,
                    stratFullyRelaxed);
            
            return new PluggableDirMut(config, collDetect,
                locopt, globConf.blowFacBondDetect, globConf.blowFacDissocDetect,
                doLocOpt, doCDFirst, provider,
                strategy);
        } else if(mutString.startsWith("directedmut:")){
            
            final String[] tokens = tokenizeThirdLevel(mutString.substring(12));
            
            CollisionDetection.CDTYPE whichCollDetect = globConf.whichCollisionEngine;
            double blowColl = globConf.blowFacBondDetect;
            double scaleFactor = 1.2;
            double rangeModFactor = 0.5;
            
            for(final String token : tokens){
                if(token.startsWith("blowcoll=")){
                    blowColl = doubleToken("blowcoll=",token);
                } else if(token.startsWith("colldetect=")){
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.startsWith("scalefactor=")){
                    scaleFactor = doubleToken("scalefactor=",token);
                } else if(token.startsWith("rangemodfactor=")){
                    rangeModFactor = doubleToken("rangemodfactor=",token);
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (directedmut).");
                }
            }
            
            return new DirectedMutation(globConf.getRefNewton().getBackend(),whichCollDetect,blowColl,scaleFactor,rangeModFactor);
        } else if(mutString.startsWith("surfacedirmut:")){
            
            final String[] tokens = tokenizeThirdLevel(mutString.substring(14));
            
            int surfMode = SurfaceDirectedMutation.SIMPLESURFMODE;
            double blowColl = globConf.blowFacBondDetect;
            CollisionDetection.CDTYPE whichCollDetect = globConf.whichCollisionEngine;
            SurfaceDetection.SURFDETECTTYPE whichSurfDetect = SurfaceDetection.DEFAUTLSURFDETECT;
            
            for(final String token : tokens){
                if(token.startsWith("blowcoll=")){
                    blowColl = doubleToken("blowcoll=",token);
                } else if(token.startsWith("surfmode=")){
                    final String s = token.substring(9).trim();
                    if(s.equalsIgnoreCase("simple")){
                        surfMode = SurfaceDirectedMutation.SIMPLESURFMODE;
                    } else if(s.equalsIgnoreCase("triang")){
                        surfMode = SurfaceDirectedMutation.TRIANGSURFMODE;
                    } else{
                        throw new RuntimeException("Unknown surfmode " + s + " in surface mutation.");
                    }
                } else if(token.startsWith("colldetect=")){
                    final String sub = stringToken("colldetect=",token);
                    whichCollDetect = CollisionDetection.parseType(sub);
                } else if(token.startsWith("surfaceengine=")){
                    final String s = token.substring(14).trim();
                    if(s.equalsIgnoreCase("hartkelonglat")){
                        whichSurfDetect = SurfaceDetection.SURFDETECTTYPE.HARTKESURFACEDETECTLONGLAT;
                    } else if(s.equalsIgnoreCase("hartketess")){
                        whichSurfDetect = SurfaceDetection.SURFDETECTTYPE.HARTKESURFACEDETECTTESS;
                    } else{
                        throw new RuntimeException("Unknown surfaceengine " + s + " in surface mutation.");
                    }
                } else {
                    throw new RuntimeException("Unknown token " + token + " in specialized mutation (surfacedirectedmut).");
                }
            }
            
            final CollisionDetectionEngine collDetect = new CollisionDetection(whichCollDetect);
            final SurfaceDetectionEngine engine = new SurfaceDetection(whichSurfDetect);
            
            return new SurfaceDirectedMutation(engine,surfMode,blowColl,globConf.getRefNewton().getBackend(),collDetect);
        } else if(mutString.startsWith("localheat:")){

            final String[] tokens = tokenizeThirdLevel(mutString.substring("localheat:".length()).trim());

            final LocalHeatPulses.Configuration heatConfig = new LocalHeatPulses.Configuration();
            CollisionDetection.CDTYPE cdType = globConf.getWhichCollisionEngine();
            double blowBonds = globConf.blowFacBondDetect;
            double blowDissoc = globConf.blowFacDissocDetect;

            for(int x = 0; x < tokens.length; x++){
                final String token = tokens[x].trim();
                if(token.startsWith("choosemode=")){
                    final String sub = token.substring(11).trim();
                    switch(sub){
                        case "pick5": heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.PICK5; break;
                        case "10percent": heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.PERCENT10; break;
                        case "upto5": heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.UPTO5; break;
                        case "upto10percent": heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.UPTO10PERCENT; break;
                        case "onsphere": heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.ONSPHERE; break;
                        case "insphere": heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.INSPHERE; break;
                        case "incenter": heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.INCENTER; break;
                        default: throw new RuntimeException("Unknown move mode " + sub + " for local heat pulses.");
                    }
                } else if(token.startsWith("docddd=")){
                    heatConfig.doCDDD = Boolean.parseBoolean(token.substring(7).trim());
                } else if(token.startsWith("eqiter=")){
                    heatConfig.eqIter = Integer.parseInt(token.substring(7).trim());
                } else if(token.startsWith("iters=")){
                    heatConfig.iters = Integer.parseInt(token.substring(6).trim());
                } else if(token.startsWith("movemode=")){
                    final String s2 = token.substring("movemode=".length()).trim();
                    if(s2.equalsIgnoreCase("coms")){
                        heatConfig.moveMode = LocalHeatPulses.Configuration.MOVEMODE.MOVECOMS;
                    } else if(s2.equalsIgnoreCase("cartesian")){
                        heatConfig.moveMode = LocalHeatPulses.Configuration.MOVEMODE.MOVECARTESIAN;
                    } else {
                        throw new RuntimeException("Illegal move mode " + s2 + " for local heat.");
                    }
                } else if(token.startsWith("scalefac=")){
                    heatConfig.scaleFac = Double.parseDouble(token.substring(9).trim());
                } else if(token.startsWith("startampl=")){
                    heatConfig.startAmplitude = Double.parseDouble(token.substring(10).trim());
                } else if(token.startsWith("temperature=")){
                    heatConfig.temperature = Double.parseDouble(token.substring(12).trim());
                } else if(token.startsWith("usemetropolis=")){
                    heatConfig.useMetropolis = Boolean.parseBoolean(token.substring(14).trim());
                } else if(token.startsWith("starteuler=")){
                    heatConfig.startEulerStrength = Double.parseDouble(token.substring(11).trim());
                } else if(token.startsWith("sigmax=")){
                    heatConfig.sigmaX = Double.parseDouble(token.substring(7).trim());
                } else if(token.startsWith("sigmay=")){
                    heatConfig.sigmaY = Double.parseDouble(token.substring(7).trim());
                } else if(token.startsWith("sigmaz=")){
                    heatConfig.sigmaZ = Double.parseDouble(token.substring(7).trim());
                } else if(token.startsWith("sigma=")){
                    final double d = Double.parseDouble(token.substring(6).trim());
                    heatConfig.sigmaX = d;
                    heatConfig.sigmaY = d;
                    heatConfig.sigmaZ = d;
                } else if(token.equalsIgnoreCase("verbose")){
                    heatConfig.printVerbose = true;
                } else {
                    throw new RuntimeException("Unknown option token " + token + " for local heat pulses.");
                }
            }

            return new LocalHeatGeometryMutation(fitness,heatConfig,cdType,blowBonds,blowDissoc,globConf.acceptableFitness);
        }
                
        // no specialized mutation found that would qualify
        return null;
    }
    
    
    private static class MolGlobOptAlgoFactory extends GenericGlobalOptimizationFactory<Double,Molecule>{
        
        private static final long serialVersionUID = (long) 20140727;
        private MolGlobOptAlgoFactory(final long totalSteps){
            super(totalSteps);
        }
        
        @Override
        protected GenericCrossover<Double, Molecule> specializedXOver(final String xOverString) throws Exception {
            
            if(xOverString.startsWith("molgermany:")){
                
                final String[] tokens = tokenizeThirdLevel(xOverString.substring(11));
                double gaussWidth = 0.3;
                boolean doOld = false;
                for(final String token : tokens){
                    if(token.startsWith("gausswidth=")){
                        gaussWidth = doubleToken("gausswidth=",token);
                    } else if(token.startsWith("dooldstyle")){
                        doOld = true;
                    } else{
                        throw new RuntimeException("Unknown token " + token + " in specialized xover (germany).");
                    }
                }
                
                return new GermanyMoleculeXOver(doOld,gaussWidth);
            } else if(xOverString.startsWith("molportugal:")){
                
                final String[] tokens = tokenizeThirdLevel(xOverString.substring(12));
                int noCuts = 1;
                for(final String token : tokens){
                    if(token.startsWith("nocuts=")){
                        noCuts = integerToken("nocuts=",token);
                    } else{
                        throw new RuntimeException("Unknown token " + token + " in specialized xover (portugal).");
                    }
                }
                
                return new PortugalMoleculeXOver(noCuts);
            }

            // nothing specialized found
            return null;
        }

        @Override
        protected GenericMutation<Double, Molecule> specializedMutation(final String mutString) throws Exception {
            
            if(mutString.startsWith("germany:")){
                
                final String[] tokens = tokenizeThirdLevel(mutString.substring(8));
                double gaussWidth = 0.3;
                boolean doOld = false;
                for(final String token : tokens){
                    if(token.startsWith("gausswidth=")){
                        gaussWidth = doubleToken("gausswidth=",token);
                    } else if(token.startsWith("dooldstyle")){
                        doOld = true;
                    } else{
                        throw new RuntimeException("Unknown token " + token + " in specialized mutation (germany).");
                    }
                }
                
                return new GermanyMoleculeMutation(doOld,gaussWidth);
            } else if(mutString.startsWith("montecarlo:")){
                
                final String[] tokens = tokenizeThirdLevel(mutString.substring(11));
                MonteCarloMutation.MOVEMODE mode = MonteCarloMutation.MOVEMODE.ONE;
                double maxMove = 0.2;
                int gaussMax = -1;
                double gaussWidth = 1.0;

                for(final String token : tokens){
                    if(token.startsWith("mode=")){
                        final String s = token.substring(5).trim();
                        if(s.equalsIgnoreCase("all")){
                            mode = MonteCarloMutation.MOVEMODE.ALL;
                        } else if(s.equalsIgnoreCase("one")){
                            mode = MonteCarloMutation.MOVEMODE.ONE;
                        } else if(s.equalsIgnoreCase("some")){
                            mode = MonteCarloMutation.MOVEMODE.SOME;
                        } else if(s.equalsIgnoreCase("some")){
                            mode = MonteCarloMutation.MOVEMODE.GAUSSIAN;
                        } else {
                            throw new RuntimeException("Unknown move mode " + mode + " in MC mutation for molecules.");
                        }
                    } else if(token.startsWith("maxmove=")){
                        maxMove = doubleToken("maxmove=",token);
                    } else if (token.startsWith("gaussmax=")){
                        gaussMax = integerToken("gaussmax=", token);
                    } else if (token.startsWith("gausswidth=")){
                        gaussWidth = doubleToken("gausswidth=", token);
                    } else{
                        throw new RuntimeException("Unknown token " + token + " in specialized mutation (MonteCarlo).");
                    }
                }

                if(gaussMax < 0 && mode == MonteCarloMutation.MOVEMODE.GAUSSIAN){
                    throw new RuntimeException("Must set gaussmax= to positive integer if Gaussian move mode is wanted.");
                }
 
                return new MonteCarloMoleculeMutation(mode,maxMove, gaussMax, gaussWidth);
            }

            // nothing specialized found
            return null;
        }

        @Override
        public GenericGlobalOptimization<Double, Molecule> translateToGlobOpt(String globOptString) throws Exception {
            // intentionally not supported (yet?)
            throw new UnsupportedOperationException("Not supported yet.");
        }
    }
    
    private Tuple<String,String[]> tokenizeThirdLevelWithBraces(final String options){
        
        if(options == null || options.trim().isEmpty()){
            return new Tuple<>(null,new String[0]);
        }
        
        final String[] beforeStrip = tokenizeThirdLevel(options);
        
        // now see if there are braces being used and if so strip them and return that
        if(!options.contains("[") && !options.contains("]")){return new Tuple<>(null,beforeStrip);}

        final int bogusStartInd = options.indexOf("[");
        final String helper = options.substring(0,bogusStartInd);
        final int realStartInd = helper.lastIndexOf(",");
        final int endInd = options.indexOf("]");
        
        final List<String> optsWithoutBraces = new ArrayList<>();
        boolean foundOp = false;
        boolean foundEnd = false;
        for(final String opt : beforeStrip){
            if(opt.contains("[")){
                foundOp = true;
                continue;
            } else if(foundOp && !foundEnd && opt.contains("]")){
                foundEnd = true;
            } else {
                optsWithoutBraces.add(opt);
            }
        }
        
        final String[] strippedOpts = new String[optsWithoutBraces.size()];
        for(int i = 0; i < optsWithoutBraces.size(); i++){
            strippedOpts[i] = optsWithoutBraces.get(i);
        }
        
        final String bracedString = options.substring(realStartInd+1,endInd+1);
        
        return new Tuple<>(bracedString,strippedOpts);
    }
}
