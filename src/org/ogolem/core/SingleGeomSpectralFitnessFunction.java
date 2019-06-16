/**
Copyright (c) 2015, J. M. Dieterich and B. Hartke
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
import org.ogolem.corrfunc.DipoleMomentCorrelator;
import org.ogolem.corrfunc.Folder;
import org.ogolem.corrfunc.FolderFactory;
import org.ogolem.corrfunc.SimplePropertyCorrelator;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.io.FileFilter;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.locopt.BOBYQALocOpt;
import org.ogolem.md.InMemoryTrajectory;
import org.ogolem.md.Initialization;
import org.ogolem.md.MDAlgorithm;
import org.ogolem.md.MDConfig;
import org.ogolem.md.MDSanityCheck;
import org.ogolem.md.MDSystem;
import org.ogolem.md.MDSystemRest;
import org.ogolem.properties.DipoleMoment;
import org.ogolem.properties.Energy;
import org.ogolem.spectral.FullSpectrum;

/**
 * A fitness function to locopt a geometry and subsequently obtain the fitness
 * value from a spectral fit against a reference.
 * @author Johannes Dieterich
 * @version 2015-04-28
 */
public class SingleGeomSpectralFitnessFunction implements GenericFitnessFunction<Molecule,Geometry>{

    /**
     * Prefactor to convert cm-1 for the Fourier transform. Last multiplier is
     * the speed of light in cm per picosecond.
     */
    private final static double PREFAC = Math.PI*2*2.99792458E-2;
    
    private static final long serialVersionUID = (long) 20150227;
    private static final boolean DEBUG = false;
    
    private final GenericFitnessFunction<Molecule,Geometry> opter;
    private final DipoleMomentCorrelator dipCorr;
    private final int corrFrames;
    private final int corrFramesPerBlock;
    private final int corrGap;
    private final FullSpectrum refSpec;
    private final double mdStepLength; // in fs
    private final int maxWaveNumbers; // in cm-1
    private final double fullThresh;
    private final boolean toBaseLine;
    private final boolean locOpt;
    private final int start;
    private final int end;
    private final Folder folder;
    private final MDConfig config;
    private final int dipoleStepsToSnap;
    private final double blowFacBonds;
    private final double blowFacDD;
    
    private InMemoryTrajectory<DipoleMoment> traj;
    
    private static final String FOLDERINFO = "folder.info";
    
    private static final int maxIter = 100;
    private static final double initialTrust = 2e-1;
    private static final double convergedTrust = 1e-5;
    
    
    private final double[] vCorr;
    private final double[] vCorrNorm;
    private final int[] iCorr;
    
            
    SingleGeomSpectralFitnessFunction(final GenericFitnessFunction<Molecule,Geometry> opter,
            final int corrFrames, final int corrFramesPerBlock, final int corrGap,
            final String refSpecFile, final int maxWaveNumbers,
            final double fullThresh, final boolean toBaseLine,
            final boolean locOpt, final int start, final int end,
            final int dipoleStepsToSnap, final MDConfig config, final double blowFacBonds, final double blowFacDD) throws Exception {
        System.err.println("#################");
        System.err.println(" PLEASE NOTE: THIS (SingleGeomSpectralFitnessFunction) IS VERY HACKY!");
        System.err.println("#################");
        
        this.opter = opter;
        this.corrFrames = corrFrames;
        this.corrFramesPerBlock = corrFramesPerBlock;
        this.corrGap = corrGap;
        this.fullThresh = fullThresh;
        this.maxWaveNumbers = maxWaveNumbers;
        this.toBaseLine = toBaseLine;
        this.locOpt = locOpt;
        this.start = start;
        this.end = end;
        this.dipoleStepsToSnap = dipoleStepsToSnap;
        this.blowFacBonds = blowFacBonds;
        this.blowFacDD = blowFacDD;
        
        final String[] refData = InputPrimitives.readFileIn(refSpecFile);
        final FullSpectrum fullSpec = new FullSpectrum(refData);
        this.refSpec = fullSpec.getSubSpectrumNormed(start, end, toBaseLine);
        
        this.config = config;
        this.mdStepLength = config.getStepLengthInFS();
        System.err.println("WARNING: UNCONDITIONALLY ADDING DIPOLEHACK");
        this.config.SUPERDIRTYDIPOLEHACK = true;
        this.config.superSilent = true;
        
        // the folder configuration
        final String[] confData = InputPrimitives.readFileIn(FOLDERINFO);
        String fStr = null;
        for (final String confData1 : confData) {
            fStr = confData1.trim();
            if(!(fStr.startsWith("#") || fStr.startsWith("//"))){
                break;
            }
        }
        this.folder = FolderFactory.getFolder(fStr);
       
        this.vCorr = new double[corrGap];
        this.vCorrNorm = new double[corrGap];
        this.iCorr = new int[corrGap];
        this.dipCorr = new DipoleMomentCorrelator();
    }
    
    SingleGeomSpectralFitnessFunction(final SingleGeomSpectralFitnessFunction orig){
        this.opter = orig.opter.clone();
        this.corrFrames = orig.corrFrames;
        this.corrFramesPerBlock = orig.corrFramesPerBlock;
        this.corrGap = orig.corrGap;
        this.fullThresh = orig.fullThresh;
        this.mdStepLength = orig.mdStepLength;
        this.maxWaveNumbers = orig.maxWaveNumbers;
        this.refSpec = orig.refSpec.clone(); // not really needed, but better be safe
        this.toBaseLine = orig.toBaseLine;
        this.locOpt = orig.locOpt;
        this.start = orig.start;
        this.end = orig.end;
        this.folder = orig.folder.clone();
        this.config = orig.config.clone();
        this.dipoleStepsToSnap = orig.dipoleStepsToSnap;
        this.blowFacBonds = orig.blowFacBonds;
        this.blowFacDD = orig.blowFacDD;
        
        this.vCorr = new double[corrGap];
        this.vCorrNorm = new double[corrGap];
        this.iCorr = new int[corrGap];
        this.dipCorr = orig.dipCorr.clone();
    }
    
    @Override
    public String getMyID() {
        return "SINGLE GEOMETRY SPECTRAL FITNESS FUNCTION:\n\t" + opter.getMyID();
    }
    
    @Override
    public GenericFitnessFunction<Molecule, Geometry> clone() {
        return new SingleGeomSpectralFitnessFunction(this);
    }

    @Override
    public Geometry fitness(final Geometry individual, final boolean forceOneEval) {
        
        // first optimize the structure
        final Geometry opted = opter.fitness(individual, forceOneEval);
        if(opted.getFitness() >= FixedValues.NONCONVERGEDENERGY){
            System.err.println("WARNING: Local optimization failed for " + opted.getID() + " no point in doing the spectral items.");
            return opted;
        }
        
        final Energy en = new Energy(opted.getFitness());
        opted.addProperty(en);
        
        final String stubName = "spectralfit" + individual.hashCode() + "_" + individual.getID();
        
        opted.setFitness(FixedValues.NONCONVERGEDENERGY);
        try{
            // now run the MD
            runMDTrajectory(opted.getCartesians(), opted.getBondInfo(),
                blowFacDD, blowFacBonds);
                        
            // then autocorrelate the resulting trajectory

            final int first = 0;
        
            SimplePropertyCorrelator.<DipoleMoment>autoCorrelate(first, corrFrames, corrGap, traj, vCorr, iCorr, dipCorr);
                        
            // compute the average correlation function values
            for (int i = 0; i < corrFrames; i++) {
                if (iCorr[i] > 0) {
                    vCorr[i] /= iCorr[i];
                }
            }
            
            for (int i = 0; i < corrFrames; i++) {
                vCorrNorm[i] = vCorr[i]/vCorr[0];
            }
            
            // then fold the autocorrelation function            
            final double[] extraInfo = folder.foldData(vCorr, vCorrNorm);

            final double[] fftIntensities = doFFT(vCorrNorm);
                        
            // difference integral to get the actual fitness value
            final FullSpectrum fullThisSpec = new FullSpectrum(1,1,fftIntensities);
            if(DEBUG){
                final String specCalc = fullThisSpec.getFormattedSpectrum();
                System.out.println("DEBUG: spectrum coming...");
                System.out.println(specCalc);
                OutputPrimitives.writeOut("onthefly.spec", specCalc, true);
            }
            opted.addProperty(new org.ogolem.properties.Spectrum(fullThisSpec));
            final FullSpectrum thisSpec = fullThisSpec.getSubSpectrumNormed(start, end, toBaseLine);
                        
            final Coefficient coeff = new Coefficient();
            final FitnessFunction func = new FitnessFunction(refSpec, thisSpec, fullThresh);
            double fitness;
            if(locOpt){
                // allow to optimize between 0.0 and 10.0
                final BOBYQALocOpt<Double,Coefficient> bobyqa = new BOBYQALocOpt<>(func,maxIter,initialTrust,convergedTrust);
                fitness = bobyqa.fitness(coeff, false).getFitness();
            } else {
                fitness = func.fitness(coeff, true).getFitness();
            }
            
            // set fitness
            if(DEBUG){System.out.println("DEBUG: Fitness computed is " + fitness + " with locotpt? " + locOpt);}
            opted.setFitness(fitness);
            
        } catch(Exception e){
            System.err.println("WARNING: Exception when running spectral fitness function for " + individual.getID());
            e.printStackTrace(System.err);
            opted.setFitness(FixedValues.NONCONVERGEDENERGY);
        } catch (Error err){
            System.err.println("WARNING: Error when running spectral fitness function for " + individual.getID());
            err.printStackTrace(System.err);
            opted.setFitness(FixedValues.NONCONVERGEDENERGY);
        }
                
        // cleanup stuff
        if(!DEBUG){
            try{
                cleanUp(stubName);
            } catch (Exception e){
                // System.err.println("WARNING: Couldn't remove temporary working directory " + foldName);
                e.printStackTrace(System.err);
            }
        }
        
        return opted;
    }
    
    private void runMDTrajectory(final CartesianCoordinates cartes, final BondInfo bonds,
            final double blowFacDD, final double blowFacBonds){
        
        // initialize the system
        final MDSystem system = new MDSystem(cartes,42,bonds,config);        
        final MDSystemRest rester = new MDSystemRest();
        
        Initialization.assignInitialVelocities(system, bonds, config);
        
        // always remove a full system translation/rotation
        rester.rest(system);
        
        if(config.doInitialRotTransRemoval()){
            // also remove the initial components in all molecules
            final int[] atsPerMol = cartes.getAllAtomsPerMol();
            final int noMols = atsPerMol.length;
            int offset = 0;
            int endset = 0;
            for(int mol = 0; mol < noMols; mol++){
                endset += atsPerMol[mol];
                rester.rest(system, offset, endset);
                offset += atsPerMol[mol];
            }
        }
        
        final CartesianFullBackend backend = config.getCartesianBackend();
        
        final MDSanityCheck sanity = new MDSanityCheck(config,system,blowFacDD,blowFacBonds);
        
        // propagate
        final MDAlgorithm<DipoleMoment> propagator = config.getPropagator();
        if(traj == null){
            // lazy init
            traj = new InMemoryTrajectory<>(dipoleStepsToSnap, config.getNoOfMDSteps());
        } else {
            traj.reuse();
        }
        propagator.shakeIt(system, backend, config, sanity, traj);
    }
    
    private double[] doFFT(final double[] vCorrNorm){
        
        // Fourier transforming it for every individual wavelength to get the (power) spectrum
        final double[] fftIntensities = new double[maxWaveNumbers];
        //final double[] scr = new double[3];
        for(int wave = 1; wave < maxWaveNumbers; wave++){
            final double freq = PREFAC*wave;
            double d = 0.0;
            for(int k = 0; k < corrGap; k++){
                final double time = dipoleStepsToSnap*mdStepLength/1000*k;
                //JMD: if I ever get back here, FastMath.cos was tested, checked and considered to be slower than Math.cos. End of sentence.
                d += vCorrNorm[k] * Math.cos(freq*time);
                //d += vCorrNorm[k] * FastMath.cos(freq*time, scr);
            }
            fftIntensities[wave] = d*mdStepLength*dipoleStepsToSnap;
        }
        
        return fftIntensities;
    }
    
    private static void cleanUp(final String stub) throws Exception {
        
        final File folder = new File(System.getProperty("user.dir"));
        final String[] files = folder.list(new FileFilter(stub,true,false));
        
        for(final String file : files){
            ManipulationPrimitives.remove(file);
        }
    }
    
    private static class Coefficient extends ContinuousProblem<Double>{

        private static final long serialVersionUID = (long) 20150227;
        private double coefficient;
        
        Coefficient(){
            this.coefficient = 1.0;
        }
        
        Coefficient(final Coefficient orig){
            this.coefficient = orig.coefficient;
        }
        
        @Override
        public ContinuousProblem<Double> clone() {
            return new Coefficient(this);
        }

        @Override
        public double[] getGenomeAsDouble() {
            return new double[]{coefficient};
        }

        @Override
        public Double[] getGenomeCopy() {
            return new Double[]{coefficient};
        }

        @Override
        public void setGenome(final Double[] genome) {
            this.coefficient = genome[0];
        }

        private void setCoeff(final double d){
            this.coefficient = d;
        }
    }
    
    private static class FitnessFunction implements GenericBackend<Double,Coefficient>{
        
        private static final long serialVersionUID = (long) 20150227;
        
        private final FullSpectrum ref;
        private final FullSpectrum fit;
        private final double fullThresh;
        
        FitnessFunction(final FullSpectrum ref, final FullSpectrum fit, final double fullThresh){
            this.ref = ref;
            this.fit =fit;
            this.fullThresh = fullThresh;
        }
        
        FitnessFunction(final FitnessFunction orig){
            this.ref = orig.ref.clone();
            this.fit = orig.fit.clone();
            this.fullThresh = orig.fullThresh;
        }

        @Override
        public GenericBackend<Double, Coefficient> clone() {
            return new FitnessFunction(this);
        }

        @Override
        public double gradient(final double[] currCoords, final double[] gradient, final int iteration) {
            throw new UnsupportedOperationException("Not supported yet.");
        }

        @Override
        public int numberOfActiveCoordinates(final Coefficient individual) {
            return 1;
        }

        @Override
        public double[] getActiveCoordinates(final Coefficient individual) {
            return individual.getGenomeAsDouble();
        }

        @Override
        public void resetToStable(final double[] coordinates) {
            // nothing
        }

        @Override
        public void updateActiveCoordinates(final Coefficient individual, final double[] coordinates) {
            if(DEBUG){System.out.println("DEBUG: Coeff in update " + coordinates[0]);}
            individual.setCoeff(coordinates[0]);
        }

        @Override
        public double fitness(final double[] currCoords, final int iteration) {
            
            final double[] intesRef = ref.getIntensities();
            final double[] intesComb = fit.getIntensities();
            assert(intesRef.length == intesComb.length): "intesRef " + intesRef.length + " intesComb " + intesComb.length ;
            final double dH1 = ref.getIncr();
            final double dH2 = fit.getIncr();
            assert(dH1 == dH2) : "incr 1 " + dH1 + " incr 2 " + dH2;
            
            final double coeff = currCoords[0];
            
            double fitness = 0.0;
            for(int i = 0; i < Math.min(intesRef.length,intesComb.length); i++){
                final double diff = Math.abs(intesRef[i]-coeff*intesComb[i]);
                if(DEBUG){System.out.println("DEBUG: Diff is " + diff + " in " + i);}
                if(diff <= fullThresh){continue;}
                fitness += diff*dH1;
                assert(!Double.isNaN(fitness)) : "at point " + i + " diff " + diff + ", dH1 " + dH1;
            }
        
            assert(!Double.isNaN(fitness));
            
            if(DEBUG){System.out.println("DEBUG: At iteration " + iteration + " fitness is " + fitness + " with coefficient " + coeff);}
            
            return fitness;
        }

        @Override
        public String getMyID() {
            return "Adapter for coefficient to locopt";
        }

        @Override
        public Coefficient fitness(final Coefficient individual, final boolean forceOneEval) {
            
            final double[] currCoords = individual.getGenomeAsDouble();
            final double fitness = fitness(currCoords,0);
            
            individual.setFitness(fitness);
            
            return individual;
        }

        @Override
        public BOUNDSTYPE boundariesInRepresentation(final Coefficient individual) {
            return BOUNDSTYPE.ALL;
        }

        @Override
        public void bestEstimateBoundaries(final double[] currCoords, final double[] low, final double[] high) {
            // allow to optimize between 0.0 and 10.0
            for(int i = 0; i < low.length; i++){
                low[i] = 0.0;
                high[i] = 10.0;
            }
        }
    }
}
