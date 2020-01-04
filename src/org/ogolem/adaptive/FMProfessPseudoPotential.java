/**
Copyright (c) 2015-2018, J. M. Dieterich and B. Hartke
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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.FastMath;
import org.ogolem.adaptive.genericfitness.EnergyCalculator;
import org.ogolem.adaptive.genericfitness.GenericReferencePoint;
import org.ogolem.adaptive.genericfitness.PropertyCalculator;
import org.ogolem.adaptive.genericfitness.RefBulkModulusData;
import org.ogolem.adaptive.genericfitness.RefCellVolumeData;
import org.ogolem.adaptive.genericfitness.ReferenceDeltaGaugeData;
import org.ogolem.adaptive.genericfitness.ReferenceDensityData;
import org.ogolem.adaptive.genericfitness.ReferenceEnergyOrderData;
import org.ogolem.adaptive.genericfitness.ReferenceForcesData;
import org.ogolem.adaptive.genericfitness.ReferenceGeomData;
import org.ogolem.adaptive.genericfitness.ReferenceInputData;
import org.ogolem.adaptive.genericfitness.ReferenceStressTensorData;
import org.ogolem.core.AtomicProperties;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Constants;
import org.ogolem.core.StreamGobbler;
import org.ogolem.helpers.Tuple;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.properties.BulkModulus;
import org.ogolem.properties.CellVolume;
import org.ogolem.properties.DeltaGauge;
import org.ogolem.properties.Density;
import org.ogolem.properties.Energy;
import org.ogolem.properties.EnergyOrder;
import org.ogolem.properties.Forces;
import org.ogolem.properties.Property;
import org.ogolem.properties.StressTensor;

/**
 * An interface to re-optimize pseudo-potentials using Beatriz Gonzalez Del Rio's force matching approach.
 * @author Johannes Dieterich
 * @version 2018-03-22
 */
public class FMProfessPseudoPotential implements Adaptivable {
    
    /*
     * XXX we need to sort out how to handle a combined calc of ref data that is mixed. i.e., we need to append geom data, maybe add more ref types?
     */
    
    private static final long serialVersionUID = (long) 20170925;
    
    public static final double DEFAULTLOWERALPHA = -100.0;
    public static final double DEFAULTUPPERALPHA = +100.0;
    public static final double DEFAULTLOWERBETA = 0.0000001;
    public static final double DEFAULTUPPERBETA = +100.0;
    public static final double DEFAULTLOWERPOSITION = 0.0;
    public static final double DEFAULTUPPERPOSITION = +20.0;
    
    private static final double MEV_ANG3 = 13.6058*2*Constants.BOHRTOANG*Constants.BOHRTOANG*Constants.BOHRTOANG;
    
    public static final String RESULTFILE = "energy.dat";
    public static final String DEFAULTENERGYSCRIPT = "profess_energy_calc.sh";
    
    private final int noGaussians;
    private final double potentialCutoff;
    private final double[][] origPseudoPotential;
    private final String magicEnergyScript;
    private final boolean doMirrorGauss;
    private final boolean doOverallGauss;
    private final double overAllBeta;
    private final HashMap<Short,Short> pseudized;
    private final List<Short> nonOptAtoms;
    
    private final boolean deOscillate;
    private final double wignerRadius;
    
    public FMProfessPseudoPotential(final int noGaussians, final double potCutoff, final String origPPFile,
            final boolean doMirrorGauss, final boolean doOverallGauss, final double overAllBeta, final List<Short> nonOptAtoms,
            final boolean deOscillate, final double wignerRadius) throws Exception {
        
        assert(!deOscillate || wignerRadius > 0.0);
        
        if(deOscillate){
            throw new RuntimeException("ERROR: Deoscillation of goLPS on the fly currently not supported.");
        }
        
        this.noGaussians = noGaussians;
        this.potentialCutoff = potCutoff;
        
        // read previous pseudopotential in for caching purposes
        final String[] data = InputPrimitives.readFileIn(origPPFile);
        origPseudoPotential = new double[data.length][2]; // first the position, then the potential
        for(int i = 0; i < data.length; i++){
            final String[] sa = data[i].trim().split("\\s+");
            origPseudoPotential[i][0] = Double.parseDouble(sa[0]);
            origPseudoPotential[i][1] = Double.parseDouble(sa[1]);
        }
        
        final String tmp = System.getenv("OGO_PROFESSENERGYSCRIPT");
        if(tmp == null){
            magicEnergyScript = DEFAULTENERGYSCRIPT;
        } else {
            magicEnergyScript = tmp;
        }
        this.doMirrorGauss = doMirrorGauss;
        this.doOverallGauss = doOverallGauss;
        this.overAllBeta = overAllBeta;
        
        // setup pseudized electrons information
        this.pseudized = new HashMap<>();
        final String[] psDat = InputPrimitives.readFileIn("pseudized.aux");
        for(final String s : psDat){
            final String[] sa = s.trim().split("\\s+");
            final short atomNo = AtomicProperties.giveAtomicNumber(sa[0]);
            final short pseudE = Short.parseShort(sa[1]);
            pseudized.put(atomNo, pseudE);
        }
        
        this.nonOptAtoms = nonOptAtoms;
        this.deOscillate = deOscillate;
        this.wignerRadius = wignerRadius;
    }
    
    private FMProfessPseudoPotential(final FMProfessPseudoPotential orig){
        this.noGaussians = orig.noGaussians;
        this.origPseudoPotential = orig.origPseudoPotential.clone(); // shallow copy sufficient
        this.potentialCutoff = orig.potentialCutoff;
        this.magicEnergyScript = orig.magicEnergyScript;
        this.doMirrorGauss = orig.doMirrorGauss;
        this.doOverallGauss = orig.doOverallGauss;
        this.overAllBeta = orig.overAllBeta;
        this.pseudized = orig.pseudized; // no copy needed
        this.nonOptAtoms = orig.nonOptAtoms; // no copy needed
        this.deOscillate = orig.deOscillate;
        this.wignerRadius = orig.wignerRadius;
    }

    @Override
    public FMProfessPseudoPotential clone() {
        return new FMProfessPseudoPotential(this);
    }

    @Override
    public double energyOfStructWithParams(final CartesianCoordinates cartes, final AdaptiveParameters params, final int geomID, final BondInfo bonds) {
        
        final long paramID = params.getID();
        final String folderName = "professenergycalc_" + paramID + "_" + System.currentTimeMillis();
            
        // create folder, copy script and place pseudopotential as well as the set of cartesian coordinates
        try{
            final List<Short> nos = new ArrayList<>();
            OutputPrimitives.createAFolder(folderName);
            OutputPrimitives.copyFile(magicEnergyScript, folderName + File.separator + magicEnergyScript);
            ManipulationPrimitives.markExecutable(folderName + File.separator + magicEnergyScript);
            final String[] formXYZ = cartes.createPrintableCartesians();
            final String xyzFile = folderName + File.separator + "atompos.xyz";
            OutputPrimitives.writeOut(xyzFile, formXYZ, false);
            final short[] allAtomNos = cartes.getAllAtomNumbers();
            for(int x = 0; x < allAtomNos.length; x++){
                if(!nos.contains(allAtomNos[x])){
                    nos.add(allAtomNos[x]);
                }
            }
                
            for(final short no : nos){
                if(nonOptAtoms.contains(no)){
                    final String atom = AtomicProperties.giveAtomSymbol(no);
                    final String ppFile = folderName + File.separator + "pseudopot_" + atom + ".recpot";
                    final String oldFile = "pseudopot_" + atom + ".recpot";
                    OutputPrimitives.createLink(oldFile, ppFile);
                } else {
                    final String ppFile = folderName + File.separator + "pseudopot_" + AtomicProperties.giveAtomSymbol(no) + ".recpot";
                    printModifiedPPOut(params,ppFile,no, pseudized.get(no), true);
                }
            }
        } catch(Exception e){
            System.err.println("ERROR: Could not create folder " + folderName + " or copy script " + magicEnergyScript + " into it.");
            e.printStackTrace(System.err);
            return FixedValues.NONCONVERGEDENERGY;
        }
            
        // execute the script
        Process proc;
        try {
            final Runtime rt = Runtime.getRuntime();
            final String[] cmd = new String[]{"./" + magicEnergyScript, "XXX"};
            final File dir = new File(folderName);
            proc = rt.exec(cmd, null, dir);

            // any error message?
            final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

            // any output?
            final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

            // kick them off
            errorGobbler.start();
            outputGobbler.start();

            // any error???
            if (proc.waitFor() != 0) {
                return FixedValues.NONCONVERGEDENERGY;
            }

        } catch (Exception e) {
            System.err.println("PROFESS has a problem (energy calculation). ");
            e.printStackTrace(System.err);
            try {
                ManipulationPrimitives.remove(folderName);
            } catch (Exception e1) {
                System.err.println("WARNING: Failed to clean up PROFESS forces files. ");
                e1.printStackTrace(System.err);
            }

            return FixedValues.NONCONVERGEDENERGY;
        }

        // gather the results
        double energy;
        try {
            final String resFile = folderName + File.separator + RESULTFILE;
            final String[] dat = InputPrimitives.readFileIn(resFile);
            energy = Double.parseDouble(dat[0].trim());
        } catch (Exception e) {
            System.err.println("Failure to read the energy for " + folderName);
            e.printStackTrace(System.err);

            return FixedValues.NONCONVERGEDENERGY;
        }

        // delete the folder
        try {
            ManipulationPrimitives.remove(folderName);
        } catch (Exception e) {
            System.err.println("WARNING: Couldn't remove " + folderName + ". Continuing...");
        }

        return energy;
    }

    @Override
    public double gradientOfStructWithParams(final CartesianCoordinates cartes, final AdaptiveParameters params, final int geomID, final BondInfo bonds, final double[] grad) {
        // XXX numerical gradient for the time being
        return NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, grad);
    }

    @Override
    public double[][] minMaxBordersForParams(final AdaptiveParameters params) {
        
        // according to Beatriz, alpha has no bounds, beta must be positive and the position is flexible
        // hence, come up with wide default bounds
        final int noParams = params.getNumberOfParamters();
        final double[][] borders = new double[2][noParams];
        for(int i = 0; i < noParams; i+=3){
            borders[0][i  ] = DEFAULTLOWERALPHA;
            borders[1][i  ] = DEFAULTUPPERALPHA;
            borders[0][i+1] = DEFAULTLOWERBETA;
            borders[1][i+1] = DEFAULTUPPERBETA;
            borders[0][i+2] = DEFAULTLOWERPOSITION;
            borders[1][i+2] = DEFAULTUPPERPOSITION;
        }
        
        return borders;
    }

    @Override
    public AdaptiveParameters createInitialParameterStub(final ArrayList<CartesianCoordinates> refCartes, final String sMethod) {
        
        // per element, we have n-Gaussians
        final List<Short> knownAtomicNos = new ArrayList<>();
        for(final CartesianCoordinates c : refCartes){
            final short[] nos = c.getAllAtomNumbers();
            for(int x = 0; x < nos.length; x++){
                if(!knownAtomicNos.contains(nos[x]) && !nonOptAtoms.contains(nos[x])){ // ignore atoms that we know already or are marked to ignore
                    knownAtomicNos.add(nos[x]);
                }
            }
        }
        
        final int noElements = knownAtomicNos.size();
        if(noElements == 0){
            throw new RuntimeException("No atoms specified for force matching pseudo potential fitting.");
        } else if(noElements > 1){
            System.err.println("########################################################################");
            System.err.println("########################################################################");
            System.err.println("########################################################################");
            System.err.println("");
            System.err.println("More than one element specified. This is untested. Proceed with caution!");
            System.err.println("");
            System.err.println("########################################################################");
            System.err.println("########################################################################");
            System.err.println("########################################################################");
        }
        
        // per Gaussian, we have an alpha, a beta and a position
        final String[] atoms = new String[noElements];
        final int[] paramsPerAt = new int[noElements];
        for(int x = 0; x < noElements; x++){
            atoms[x] = AtomicProperties.giveAtomSymbol(knownAtomicNos.get(x));
            paramsPerAt[x] = 3*noGaussians;
        }
        
        final int totNoParams = 3*noElements*noGaussians;
        
        final AdaptiveParameters paramStub = new AdaptiveParameters(totNoParams, -1, atoms, paramsPerAt, sMethod);
        
        // add parameter descriptions
        final String[] descr = paramStub.getAllDescriptions();
        int c = 0;
        for(int i = 0; i < noElements; i++){
            for(int gauss = 0; gauss < noGaussians; gauss++){
                descr[c  ] = "atom " + atoms[i] + " gaussian " + gauss + " alpha";
                descr[c+1] = "atom " + atoms[i] + " gaussian " + gauss + " beta";
                descr[c+2] = "atom " + atoms[i] + " gaussian " + gauss + " position";
                c+=3;
            }
        }
        
        return paramStub;
    }

    @SuppressWarnings("unchecked")
    @Override
    public <T extends Property, V extends ReferenceInputData<T>> PropertyCalculator<T,V> getCalculatorForProperty(final T property, final V data) {
        
        // instantiate nested classes
        if(property instanceof Energy){
            if(data instanceof ReferenceGeomData){
                return (PropertyCalculator<T,V>) new EnergyCalculator(this.clone());
            }
        } else if(property instanceof BulkModulus){
            if(data instanceof RefBulkModulusData){
                return (PropertyCalculator<T,V>) this.new BulkModulusProfessCalculator();
            }
        } else if(property instanceof Forces){
            if(data instanceof ReferenceForcesData){
                return (PropertyCalculator<T,V>) this.new ForcesProfessCalculator();
            }
        } else if(property instanceof CellVolume){
            if(data instanceof RefCellVolumeData){
                return (PropertyCalculator<T,V>) this.new CellVolumeProfessCalculator();
            }
        } else if(property instanceof EnergyOrder){
            if(data instanceof ReferenceEnergyOrderData){
                final EnergyOrder eo = (EnergyOrder) property;
                return (PropertyCalculator<T,V>) this.new EnergyOrderProfessCalculator(eo.shouldConsiderRelEnergies());
            }
        } else if(property instanceof DeltaGauge){
            if(data instanceof ReferenceDeltaGaugeData){
                final DeltaGauge eo = (DeltaGauge) property;
                return (PropertyCalculator<T,V>) this.new DeltaGaugeProfessCalculator();
            }
        } else if(property instanceof Density){
            if(data instanceof ReferenceDensityData){
                final Density de = (Density) property;
                return (PropertyCalculator<T,V>) this.new DensityProfessCalculator();
            }
        } else if(property instanceof StressTensor){
            if(data instanceof ReferenceStressTensorData){
                final StressTensor st = (StressTensor) property;
                return (PropertyCalculator<T,V>) this.new StressTensorProfessCalculator();
            }
        }
        
        System.err.println("No calculator for property " + property.name() + " with data " + data.getClass().getName() + " in FMProfessPseudoPotential.");
        
        return null;
    }
    
    public void printModifiedPPOut(final AdaptiveParameters p, final String toFile,
            final short atomNo, final short pseudizedE, final boolean withCoulomb) throws Exception{
        
        if(deOscillate){
            printModifiedDeoscillatedPPOut(p, toFile, atomNo, pseudizedE, withCoulomb);
            return;
        }
        
        final NumberFormat formatter = new DecimalFormat("0.#########E0");
        
        final double[] params = p.getAllParamters();
        final String sep = System.lineSeparator();
        try(final BufferedWriter buffwriter = new BufferedWriter(new FileWriter(toFile, false))) {
            
            // oh the beauty of file formats in our domain. this one proudly brought to you by CASTEP, it seems.
            
            // some header (comment block)
            buffwriter.write("START COMMENT" + sep);
            final DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
            final Date date = new Date();
            buffwriter.write("BLPS: generated by OGOLEM on " + dateFormat.format(date) + ", developed by Beatriz" + sep);
            buffwriter.write("NONE" + sep);
            buffwriter.write("END COMMENT" + sep);
            // versioning number, of course OFDFT ignores it
            buffwriter.write("3     5" + sep);
            // cutoff for the potential
            buffwriter.write(formatter.format(potentialCutoff*Constants.ANGTOBOHR) + sep);
            
            // for each stored point
            int countVals = 0;
            for(int i = 0; i < origPseudoPotential.length; i++){
                
                final double dist = origPseudoPotential[i][0];
                
                // add the gaussian contribution
                final double gaussContr = gaussContr(params,dist);
                
                // add the electrostatic contribution
                double coulContr = 0.0;
                double convFac = 1.0;
                if(withCoulomb){
                    // assuming atomNo to be the valence electron number
                    coulContr = (i == 0) ? 0.0 : -4*Math.PI*(atomNo-pseudizedE)/(dist*dist);
                    convFac = MEV_ANG3;
                }
                
                double overAllGauss  = 1.0;
                if(doOverallGauss){
                    overAllGauss *= FastMath.exp(-overAllBeta*dist*dist);
                }
                
                final double ppVal = (origPseudoPotential[i][1] + overAllGauss*gaussContr + coulContr)*convFac;
                
                // always three pseudopotential values in one line. no distance. cause potatoe.                
                final String outp = formatter.format(ppVal) + "  ";
                buffwriter.write(outp);
                
                countVals++;
                if(countVals == 3 || i == origPseudoPotential.length-1){
                    buffwriter.write(sep);
                    countVals = 0;
                }
            }
            
            // and some footer (the magic number is 1000)
            buffwriter.write("1000" + sep);
        } catch (IOException e) {
            throw e;
        }
        
        if(withCoulomb){
            
            // also print out the Columb free LPS
            try(final BufferedWriter buffwriter = new BufferedWriter(new FileWriter(toFile + "-nocoulomb", false))) {
                for (int i = 0; i < origPseudoPotential.length; i++) {

                    final double dist = origPseudoPotential[i][0];

                    // add the gaussian contribution
                    final double gaussContr = gaussContr(params, dist);
                    
                    double overAllGauss = 1.0;
                    if (doOverallGauss) {
                        overAllGauss *= FastMath.exp(-overAllBeta * dist * dist);
                    }

                    final double ppVal = (origPseudoPotential[i][1] + overAllGauss * gaussContr);

                    // always three pseudopotential values in one line. no distance. cause potatoe.
                    final String outp = formatter.format(dist) + "    " + formatter.format(ppVal) + "\n";
                    buffwriter.write(outp);
                }
            } catch (IOException e) {
                throw e;
            }
        }
    }
    
    public void printModifiedDeoscillatedPPOut(final AdaptiveParameters p, final String toFile,
            final short atomNo, final short pseudizedE, final boolean withCoulomb) throws Exception{
        
        final double[] params = p.getAllParamters();
        
        /*
         * BIG NOTE HERE: WE ARE ASSUMING THAT THE PSEUDOPOTENTIAL IS EQUIDISTANTLY
         * SPACED. OTHERWISE THE FFTs WILL JUST BE GARBAGE IN, GARBAGE OUT!
         */
        
        // it also turns out that apache commons wants the input to be a power of two.
        // hence, round up and pad
        final int next = (int) Math.pow(2, Math.ceil(Math.log(origPseudoPotential.length)/Math.log(2)));
        
        // we are in reciprocal space but all values are assumed to be real
        final Complex[] reci = new Complex[next];
        for(int i = 0; i < origPseudoPotential.length; i++){
            
            // apply the Gaussian contribution
            final double dist = origPseudoPotential[i][0];
                
            // add the gaussian contribution
            final double gaussContr = gaussContr(params,dist);
            
            double overAllGauss  = 1.0;
            if(doOverallGauss){
                overAllGauss *= FastMath.exp(-overAllBeta*dist*dist);
            }
            
            final double ppVal = origPseudoPotential[i][1] + overAllGauss*gaussContr;
            reci[i] = new Complex(ppVal,0.0);
        }
        for(int i = origPseudoPotential.length; i < next; i++){
            reci[i] = new Complex(0.0,0.0);
        }
        
        // now get it into real space
        final FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        final Complex[] real = fft.transform(reci, TransformType.INVERSE);
        
        assert(real.length == origPseudoPotential.length);
        
        // filter the oscillations beyond the Wigner-Seitz radius out
        for(int i = 0; i < origPseudoPotential.length; i++){
                
            final double dist = origPseudoPotential[i][0];
            final double ri = wignerRadius - 1.0;
            if(dist <= ri){
                // nothing, we keep the current
            } else if(dist > ri && dist <= wignerRadius){
                real[i] = new Complex(real[i].getReal()*(0.5+0.5*Math.cos((dist-ri)*Math.PI)));
            } else {
                // outside the sphere
                real[i] = new Complex(0.0,0.0);
            }
        }
        
        // get back into reciprocal space
        final Complex[] reciDone = fft.transform(real, TransformType.FORWARD);
        
        assert(reciDone.length == origPseudoPotential.length);
        
        // apply Coulomb, if applicable, and print it out
        
        final NumberFormat formatter = new DecimalFormat("0.#########E0");
        
        final String sep = System.lineSeparator();
        try(final BufferedWriter buffwriter = new BufferedWriter(new FileWriter(toFile, false))) {
            
            // oh the beauty of file formats in our domain. this one proudly brought to you by CASTEP, it seems.
            
            // some header (comment block)
            buffwriter.write("START COMMENT" + sep);
            final DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
            final Date date = new Date();
            buffwriter.write("BLPS: generated by OGOLEM on " + dateFormat.format(date) + ", developed by Beatriz" + sep);
            buffwriter.write("NONE" + sep);
            buffwriter.write("END COMMENT" + sep);
            // versioning number, of course OFDFT ignores it
            buffwriter.write("3     5" + sep);
            // cutoff for the potential
            buffwriter.write(formatter.format(potentialCutoff*Constants.ANGTOBOHR) + sep);
            
            // for each stored point
            int countVals = 0;
            for(int i = 0; i < origPseudoPotential.length; i++){
                
                final double dist = origPseudoPotential[i][0];
                
                // add the electrostatic contribution
                double coulContr = 0.0;
                double convFac = 1.0;
                if(withCoulomb){
                    // assuming atomNo to be the valence electron number
                    coulContr = (i == 0) ? 0.0 : -4*Math.PI*(atomNo-pseudizedE)/(dist*dist);
                    convFac = MEV_ANG3;
                }
                
                final double ppVal = (reciDone[i].getReal() + coulContr)*convFac;
                
                // always three pseudopotential values in one line. no distance. cause potatoe.                
                final String outp = formatter.format(ppVal) + "  ";
                buffwriter.write(outp);
                
                countVals++;
                if(countVals == 3 || i == origPseudoPotential.length-1){
                    buffwriter.write(sep);
                    countVals = 0;
                }
            }
            
            // and some footer (the magic number is 1000)
            buffwriter.write("1000" + sep);
        } catch (IOException e) {
            throw e;
        }
    }
    
    private double gaussContr(final double[] params, final double point){
        
        double contr = 0;
        for(int i = 0; i < params.length; i+=3){
            final double distTerm = (point-params[i+2]);
            contr += params[i]*FastMath.exp(-params[i+1]*distTerm*distTerm);
        }
        
        if(doMirrorGauss){
            for(int i = 0; i < params.length; i+=3){
                final double distTerm = (-point-params[i+2]);
                contr += params[i]*FastMath.exp(-params[i+1]*distTerm*distTerm);
            }
        }
        
        return contr;
    }

    @SuppressWarnings("unchecked")
    @Override
    public List<? extends Property> runAllPropertyCalcs(final AdaptiveParameters params, final List<GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>>> referencePoints) {
        
        // we can calculate those:
        int enTerm = -1;
        int foTerm = -1;
        int buTerm = -1;
        int ceTerm = -1;
        int eoTerm = -1;
        int dgTerm = -1;
        int deTerm = -1;
        int stTerm = -1;

        int c = 0;
        for(final GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>> refPoint : referencePoints){
            final Property p = refPoint.getReferenceProperty();
            if(p instanceof Energy){
                enTerm = c;
            } else if(p instanceof Forces){
                foTerm = c;
            } else if(p instanceof BulkModulus){
                buTerm = c;
            } else if(p instanceof CellVolume){
                ceTerm = c;
            } else if(p instanceof EnergyOrder){
                eoTerm = c;
            } else if(p instanceof DeltaGauge){
                dgTerm = c;
            } else if(p instanceof Density){
                deTerm = c;
            } else if(p instanceof StressTensor){
                stTerm = c;
            } else {
                // error
                throw new RuntimeException("Unknown property type " + p.printableProperty() + ".");
            }
            c++;
        }
        
        final List<Property> allProps = new ArrayList<>();
        
        // we know that we can calculate energy, bulkmod and cellvol together (if the later are true)
        // forces must be calculated individually and energy w/o either bulkmod or cellvol as well
        // same applies to energy order, we are currently not caching results
        if(foTerm >= 0){
            final Forces ref = (Forces) referencePoints.get(foTerm).getReferenceProperty();
            final ReferenceForcesData<CartesianCoordinates> dat = (ReferenceForcesData<CartesianCoordinates>) referencePoints.get(foTerm).getReferenceInputData();
            final PropertyCalculator<Forces,ReferenceForcesData<CartesianCoordinates>> forceCalc = getCalculatorForProperty(ref, dat);
            
            final Forces forces = forceCalc.calculateProperty(params, dat);
            allProps.add(forces);
        }
        
        if(eoTerm >= 0){
            final EnergyOrder ref = (EnergyOrder) referencePoints.get(eoTerm).getReferenceProperty();
            final ReferenceEnergyOrderData<CartesianCoordinates> dat = (ReferenceEnergyOrderData<CartesianCoordinates>) referencePoints.get(eoTerm).getReferenceInputData();
            final PropertyCalculator<EnergyOrder,ReferenceEnergyOrderData<CartesianCoordinates>> eOrderCalc = getCalculatorForProperty(ref, dat);
            
            final EnergyOrder eOrder = eOrderCalc.calculateProperty(params, dat);
            allProps.add(eOrder);
        }
        
        if(enTerm >= 0 && !(buTerm >= 0 && ceTerm >= 0)){
            final Energy ref = (Energy) referencePoints.get(enTerm).getReferenceProperty();
            final ReferenceGeomData<Energy,CartesianCoordinates> dat = (ReferenceGeomData<Energy,CartesianCoordinates>) referencePoints.get(enTerm).getReferenceInputData();
            final PropertyCalculator<Energy,ReferenceGeomData<Energy,CartesianCoordinates>> enCalc = getCalculatorForProperty(ref, dat);
            
            final Energy energy = enCalc.calculateProperty(params, dat);
            allProps.add(energy);
        }
        
        if(dgTerm >= 0){
            final DeltaGauge ref = (DeltaGauge) referencePoints.get(dgTerm).getReferenceProperty();
            final ReferenceDeltaGaugeData<CartesianCoordinates> dat = (ReferenceDeltaGaugeData<CartesianCoordinates>) referencePoints.get(dgTerm).getReferenceInputData();
            final PropertyCalculator<DeltaGauge,ReferenceDeltaGaugeData<CartesianCoordinates>> deltaGaugeCalc = getCalculatorForProperty(ref, dat);
            
            final DeltaGauge deltaGauge = deltaGaugeCalc.calculateProperty(params, dat);
            allProps.add(deltaGauge);
        }
        
        if(deTerm >= 0){
            final Density ref = (Density) referencePoints.get(deTerm).getReferenceProperty();
            final ReferenceDensityData<CartesianCoordinates> dat = (ReferenceDensityData<CartesianCoordinates>) referencePoints.get(deTerm).getReferenceInputData();
            final PropertyCalculator<Density,ReferenceDensityData<CartesianCoordinates>> densityCalc = getCalculatorForProperty(ref, dat);
            
            final Density density = densityCalc.calculateProperty(params, dat);
            allProps.add(density);
        }
        
        if(stTerm >= 0){
            final StressTensor ref = (StressTensor) referencePoints.get(stTerm).getReferenceProperty();
            final ReferenceStressTensorData<CartesianCoordinates> dat = (ReferenceStressTensorData<CartesianCoordinates>) referencePoints.get(stTerm).getReferenceInputData();
            final PropertyCalculator<StressTensor,ReferenceStressTensorData<CartesianCoordinates>> stressCalc = getCalculatorForProperty(ref, dat);
            
            final StressTensor stress = stressCalc.calculateProperty(params, dat);
            allProps.add(stress);
        }
        
        if(buTerm >= 0 || ceTerm >= 0){
            // the combined calculation: INTRINSICALLY ASSUMING THAT THE REFERENCE POINTS ARE INDEED EQUIVALENT!
            String sym = "N/A";
            short atomNo = -1;
            if(buTerm >= 0){
                final RefBulkModulusData<CartesianCoordinates> refBM = (RefBulkModulusData<CartesianCoordinates>) referencePoints.get(buTerm).getReferenceInputData();
                sym = refBM.getSymmetry();
                atomNo = refBM.getAtomNo();
            } else {
                final RefCellVolumeData<CartesianCoordinates> refCV = (RefCellVolumeData<CartesianCoordinates>) referencePoints.get(ceTerm).getReferenceInputData();
                sym = refCV.getSymmetry();
                atomNo = refCV.getAtomNo();
            }
            
            final Tuple3D<Energy,BulkModulus,CellVolume> tup = runCombinedCalc(sym,atomNo,params);
            if(enTerm >= 0){
                allProps.add(tup.getObject1());
            }
            
            if(buTerm >= 0){
                allProps.add(tup.getObject2());
            }
            
            if(ceTerm >= 0){
                allProps.add(tup.getObject3());
            }
        }
        
        return allProps;
    }
    
    private static final String MAGICCOMBINEDSCRIPT = "profess_combined_calc.sh";
    private static final String RESULTFILEBULK = "bulk.modulus";
    private static final String RESULTFILECELLVOL = "cell.volume";
    
    private Tuple3D<Energy,BulkModulus,CellVolume> runCombinedCalc(final String crystalLattice, final short atomNo, final AdaptiveParameters p){
        
        // the idea here is simple: we create a directory, copy a script there (with a magic name, sorry)
        // and execute it with the crystal structure name as the argument. We assume the script
        // to return the bulk modulus (in the same units as the reference) into a file called bulk.modulus
        // the cell volume (in the same units as the reference) into a file called cell.colume
        // and the energy (in the same units as the reference) into a file called energy.dat
        final long paramID = p.getID();
        final String folderName = "professcombcalc_" + paramID + "_" + System.currentTimeMillis();

        // create folder, copy script and place pseudopotential
        try {
            OutputPrimitives.createAFolder(folderName);
            OutputPrimitives.copyFile(MAGICCOMBINEDSCRIPT, folderName + File.separator + MAGICCOMBINEDSCRIPT);
            ManipulationPrimitives.markExecutable(folderName + File.separator + MAGICCOMBINEDSCRIPT);
            final String ppFile = folderName + File.separator + "pseudopot.recpot";
            printModifiedPPOut(p, ppFile, atomNo, pseudized.get(atomNo), true);
        } catch (Exception e) {
            System.err.println("ERROR: Could not create folder " + folderName + " or copy script " + MAGICCOMBINEDSCRIPT + " into it.");
            e.printStackTrace(System.err);
            return new Tuple3D<>(new Energy(FixedValues.NONCONVERGEDENERGY),
                    new BulkModulus(BulkModulus.DEFAULTBULKMODULUS),
                    new CellVolume(0.0));
        }
            
        // execute the script
        Process proc;
        try {
            final Runtime rt = Runtime.getRuntime();
            final String[] cmd = new String[]{"./" + MAGICCOMBINEDSCRIPT, crystalLattice};
            final File dir = new File(folderName);
            proc = rt.exec(cmd, null, dir);

            // any error message?
            final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

            // any output?
            final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

            // kick them off
            errorGobbler.start();
            outputGobbler.start();

            // any error???
            if (proc.waitFor() != 0) {
                return new Tuple3D<>(new Energy(FixedValues.NONCONVERGEDENERGY),
                        new BulkModulus(BulkModulus.DEFAULTBULKMODULUS),
                        new CellVolume(0.0));
            }

        } catch (Exception e) {
            System.err.println("PROFESS has a problem (bulk modulus calculation). ");
            e.printStackTrace(System.err);
            try {
                ManipulationPrimitives.remove(folderName);
            } catch (Exception e1) {
                System.err.println("WARNING: Failed to clean up PROFESS bulk modulus files. ");
                e1.printStackTrace(System.err);
            }

            return new Tuple3D<>(new Energy(FixedValues.NONCONVERGEDENERGY),
                    new BulkModulus(BulkModulus.DEFAULTBULKMODULUS),
                    new CellVolume(0.0));
        }
            
        // gather the results
        BulkModulus mod;
        CellVolume cell;
        Energy en;
        try {
            final String resFile = folderName + File.separator + RESULTFILEBULK;
            final String[] dat = InputPrimitives.readFileIn(resFile);
            final String resLine = dat[0].trim();
            final double rawVal = Double.parseDouble(resLine);
            mod = new BulkModulus(rawVal);
        } catch (Exception e) {
            System.err.println("Failure to get the calculated bulk modulus for " + folderName);
            e.printStackTrace(System.err);

            return new Tuple3D<>(new Energy(FixedValues.NONCONVERGEDENERGY),
                    new BulkModulus(BulkModulus.DEFAULTBULKMODULUS),
                    new CellVolume(0.0));
        }
        
        try {
            final String resFile = folderName + File.separator + RESULTFILECELLVOL;
            final String[] dat = InputPrimitives.readFileIn(resFile);
            final String resLine = dat[0].trim();
            final double rawVal = Double.parseDouble(resLine);
            cell = new CellVolume(rawVal);
        } catch (Exception e) {
            System.err.println("Failure to get the calculated cell volume for " + folderName);
            e.printStackTrace(System.err);

            return new Tuple3D<>(new Energy(FixedValues.NONCONVERGEDENERGY),
                    new BulkModulus(BulkModulus.DEFAULTBULKMODULUS),
                    new CellVolume(0.0));
        }
        
        try {
            final String resFile = folderName + File.separator + RESULTFILE;
            final String[] dat = InputPrimitives.readFileIn(resFile);
            final String resLine = dat[0].trim();
            final double rawVal = Double.parseDouble(resLine);
            en = new Energy(rawVal);
        } catch (Exception e) {
            System.err.println("Failure to get the calculated energy for " + folderName);
            e.printStackTrace(System.err);

            return new Tuple3D<>(new Energy(FixedValues.NONCONVERGEDENERGY),
                    new BulkModulus(BulkModulus.DEFAULTBULKMODULUS),
                    new CellVolume(0.0));
        }
            
        // delete the folder
        try {
            ManipulationPrimitives.remove(folderName);
        } catch (Exception e) {
            System.err.println("WARNING: Couldn't remove " + folderName + ". Continuing...");
        }

        return new Tuple3D<>(en, mod, cell);
    }
    
    class BulkModulusProfessCalculator implements PropertyCalculator<BulkModulus,RefBulkModulusData<CartesianCoordinates>> {

        private static final long serialVersionUID = (long) 20150801;
        
        public static final String RESULTFILE = "bulk.modulus";
        public static final String DEFAULTBULKSCRIPT = "profess_bulk_calc.sh";
        
        private final String magicBulkScript;
        
        BulkModulusProfessCalculator(){
            
            // default: profess_bulk_calc.sh
            final String tmp = System.getenv("OGO_PROFESSBULKSCRIPT");
            if(tmp == null){
                magicBulkScript = DEFAULTBULKSCRIPT;
            } else {
                magicBulkScript = tmp;
            }
        }
        
        private BulkModulusProfessCalculator(final BulkModulusProfessCalculator orig){
            this.magicBulkScript = orig.magicBulkScript;
        }
        
        @Override
        public BulkModulusProfessCalculator clone() {
            return new BulkModulusProfessCalculator(this);
        }

        @Override
        public BulkModulus calculateProperty(final AdaptiveParameters p, final RefBulkModulusData<CartesianCoordinates> data) {
            
            // the idea here is simple: we create a directory, copy a script there (with a magic name, sorry)
            // and execute it with the crystal structure name as the argument. We assume the script
            // to return the bulk modulus (in the same units as the reference) into a file called bulk.modulus
            
            final long paramID = p.getID();
            final String folderName = "professbulkcalc_" + paramID + "_" + System.currentTimeMillis();
            
            // create folder, copy script and place pseudopotential
            try{
                OutputPrimitives.createAFolder(folderName);
                OutputPrimitives.copyFile(magicBulkScript, folderName + File.separator + magicBulkScript);
                ManipulationPrimitives.markExecutable(folderName + File.separator + magicBulkScript);
                final String ppFile = folderName + File.separator + "pseudopot.recpot";
                printModifiedPPOut(p,ppFile,data.getAtomNo(), pseudized.get(data.getAtomNo()), true);
            } catch(Exception e){
                System.err.println("ERROR: Could not create folder " + folderName + " or copy script " + magicBulkScript + " into it.");
                e.printStackTrace(System.err);
                return new BulkModulus(BulkModulus.DEFAULTBULKMODULUS);
            }
            
            // execute the script
            Process proc;
            try{
                final Runtime rt = Runtime.getRuntime();
                final String[] cmd = new String[]{"./" + magicBulkScript, data.getSymmetry()};
                final File dir = new File(folderName);
                proc = rt.exec(cmd,null,dir);

                // any error message?
                final StreamGobbler errorGobbler = new
                    StreamGobbler(proc.getErrorStream(), "ERROR");

                // any output?
                final StreamGobbler outputGobbler = new
                    StreamGobbler(proc.getInputStream(), "OUTPUT");

                // kick them off
                errorGobbler.start();
                outputGobbler.start();

                // any error???
                if(proc.waitFor() != 0){
                    return new BulkModulus(BulkModulus.DEFAULTBULKMODULUS);
                }

            } catch (Exception e) {
                System.err.println("PROFESS has a problem (bulk modulus calculation). ");
                e.printStackTrace(System.err);
                try{
                    ManipulationPrimitives.remove(folderName);
                } catch(Exception e1){
                    System.err.println("WARNING: Failed to clean up PROFESS bulk modulus files. ");
                    e1.printStackTrace(System.err);
                }
                
                return new BulkModulus(BulkModulus.DEFAULTBULKMODULUS);
            }
            
            // gather the results
            BulkModulus mod;
            try{
                final String resFile = folderName + File.separator + RESULTFILE;
                final String[] dat = InputPrimitives.readFileIn(resFile);
                final String resLine = dat[0].trim();
                final double rawVal = Double.parseDouble(resLine);
                mod = new BulkModulus(rawVal);
            } catch(Exception e){
                System.err.println("Failure to get the calculated bulk modulus for " + folderName);
                e.printStackTrace(System.err);
                
                return new BulkModulus(BulkModulus.DEFAULTBULKMODULUS);
            }
            
            // delete the folder
            try{
                ManipulationPrimitives.remove(folderName);
            } catch(Exception e){
                System.err.println("WARNING: Couldn't remove " + folderName + ". Continuing...");
            }
            
            return mod;
        }

        @Override
        public BulkModulus calculatePropertyGradient(final AdaptiveParameters p, final RefBulkModulusData<CartesianCoordinates> data,
                final double[] grad) {
            return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
        }
    }
    
    class CellVolumeProfessCalculator implements PropertyCalculator<CellVolume,RefCellVolumeData<CartesianCoordinates>> {

        private static final long serialVersionUID = (long) 20150920;
        
        public static final String RESULTFILE = "cell.volume";
        public static final String DEFAULTCELLVOLSCRIPT = "profess_cellvol_calc.sh";
        
        private final String magicCellVolScript;
        
        CellVolumeProfessCalculator(){
            
            // default: profess_cellvol_calc.sh
            final String tmp = System.getenv("OGO_PROFESSCELLVOLSCRIPT");
            if(tmp == null){
                magicCellVolScript = DEFAULTCELLVOLSCRIPT;
            } else {
                magicCellVolScript = tmp;
            }
        }
        
        private CellVolumeProfessCalculator(final CellVolumeProfessCalculator orig){
            this.magicCellVolScript = orig.magicCellVolScript;
        }
        
        @Override
        public CellVolumeProfessCalculator clone() {
            return new CellVolumeProfessCalculator(this);
        }

        @Override
        public CellVolume calculateProperty(final AdaptiveParameters p, final RefCellVolumeData<CartesianCoordinates> data) {
            
            // the idea here is simple: we create a directory, copy a script there (with a magic name, sorry)
            // and execute it with the crystal structure name as the argument. We assume the script
            // to return the cell volume (in the same units as the reference) into a file called cell.volume
            
            final long paramID = p.getID();
            final String folderName = "professcellvolcalc_" + paramID + "_" + System.currentTimeMillis();
            
            // create folder, copy script and place pseudopotential
            try{
                OutputPrimitives.createAFolder(folderName);
                OutputPrimitives.copyFile(magicCellVolScript, folderName + File.separator + magicCellVolScript);
                ManipulationPrimitives.markExecutable(folderName + File.separator + magicCellVolScript);
                final String ppFile = folderName + File.separator + "pseudopot.recpot";
                printModifiedPPOut(p,ppFile,data.getAtomNo(), pseudized.get(data.getAtomNo()), true);
            } catch(Exception e){
                System.err.println("ERROR: Could not create folder " + folderName + " or copy script " + magicCellVolScript + " into it.");
                e.printStackTrace(System.err);
                return new CellVolume(0.0);
            }
            
            // execute the script
            Process proc;
            try{
                final Runtime rt = Runtime.getRuntime();
                final String[] cmd =  new String[]{"./" + magicCellVolScript, data.getSymmetry()};
                final File dir = new File(folderName);
                proc = rt.exec(cmd,null,dir);

                // any error message?
                final StreamGobbler errorGobbler = new
                    StreamGobbler(proc.getErrorStream(), "ERROR");

                // any output?
                final StreamGobbler outputGobbler = new
                    StreamGobbler(proc.getInputStream(), "OUTPUT");

                // kick them off
                errorGobbler.start();
                outputGobbler.start();

                // any error???
                if(proc.waitFor() != 0){
                    return new CellVolume(0.0);
                }

            } catch (Exception e) {
                System.err.println("PROFESS has a problem (cell volume calculation). ");
                e.printStackTrace(System.err);
                try{
                    ManipulationPrimitives.remove(folderName);
                } catch(Exception e1){
                    System.err.println("WARNING: Failed to clean up PROFESS cell volume files. ");
                    e1.printStackTrace(System.err);
                }
                
               return new CellVolume(0.0);
            }
            
            // gather the results
            CellVolume cellVol;
            try{
                final String resFile = folderName + File.separator + RESULTFILE;
                final String[] dat = InputPrimitives.readFileIn(resFile);
                final String resLine = dat[0].trim();
                final double rawVal = Double.parseDouble(resLine);
                cellVol = new CellVolume(rawVal);
            } catch(Exception e){
                System.err.println("Failure to get the calculated bulk modulus for " + folderName);
                e.printStackTrace(System.err);
                
                return new CellVolume(0.0);
            }
            
            // delete the folder
            try{
                ManipulationPrimitives.remove(folderName);
            } catch(Exception e){
                System.err.println("WARNING: Couldn't remove " + folderName + ". Continuing...");
            }
            
            return cellVol;
        }

        @Override
        public CellVolume calculatePropertyGradient(final AdaptiveParameters p, final RefCellVolumeData<CartesianCoordinates> data,
                final double[] grad) {
            return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
        }
    }
    
    class ForcesProfessCalculator implements PropertyCalculator<Forces,ReferenceForcesData<CartesianCoordinates>>{

        private static final long serialVersionUID = (long) 20151202;
        
        public static final String RESULTFILE = "forces.out";
        public static final String DEFAULTFORCESSCRIPT = "profess_forces_calc.sh";
        
        private final String magicForcesScript;
        
        ForcesProfessCalculator(){
            // default: profess_bulk_calc.sh
            final String tmp = System.getenv("OGO_PROFESSFORCESSCRIPT");
            if(tmp == null){
                magicForcesScript = DEFAULTFORCESSCRIPT;
            } else {
                magicForcesScript = tmp;
            }
        }
        
        private ForcesProfessCalculator(final ForcesProfessCalculator orig){
            this.magicForcesScript = orig.magicForcesScript;
        }
        
        @Override
        public ForcesProfessCalculator clone() {
            return new ForcesProfessCalculator(this);
        }

        @Override
        public Forces calculateProperty(final AdaptiveParameters p, final ReferenceForcesData<CartesianCoordinates> data) {
            
            final long paramID = p.getID();
            final String folderName = "professforcescalc_" + paramID + "_" + System.currentTimeMillis();
            
            // create folder, copy script and place pseudopotential as well as the set of cartesian coordinates
            final int noAtoms = data.getGeomData().c.getNoOfAtoms();
            try{
                final List<Short> nos = new ArrayList<>();
                OutputPrimitives.createAFolder(folderName);
                OutputPrimitives.copyFile(magicForcesScript, folderName + File.separator + magicForcesScript);
                ManipulationPrimitives.markExecutable(folderName + File.separator + magicForcesScript);
                final String[] formXYZ = data.getGeomData().c.getCartesianCoordinates().createPrintableCartesians();
                final String xyzFile = folderName + File.separator + "atompos.xyz";
                OutputPrimitives.writeOut(xyzFile, formXYZ, false);
                final short[] allAtomNos = data.getGeomData().c.getAllAtomNumbers();
                for(int x = 0; x < allAtomNos.length; x++){
                    if(!nos.contains(allAtomNos[x])){
                        nos.add(allAtomNos[x]);
                    }
                }
                
                for(final short no : nos){
                    if(nonOptAtoms.contains(no)){
                        final String atom = AtomicProperties.giveAtomSymbol(no);
                        final String ppFile = folderName + File.separator + "pseudopot_" + atom + ".recpot";
                        final String oldFile = "pseudopot_" + atom + ".recpot";
                        OutputPrimitives.createLink(oldFile, ppFile);
                    } else {
                        final String ppFile = folderName + File.separator + "pseudopot_" + AtomicProperties.giveAtomSymbol(no) + ".recpot";
                        printModifiedPPOut(p,ppFile,no, pseudized.get(no), true);
                    }
                }
            } catch(Exception e){
                System.err.println("ERROR: Could not create folder " + folderName + " or copy script " + magicForcesScript + " into it.");
                e.printStackTrace(System.err);
                return Forces.getDefaultForces(noAtoms);
            }
            
            // execute the script
            Process proc;
            try{
                final Runtime rt = Runtime.getRuntime();
                final String[] cmd = new String[]{"./" + magicForcesScript, "XXX"};
                final File dir = new File(folderName);
                proc = rt.exec(cmd,null,dir);

                // any error message?
                final StreamGobbler errorGobbler = new
                    StreamGobbler(proc.getErrorStream(), "ERROR");

                // any output?
                final StreamGobbler outputGobbler = new
                    StreamGobbler(proc.getInputStream(), "OUTPUT");

                // kick them off
                errorGobbler.start();
                outputGobbler.start();

                // any error???
                if(proc.waitFor() != 0){
                    return Forces.getDefaultForces(noAtoms);
                }

            } catch (Exception e) {
                System.err.println("PROFESS has a problem (forces calculation). ");
                e.printStackTrace(System.err);
                try{
                    ManipulationPrimitives.remove(folderName);
                } catch(Exception e1){
                    System.err.println("WARNING: Failed to clean up PROFESS forces files. ");
                    e1.printStackTrace(System.err);
                }
                
                return Forces.getDefaultForces(noAtoms);
            }
            
            // gather the results
            final Forces forces;
            try{
                final String resFile = folderName + File.separator + RESULTFILE;
                final String[] dat = InputPrimitives.readFileIn(resFile);
                final double[][] forceVals = new double[3][noAtoms];
                for(int at = 0; at < noAtoms; at++){
                    final String[] line = dat[at+5].trim().split("\\s+"); // ignore the first five lines
                    forceVals[0][at] = Double.parseDouble(line[2]);
                    forceVals[1][at] = Double.parseDouble(line[3]);
                    forceVals[2][at] = Double.parseDouble(line[4]);
                }
                
                forces = new Forces(forceVals);
            } catch(Exception e){
                System.err.println("Failure to read the forces for " + folderName);
                e.printStackTrace(System.err);
                
                return Forces.getDefaultForces(noAtoms);
            }
            
            // delete the folder
            try{
                ManipulationPrimitives.remove(folderName);
            } catch(Exception e){
                System.err.println("WARNING: Couldn't remove " + folderName + ". Continuing...");
            }
            
            return forces;
        }

        @Override
        public Forces calculatePropertyGradient(final AdaptiveParameters p, final ReferenceForcesData<CartesianCoordinates> data,
                final double[] grad) {
            return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
        }
    }
    
    class EnergyOrderProfessCalculator implements PropertyCalculator<EnergyOrder,ReferenceEnergyOrderData<CartesianCoordinates>>{

        private static final long serialVersionUID = (long) 20160723;
        private static final boolean DEBUG = false;
        
        public static final String RESULTFILE = "energy.dat";
        public static final String DEFAULTENERGYSCRIPT = "profess_energy_calc.sh";
        
        private final String magicEnergyScript;
        
        private final boolean considerRelEnergies;
        
        EnergyOrderProfessCalculator(final boolean considerRealEnergies){
            // default: profess_bulk_calc.sh
            final String tmp = System.getenv("OGO_PROFESSENERGYSCRIPT");
            if(tmp == null){
                magicEnergyScript = DEFAULTENERGYSCRIPT;
            } else {
                magicEnergyScript = tmp;
            }
            this.considerRelEnergies = considerRealEnergies;
        }
        
        private EnergyOrderProfessCalculator(final EnergyOrderProfessCalculator orig){
            this.magicEnergyScript = orig.magicEnergyScript;
            this.considerRelEnergies = orig.considerRelEnergies;
        }
        
        @Override
        public EnergyOrderProfessCalculator clone() {
            return new EnergyOrderProfessCalculator(this);
        }

        @Override
        public EnergyOrder calculateProperty(final AdaptiveParameters p, final ReferenceEnergyOrderData<CartesianCoordinates> data) {
            
            final long paramID = p.getID();
            final String folderNamePrefix = "professenergycalc_" + paramID + "_" + System.currentTimeMillis() + "_";
            
            final List<Tuple<String,ReferenceGeomData<EnergyOrder,CartesianCoordinates>>> allGeoms = data.getAllGeomData();
            
            final List<Tuple<String,Double>> energies = new ArrayList<>();
            for(final Tuple<String,ReferenceGeomData<EnergyOrder,CartesianCoordinates>> tup : allGeoms){
                
                final String key = tup.getObject1();
                final ReferenceGeomData<EnergyOrder,CartesianCoordinates> geom = tup.getObject2();
            
                final String folderName = folderNamePrefix + key;
                final CartesianCoordinates cartes = geom.c;
                
                // create folder, copy script and place pseudopotential as well as the set of cartesian coordinates
                try{
                    final List<Short> nos = new ArrayList<>();
                    OutputPrimitives.createAFolder(folderName);
                    OutputPrimitives.copyFile(magicEnergyScript, folderName + File.separator + magicEnergyScript);
                    ManipulationPrimitives.markExecutable(folderName + File.separator + magicEnergyScript);
                    final String[] formXYZ = cartes.createPrintableCartesians();
                    final String xyzFile = folderName + File.separator + "atompos.xyz";
                    OutputPrimitives.writeOut(xyzFile, formXYZ, false);
                    final short[] allAtomNos = cartes.getAllAtomNumbers();
                    for(int x = 0; x < allAtomNos.length; x++){
                        if(!nos.contains(allAtomNos[x])){
                            nos.add(allAtomNos[x]);
                        }
                    }
                
                    for(final short no : nos){
                        if(nonOptAtoms.contains(no)){
                            final String atom = AtomicProperties.giveAtomSymbol(no);
                            final String ppFile = folderName + File.separator + "pseudopot_" + atom + ".recpot";
                            final String oldFile = "pseudopot_" + atom + ".recpot";
                            OutputPrimitives.createLink(oldFile, ppFile);
                        } else {
                            final String ppFile = folderName + File.separator + "pseudopot_" + AtomicProperties.giveAtomSymbol(no) + ".recpot";
                            printModifiedPPOut(p,ppFile,no, pseudized.get(no), true);
                        }
                    }
                } catch(Exception e){
                    System.err.println("ERROR: Could not create folder " + folderName + " or copy script " + magicEnergyScript + " into it.");
                    e.printStackTrace(System.err);
                    final Tuple<String,Double> et = new Tuple<>(key,FixedValues.NONCONVERGEDENERGY);
                    energies.add(et);
                    continue;
                }
            
                // execute the script
                Process proc;
                try {
                    final Runtime rt = Runtime.getRuntime();
                    final String[] cmd = new String[]{"./" + magicEnergyScript, key};
                    final File dir = new File(folderName);
                    proc = rt.exec(cmd, null, dir);

                    // any error message?
                    final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

                    // any output?
                    final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

                    // kick them off
                    errorGobbler.start();
                    outputGobbler.start();

                    // any error???
                    if (proc.waitFor() != 0) {
                        System.err.println("PROFESS has a non-zero return value (energy calculation). ");
                        if(DEBUG){
                            final String[] out = outputGobbler.getData();
                            final String[] err = errorGobbler.getData();
                            System.out.println("DEBUG: system output...");
                            for(final String s : out){
                                System.out.println(s);
                            }
                            System.out.println("DEBUG: system error...");
                            for(final String s : err){
                                System.out.println(s);
                            }
                        }
                        final Tuple<String,Double> et = new Tuple<>(key,FixedValues.NONCONVERGEDENERGY);
                        energies.add(et);
                        continue;
                    }

                } catch (Exception e) {
                    System.err.println("PROFESS has a problem (energy calculation). ");
                    e.printStackTrace(System.err);
                    try {
                        ManipulationPrimitives.remove(folderName);
                    } catch (Exception e1) {
                        System.err.println("WARNING: Failed to clean up PROFESS forces files. ");
                        e1.printStackTrace(System.err);
                    }

                    final Tuple<String,Double> et = new Tuple<>(key,FixedValues.NONCONVERGEDENERGY);
                    energies.add(et);
                    continue;
                }

                // gather the results
                try {
                    final String resFile = folderName + File.separator + RESULTFILE;
                    final String[] dat = InputPrimitives.readFileIn(resFile);
                    final double energy = Double.parseDouble(dat[0].trim());
                    
                    if(DEBUG){
                        System.out.println("DEBUG: Energy read: " + energy);
                    }
                    final Tuple<String,Double> et = new Tuple<>(key,energy);
                    energies.add(et);
                } catch (Exception e) {
                    System.err.println("Failure to read the energy for " + folderName);
                    e.printStackTrace(System.err);

                    final Tuple<String,Double> et = new Tuple<>(key,FixedValues.NONCONVERGEDENERGY);
                    energies.add(et);
                    continue;
                }
                
                // delete the folder
                try {
                    ManipulationPrimitives.remove(folderName);
                } catch (Exception e) {
                    System.err.println("WARNING: Couldn't remove " + folderName + ". Continuing...");
                }
            }
            
            final EnergyOrder order = new EnergyOrder(energies, this.considerRelEnergies);
            
            if(DEBUG){
                System.out.println("DEBUG: returning energy order: ");
                energies.forEach((et) -> {
                    System.out.println("DEBUG: " + et.getObject1() + " with " + et.getObject2());
                });
            }
            
            return order;
        }

        @Override
        public EnergyOrder calculatePropertyGradient(final AdaptiveParameters p, final ReferenceEnergyOrderData<CartesianCoordinates> data,
                final double[] grad) {
            return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
        }
    }
    
    class DeltaGaugeProfessCalculator implements PropertyCalculator<DeltaGauge,ReferenceDeltaGaugeData<CartesianCoordinates>> {

        private static final long serialVersionUID = (long) 20150801;
        
        public static final String RESULTFILE = "Delta.out";
        public static final String DEFAULTDELTAGAUGESCRIPT = "profess_deltagauge_calc.sh";
        
        private final String magicDeltaGaugeScript;
        
        DeltaGaugeProfessCalculator(){
            
            // default: profess_bulk_calc.sh
            final String tmp = System.getenv("OGO_PROFESSDELTAGAUGESCRIPT");
            if(tmp == null){
                magicDeltaGaugeScript = DEFAULTDELTAGAUGESCRIPT;
            } else {
                magicDeltaGaugeScript = tmp;
            }
        }
        
        private DeltaGaugeProfessCalculator(final DeltaGaugeProfessCalculator orig){
            this.magicDeltaGaugeScript = orig.magicDeltaGaugeScript;
        }
        
        @Override
        public DeltaGaugeProfessCalculator clone() {
            return new DeltaGaugeProfessCalculator(this);
        }

        @Override
        public DeltaGauge calculateProperty(final AdaptiveParameters p, final ReferenceDeltaGaugeData<CartesianCoordinates> data) {
            
            // the idea here is simple: we create a directory, copy a script there (with a magic name, sorry)
            // and execute it with the crystal structure name as the argument. We assume the script
            // to return the delta gauge (in the same units as the reference) into a file called Delta.out
            
            final long paramID = p.getID();
            final String folderName = "professdeltagaugecalc_" + paramID + "_" + System.currentTimeMillis();
            
            // create folder, copy script and place pseudopotential
            try{
                OutputPrimitives.createAFolder(folderName);
                OutputPrimitives.copyFile(magicDeltaGaugeScript, folderName + File.separator + magicDeltaGaugeScript);
                ManipulationPrimitives.markExecutable(folderName + File.separator + magicDeltaGaugeScript);
                final String ppFile = folderName + File.separator + "pseudopot.recpot";
                printModifiedPPOut(p,ppFile,data.getAtomNo(), pseudized.get(data.getAtomNo()), true);
            } catch(Exception e){
                System.err.println("ERROR: Could not create folder " + folderName + " or copy script " + magicDeltaGaugeScript + " into it.");
                e.printStackTrace(System.err);
                return new DeltaGauge(1000.0);
            }
            
            // execute the script
            Process proc;
            try{
                final Runtime rt = Runtime.getRuntime();
                final String[] cmd = new String[]{"./" + magicDeltaGaugeScript, data.getSymmetry()};
                final File dir = new File(folderName);
                proc = rt.exec(cmd,null,dir);

                // any error message?
                final StreamGobbler errorGobbler = new
                    StreamGobbler(proc.getErrorStream(), "ERROR");

                // any output?
                final StreamGobbler outputGobbler = new
                    StreamGobbler(proc.getInputStream(), "OUTPUT");

                // kick them off
                errorGobbler.start();
                outputGobbler.start();

                // any error???
                if(proc.waitFor() != 0){
                    return new DeltaGauge(100.0);
                }

            } catch (Exception e) {
                System.err.println("PROFESS has a problem (delta gauge calculation). ");
                e.printStackTrace(System.err);
                try{
                    ManipulationPrimitives.remove(folderName);
                } catch(Exception e1){
                    System.err.println("WARNING: Failed to clean up PROFESS delta gauge files. ");
                    e1.printStackTrace(System.err);
                }
                
                return new DeltaGauge(1000.0);
            }
            
            // gather the results
            DeltaGauge mod;
            try{
                final String resFile = folderName + File.separator + RESULTFILE;
                final String[] dat = InputPrimitives.readFileIn(resFile);
                final String resLine = dat[0].trim();
                final double rawVal = Double.parseDouble(resLine);
                mod = new DeltaGauge(rawVal);
            } catch(Exception e){
                System.err.println("Failure to get the calculated delta gauge for " + folderName);
                e.printStackTrace(System.err);
                
                return new DeltaGauge(1000.0);
            }
            
            // delete the folder
            try{
                ManipulationPrimitives.remove(folderName);
            } catch(Exception e){
                System.err.println("WARNING: Couldn't remove " + folderName + ". Continuing...");
            }
            
            return mod;
        }

        @Override
        public DeltaGauge calculatePropertyGradient(final AdaptiveParameters p, final ReferenceDeltaGaugeData<CartesianCoordinates> data,
                final double[] grad) {
            return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
        }
    }
    
    class DensityProfessCalculator implements PropertyCalculator<Density,ReferenceDensityData<CartesianCoordinates>>{

        private static final long serialVersionUID = (long) 20151202;
        
        public static final String RESULTFILE = "density.out";
        public static final String DEFAULTDENSITYSCRIPT = "profess_density_calc.sh";
        
        private final String magicDensityScript;
        
        DensityProfessCalculator(){
            // default: profess_density_calc.sh
            final String tmp = System.getenv("OGO_PROFESSDENSITYSCRIPT");
            if(tmp == null){
                magicDensityScript = DEFAULTDENSITYSCRIPT;
            } else {
                magicDensityScript = tmp;
            }
        }
        
        private DensityProfessCalculator(final DensityProfessCalculator orig){
            this.magicDensityScript = orig.magicDensityScript;
        }
        
        @Override
        public DensityProfessCalculator clone() {
            return new DensityProfessCalculator(this);
        }

        @Override
        public Density calculateProperty(final AdaptiveParameters p, final ReferenceDensityData<CartesianCoordinates> data) {
            
            final long paramID = p.getID();
            final String folderName = "professdensitycalc_" + paramID + "_" + System.currentTimeMillis();
            
            // create folder, copy script and place pseudopotential as well as the set of cartesian coordinates
            final int noAtoms = data.getGeomData().c.getNoOfAtoms();
            try{
                final List<Short> nos = new ArrayList<>();
                OutputPrimitives.createAFolder(folderName);
                OutputPrimitives.copyFile(magicDensityScript, folderName + File.separator + magicDensityScript);
                ManipulationPrimitives.markExecutable(folderName + File.separator + magicDensityScript);
                final String[] formXYZ = data.getGeomData().c.getCartesianCoordinates().createPrintableCartesians();
                final String xyzFile = folderName + File.separator + "atompos.xyz";
                OutputPrimitives.writeOut(xyzFile, formXYZ, false);
                final short[] allAtomNos = data.getGeomData().c.getAllAtomNumbers();
                for(int x = 0; x < allAtomNos.length; x++){
                    if(!nos.contains(allAtomNos[x])){
                        nos.add(allAtomNos[x]);
                    }
                }
                
                for(final short no : nos){
                    if(nonOptAtoms.contains(no)){
                        final String atom = AtomicProperties.giveAtomSymbol(no);
                        final String ppFile = folderName + File.separator + "pseudopot_" + atom + ".recpot";
                        final String oldFile = "pseudopot_" + atom + ".recpot";
                        OutputPrimitives.createLink(oldFile, ppFile);
                    } else {
                        final String ppFile = folderName + File.separator + "pseudopot_" + AtomicProperties.giveAtomSymbol(no) + ".recpot";
                        printModifiedPPOut(p,ppFile,no, pseudized.get(no), true);
                    }
                }
            } catch(Exception e){
                System.err.println("ERROR: Could not create folder " + folderName + " or copy script " + magicDensityScript + " into it.");
                e.printStackTrace(System.err);
                return new Density(null);
            }
            
            // execute the script
            Process proc;
            try{
                final Runtime rt = Runtime.getRuntime();
                final String[] cmd = new String[]{"./" + magicDensityScript, data.getTag()};
                final File dir = new File(folderName);
                proc = rt.exec(cmd,null,dir);

                // any error message?
                final StreamGobbler errorGobbler = new
                    StreamGobbler(proc.getErrorStream(), "ERROR");

                // any output?
                final StreamGobbler outputGobbler = new
                    StreamGobbler(proc.getInputStream(), "OUTPUT");

                // kick them off
                errorGobbler.start();
                outputGobbler.start();

                // any error???
                if(proc.waitFor() != 0){
                    return new Density(null);
                }

            } catch (Exception e) {
                System.err.println("PROFESS has a problem (density calculation). ");
                e.printStackTrace(System.err);
                try{
                    ManipulationPrimitives.remove(folderName);
                } catch(Exception e1){
                    System.err.println("WARNING: Failed to clean up PROFESS density files. ");
                    e1.printStackTrace(System.err);
                }
                
                return new Density(null);
            }
            
            // gather the results
            Density density;
            try{
                final String resFile = folderName + File.separator + RESULTFILE;
                final String[] dat = InputPrimitives.readFileIn(resFile);
                
                // this is simply an insane density format.
                final String[] line = dat[0].trim().split("\\s+");
                
                final int xDim = Integer.parseInt(line[1]);
                final int yDim = Integer.parseInt(line[3]);
                final int zDim = Integer.parseInt(line[5]);

                final double[][][] densityValues = new double[xDim][yDim][zDim];
                
                int count = 10; // only the tenth element onwards contains density                
                for (int x = 0; x < xDim; x++) {
                    for (int y = 0; y < yDim; y++) {
                        for (int z = 0; z < zDim; z++) {
                            densityValues[x][y][z] = Double.parseDouble(line[count]);
                            count++;
                        }
                    }
                }
                
                if(false){
                    System.out.println("DEBUG: For  " + data.getTag() + " read in dimensions are " + xDim + " " + yDim + " " + zDim);
                }
                
                density = new Density(densityValues);
            } catch(Exception e){
                System.err.println("Failure to read the forces for " + folderName);
                e.printStackTrace(System.err);
                
                return new Density(null);
            }
            
            // delete the folder
            try{
                ManipulationPrimitives.remove(folderName);
            } catch(Exception e){
                System.err.println("WARNING: Couldn't remove " + folderName + ". Continuing...");
            }
            
            return density;
        }

        @Override
        public Density calculatePropertyGradient(final AdaptiveParameters p, final ReferenceDensityData<CartesianCoordinates> data,
                final double[] grad) {
            return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
        }
    }
    
    class StressTensorProfessCalculator implements PropertyCalculator<StressTensor,ReferenceStressTensorData<CartesianCoordinates>>{

        private static final long serialVersionUID = (long) 20180123;
        
        public static final String RESULTFILE = "stresstensor.out";
        public static final String DEFAULTDENSITYSCRIPT = "profess_stresstensor_calc.sh";
        
        private final String magicStressTensorScript;
        
        StressTensorProfessCalculator(){
            // default: profess_density_calc.sh
            final String tmp = System.getenv("OGO_PROFESSDENSITYSCRIPT");
            if(tmp == null){
                magicStressTensorScript = DEFAULTDENSITYSCRIPT;
            } else {
                magicStressTensorScript = tmp;
            }
        }
        
        private StressTensorProfessCalculator(final StressTensorProfessCalculator orig){
            this.magicStressTensorScript = orig.magicStressTensorScript;
        }
        
        @Override
        public StressTensorProfessCalculator clone() {
            return new StressTensorProfessCalculator(this);
        }

        @Override
        public StressTensor calculateProperty(final AdaptiveParameters p, final ReferenceStressTensorData<CartesianCoordinates> data) {
            
            final long paramID = p.getID();
            final String folderName = "professstresstensorcalc_" + paramID + "_" + System.currentTimeMillis();
            
            // create folder, copy script and place pseudopotential as well as the set of cartesian coordinates
            final int noAtoms = data.getGeomData().c.getNoOfAtoms();
            try{
                final List<Short> nos = new ArrayList<>();
                OutputPrimitives.createAFolder(folderName);
                OutputPrimitives.copyFile(magicStressTensorScript, folderName + File.separator + magicStressTensorScript);
                ManipulationPrimitives.markExecutable(folderName + File.separator + magicStressTensorScript);
                final String[] formXYZ = data.getGeomData().c.getCartesianCoordinates().createPrintableCartesians();
                final String xyzFile = folderName + File.separator + "atompos.xyz";
                OutputPrimitives.writeOut(xyzFile, formXYZ, false);
                final short[] allAtomNos = data.getGeomData().c.getAllAtomNumbers();
                for(int x = 0; x < allAtomNos.length; x++){
                    if(!nos.contains(allAtomNos[x])){
                        nos.add(allAtomNos[x]);
                    }
                }
                
                for(final short no : nos){
                    if(nonOptAtoms.contains(no)){
                        final String atom = AtomicProperties.giveAtomSymbol(no);
                        final String ppFile = folderName + File.separator + "pseudopot_" + atom + ".recpot";
                        final String oldFile = "pseudopot_" + atom + ".recpot";
                        OutputPrimitives.createLink(oldFile, ppFile);
                    } else {
                        final String ppFile = folderName + File.separator + "pseudopot_" + AtomicProperties.giveAtomSymbol(no) + ".recpot";
                        printModifiedPPOut(p,ppFile,no, pseudized.get(no), true);
                    }
                }
            } catch(Exception e){
                System.err.println("ERROR: Could not create folder " + folderName + " or copy script " + magicStressTensorScript + " into it.");
                e.printStackTrace(System.err);
                return new StressTensor(null);
            }
            
            // execute the script
            Process proc;
            try{
                final Runtime rt = Runtime.getRuntime();
                final String[] cmd = new String[]{"./" + magicStressTensorScript};
                final File dir = new File(folderName);
                proc = rt.exec(cmd,null,dir);

                // any error message?
                final StreamGobbler errorGobbler = new
                    StreamGobbler(proc.getErrorStream(), "ERROR");

                // any output?
                final StreamGobbler outputGobbler = new
                    StreamGobbler(proc.getInputStream(), "OUTPUT");

                // kick them off
                errorGobbler.start();
                outputGobbler.start();

                // any error???
                if(proc.waitFor() != 0){
                    return new StressTensor(null);
                }

            } catch (Exception e) {
                System.err.println("PROFESS has a problem (stress tensor calculation). ");
                e.printStackTrace(System.err);
                try{
                    ManipulationPrimitives.remove(folderName);
                } catch(Exception e1){
                    System.err.println("WARNING: Failed to clean up PROFESS stress tensor files. ");
                    e1.printStackTrace(System.err);
                }
                
                return new StressTensor(null);
            }
            
            // gather the results
            StressTensor stress;
            try{
                final String resFile = folderName + File.separator + RESULTFILE;
                final String[] dat = InputPrimitives.readFileIn(resFile);
                
                final double[][] stressVals = new double[3][3];
                for(int i = 0; i < 3; i++){
                    final String[] line = dat[i].trim().split("\\s+");
                    for(int j = 0; j < 3; j++){
                        stressVals[i][j] = Double.parseDouble(line[j]);
                    }
                }
                
                stress = new StressTensor(stressVals);
            } catch(Exception e){
                System.err.println("Failure to read the stress tensor for " + folderName);
                e.printStackTrace(System.err);
                
                return new StressTensor(null);
            }
            
            // delete the folder
            try{
                ManipulationPrimitives.remove(folderName);
            } catch(Exception e){
                System.err.println("WARNING: Couldn't remove " + folderName + ". Continuing...");
            }
            
            return stress;
        }

        @Override
        public StressTensor calculatePropertyGradient(final AdaptiveParameters p, final ReferenceStressTensorData<CartesianCoordinates> data,
                final double[] grad) {
            return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
        }
    }
}
