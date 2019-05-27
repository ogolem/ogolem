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
package org.ogolem.parameters;

import java.util.List;
import java.util.Locale;
import org.ogolem.adaptive.*;
import org.ogolem.adaptive.genericfitness.ReferenceGeomData;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.InitIOException;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.helpers.Fortune;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.properties.Energy;

/**
 * Post-run analyzation of parameter optimizations.
 * @author Johannes Dieterich
 * @version 2016-07-20
 */
public class MainParameters {

    /**
     * Starting point. Valid arguments are:
     * - -dongetparams
     * - -pool=
     * - -makegaussplot
     * - -calcFitness
     * - -ogoinput=
     * @param args
     */
    public static void run(String[] args) {

        boolean bGetParamsOut = true;
        boolean bMakeGaussPlot = false;
        boolean calcFitness = false;
        String sParamPoolPath = "paramsPool.bin";
        String sParamFolder = "params";
        String ogoInput = "inp.ogo";
        
        for(final String sCurrArg : args){
            if(sCurrArg.equalsIgnoreCase("-dontgetparams")){
                bGetParamsOut = false;
            } else if(sCurrArg.startsWith("-pool=")){
                sParamPoolPath = sCurrArg.substring(6);
            } else if(sCurrArg.equalsIgnoreCase("-makegaussplot")){
                bMakeGaussPlot = true;
            } else if(sCurrArg.equalsIgnoreCase("-calcFitness")){
                calcFitness = true;
            } else if(sCurrArg.startsWith("-ogoinp=")){
                ogoInput = sCurrArg.substring(8);
            } else{
                System.err.println("ERROR: I do not understand the argument "
                        + sCurrArg + ". Ignoring it.");
            }
        }

        // FIRST THING: READ THE POOL IN
        GenericPool<Double,AdaptiveParameters> pool = null;
        try{
            pool = Input.readParamPoolIn(sParamPoolPath);
        } catch(Exception e){
            System.err.println("ERROR: Couldn't read in parameter pool. Aborting. "
                    + e.toString());
            System.exit(123);
        }

        // first create a folder for the parameter sets
        try {
            OutputPrimitives.createAFolder(sParamFolder);
        } catch (Exception e) {
            System.err.println("ERROR: Couldn't create folder for parameters. "
                    + "Aborting. " + e.toString());
            System.exit(145);
        }

        if(bGetParamsOut){

            // then write the parameter sets out
            System.out.println("INFO: Pool contains " + pool.getPoolSize() + " individuals.");
            for(int i = 0; i < pool.getPoolSize(); i++){
                final AdaptiveParameters param = pool.getIndividualAtPosition(i);
                final String[] saParamCont = param.createFormattedOutput();
                final String sOutFile = sParamFolder + System.getProperty("file.separator") + "rank" + i + "param" + param.getID() + ".prm";
                try{
                    OutputPrimitives.writeOut(sOutFile, saParamCont, false);
                } catch(Exception e){
                    System.err.println("ERROR: Failure to write paramter result " + i + " out. Continuing.");
                    e.printStackTrace(System.err);
                }
            }
       }

        // JUST FOR GAUSSIAN BASED FORCE FIELDS: PLOT EVERY INDIVIDUAL GAUSSIAN
        // ADDITIONALLY: this produces gnuplot .gnu files which effectively restricts it
        //               to UNIX. That's life... ;-)
        if(bMakeGaussPlot){
            // get the best parameter set
            AdaptiveParameters params = pool.getIndividualAtPosition(0);

            // get all keys
            String[] saKeys = params.getAllKeysCopy();

            for(String key : saKeys){
                // output
                String sOutput = sParamFolder + System.getProperty("file.separator")
                        + key + "_plot.gnu";

                // get the actual params
                double[] daParameters = params.getParametersForKey(key);

                String sTotalGauss = "plot ";
                String sOverall = "";
                for(int i = 0; i < daParameters.length; i=i+3){
                    String sThisGauss = ""+daParameters[i]+ "*exp(-"
                            + daParameters[i+1] + "*(x-" + daParameters[i+2]
                            + ")**2)";
                    sTotalGauss += " + " + sThisGauss;
                    sThisGauss += "lw 2 notitle";
                    
                    if(i < daParameters.length-3){
                        // it's not the last one
                        sThisGauss += ",\\";
                    }
                    sOverall += "\n" + sThisGauss;
                }

                sTotalGauss += "lw 2 notitle,\\" + sOverall + "\n ";

                String [] saPrint = {
                    "set terminal post eps enhanced color 'Helvetica,16'",
                    "set size 0.6",
                    "set output '"+ key +".eps'",
                    "set label font 'Helvetica,16'",
                    "set xtics nomirror border out font 'Helvetica,16'",
                    "set ytics nomirror font 'Helvetica,16'",
                    "set xlabel font 'Helvetica,20'",
                    "set ylabel font 'Helvetica,20'",
                    // better too much than too little.. ;-)
                    "set samples 5000",
                    sTotalGauss};
                try{
                    OutputPrimitives.writeOut(sOutput, saPrint, true);
                }catch(Exception e){
                    e.printStackTrace(System.err);
                }
            }
        }
        
        if(calcFitness){
            
            // read input in
            GlobalConfig globConf = null;
            try {
                globConf = org.ogolem.core.Input.ConfigureMe(ogoInput);
            } catch (InitIOException e) {
                System.err.println("ERROR: Failure in reading the global configuration!");
                e.printStackTrace(System.err);
                System.exit(1111);
            } catch (CastException e) {
                System.err.println("ERROR: Failure in casting the global configuration!");
                e.printStackTrace(System.err);
                System.exit(1111);
            } catch (Exception e) {
                // catch also the whole rest of exceptions, out of safety
                System.err.println("ERROR: Failure in creating the global configuration!");
                e.printStackTrace(System.err);
                System.exit(1111);
            }
    
            System.out.println("INFO: Raw energies");
            System.out.println("#No\t Ref energy\t Eq ref \t actual energy");
            final AdaptiveConf adaptiveConf = globConf.getAdaptiveConf();
            if(adaptiveConf.structuralDataType() != AdaptiveConf.StructureDataType.Cartesian){
                System.err.println("ERROR: Adaptive configuration does not use a cartesian coordinates structural data type. Implementation incomplete. Please contact author(s).");
                System.exit(42);
            }
            final List<Energy> referenceEnergies = adaptiveConf.getReferenceEnergies();
            final List<ReferenceGeomData<Energy,CartesianCoordinates>> referenceData = adaptiveConf.getReferenceGeomDataForEnergy();
            final AdaptiveGateway adaptivable = new AdaptiveGateway(globConf.getAdaptiveConf());
            final double first = referenceEnergies.get(0).getValue();
            for (int i = 0; i < referenceData.size(); i++) {
                final ReferenceGeomData<Energy,CartesianCoordinates> db = referenceData.get(i);
                final CartesianCoordinates cartes = db.c;
                final double energy = adaptivable.energyOfStructWithParams(cartes, pool.getIndividualAtPosition(0), i, db.bonds);
                final double ref = referenceEnergies.get(i).getValue();
                System.out.println(" " + i + "    " + String.format(Locale.US, "%20.9f", ref) + "    " + String.format(Locale.US, "%20.9f",(ref - first)) + "    " + String.format(Locale.US, "%20.9f",energy));
            }
        }

        System.out.println("Thank you for running OGOLEM. Hope it was satisfying. :-)");
        System.out.println(Fortune.randomFortune());
    }
}
