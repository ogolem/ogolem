/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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
package org.ogolem.switches;

import java.io.Serializable;
import org.ogolem.generic.genericpool.*;

/**
 * The central place for the runtime configuration of the program.
 * @author Johannes Dieterich
 * @version 2014-07-18
 */
public final class SwitchesConfig implements Serializable{

    private static final long serialVersionUID = (long) 20140718;

    /**
     * gets read in and constructed by the input subroutine.
     * objects includes connectivity and bonding information
     */
    public static Backbone backbone = null;

    /**
     * The size of our population.
     * PoolSize=
     */
    public static int iPoolSize = 100;

    /**
     * how many global iterations we will do.
     * NoOfGlobIters=
     */
    public static int iNoOfGlobIters = 900;

    /**
     * should we even check the diversity or trust taboo search
     * AnyDiversityCheck=
     */
    static boolean bAnyDivCheck = true;

    /**
     * the fitness diversity that needs to be maintained
     * FitnessDiversity=
     */
    static double dFitnessDiversity = 1E-6;

    /**
     * after how many steps the pool gets automatically serialized
     * StepsToSerialization=
     */
    static int iSwitchesToSerial = 100;
    
    /**
     * TODO keyword and doc
     */
    static String switchDumpFolder = "allswitches";

    /**
     * Which global optimization to be used. No input option yet since we
     * only provide Germany ATM.
     */
    int iWhichGlobAlgo = 0;

    /**
     * Which local optimization to be used.
     * LocOptAlgo=
     */
    int iWhichLocAlgo = 100;

    /**
     * how many iterations the locopt is allowed to do
     * MaxIterLocOpt=
     */
    int iMaxIterLocOpt = 1000;

    /**
     * Which excitor to be used.
     * ExcitorMethod=
     */
    int iWhichExcitor = 0;

    /**
     * Our blowing factor for detection of bonds.
     * BlowBondDetect=
     */
    public double dBlowBondsFac = 1.2;

    /**
     * more mutation or less
     * MoreMutation=
     */
    boolean bMoreMutation = true;

    /**
     * Which atoms to be used for the dihedral check
     * AtomsForDihedralCis=
     * AtomsForDihedralTrans=
     */
    int[] iaWhichAtsForDihedralCis = null;
    int[] iaWhichAtsForDihedralTrans = null;

    /*
     * some factors for the fitness function
     */
    //MinimumEnergyDiff=
    double dMinEnergyDiff = 0.4;
    //ExactEnergyCis=
    double dExactEnergyCis = 1.0;
    //ExactEnergyTrans=
    double dExactEnergyTrans = 1.4;
    //ExactEnergyFactor=
    double dExactEnergyFactor = 10.0;
    //SigmoidFactorA=
    double dSigFactorA = 1000.0;
    //SigmoidFactorB=
    double dSigFactorB = 10.0;
    //DihedralThreshCis=
    double dThreshDihedralCis = 0.05;
    //DihedralThreshTrans=
    double dThreshDihedralTrans = 0.05;
    //OptimalCisDihedral=
    double dOptimalCisVal = 0.0;
    //OptimalTransDihedral
    double dOptimalTransVal = 0.5 * Math.PI;

    /*
     * S0->S1 related stuff
     */

    // NumberOfActiveOcc=
    int iActiveOcc = 7;
    // NumberOfNonActiveOcc=
    int iActiveNonOcc = 6;
    // NumberOfOccPi=
    int iOccPi = 3;
    // NumberOfNonOccPi=
    int iNonOccPi = 0;
    // GUGADefOfOrbs=
    int iOrbDef = -3;
    // NumberOfRefOcc=
    int iNoOfRefOcc = 1;
    // DefOfRefOcc=
    int iDefOfRefOcc = 0;
    // MaxExcitationLevel=
    int iMaxExcitLevel = 4;
    // NumberOfLowestStates=
    int iNoOfLowestStates = 4;
    // WantedCIState=
    int iWantedCiState = 4;

    /**
     * Debugging on/off.
     * PhonyDebug=
     */
    boolean bDebug = false;
    
    /**
     * TODO keyword and doc
     */
    String whichParentsChoice = "fitnessrankbased:gausswidth=0.05";

    
    static int parseIntLocOptFromString(final String sLocOpt){
        if(sLocOpt.startsWith("tinker:")){
            final String sTemp = sLocOpt.substring(7).trim();
            if(sTemp.equalsIgnoreCase("minimize")){
                return 0;
            } else if(sTemp.equalsIgnoreCase("newton")){
                return 1;
            } else if(sTemp.equalsIgnoreCase("optimize")){
                return 2;
            } else {
                System.err.println("WARNING: Not aware of locopt method " + sLocOpt + " using tinker:minimize now.");
                return 0;
            }
        } else if(sLocOpt.startsWith("openbabel:")){
            final String sTemp3 = sLocOpt.substring(10).trim();
            if (sTemp3.equalsIgnoreCase("ghemical")) {
                return 100;
            } else if (sTemp3.equalsIgnoreCase("mmff94")) {
                return 101;
            } else if (sTemp3.equalsIgnoreCase("mmff94s")) {
                return 102;
            } else if (sTemp3.equalsIgnoreCase("uff")) {
                return 103;
            } else {
                System.err.println("WARNING: Wrong input to configure OpenBabel: " + sTemp3 + " using the GHEMICAL force field.");
                return 100;
            }
        } else if (sLocOpt.startsWith("mopac:")) {
            String sTemp3 = sLocOpt.substring(6).trim();
            if (sTemp3.equalsIgnoreCase("mndo")) {
                return 200;
            } else if (sTemp3.equalsIgnoreCase("am1")) {
                return 201;
            } else if (sTemp3.equalsIgnoreCase("pm3")) {
                return 202;
            } else if (sTemp3.equalsIgnoreCase("pm5")) {
                return 203;
            } else if (sTemp3.equalsIgnoreCase("pm6")) {
                return 204;
            } else {
                System.err.println("Wrong input to configure MOPAC: " + sTemp3 + " using pm3.");
                return 202;
            }
        } else {
            System.err.println("WARNING: Not aware of locopt method " + sLocOpt + " using openbabel:ghemical now.");
            return 100;
        }
    }

    static int parseIntExcitorFromString(final String sExcitor){
        if(sExcitor.startsWith("guga:")){
            final String sTemp = sExcitor.substring(5).trim();
            if(sTemp.equalsIgnoreCase("mndo/d")){
                return 0;
            } else if(sTemp.equalsIgnoreCase("om3")){
                return 1;
            } else if(sTemp.equalsIgnoreCase("pm3")){
                return 2;
            } else if(sTemp.equalsIgnoreCase("om2")){
                return 3;
            } else if(sTemp.equalsIgnoreCase("om1")){
                return 4;
            } else if(sTemp.equalsIgnoreCase("am1")){
                return 5;
            } else if(sTemp.equalsIgnoreCase("mndoc")){
                return 6;
            } else if(sTemp.equalsIgnoreCase("mndo")){
                return 7;
            } else if(sTemp.equalsIgnoreCase("mindo/3")){
                return 8;
            } else if(sTemp.equalsIgnoreCase("cndo/2")){
                return 9;
            } else if(sTemp.equalsIgnoreCase("SCC-DFTB")){
                return 10;
            } else if(sTemp.equalsIgnoreCase("SCC-DFTB-Jorgensen")){
                return 11;
            } else{
                System.err.println("WARNING: Not aware of excitor " + sExcitor + " in guga using guga:om2 now.");
                return 7;
            }
        } else if(sExcitor.startsWith("mopac:")){
            final String sTemp = sExcitor.substring(6).trim();
            if(sTemp.equalsIgnoreCase("am1,noperisico,noexternal")){
                return 100;
            } else if(sTemp.equalsIgnoreCase("pm3,nopersico,noexternal")){
                return 101;
            } else if(sTemp.equalsIgnoreCase("pm5,nopersico,noexternal")){
                return 102;
            } else if(sTemp.equalsIgnoreCase("pm6,nopersico,noexternal")){
                return 103;
            } else if(sTemp.equalsIgnoreCase("am1,perisco,external")){
                return 104;
            } else if(sTemp.equalsIgnoreCase("am1,persico,noexternal")){
                return 105;
            } else if(sTemp.equalsIgnoreCase("pm3,persico,noexternal")){
                return 106;
            } else if(sTemp.equalsIgnoreCase("pm5,persico,noexternal")){
                return 107;
            } else if(sTemp.equalsIgnoreCase("pm6,persico,noexternal")){
                return 108;
            } else{
                System.err.println("WARNING: Not aware of excitor " + sExcitor + " in mopac using mopac:pm3,nopersico,noexternal now.");
                return 101;
            }
        } else if(sExcitor.startsWith("orca:")){
            final String sTemp = sExcitor.substring(5).trim();
            if(sTemp.equalsIgnoreCase("bp,def2-tzvp")){
                return 150;
            } else if(sTemp.equalsIgnoreCase("blyp,def2-tzvp")){
                return 151;
            } else {
                System.err.println("WARNING: Not aware of excitor " + sExcitor + " in orca using orca:bp,def2-TZVP now.");
                return 150;
            }
        } else{
            System.err.println("WARNING: Not aware of excitor " + sExcitor + " using guga:om2 now.");
            return 7;
        }
    }
    //TODO consider whether charges (more no) and spins (perhaps) might be in
    // order... but I guess one can add that later.
    
    public GenericPoolConfig<Color,Switch> getGenericConfig(){
        
        final GenericPoolConfig<Color,Switch> config = new GenericPoolConfig<>();
        
        config.setDoNiching(false); //XXX
        config.setSerializeAfterNewBest(true);
        config.setAcceptableFitness(0.0);
        config.setAddsToSerial(SwitchesConfig.iSwitchesToSerial);
        config.setWriteEveryAdd(false);
        config.setPoolSize(SwitchesConfig.iPoolSize);
        config.setInterBinFile("IntermediateSwitchPool.bin");
        
        DiversityChecker<Color,Switch> diver;
        switch(0){
            case 0: diver = new GenericDiversityCheckers.FitnessDiversityChecker<>(SwitchesConfig.dFitnessDiversity); break;
                //TODO more!!!
            default: diver = new GenericDiversityCheckers.FitnessDiversityChecker<>(SwitchesConfig.dFitnessDiversity); break;
        }
        config.setDiversityChecker(diver);
        
        if(false){
            //config.setAddsToStats(this.addsToNicheStats);
            //switch(this.whichNicher){
            //    case 0: config.setNicheComp(new HydrogenBondNicheComp(2.5*Constants.ANGTOBOHR, 3.5*Constants.ANGTOBOHR, 45.0)); break;
            //    default: config.setNicheComp(new HydrogenBondNicheComp(2.5*Constants.ANGTOBOHR, 3.5*Constants.ANGTOBOHR, 45.0)); break;
            //}
            //config.setNicher(new SimpleNicher<Switch>(this.noOfIndividualsPerNicheAtMax));
        }
        
        ParentSelector<Color,Switch> selec;
        try{
            selec = GenericParentSelectors.buildSelector(whichParentsChoice);
        } catch(Exception e){
            throw new RuntimeException("Error in creating parent selector.",e);
        }
        config.setSelector(selec);
        
        config.setWriter(new SwitchWriter());
        config.setStats(new GenericStatistics("switchprogress.log",10000)); // XXX hard coded

        return config;
    }
}
