/**
Copyright (c) 2015-2020, J. M. Dieterich and B. Hartke
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
import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.generic.GenericBackend;
import org.ogolem.helpers.Tuple;
import org.ogolem.qmdff.AdditionalTerm;
import org.ogolem.qmdff.EVB;
import org.ogolem.qmdff.EVBQMDForceField;
import org.ogolem.qmdff.ElectrostaticEmbedding;
import org.ogolem.qmdff.QMDFFData;
import org.ogolem.qmdff.QMDForceField;
import org.ogolem.qmdff.QMDForceFieldParameters;
import org.ogolem.qmdff.SimpleEVB2;

/**
 * A backend constructing factory.
 * @author Johannes Dieterich
 * @version 2020-07-19
 */
public class BackendFactory {
    
    private final AdaptiveParameters params;
    private final GlobalConfig config;
    
    public BackendFactory(final GlobalConfig config, final AdaptiveParameters params){
        this.config = config;
        this.params = params;
    }
    
    public GenericBackend<Molecule, Geometry> parseBackend(final String backend) throws Exception {
        
        // we either have an xyz or rigid body backend
        if(backend.startsWith("xyz:")){
            final CartesianFullBackend cBack = mapStringToBackend(config, backend.substring(4),params,config.outputFolder);
            final FullyCartesianCoordinates cCoords = new FullyCartesianCoordinates(cBack);

            return cCoords;
        } else if(backend.startsWith("environment+xyz:")){
            final CartesianFullBackend cBack = mapStringToBackend(config, backend.substring(16),params,config.outputFolder);
            final EnvironmentCartesCoordinates ecCoords = new EnvironmentCartesCoordinates(cBack);
            
            return ecCoords;
        } else if(backend.startsWith("rigid:")){
            try{
                final RigidBodyBackend rBack = mapStringToRigidBackend(backend.substring(6),params);
                final RigidBodyCoordinates rCoords = new RigidBodyCoordinates(rBack);
            
                return rCoords;
            } catch (Exception e){
                // no need to act on this, as we try compatibility next
                System.out.println("INFO: Parsing backend as fully rigid failed. Trying parsing as Cartesian2Rigid next.");
            }
            
            // XXX think about if own key (like xyz2rigid: or something) is better?
            
            final CartesianFullBackend cBack = mapStringToBackend(config, backend.substring(6),params,config.outputFolder);
            final CartesianToRigidCoordinates rCoords = new CartesianToRigidCoordinates(cBack);
            
            return rCoords;
        }
        
        throw new RuntimeException("No backend for input string " + backend + " found.");
    }
    
    private static RigidBodyBackend mapStringToRigidBackend(final String backend, final AdaptiveParameters params) throws Exception {
        
        if(backend.equalsIgnoreCase("tip3p")){
            return new TIP3PForceField();
        } else if(backend.equalsIgnoreCase("tip4p")){
            return new TIP4PForceField();
        }
        
        throw new RuntimeException("No rigid body backend for key " + backend);
    }

    public static CartesianFullBackend mapStringToBackend(final GlobalConfig config, final String sBackend, final AdaptiveParameters params, final String outFolder) throws Exception {
        
        CartesianFullBackend back;
        
        if (sBackend.equalsIgnoreCase("lennardjones")) {
            back = new LennardJonesFF(true);
        } else if(sBackend.equalsIgnoreCase("lennardjones,nocache")){
            back = new LennardJonesFF(false);
        } else if(sBackend.equalsIgnoreCase("mixedlj")){
            back = new MixedLJForceField(true);
        } else if(sBackend.equalsIgnoreCase("mixedlj,nocache")){
            back = new MixedLJForceField(false);
        } else if(sBackend.startsWith("native:")){
            final String[] sub = sBackend.substring(7).trim().split("\\,");
            System.out.println("DEBUG " + sBackend.substring(7).trim() + "    " + sBackend);
            final int which = Integer.parseInt(sub[0].trim());
            final int howMany = Integer.parseInt(sub[1].trim());
            final double cutE = Double.parseDouble(sub[2].trim());
            back = new FFEngineWrapper(which,howMany,cutE);
        } else if (sBackend.startsWith("molpro:")) {
            String sTemp3 = sBackend.substring(7).trim();
            if (sTemp3.equalsIgnoreCase("am1/vdz") || sTemp3.equalsIgnoreCase("am1")) {
                back = new MolproCaller(MolproCaller.METHOD.AM1);
            } else if (sTemp3.equalsIgnoreCase("hf/vdz")) {
                back = new MolproCaller(MolproCaller.METHOD.HFVDZ);
            } else if (sTemp3.equalsIgnoreCase("b86/vdz")) {
                back = new MolproCaller(MolproCaller.METHOD.B86VDZ);
            } else if (sTemp3.equalsIgnoreCase("mp2/avtz")) {
                back = new MolproCaller(MolproCaller.METHOD.MP2AVTZ);
            } else {
                System.err.println("Wrong input to configure MOLPRO: " + sTemp3 + " using b86/vdz.");
                back = new MolproCaller(MolproCaller.METHOD.B86VDZ);
            }
        } else if (sBackend.startsWith("universal")) {
            if(sBackend.startsWith("universal:")){
                final String s = sBackend.substring(10).trim();
                if(s.equalsIgnoreCase("babel")) {back = new UniversalFF(1,true);}
                else if(s.equalsIgnoreCase("babel,nocache")) {back = new UniversalFF(1,false);}
                else if(s.equalsIgnoreCase("nocache")) {back = new UniversalFF(0,false);}
                else {back = new UniversalFF(0,true);}
            } else {back = new UniversalFF(0,true);}
        } else if (sBackend.startsWith("untangle:")) {
            final String s = sBackend.substring(9).trim();
            double sigma = 1.0;
            boolean moleculeWise = true;
            if (s.length() != 0) {
                final String[] sub = s.split(",");
                for (String opt: sub) {
                    if (opt.startsWith("sigma=")) {
                        try {
                            sigma = Double.parseDouble(opt.substring(6));
                            if (sigma <= 0)
                                throw new RuntimeException("ERROR: Sigma value in untangler must be > 0!");
                        } catch (Exception e) {
                            throw new RuntimeException("ERROR: Could not parse untangler sigma value!", e);
                        }
                    } else if (opt.startsWith("moleculeWise=")){
                        try {
                            moleculeWise = Boolean.parseBoolean(opt.substring(13));
                        } catch (Exception e) {
                            throw new RuntimeException("ERROR: Could not parse untangler moleculeWise flag!", e);
                        }
                    } else {
                        throw new RuntimeException("ERROR: Unknown option in untangler: " + opt);
                    }
                }
            }
            back = new MoleculeUntangler(sigma, moleculeWise);
        } else if(sBackend.equalsIgnoreCase("adaptiveUFF")){
            back = new org.ogolem.adaptive.AdaptiveUFF(false);
        } else if(sBackend.startsWith("adaptiveSWGFF")){
            double blowFacClose = 0.2;
            if(sBackend.startsWith("adaptiveSWGFF:")){
                final String s = sBackend.substring("adaptiveSWGFF:".length()).trim();
                if(s.startsWith("blowfacclose=")){
                    blowFacClose = Double.parseDouble(s.substring("blowfacclose=".length()).trim());
                } else {
                    throw new RuntimeException("Illegal option " + s + " for adaptiveSWGFF.");
                }
            }
            back = new org.ogolem.adaptive.AdaptiveSWGFF(false, params, true, blowFacClose);
        } else if (sBackend.startsWith("adaptiveLJFF")) {
            if (sBackend.endsWith("easy;6,12,6")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, true, 6, 12, 6, params, false,1.2,20.0,false,false);
            } else if (sBackend.endsWith("easy;2,12,2")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, true, 2, 12, 2, params, false,1.2,20.0,false,false);
            } else if (sBackend.endsWith("easy;6,16,2")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, true, 6, 16, 2, params, false,1.2,20.0,false,false);
            } else if (sBackend.endsWith("nomix;6,12,6")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, false, 6, 12, 6, params, false,1.2,20.0,false,false);
            } else if (sBackend.endsWith("nomix;2,12,2")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, false, 2, 12, 2, params, false,1.2,20.0,false,false);
            } else if (sBackend.endsWith("nomix;6,16,2")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, false, 6, 16, 2, params, false,1.2,20.0,false,false);
            } else if (sBackend.endsWith("easy;6,12,6;3body")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, true, 6, 12, 6, params, true,1.2,20.0,false,false);
            } else if (sBackend.endsWith("easy;2,12,2;3body")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, true, 2, 12, 2, params, true,1.2,20.0,false,false);
            } else if (sBackend.endsWith("easy;6,16,2;3body")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, true, 6, 16, 2, params, true,1.2,20.0,false,false);
            } else if (sBackend.endsWith("nomix;6,12,6;3body")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, false, 6, 12, 6, params, true,1.2,20.0,false,false);
            } else if (sBackend.endsWith("nomix;2,12,2;3body")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, false, 2, 12, 2, params, true,1.2,20.0,false,false);
            } else if (sBackend.endsWith("nomix;6,16,2;3body")) {
                back = new org.ogolem.adaptive.AdaptiveLJFF(false, false, 6, 16, 2, params, true,1.2,20.0,false,false);
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
                    final String s2 = sBackend.substring(12).trim();
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
                    back = new org.ogolem.adaptive.AdaptiveLJFF(false,mix,start,end,incr,params,threeB,blowClose,blowDist,cache,disc2body);
                } catch(Exception e){
                    throw new Exception("ERROR: Couldn't parse adaptiveLJ Input.", e);
                }
            }
        } else if(sBackend.startsWith("adaptivegupta:")){
            final String s2 = sBackend.substring(14).trim();
            final String[] sa = s2.split("\\,");
            try{
                final double distBlow = Double.parseDouble(sa[0]);
                final boolean cache = Boolean.parseBoolean(sa[1]);
                int expF = 0;
                if(sa.length == 3) {expF = Integer.parseInt(sa[2]);}
                back = new org.ogolem.adaptive.AdaptiveGUPTA(false, params, distBlow, cache, expF);
            } catch(Exception e){
                System.err.println("ERROR: Couldn't parse required information for adaptivegupta. " + e.toString());
                return null;
            }
        } else if(sBackend.startsWith("adaptivemorse:")){
            back = new org.ogolem.adaptive.AdaptiveMorse(false,sBackend,params);
        } else if(sBackend.startsWith("amberff")){
            // XXX this parsing is crap!
            boolean useCache = true;
            int whichAcos = 0;
            int whichCos = 0;
            int whichSin = 0;
            boolean totalShift = false;
            if(sBackend.contains("cos=0")) {whichCos=0;}
            if(sBackend.contains("cos=1")) {whichCos=1;}
            if(sBackend.contains("cos=2")) {whichCos=2;}
            if(sBackend.contains("cos=3")) {whichCos=3;}
            if(sBackend.contains("cos=4")) {whichCos=4;}
            if(sBackend.contains("cos=5")) {whichCos=5;}
            if(sBackend.contains("cos=10")) {whichCos=10;}
            if(sBackend.contains("cos=11")) {whichCos=11;}
            if(sBackend.contains("cos=12")) {whichCos=12;}
            if(sBackend.contains("cos=13")) {whichCos=13;}
            if(sBackend.contains("cos=14")) {whichCos=14;}
            if(sBackend.contains("cos=15")) {whichCos=15;}
            if(sBackend.contains("acos=0")) {whichAcos=0;}
            if(sBackend.contains("acos=1")) {whichAcos=1;}
            if(sBackend.contains("acos=2")) {whichAcos=2;}
            if(sBackend.contains("acos=3")) {whichAcos=3;}
            if(sBackend.contains("acos=4")) {whichAcos=4;}
            if(sBackend.contains("acos=5")) {whichAcos=5;}
            if(sBackend.contains("acos=10")) {whichAcos=10;}
            if(sBackend.contains("acos=11")) {whichAcos=11;}
            if(sBackend.contains("acos=12")) {whichAcos=12;}
            if(sBackend.contains("acos=13")) {whichAcos=13;}
            if(sBackend.contains("acos=14")) {whichAcos=14;}
            if(sBackend.contains("acos=15")) {whichAcos=15;}
            if(sBackend.contains("fast")){
                whichCos = 4;
                whichAcos = 4;
                whichSin = 4;
            }
            if(sBackend.contains("totalshift")) {totalShift = true;}
            if(sBackend.endsWith("nocache")) {useCache = false;}
            try{
                back = new org.ogolem.adaptive.AdaptiveAmberFF(false, params, useCache, whichAcos, whichCos, whichSin, totalShift);
            } catch(Exception e){
                System.err.println("ERROR: Couldn't setup adaptive AMBER FF. Stack coming...");
                e.printStackTrace(System.err);
                return null;
            }
        } else if(sBackend.startsWith("skalevala")){
            back = new SkalevalaCaller();
        } else if(sBackend.startsWith("scattm3f:")){
            final String[] s2 = sBackend.substring(9).trim().split("\\,");
            final int no = Integer.parseInt(s2[0].trim());
            final double cut = Double.parseDouble(s2[1].trim());
            final boolean getDipols = (s2.length >= 3) ? Boolean.parseBoolean(s2[2].trim()) : false;
            final boolean addCutoff = (s2.length >= 4) ? Boolean.parseBoolean(s2[3].trim()) : false;
            back = new ScaTTM3FBackend(no,cut,getDipols, outFolder, addCutoff);
        } else if(sBackend.startsWith("qmdff:")){
            final String[] s2 = sBackend.substring(6).trim().split("\\,");
            
            String paramFile = "solvent";
            List<AdditionalTerm> addTerms = new ArrayList<>();
            for(final String token : s2){
                if(token.startsWith("params=")){
                    paramFile = token.substring(7).trim();
                } else if(token.startsWith("addterm=")){
                    final String subtok = token.substring(8).trim();
                    if(subtok.startsWith("electrostatic:")){
                        final String chargeFile = subtok.substring(14).trim();
                        final ElectrostaticEmbedding emb = new ElectrostaticEmbedding(chargeFile);
                        addTerms.add(emb);
                    } else {
                        throw new RuntimeException("Unknown additional term " + subtok + " for QMDFF.");
                    }
                } else {
                    throw new RuntimeException("Unknown token " + token + " for QMDFF calculation.");
                }
            }
            
            final Tuple<QMDForceFieldParameters,QMDFFData> tup = org.ogolem.qmdff.SetupQMDFF.setup(paramFile);
            
            if(addTerms.isEmpty()){
                back = new QMDForceField(paramFile,tup.getObject2(),tup.getObject1());
            } else {
                back = new QMDForceField(paramFile,tup.getObject2(),tup.getObject1(), addTerms);
            }
        } else if(sBackend.startsWith("evb-qmdff:")){
            
            final String[] s2 = sBackend.substring(10).trim().split("\\,");
            
            final List<Tuple<QMDForceFieldParameters,QMDForceField>> ffs = new ArrayList<>();
            EVB evb = null;
            for(final String token : s2){
                if(token.startsWith("evb=")){
                    final String evbString = token.substring(4).trim();
                    if(evbString.startsWith("simpleevb2:")){
                        final String[] sa = evbString.substring(11).trim().split("\\:");
                        
                        boolean useGauss = false;
                        double offA = 0.0;
                        double offB = 0.0;
                        double eShift1 = 0.0;
                        double eShift2 = 0.0;
                        for(final String subTok : sa){
                            if(subTok.startsWith("offa=")){
                                final String s = subTok.substring(5).trim();
                                offA = Double.parseDouble(s);
                            } else if(subTok.startsWith("offb=")){
                                final String s = subTok.substring(5).trim();
                                offB = Double.parseDouble(s);
                            } else if(subTok.startsWith("eshift2=")){
                                final String s = subTok.substring(8).trim();
                                eShift2 = Double.parseDouble(s);
                            } else if(subTok.equalsIgnoreCase("gaussmode")){
                                useGauss = true;
                            } else {
                                throw new RuntimeException("");
                            }
                        }
                        
                        if(useGauss){
                            evb = new SimpleEVB2(offA,offB,eShift1,eShift2);
                        } else {
                            evb = new SimpleEVB2(offA,eShift1,eShift2);
                        }
                    } else {
                        throw new RuntimeException("Unknown EVB " + evbString + " specified.");
                    }
                } else if(token.startsWith("qmdffs=")){
                    final String[] saX = token.substring(7).trim().split("\\:");
                    
                    for(final String subTok : saX){
                        final String[] sx = subTok.trim().split("\\,");
            
                        String paramFile = "solvent";
                        List<AdditionalTerm> addTerms = new ArrayList<>();
                        for (final String tokenQMDFF : sx) {
                            if (tokenQMDFF.startsWith("params=")) {
                                paramFile = tokenQMDFF.substring(7).trim();
                            } else if (tokenQMDFF.startsWith("addterm=")) {
                                final String subtok = tokenQMDFF.substring(8).trim();
                                if (subtok.startsWith("electrostatic:")) {
                                    final String chargeFile = subtok.substring(14).trim();
                                    final ElectrostaticEmbedding emb = new ElectrostaticEmbedding(chargeFile);
                                    addTerms.add(emb);
                                } else {
                                    throw new RuntimeException("Unknown additional term " + subtok + " for QMDFF.");
                                }
                            } else {
                                throw new RuntimeException("Unknown token " + tokenQMDFF + " for QMDFF calculation.");
                            }
                        }
            
                        final Tuple<QMDForceFieldParameters, QMDFFData> tup = org.ogolem.qmdff.SetupQMDFF.setup(paramFile);

                        QMDForceField ff = null;
                        if (addTerms.isEmpty()) {
                            ff = new QMDForceField(paramFile, tup.getObject2(), tup.getObject1());
                        } else {
                            ff = new QMDForceField(paramFile, tup.getObject2(), tup.getObject1(), addTerms);
                        }
                        
                        final Tuple<QMDForceFieldParameters,QMDForceField> tupX = new Tuple<>(tup.getObject1(),ff);
                        ffs.add(tupX);
                    }
                } else {
                    throw new RuntimeException("Unknown token " + token + " for EVB-QMDFF calculation.");
                }
            }
            
            if(evb == null){throw new RuntimeException("EVB for EVB-QMDFF should really be set...");}
            if(ffs.isEmpty()){throw new RuntimeException("No force fields for EVB-QMDFF specified.");}
            
            back = new EVBQMDForceField(evb,ffs);
        } else if(sBackend.startsWith("surfacetopadsorptionbiaser:")) {

            double eIncr = FixedValues.NONCONVERGEDENERGY;
            double gradIncr = FixedValues.NONCONVERGEDGRADIENT;

            final String[] s2 = sBackend.substring("surfacetopadsorptionbiaser:".length()).trim().split("\\,");
            for(final String token : s2){
                if(token.isEmpty()){
                    continue;
                } else if(token.startsWith("eincr=")){
                    eIncr = Double.parseDouble(token.substring(6).trim());
                } else if(token.startsWith("gradincr=")){
                    gradIncr = Double.parseDouble(token.substring(9).trim());
                } else {
                    throw new RuntimeException("Unknown token " + token + " for surface top adsorption biaser calculation.");
                }
            }

            back = new HackySurfaceAdsorptionBiaser(eIncr, gradIncr);
        } else if(sBackend.startsWith("chainedfullcartesian:")) {

            final String[] s2 = sBackend.substring("chainedfullcartesian:".length()).trim().split("\\|");

            final List<CartesianFullBackend> backends = new ArrayList<>();
            for(final String token : s2){
                final CartesianFullBackend singleBack = mapStringToBackend(config, token.trim(), params, outFolder);
                backends.add(singleBack);
            }

            back = new ChainedCartesianFullBackend(backends);
        } else if(sBackend.startsWith("tinker:")){
            String sTemp3 = sBackend.substring(7).trim();
            if (sTemp3.equalsIgnoreCase("custom")){
                back = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, false, false, false);
            } else if (sTemp3.equalsIgnoreCase("custom,solvate")){
                back = new TinkerCaller(config, TinkerCaller.METHOD.custom, true, false, false, false);
            } else if (sTemp3.equalsIgnoreCase("custom,fileout")){
                back = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, true, false, false);
            } else if (sTemp3.equalsIgnoreCase("custom,fileout,solvate")){
                back = new TinkerCaller(config, TinkerCaller.METHOD.custom, true, true, false, false);
            } else if (sTemp3.equalsIgnoreCase("ubercustom")){
                back = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, false, true, false);
            } else if (sTemp3.equalsIgnoreCase("ubercustom,fileout")){
                back = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, true, true, false);
            } else if (sTemp3.equalsIgnoreCase("ubercustom2")){
                back = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, false, true, true);
            } else if (sTemp3.equalsIgnoreCase("ubercustom2,fileout")){
                back = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, true, true, true);
            } else {
                System.err.println("Wrong input to configure tinker: " + sTemp3 + ".");
                throw new RuntimeException("Illlegal input to configure tinker as a backend.");
            }
        } else {
            System.err.println("Wrong input for EnergyBackend: " + sBackend + " returning null.");
            back = null;
        }

        return back;
    }
}
