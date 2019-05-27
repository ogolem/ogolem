/**
Copyright (c) 2015-2017, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive.genericfitness;

import java.io.Serializable;
import org.ogolem.properties.Property;

/**
 * A configuration object for each property term.
 * @author Johannes Dieterich
 * @version 2017-10-17
 */
public class FitnessTermConfig<T extends Property> implements Serializable, Cloneable {
    
    private static final long serialVersionUID = (long) 20171017;
    
    private double termWeight = 1.0;
    private boolean useMaxAllowedDiffs = false;
    private boolean useOnlyMaxAllowedDiffs = false;
    private boolean exactEnergyComp = true;
    private boolean referenceToFirst = false;
    private double penaltyPrefactor = 10.0;
    private double penaltyPower = 2.0;
    private boolean dontNormalizePenalty = false;
    
    public FitnessTermConfig(final String[] input) throws Exception {
        
        for(final String s : input){
            if(s.startsWith("//") || s.startsWith("#")){
                continue;
            } else if (s.startsWith("ExactFitnessComp=")) {
                String sTemp2 = s.substring(17).trim();
                try {
                    boolean bExactComp = Boolean.parseBoolean(sTemp2);
                    exactEnergyComp = bExactComp;
                } catch (Exception e) {
                    System.err.println("Wrong input in ExactFitnessComp= . Using default. " + e.toString());
                }
            } else if (s.startsWith("ReferenceEnergiesToFirst=")) {
                String sTemp2 = s.substring(25).trim();
                try {
                    boolean b = Boolean.parseBoolean(sTemp2);
                    referenceToFirst = b;
                } catch (Exception e) {
                    System.err.println("Wrong input in ReferenceEnergiesToFirst= . Using default. " + e.toString());
                }
            } else if (s.startsWith("UseMaxAllowedDiffs=")) {
                final String sTemp2 = s.substring(19).trim();
                try {
                    final boolean bUseDiffs = Boolean.parseBoolean(sTemp2);
                    useMaxAllowedDiffs = bUseDiffs;
                } catch (Exception e) {
                    System.err.println("Wrong input in UseMaxAllowedDiffs= . Using default. " + e.toString());
                }
            } else if (s.startsWith("UseOnlyMaxAllowedDiffs=")) {
                final String sTemp2 = s.substring(23).trim();
                try {
                    final boolean b = Boolean.parseBoolean(sTemp2);
                    useOnlyMaxAllowedDiffs = b;
                } catch (Exception e) {
                    System.err.println("Wrong input in UseOnlyMaxAllowedDiffs= . Using default. " + e.toString());
                }
            } else if (s.startsWith("NormalizePenality=")) {
                final String sTemp2 = s.substring("NormalizePenality=".length()).trim();
                try {
                    final boolean b = Boolean.parseBoolean(sTemp2);
                    dontNormalizePenalty = !b; // b/c it's DONT
                } catch (Exception e) {
                    System.err.println("Wrong input in NormalizePenality= . Using default. " + e.toString());
                }
            } else if (s.startsWith("PenaltyConst=")) {
                final String sTemp2 = s.substring(13).trim();
                try {
                    final double dIncrease = Double.parseDouble(sTemp2);
                    penaltyPrefactor = dIncrease;
                } catch (Exception e) {
                    System.err.println("Wrong input in PenaltyConst= . Using default. " + e.toString());
                }
            } else if (s.startsWith("PenaltyPow=")) {
                final String sTemp2 = s.substring(11).trim();
                try {
                    final double dPow = Double.parseDouble(sTemp2);
                    penaltyPower = dPow;
                } catch (Exception e) {
                    System.err.println("Wrong input in PenaltyPow= . Using default. " + e.toString());
                }
            } else if (s.startsWith("TermWeight=")) {
                final String sTemp2 = s.substring(11).trim();
                try {
                    final double d = Double.parseDouble(sTemp2);
                    termWeight = d;
                } catch (Exception e) {
                    System.err.println("Wrong input in TermWeight= . Using default. " + e.toString());
                }
            } else {
                throw new RuntimeException("Unknown option: " + s);
            }
        }
    }
    
    private FitnessTermConfig(final FitnessTermConfig<T> orig){
        this.exactEnergyComp = orig.exactEnergyComp;
        this.penaltyPower = orig.penaltyPower;
        this.penaltyPrefactor = orig.penaltyPrefactor;
        this.referenceToFirst = orig.referenceToFirst;
        this.termWeight = orig.termWeight;
        this.useMaxAllowedDiffs = orig.useMaxAllowedDiffs;
        this.useOnlyMaxAllowedDiffs = orig.useOnlyMaxAllowedDiffs;
        this.dontNormalizePenalty = orig.dontNormalizePenalty;
    }
    
    @Override
    public FitnessTermConfig<T> clone(){
        return new FitnessTermConfig<>(this);
    }
    
    public double getTermWeight(){
        assert(termWeight > 0.0);
        return termWeight;
    }

    public boolean useMaxAllowedDiffs() {
        return useMaxAllowedDiffs;
    }

    public boolean useOnlyMaxAllowedDiffs() {
        return useOnlyMaxAllowedDiffs;
    }

    public boolean doExactEnergyComp() {
        return exactEnergyComp;
    }

    public boolean doReferenceToFirst() {
        return referenceToFirst;
    }

    public double getPenaltyPrefactor() {
        return penaltyPrefactor;
    }

    public double getPenaltyPower() {
        return penaltyPower;
    }
    
    public boolean normalizePenality(){
        return !dontNormalizePenalty;
    }
    
    public static String[] findBlockFor(final String[] all, final String tag) throws Exception {
        
        int begin = -1;
        int end = -1;
        
        final String begTag = "<" + tag + ">";
        final String endTag = "</" + tag + ">";
        
        // return an empty block if there is nothing to be found
        for(int i = 0; i < all.length; i++){
            if(all[i].trim().startsWith(begTag)){
                begin = i;
            } else if(all[i].trim().startsWith(endTag)){
                end = i;
                break;
            }
        }
        
        if(begin < 0 && end < 0){
            // nothing found
            return new String[0];
        } else if((begin < 0 && end >= 0) || begin >= 0 && end < 0){
            throw new RuntimeException("Tags for" + tag + " are either not properly closed or where never opened.");
        } else if(end < begin){
            throw new RuntimeException("Opening tag MUST come before closing tag.");
        }
        
        final String[] ret = new String[end-begin-1];
        int c = 0;
        for(int x = begin+1; x < end; x++){
            ret[c] = all[x];
            c++;
        }
        
        return ret;
    }
}
