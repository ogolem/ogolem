/**
Copyright (c) 2016-2017, J. M. Dieterich and B. Hartke
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
package org.ogolem.properties;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import org.ogolem.core.FixedValues;
import org.ogolem.helpers.Tuple;

/**
 * The energy ordering of a few structures (including their gaps).
 * @author Johannes Dieterich
 * @version 2017-10-09
 */
public class EnergyOrder implements Property {

    private static final long serialVersionUID = (long) 20171009;
    private static final boolean DEBUG = false;
    
    public static final double WRONGORDERPENALTY = 10.0;
    public static final double RELENERGYWEIGHT = 1.0;
    
    private final List<Tuple<String,Double>> energies;
    private final boolean considerRelEnergies;
    
    public EnergyOrder(final List<Tuple<String,Double>> energies, final boolean considerRelEnergies){
        assert(energies != null);
        assert(!energies.isEmpty());
        this.energies = energies;
        this.considerRelEnergies = considerRelEnergies;
        
        final Comparator<Tuple<String,Double>> comp = new EnergyComparator();
        Collections.sort(energies, comp);
    }
    
    private EnergyOrder(final EnergyOrder orig){
        this.energies = new ArrayList<>();
        this.considerRelEnergies = orig.considerRelEnergies;
        for(final Tuple<String,Double> tup : orig.energies){
            final String key = tup.getObject1();
            final double e = tup.getObject2();
            final Tuple<String,Double> t = new Tuple<>(key,e);
            energies.add(t);
        }
        
        final Comparator<Tuple<String,Double>> comp = new EnergyComparator();
        Collections.sort(energies, comp);
    }
    
    @Override
    public EnergyOrder clone() {
        return new EnergyOrder(this);
    }

    @Override
    public double getValue() {
        System.out.println("WARNING: EnergyOrder does not support getValue();");
        return 0.0;
    }

    @Override
    public double signedDifference(Property p) {
        
        return absoluteDifference(p);
    }
    
    private static int location(final List<Tuple<String,Double>> ll, final String key){
        
        int loc = -1;
        for(int i = 0; i < ll.size(); i++){
            final Tuple<String,Double> ent = ll.get(i);
            if(ent.getObject1().equalsIgnoreCase(key)){
                return i;
            }
        }
        
        return loc;
    }

    @Override
    public double absoluteDifference(Property p) {
        
        if(!(p instanceof EnergyOrder)) {throw new IllegalArgumentException("Property should be an instance of EnergyOrder!");}
        
        if(DEBUG){System.out.println("DEBUG: Calculating difference of energy orders...");}
        
        double diff = 0.0;
        
        final EnergyOrder eo = (EnergyOrder) p;
        for(int i = 0; i < energies.size(); i++){
            final Tuple<String,Double> ent = energies.get(i);
            final String key = ent.getObject1();
            final double e = ent.getObject2();
            final int locO = location(eo.energies, key);
            if(locO < 0){
                // ok, this spells trouble
                System.err.println("ERROR: Couldn't find key " + key + " in other set of energies. This is NOT supported.");
                diff += FixedValues.NONCONVERGEDENERGY;
                continue;
            }
            
            final double eO = eo.energies.get(locO).getObject2();
            
            if(eO >= FixedValues.NONCONVERGEDENERGY || e >= FixedValues.NONCONVERGEDENERGY){
                // just as a safety net
                diff += FixedValues.NONCONVERGEDENERGY;
                if(DEBUG){
                    System.out.println("DEBUG: setting " + eO + " from " + i + " to NONCONVERGED.");
                }
                continue;
            }
            
            if(locO != i){
                if(DEBUG){
                    System.out.println("DEBUG: energy order wrong " + locO + " vs " + i + ", setting WRONGORDERPENALTY.");
                }
                diff += WRONGORDERPENALTY;
                if(considerRelEnergies){
                
                    if(i < energies.size() - 1){
                    
                        final double myDiff = Math.abs(e - energies.get(i+1).getObject2());
                        final double theirDiff = Math.abs(eo.energies.get(i).getObject2() - eo.energies.get(i+1).getObject2());
                        final double comp = Math.abs(myDiff-theirDiff);
                    
                        diff += RELENERGYWEIGHT*comp;
                    }
                }
            } else if(considerRelEnergies){
                
                if(i < energies.size() - 1){
                    
                    final double myDiff = Math.abs(e - energies.get(i+1).getObject2());
                    final double theirDiff = Math.abs(eo.energies.get(i).getObject2() - eo.energies.get(i+1).getObject2());
                    final double comp = Math.abs(myDiff-theirDiff);
                    
                    diff += RELENERGYWEIGHT*comp;
                }
            }
        }
        
        if(DEBUG){
            System.out.println("DEBUG: energy order difference: " + diff);
        }
        
        return diff;
    }

    @Override
    public boolean makeSensible() {
        
        boolean allSane = true;
        for(final Tuple<String,Double> tup : energies){
            final double e = tup.getObject2();
            if(Double.isNaN(e) || Double.isInfinite(e)){
                allSane = false;
                if(DEBUG){
                    System.out.println("DEBUG: setting " + e + " to NONCONVERGED.");
                }
                tup.setObject2(FixedValues.NONCONVERGEDENERGY);
            }
        }
        
        return !allSane;
    }

    @Override
    public String printableProperty() {
        
        String st = "";
        for(final Tuple<String,Double> tup : energies){
            final double e = tup.getObject2();
            st += "" + tup.getObject1() + ":\t" + e + "\n";
        }
        
        return st;
    }

    @Override
    public String name() {
        return "ENERGY ORDER";
    }
    
    private static class EnergyComparator implements Comparator<Tuple<String,Double>> {

        @Override
        public int compare(Tuple<String, Double> o1, Tuple<String, Double> o2) {
            
            final double e1 = o1.getObject2();
            final double e2 = o2.getObject2();
            if(e1 < e2){
                return -1;
            } else if(e1 == e2){
                // unlikely...
                return 0;
            } else {
                return 1;
            }
        }
    }
    
    public boolean shouldConsiderRelEnergies(){
        return this.considerRelEnergies;
    }
}
