/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2012, J. M. Dieterich
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

import java.util.ArrayList;
import java.util.Iterator;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.generic.DiscreteProblem;

/**
 * The actual switch.
 * @author Johannes Dieterich
 * @version 2013-11-22
 */
public class Switch extends DiscreteProblem<Color> {

    private static final long serialVersionUID = (long) 20131122;

    private final Backbone back;

    private final ArrayList<Sidechain> alSides;

    private ArrayList<Color> alColorCode;

    private ArrayList<CartesianCoordinates> alCartes = null;

    private double dFitness;

    private double dS0S1EnergyCis;

    private double dS0S1EnergyTrans;

    private long id;
    private long fatherID; // TODO enable
    private long motherID; // TODO enable

    public Switch(final Backbone backbone, final ArrayList<Color> alColors){
        this.back = backbone;
        this.alSides = new ArrayList<>(alColors.size());
        this.alColorCode = alColors;

        // fill the sidechains in
        final Iterator<Color> itColors = alColors.iterator();
        while(itColors.hasNext()){
            alSides.add(new Sidechain(itColors.next()));
        }
        this.dS0S1EnergyCis = Double.NaN;
        this.dS0S1EnergyTrans = Double.NaN;
    }

    public Switch(final Switch refSwitch){
        this.back = new Backbone(refSwitch.getBackbone());
        this.alColorCode = refSwitch.getCopyOfColors();
        this.alSides = new ArrayList<>(alColorCode.size());

        // fill the sidechains in
        final Iterator<Color> itColors = alColorCode.iterator();
        while(itColors.hasNext()){
            alSides.add(new Sidechain(itColors.next()));
        }

        this.dFitness = refSwitch.getFitness();
        this.dS0S1EnergyCis = refSwitch.getS0S1EnergyCis();
        this.dS0S1EnergyTrans = refSwitch.getS0S1EnergyTrans();
        this.id = refSwitch.id;
        this.fatherID = refSwitch.fatherID;
        this.motherID = refSwitch.motherID;
    }
    
    @Override
    public Switch clone(){
        return new Switch(this);
    }
    
    /*
     * Getters/Setters
     */
    @Override
    public double getFitness(){
        return dFitness;
    }

    double getS0S1EnergyCis(){
        return dS0S1EnergyCis;
    }

    double getS0S1EnergyTrans(){
        return dS0S1EnergyTrans;
    }

    @Override
    public void setFitness(final double fitness){
        this.dFitness = fitness;
    }

    public void setS0S1EnergyCis(final double s0s1Cis){
        this.dS0S1EnergyCis = s0s1Cis;
    }

    public void setS0S1EnergyTrans(final double s0s1Trans){
        this.dS0S1EnergyTrans = s0s1Trans;
    }

    @Override
    public long getID(){
        return id;
    }
    
    @Override
    public long getMotherID(){
        return motherID;
    }
    
    @Override
    public long getFatherID(){
        return fatherID;
    }

    @Override
    public void setID(final long id){
        this.id = id;
    }

    void setAllColors(final ArrayList<Color> alColors){
        alColorCode = alColors;
    }

    Backbone getBackbone(){
        return back;
    }

    ArrayList<Sidechain> getSidechains(){
        return alSides;
    }

    void setCartesians(ArrayList<CartesianCoordinates> alCartesians){
        alCartes = alCartesians;
    }

    /*
     * Methods
     */

    /**
     * Returns the colorcode of this Switch as a String.
     * @return the color code
     */
    String myColorCode(){
        final Iterator<Color> itCode = alColorCode.iterator();

        final StringBuffer sBuff = new StringBuffer();

        while(itCode.hasNext()){
            sBuff.append(itCode.next().getThisColor());
        }

        final String sColorCode = sBuff.toString();

        return sColorCode;
    }

    /**
     * Creates a deep copy of the color code.
     * @return a deep copies ArrayList
     */
    ArrayList<Color> getCopyOfColors(){

        final ArrayList<Color> alColorCopy = new ArrayList<>(alColorCode.size());

        final Iterator<Color> itColor = alColorCode.iterator();
        while(itColor.hasNext()){
            alColorCopy.add(new Color(itColor.next()));
        }

        return alColorCopy;
    }

    public void changeToRandomColors(){
        // assumes that there is a sufficient amount of colors already in the colorcode
        final int iNoOfColors = alColorCode.size();

        final ColorPalette palette = ColorPalette.getReference();

        for(int i = 0; i < iNoOfColors; i++){
            // get a new color
            final Color col = palette.getRandomColor();

            // translate it to sidechain
            final Sidechain sidech = new Sidechain(col);

            // put the sidechain and color in the arays
            alColorCode.set(i, col);
            alSides.set(i, sidech);
        }

    }

    /**
     * Creates a set of printable cartesians.
     * @return An array of Strings. The first is the cis, the second trans the isomer.
     */
    public String[][] createPrintableCisTrans(final double dBlowBondDetect){

        if(alCartes == null){

            // first try to glue it together
            alCartes = MagicGlue.glueTogether(back, alSides, dBlowBondDetect);
            
            if(alCartes == null){
                // this is really not glueable
                System.out.println("INFO: This is really not glueable.");
                final String[][] saRet = new String[2][1];
                saRet[0][0] = "Cartesians were null, fitness is " + dFitness;
                saRet[1][0] = "Cartesians were null, fitness is " + dFitness;
                return saRet;
            }
        }

        if (alCartes.get(0) == null) {
            final String[][] saRet = new String[2][1];
            saRet[0][0] = "Cartesians were null, fitness is " + dFitness;
            saRet[1][0] = "Cartesians were null, fitness is " + dFitness;
            return saRet;
        } else if (alCartes.get(1) == null) {

            final String[][] saRet = new String[2][1];
            saRet[0][0] = "Cartesians were null, fitness is " + dFitness;
            saRet[1][0] = "Cartesians were null, fitness is " + dFitness;
            return saRet;
        }

        // at least we have some useful cartesian coordinates
        final CartesianCoordinates cartesCis = alCartes.get(0);
        final CartesianCoordinates cartesTrans = alCartes.get(1);
        cartesCis.setEnergy(dFitness);
        cartesTrans.setEnergy(dFitness);

        final String[][] saOutput = new String[2][];
        saOutput[0] = cartesCis.createPrintableCartesians();
        saOutput[1] = cartesTrans.createPrintableCartesians();

        /*
         * actually, since we in this case really have a fitness, we change the lines
         */
        saOutput[0][1] = "Fitness: " + dFitness + "\tS0S1Excitation(cis): " + dS0S1EnergyCis + "\tS0S1Excitation(trans): " + dS0S1EnergyTrans;
        saOutput[1][1] = "Fitness: " + dFitness + "\tS0S1Excitation(cis): " + dS0S1EnergyCis + "\tS0S1Excitation(trans): " + dS0S1EnergyTrans;

        return saOutput;
    }

    public String[] createPrintableColors(){

        final String[] saColors = new String[alColorCode.size()+1];

        saColors[0] = "Fitness: " + dFitness + "\tS0S1Excitation(cis): " + dS0S1EnergyCis + "\tS0S1Excitation(trans): " + dS0S1EnergyTrans;

        for(int i = 1; i < alColorCode.size()+1; i++){
            saColors[i] = alColorCode.get(i-1).getThisColorString();
        }

        return saColors;
    }
    
    @Override
    public double[] getGenomeAsDouble(){
        return null;
    }
    
    @Override
    public Color[] getGenomeCopy(){
        
        final Color[] genome = new Color[alColorCode.size()];
        this.alColorCode.toArray(genome);
        
        return genome;
    }
    
    @Override
    public void setGenome(final Color[] genome){
        
        assert(genome != null);
        assert(genome.length == alColorCode.size());
        
        // assumes that there is a sufficient amount of colors already in the colorcode
        for(int i = 0; i < genome.length; i++){
            final Color col = genome[i];
            // translate it to sidechain
            final Sidechain sidech = new Sidechain(col);

            // put the sidechain and color in the arays
            alColorCode.set(i, col);
            alSides.set(i, sidech);
        }

    }

    @Override
    public void setFatherID(long fatherID) {
        this.fatherID = fatherID;
    }

    @Override
    public void setMotherID(long motherID) {
        this.motherID = motherID;
    }
}
