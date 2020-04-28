/**
Copyright (c) 2014, J. M. Dieterich
              2016-2020, J. M. Dieterich and B. Hartke
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
import java.util.Iterator;
import java.util.List;
import static org.ogolem.core.GlobOptAtomics.findOptimalPlaneHeight;
import static org.ogolem.core.GlobOptAtomics.randomCuttingPlane;
import static org.ogolem.core.GlobOptAtomics.randomRotation;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.Lottery;

/**
 * Phenotype algorithm with implicit exchange of atom types. Like Sweden, just different. ;-)
 * @author Johannes Dieterich
 * @version 2020-04-25
 */
public class LaplandGeometryXOver implements GenericCrossover<Molecule,Geometry>{
    
    private static final long serialVersionUID = (long) 20200425;
    private static final boolean DEBUG = false;
    private final Lottery random = Lottery.getInstance();
    private final int whichGlobOpt;

    LaplandGeometryXOver(final int whichGlobOpt){
        this.whichGlobOpt = whichGlobOpt;
    }
    
    LaplandGeometryXOver(final LaplandGeometryXOver orig){
        this.whichGlobOpt = orig.whichGlobOpt;
    }
    
    @Override
    public LaplandGeometryXOver clone() {
        return new LaplandGeometryXOver(this);
    }

    @Override
    public String getMyID() {
        return "LAPLAND GEOMETRY CROSSOVER\n\tcutting type: " + whichGlobOpt;
    }

    @Override
    public Tuple<Geometry, Geometry> crossover(final Geometry mother, final Geometry father, final long futureID) {
        
        // create the children
        Geometry gChild1 = new Geometry(father);
        gChild1.setID(futureID);
        Geometry gChild2 = new Geometry(mother);
        gChild2.setID(futureID);

        CartesianCoordinates cartesMother = gChild2.getCartesians();
        CartesianCoordinates cartesFather = gChild1.getCartesians();

        cartesMother.moveCoordsToCOM();
        cartesFather.moveCoordsToCOM();

        // rotate the two
        randomRotation(cartesMother);
        randomRotation(cartesFather);


         // translate it to geometries again
         gChild2 = new Geometry(cartesMother, gChild2.getID(), gChild2.getNumberOfIndieParticles(),
                cartesMother.getAllAtomsPerMol(), gChild2.getAllFlexies(), gChild2.getExplicitDoFs(), gChild2.getAllConstraints(false),
                gChild2.getAllConstraintsXYZ(false), gChild2.getSIDs(), gChild2.getBondInfo());
         gChild1 = new Geometry(cartesFather, gChild1.getID(), gChild1.getNumberOfIndieParticles(),
                cartesFather.getAllAtomsPerMol(), gChild1.getAllFlexies(), gChild2.getExplicitDoFs(), gChild1.getAllConstraints(false),
                gChild1.getAllConstraintsXYZ(false), gChild1.getSIDs(), gChild1.getBondInfo());

         // in case of environments, add it in at this point
        if(mother.containsEnvironment()){
            // the environment
            final List<Environment> envChilds = mother.getEnvironment().createOffspring(father.getEnvironment());
            gChild1.setEnvironment(envChilds.get(0));
            gChild2.setEnvironment(envChilds.get(1));
        }

        // first get all COMs and put them into matrices/double arrays
        double[][] daMotherCOMs = new double[3][mother.getNumberOfIndieParticles()];
        double[][] daFatherCOMs = new double[3][mother.getNumberOfIndieParticles()];

        for(int i = 0; i < mother.getNumberOfIndieParticles(); i++){
            final double[] daTemp1 = gChild2.getMoleculeAtPosition(i).getExternalCenterOfMass();
            final double[] daTemp2 = gChild1.getMoleculeAtPosition(i).getExternalCenterOfMass();
            for(int j = 0; j < 3; j++){
                daMotherCOMs[j][i] = daTemp1[j];
                daFatherCOMs[j][i] = daTemp2[j];
            }
        }

         /*
          * now we want to take a plane cutting through both
          * 1) which height on the z axis?
          */
        final int iNoOfMolecules = father.getNumberOfIndieParticles();
        final ArrayList<Integer> alFatherAbove = new ArrayList<>(iNoOfMolecules);
        final ArrayList<Integer> alFatherUnder = new ArrayList<>(iNoOfMolecules);
        
        double dPlaneHeightFather;
        int ec = 0;
        boolean cont;
        int iFatherCOMs = -1;
        do{
            dPlaneHeightFather = randomCuttingPlane(whichGlobOpt, daFatherCOMs, random);
            // check which COMs are above and underneath the plane (for the father)
            for (int i = 0; i < iNoOfMolecules; i++) {
                if (daFatherCOMs[2][i] <= dPlaneHeightFather) {
                    // underneath
                    alFatherUnder.add(i);
                } else {
                    // above
                    alFatherAbove.add(i);
                }
            }

            iFatherCOMs = alFatherUnder.size();
            if (iFatherCOMs == 0) {
                // then there was none underneath the plane, this is NOT ok.
                cont = true;
                alFatherUnder.clear();
                alFatherAbove.clear();
            } else if (iFatherCOMs >= iNoOfMolecules) {
                // all are above the plane, also not OK.
                cont = true;
                alFatherUnder.clear();
                alFatherAbove.clear();
            } else{
                cont = false;
            }
            
            ec++;
        } while(cont && ec < FixedValues.MAXTOEMERGENCY);
        
        if(cont){
            System.err.print("Too many attemps in finding a first plane height!");
            return new Tuple<>(null,null);
        }
        
        // now move the plane in the mother geometry till enough COMs are above/underneath
        final double dPlaneHeightMother = findOptimalPlaneHeight(daMotherCOMs, iFatherCOMs);

        if(Double.isNaN(dPlaneHeightMother)){
            // two COMs were too close together, return null'd geometries
            return new Tuple<>(null,null);
        }

        final List<Integer> alMotherUnder = new ArrayList<>(iNoOfMolecules);
        final List<Integer> alMotherAbove = new ArrayList<>(iNoOfMolecules);

        for (int ii = 0; ii < iNoOfMolecules; ii++) {
            if (daMotherCOMs[2][ii] <= dPlaneHeightMother) {
                // underneath
                alMotherUnder.add(ii);
            } else {
                // above
                alMotherAbove.add(ii);
            }
        }

        /*
         * still we need to exchange parts though
         * remember that we need to adjust with the cutting height (so
         * substracting it from the z component of the mother COMs)
         */

        int iMotherCOMs = alMotherUnder.size();

        if (iMotherCOMs != iFatherCOMs) {
            // we f***** it up
            System.err.println("DEBUG: We f***** it up in Lapland. Mother: " + iMotherCOMs + " father: " + iFatherCOMs);
            return new Tuple<>(null,null);
        }

        // exchange the molecules but keep in mind that the order needs to
        // be conserved!
        
        final Iterator<Integer> itFatherAbove = alFatherAbove.iterator();
        final Iterator<Integer> itMotherAbove = alMotherAbove.iterator();

        while (itFatherAbove.hasNext()){
            
            // get infos
            final int iMol1 = itFatherAbove.next();
            final double[] daCOM1 = gChild1.getMoleculeAtPosition(iMol1).getExternalCenterOfMass().clone();
            final double[] daEulers1 = gChild1.getMoleculeAtPosition(iMol1).getOrientation().clone();

            final int iMol2 = itMotherAbove.next();
            final double[] daCOM2 = gChild2.getMoleculeAtPosition(iMol2).getExternalCenterOfMass().clone();
            final double[] daEulers2 = gChild2.getMoleculeAtPosition(iMol2).getOrientation().clone();

            // change
            gChild1.getMoleculeAtPosition(iMol1).setExternalCenterOfMass(daCOM2);
            gChild1.getMoleculeAtPosition(iMol1).setOrientation(daEulers2);

            gChild2.getMoleculeAtPosition(iMol2).setExternalCenterOfMass(daCOM1);
            gChild2.getMoleculeAtPosition(iMol2).setOrientation(daEulers1);
        }

        // adding huge fitness values, just in case something goes really wrong
        gChild1.setFitness(Double.MAX_VALUE);
        gChild2.setFitness(Double.MAX_VALUE);

        if(DEBUG){
            final String[] saCoords1 = gChild1.makePrintableAbsoluteCoord(true);
            final String[] saCoords2= gChild2.makePrintableAbsoluteCoord(true);

            System.out.println("DEBUG: First child coming:");
            for(final String s: saCoords1){
                System.out.println(s);
            }

            System.out.println("DEBUG: Second child coming:");
            for(final String s: saCoords2){
                System.out.println(s);
            }
        }
        
        return new Tuple<>(gChild1,gChild2);
    }

    @Override
    public short hasPriority() {
        return -1;
    }
}
