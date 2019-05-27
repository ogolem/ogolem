/**
Copyright (c) 2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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
import java.util.Random;
import static org.ogolem.core.GlobOptAtomics.randomCuttingPlane;
import org.ogolem.generic.GenericMutation;
import org.ogolem.helpers.RandomUtils;

/**
 * This is a phenotype-like algorithm. It cuts the clusters and puts them back
 * together after tilting, moving and rotating.
 * Or, to sum it up, it behaves like a typical Finn on a typical Saturday evening.
 * ;-)
 * The molecular global optimization algorithm is a straightforward genotype. Keep
 * it plain. And Finnish. ;-)
 * The molecular mutation is done following to a sauna algorithm stretching ALL bond
 * lengths of flexible molecules. (yeah... not yet...)
 * @author Johannes Dieterich
 * @version 2015-07-20
 */
public class FinlandGeometryMutation implements GenericMutation<Molecule,Geometry> {
    
    private static final long serialVersionUID = (long) 20140401;
    private static final boolean DEBUG = false;
    private static final double MOVETOUNIVERSE = 100000;
    private static final double ALMOSTUNIVERSE = 90000;
    private static final double FUDGEBLOWFACTOR = 0.95;
    
    private final Random random = new Random();
    private final CollisionDetectionEngine colldetect = new AdvancedPairWise(false, new DummyCollisionStrengthComputer());
    private final double blowDiss;
    private final double blowColl;
    private final int whichGlobOpt;
    
    FinlandGeometryMutation(final int whichCutType, final double blowDiss, final double blowColl){
        this.blowColl = blowColl;
        this.blowDiss = blowDiss;
        this.whichGlobOpt = whichCutType;
    }
    
    FinlandGeometryMutation(final FinlandGeometryMutation orig){
        this.blowColl = orig.blowColl;
        this.blowDiss = orig.blowDiss;
        this.whichGlobOpt = orig.whichGlobOpt;
    }
    
    @Override
    public FinlandGeometryMutation clone() {
        return new FinlandGeometryMutation(this);
    }

    @Override
    public String getMyID() {
        return "FINLAND GEOMETRY MUTATION\n\tcutting type: " + whichGlobOpt;
    }

    @Override
    public Geometry mutate(final Geometry individual) {

        // first get all COMs and put them into matrices/double arrays
        final double[][] allCOMs = new double[3][individual.getNumberOfIndieParticles()];

        for (int i = 0; i < individual.getNumberOfIndieParticles(); i++) {
            final double[] tmp1 = individual.getMoleculeAtPosition(i).getExternalCenterOfMass();
            for (int j = 0; j < 3; j++) {
                allCOMs[j][i] = tmp1[j];
            }
        }

        // get the total weights of all molecules
        final double[] totalWeights = new double[individual.getNumberOfIndieParticles()];
        for (int i = 0; i < individual.getNumberOfIndieParticles(); i++) {
            totalWeights[i] = individual.getMoleculeAtPosition(i).getMyTotalWeight();
        }

        // calculate the "fake" COM
        final double[] myCOM = calculateCOM(allCOMs, totalWeights);

        // move the COMs to the overall COM
        moveCOMToZero(allCOMs, myCOM);
        
        // rotate both sets of COMs randomly
        final double[][] rotated = RandomUtils.randomRotation(allCOMs);

        /*
         * now we want to take a plane cutting through both
         * 1) which heigth on the z axis?
         */
        final double planeHeigth = randomCuttingPlane(whichGlobOpt, rotated, random);

        // check which COMs are above and underneath the plane
        final List<Integer> above = new ArrayList<>(rotated[0].length);
        final List<Integer> under = new ArrayList<>(rotated[0].length);

        for (int i = 0; i < rotated[0].length; i++) {
            if (rotated[2][i] <= planeHeigth) {
                // underneath
                under.add(i);
            } else {
                // above
                above.add(i);
            }
        }

        // create two new seperate matrices for each parent
        final double[][] aboveMat = new double[3][above.size()];
        final double[][] underMat = new double[3][under.size()];

        // fill them
        for (int i = 0; i < aboveMat[0].length; i++) {
            final int which = above.get(i);
            aboveMat[0][i] = rotated[0][which];
            aboveMat[1][i] = rotated[1][which];
            aboveMat[2][i] = rotated[2][which];
        }

        for (int i = 0; i < underMat[0].length; i++) {
            final int which = under.get(i);
            underMat[0][i] = rotated[0][which];
            underMat[1][i] = rotated[1][which];
            underMat[2][i] = rotated[2][which];
        }

        /*
         * now that we have two distinct matrices representing the parts of the geometry,
         * we can start manipulating them.
         * I) rotate both upper parts randomly
         */
        
        // get the total weights
        final double[] weightsAbove = new double[above.size()];

        // fill them
        for (int i = 0; i < above.size(); i++) {
            final int which = above.get(i);
            weightsAbove[i] = totalWeights[which];
        }
        
        // calculate the fake COMs
        final double[] partCOM = calculateCOM(aboveMat, weightsAbove);

        // move them to that COM
        moveCOMToZero(aboveMat, partCOM);
        
        final double[][] rotAboveMat = RandomUtils.randomRotation(aboveMat);

        /*
         * III) translate the above parts to somewhere in space
         */

        // draw some BIG random numbers
        final double[] translation = new double[3];
        RandomUtils.randomVector(translation, MOVETOUNIVERSE);

        /*
         * IV) put the parts back together
         */
        // first move the above parts to the right spot
        addValueToCOMs(rotAboveMat, translation);

        // generate mutated geometry
        final Geometry mutated = new Geometry(individual);

        // put the new COMs into the correct positions
        for(int i = 0; i < above.size(); i++){
            final int which = above.get(i);
            final double[] com = new double[3];
            com[0] = rotAboveMat[0][i];
            com[1] = rotAboveMat[1][i];
            com[2] = rotAboveMat[2][i];
            mutated.getMoleculeAtPosition(which).setExternalCenterOfMass(com);
        }

        for(int i = 0; i < under.size(); i++){
            final int which = under.get(i);
            final double[] com = new double[3];
            com[0] = underMat[0][i];
            com[1] = underMat[1][i];
            com[2] = underMat[2][i];
            mutated.getMoleculeAtPosition(which).setExternalCenterOfMass(com);
        }

        // now we have geometries that are FOR SURE dissociated into two parts
        final CartesianCoordinates mutatedCartes = mutated.getCartesians();

        // get the collision infos anyway, we need the distances
        final CollisionInfo childColl = colldetect.checkForCollision(mutatedCartes, blowColl, mutated.getBondInfo());

        // now we want to figure out what the smallest distance of the BIG distances is
        final double[][] collInfo = childColl.getPairWiseDistances();

        final int[] whichAts = new int[2];
        double smallesDist = Double.MAX_VALUE;
        for(int i = 0; i < collInfo.length-1; i++){
            for(int j = i; j < collInfo.length; j++){
                // > 90 000 sounds like a reasonable cutoff
                if(collInfo[i][j] > ALMOSTUNIVERSE && collInfo[i][j] < smallesDist){
                    whichAts[0] = i;
                    whichAts[1] = j;
                    smallesDist = collInfo[i][j];
                }
            }
        }

        final double[] posConn1 = mutatedCartes.getXYZCoordinatesOfAtom(whichAts[0]);
        final double[] posConn2 = mutatedCartes.getXYZCoordinatesOfAtom(whichAts[1]);

        final double atRad1 = AtomicProperties.giveRadius(mutatedCartes.getAtomNumberOfAtom(whichAts[0]));
        final double atRad2 = AtomicProperties.giveRadius(mutatedCartes.getAtomNumberOfAtom(whichAts[1]));

        // Using 0.95 as a fudge factor should just work (TM) ;-)
        final double wantedLength = FUDGEBLOWFACTOR * blowDiss * (atRad1 + atRad2);

        // we want to calculate the ratios of delta x, delta y, delta z
        final double diffX = posConn2[0] - posConn1[0];
        final double diffY = posConn2[1] - posConn1[1];
        final double diffZ = posConn2[2] - posConn1[2];

        final double ratioXY = diffY/diffX;
        final double ratioXZ = diffZ/diffX;

        // now we can calculate the new deltas
        final double deltaXNew = wantedLength*wantedLength / (1 + ratioXY*ratioXY + ratioXZ*ratioXZ);
        final double deltaYNew = deltaXNew * ratioXY;
        final double deltaZNew = deltaXNew * ratioXZ;

        // now we calculate the new position of atom two
        final double[] newAt2 = new double[3];
        newAt2[0] = deltaXNew + posConn1[0];
        newAt2[1] = deltaYNew + posConn1[1];
        newAt2[2] = deltaZNew + posConn1[2];

        // now we can calculate how much we need to move the part
        final double[] movement = new double[3];
        for(int i = 0; i < 3; i++){
            movement[i] = Math.abs(newAt2[i] - posConn2[i]);
        }

        for (final int which : above) {
            final double[] com = mutated.getMoleculeAtPosition(which).getExternalCenterOfMass();
            com[0] -= movement[0];
            com[1] -= movement[1];
            com[2] -= movement[2];
            mutated.getMoleculeAtPosition(which).setExternalCenterOfMass(com);
        }

        // lets move it completely to the total/fake COM, to make sure it is not somewhere in space
        final double[][] newCOMs = new double[3][mutated.getNumberOfIndieParticles()];
        for(int i = 0; i < mutated.getNumberOfIndieParticles(); i++){
            final double[] com = mutated.getMoleculeAtPosition(i).getExternalCenterOfMass();
            for(int j = 0; j < 3; j++){
                newCOMs[j][i] = com[j];
            }
        }

        final double[] fakeCOM = calculateCOM(newCOMs,totalWeights);
        moveCOMToZero(newCOMs, fakeCOM);

        for(int i = 0; i < mutated.getNumberOfIndieParticles(); i++){
            final double[] com = new double[3];
            com[0] = newCOMs[0][i];
            com[1] = newCOMs[1][i];
            com[2] = newCOMs[2][i];
            mutated.getMoleculeAtPosition(i).setExternalCenterOfMass(com);
        }


        if(DEBUG){
            // dump resulting geom
            System.out.println("DEBUG: resulting finland geom is: ");
            final String[] g1 = mutated.makePrintableAbsoluteCoord(true);
            for(final String s : g1){
                System.out.println(s);
            }
        }
        
        return mutated;
    }
    
    private static double[] calculateCOM(final double[][] allCOMs, final double[] allWeights){
        
        final double[] com = new double[3];
        // formula: \mathbf{R}=\dfrac{\sum_i mi_i\cdot \mathbf{r}_i}{\sum_i m_i}
        final int noCOMs = allWeights.length;
        final double[] numerator = {0.0, 0.0, 0.0};

        double denominator = 0.0;
        for (int i = 0; i < noCOMs; i++) {
            // denominator: simply add up all weights for all molecules
            denominator += allWeights[i];

            // numerator: multiply position with weight and add up
            for (int j = 0; j < 3; j++) {
                numerator[j] += allWeights[i] * allCOMs[j][i];
            }
        }
        
        for (int k = 0; k < 3; k++) {
            // divide all numerators by denominator
            com[k] = numerator[k] / denominator;
        }

        return com;
    }

    private static void moveCOMToZero(final double[][] daCOMs, final double[] daTotalCOM){

        for(int j = 0; j < 3; j++){
            for(int i = 0; i < daCOMs[0].length; i++){
                daCOMs[j][i] -= daTotalCOM[j];
            }
        }
    }

    private static void addValueToCOMs(final double[][] daCOMs, final double[] daValue){

        for(int j = 0; j < 3; j++){
            for(int i = 0 ; i < daCOMs[0].length; i++){
                daCOMs[j][i] += daValue[j];
            }
        }
    }
}
