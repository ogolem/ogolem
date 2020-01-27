/**
Copyright (c) 2014-2015, J. M. Dieterich
              2016, J. M. Dieterich and B. Hartke
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
import org.junit.Test;
import static org.junit.Assert.*;
import static org.ogolem.core.Constants.ANGTOBOHR;

/**
 *
 * @author Johannes Dieterich
 * @version 2016-09-03
 */
public class RigidBodyCoordinatesTest {
    
    private static final double NUMACC = 1e-6;
    private static final double NUMACCLOOSE = 1e-4;
    private final Geometry gWaterWrong;
    private final Geometry gWater2;
    private final Geometry gWater2Distorted;
    private final Geometry gWater20;
    private final TIP3PForceField tip3p = new TIP3PForceField();
    private final TIP4PForceField tip4p = new TIP4PForceField();
    private final RigidBodyCoordinates rigidEmpt;
    private final RigidBodyCoordinates rigidW2TIP3P;
    private final RigidBodyCoordinates rigidW2DisTIP3P;
    private final RigidBodyCoordinates rigidW20TIP3P;
    private final RigidBodyCoordinates rigidW2TIP4P;
    private final RigidBodyCoordinates rigidW2DisTIP4P;
    private final RigidBodyCoordinates rigidW20TIP4P;
    private final double[] gradW2Phenix = new double[]{
            -2.49299292E-09,
            -1.74202093E-08,
            1.09019958E-08,
            -1.00777087E-09, // ATTENTION, PHENIX Euler definition different!
            3.99753608E-09, // ATTENTION, PHENIX Euler definition different!
            -1.47489343E-09, // ATTENTION, PHENIX Euler definition different!
            2.49299292E-09,
            1.74202093E-08,
            -1.09019958E-08,
            -4.38705738E-10, // ATTENTION, PHENIX Euler definition different!
            2.52028243E-09, // ATTENTION, PHENIX Euler definition different!
            4.31826158E-10 // ATTENTION, PHENIX Euler definition different!
        };
    private final double[] gradW2DisPhenix = new double[]{
            -4.61163116E-04,
            -1.06346468E-03,
             3.18424893E-04,
             2.27286480E-03, // ATTENTION, PHENIX Euler definition different!
            -1.51189405E-03, // ATTENTION, PHENIX Euler definition different!
            -1.47368794E-03, // ATTENTION, PHENIX Euler definition different!
            4.61163116E-04,
            1.06346468E-03,
            -3.18424893E-04,
            -9.05349967E-04, // ATTENTION, PHENIX Euler definition different!
            5.10955811E-04, // ATTENTION, PHENIX Euler definition different!
            -1.25210837E-03 // ATTENTION, PHENIX Euler definition different!
        };
    

    
    public RigidBodyCoordinatesTest(){
        // build a few test geometries...
        this.gWater2 = getWater2TIP4P();
        this.gWater2Distorted = getWater2TIP4PDistorted();
        this.gWater20 = getWater20TIP4P();
        this.gWaterWrong = getNoWaterTIP4P();
        this.rigidEmpt = new RigidBodyCoordinates(tip4p);
        this.rigidW2TIP3P = new RigidBodyCoordinates(tip3p);
        this.rigidW2DisTIP3P = new RigidBodyCoordinates(tip3p);
        this.rigidW20TIP3P = new RigidBodyCoordinates(tip3p);
        this.rigidW2TIP4P = new RigidBodyCoordinates(tip4p);
        this.rigidW2DisTIP4P = new RigidBodyCoordinates(tip4p);
        this.rigidW20TIP4P = new RigidBodyCoordinates(tip4p);
    }
    
    /**
     * Test of numberOfActiveCoordinates method, of class RigidBodyCoordinates.
     */
    @Test
    public void testNumberOfActiveCoordinates() {
        System.out.println("numberOfActiveCoordinates");

        // should be empty
        final int expResEmpty = 0;
        final int resEmpty = rigidEmpt.numberOfActiveCoordinates(gWaterWrong);
        assertEquals(expResEmpty, resEmpty);
        
        // should be 12 and 120, respectively
        final int expResW2 = 12;
        final int expResW20 = 120;
        final int resW2 = rigidW2TIP4P.numberOfActiveCoordinates(gWater2);
        final int resW20 = rigidW2TIP4P.numberOfActiveCoordinates(gWater20);
        assertEquals(expResW2, resW2);
        assertEquals(expResW20, resW20);
    }

    /**
     * Test of getActiveCoordinates method, of class RigidBodyCoordinates.
     */
    @Test
    public void testGetActiveCoordinates() {
        System.out.println("getActiveCoordinates");
        // basically just checks if we get arrays of proper length back, can certainly be extended
        // should be 12 and 120, respectively
        final int expResW2 = 12;
        final int expResW20 = 120;
        final int resW2 = rigidW2TIP4P.getActiveCoordinates(gWater2).length;
        final int resW20 = rigidW2TIP4P.getActiveCoordinates(gWater20).length;
        assertEquals(expResW2, resW2);
        assertEquals(expResW20, resW20);
    }

    /**
     * Test of updateActiveCoordinates method, of class RigidBodyCoordinates.
     */
    /*@Test
    public void testUpdateActiveCoordinates() {
        System.out.println("updateActiveCoordinates");
        
        // required to initialize..
        rigidW2TIP4P.getActiveCoordinates(gWater2.clone());

        // again, very basic...
        final Geometry g2Work = gWater2.clone();
        final double[] comEuler = new double[12];
        rigidW2TIP4P.updateActiveCoordinates(g2Work, comEuler);
        
        // check if they are all zero
        final double[] extCOM = g2Work.getAllExtCoords();
        for(int i = 0; i < 12; i++){
            assertEquals(0.0,extCOM[i],1e-10);
        }
        
        // put some random ones in
        final Random r = new Random();
        for(int i = 0; i < 12; i++){
            comEuler[i] = (r.nextBoolean()) ? r.nextDouble() : -r.nextDouble(); // will work even with the Eulers
        }
        
        rigidW2TIP4P.updateActiveCoordinates(g2Work, comEuler);
        
        // check if they are all the same
        final double[] extCOMRand = g2Work.getAllExtCoords();
        for(int i = 0; i < 12; i++){
            
            //TODO since the rigidify fix this will need to be compared more... intricately...
            assertEquals(comEuler[i],extCOMRand[i],1e-10);
        }
    }*/

    /**
     * Test of fitness method, of class RigidBodyCoordinates.
     */
    @Test
    public void testFitness_doubleArr_int() {
        System.out.println("fitness");

        // make sure that w2 and w20 have the right energy
        final double[] extCoord2 = rigidW2TIP4P.getActiveCoordinates(gWater2.clone());
        final double[] extCoord20 = rigidW20TIP4P.getActiveCoordinates(gWater20.clone());
        
        final double e2 = rigidW2TIP4P.fitness(extCoord2, 0);
        final double e20 = rigidW20TIP4P.fitness(extCoord20, 0);
        
        // from DJW
        final double expResult2 = -0.009936231756023471;
        final double expResult20 = -0.3325038695132564;
        
        assertEquals(expResult2, e2, NUMACC);
        assertEquals(expResult20, e20, NUMACC);
    }

    /**
     * Test of fitness method, of class RigidBodyCoordinates.
     */
    @Test
    public void testFitness_Geometry_boolean() {
        System.out.println("fitness");
        
        // make sure that w2 and w20 have the right energy

        final double e2 = rigidW2TIP4P.fitness(gWater2.clone(), true).getFitness();
        final double e20 = rigidW20TIP4P.fitness(gWater20.clone(), true).getFitness();
        
        // from DJW
        final double expResult2 = -0.009936231756023471;
        final double expResult20 = -0.3325038695132564;
        
        assertEquals(expResult2, e2, NUMACC);
        assertEquals(expResult20, e20, NUMACC);
    }

    /**
     * Test of gradient method, of class RigidBodyCoordinates.
     */
    @Test
    public void testGradient() {
        System.out.println("gradient");

        // make sure that w2 and w20 have the right energy and gradient
        final double[] extCoord2 = rigidW2TIP4P.getActiveCoordinates(gWater2.clone());
        final double[] extCoord2Dis = rigidW2DisTIP4P.getActiveCoordinates(gWater2Distorted.clone());
        final double[] extCoord20 = rigidW20TIP4P.getActiveCoordinates(gWater20.clone());
        
        final double[] grad2W = new double[6*2];
        final double[] grad2WDis = new double[6*2];
        final double[] grad20W = new double[6*20];
        
        final double e2 = rigidW2TIP4P.gradient(extCoord2, grad2W, 0);
        final double e2Dis = rigidW2DisTIP4P.gradient(extCoord2Dis, grad2WDis, 0);
        final double e20 = rigidW20TIP4P.gradient(extCoord20, grad20W, 0);
        
        // from DJW
        final double expResult2 = -0.009936231756023471;
        final double expResult20 = -0.3325038695132564;
        // from Phenix
        final double expResult2Dis = -0.002495066061667801;
        
        assertEquals(expResult2, e2, NUMACC);
        assertEquals(expResult2Dis, e2Dis, NUMACC);
        assertEquals(expResult20, e20, NUMACC);
        
        // check the gradients vs Bernd's Phenix implementation of TIP4P
        // and yes, since these are minimum geometries, the gradient elements SHOULD BE very small
        
        for(int i = 0; i < 12; i++){
            assertEquals("Failure in gradient element " + i + " reference element "
                    + gradW2Phenix[i] + " ours " + grad2W[i] + " diff " + Math.abs(gradW2Phenix[i]-grad2W[i])
                    ,gradW2Phenix[i],grad2W[i],NUMACCLOOSE);
        }
        
        // distorted gradient
        // only check the translational elements as the rotational ones are different
        for(int i = 0; i < 3; i++){
            assertEquals("Failure in gradient element " + i + " reference element "
                    + gradW2DisPhenix[i] + " ours " + grad2WDis[i] + " diff " + Math.abs(gradW2DisPhenix[i]-grad2WDis[i]),
                    gradW2DisPhenix[i],grad2WDis[i],NUMACC);
        }
        for(int i = 6; i < 9; i++){
            assertEquals("Failure in gradient element " + i + " reference element "
                    + gradW2DisPhenix[i] + " ours " + grad2WDis[i] + " diff " + Math.abs(gradW2DisPhenix[i]-grad2WDis[i]),
                    gradW2DisPhenix[i],grad2WDis[i],NUMACC);
        }
        
        
        // don't bother to put in individual elements...
        for(int i = 0; i < 120; i++){
            assertEquals("Failure in gradient element " + i + " reference element "
                    + 0.0 + " ours " + grad20W[i] + " diff " + Math.abs(grad20W[i]),
                    0.0,grad20W[i],NUMACC);
        }
    }
    
    @Test
    public void testNumericalGradient() {
        
        System.out.println("numerical gradients");

        // make sure that w2 and w20 have the right energy and gradient
        final double[] extCoord2 = rigidW2TIP4P.getActiveCoordinates(gWater2.clone());
        final double[] extCoord2Dis = rigidW2DisTIP4P.getActiveCoordinates(gWater2Distorted.clone());
        final double[] extCoord20 = rigidW20TIP4P.getActiveCoordinates(gWater20.clone());
        
        final double[] grad2W = new double[6*2];
        final double[] grad2WDis = new double[6*2];
        final double[] grad20W = new double[6*20];
        final double[] grad2WAna = new double[6*2];
        final double[] grad2WDisAna = new double[6*2];
        final double[] grad20WAna = new double[6*20];
        
        // this is kinda weird, so let me explain: we need to calculate the gradient analytically to make sure
        // that the internal caches are properly filled. yes, this is a test that relies on a specific implementation
        // behavior but that is life
        
        rigidW2TIP4P.gradient(extCoord2, grad2WAna, 0);
        final double e2 = NumericalGradients.gradientForRigidBody(gWater2,
                rigidW2TIP4P.getCartesBackupCopy(),
                RigidBodyCoordinates.getCOMs(extCoord2, 2),
                RigidBodyCoordinates.getEulers(extCoord2, 2),
                grad2W, tip4p, 0);
        
        rigidW2DisTIP4P.gradient(extCoord2Dis, grad2WDisAna, 0);
        final double e2Dis = NumericalGradients.gradientForRigidBody(gWater2Distorted,
                rigidW2DisTIP4P.getCartesBackupCopy(),
                RigidBodyCoordinates.getCOMs(extCoord2Dis, 2),
                RigidBodyCoordinates.getEulers(extCoord2Dis, 2),
                grad2WDis, tip4p, 0);
        
        rigidW20TIP4P.gradient(extCoord20, grad20WAna, 0);
        final double e20 = NumericalGradients.gradientForRigidBody(gWater20,
                rigidW20TIP4P.getCartesBackupCopy(),
                RigidBodyCoordinates.getCOMs(extCoord20, 20),
                RigidBodyCoordinates.getEulers(extCoord20, 20),
                grad20W, tip4p, 0);
        
        // from DJW
        final double expResult2 = -0.009936231756023471;
        final double expResult20 = -0.3325038695132564;
        // from Phenix
        final double expResult2Dis = -0.002495066061667801;
        
        assertEquals(expResult2, e2, NUMACC);
        assertEquals(expResult2Dis, e2Dis, NUMACC);
        assertEquals(expResult20, e20, NUMACC);
        
        // check numerical vs analytical gradient
        for(int i = 0; i < 12; i++){
            assertEquals("Failure in numerical gradient element " + i + " analytical element "
                    + grad2WAna[i] + " ours " + grad2W[i] + " diff " + Math.abs(grad2WAna[i]-grad2W[i]),
                    grad2WAna[i],grad2W[i],NUMACCLOOSE);
        }
        for(int i = 0; i < 12; i++){
            assertEquals("Failure in numerical gradient element " + i + " analytical element "
                    + grad2WDisAna[i] + " ours " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna[i]-grad2WDis[i]),
                    grad2WDisAna[i],grad2WDis[i],NUMACC);
        }
        for(int i = 0; i < 120; i++){
            assertEquals("Failure in numerical gradient element " + i + " analytical element "
                    + grad20WAna[i] + " ours " + grad20W[i] + " diff " + Math.abs(grad20WAna[i]-grad20W[i]),
                    grad20WAna[i],grad20W[i],NUMACC);
        }
        
        
        // check the gradients vs Bernd's Phenix implementation of TIP4P
        // and yes, since these are minimum geometries, the gradient elements SHOULD BE very small
        
        for(int i = 0; i < 12; i++){
            assertEquals("Failure in gradient element " + i + " reference element "
                    + gradW2Phenix[i] + " ours " + grad2W[i] + " diff " + Math.abs(gradW2Phenix[i]-grad2W[i])
                    ,gradW2Phenix[i],grad2W[i],NUMACCLOOSE);
        }
                
        // only check the translational elements as the rotational ones are different
        for(int i = 0; i < 3; i++){
            assertEquals("Failure in gradient element " + i + " reference element "
                    + gradW2DisPhenix[i] + " ours " + grad2WDis[i] + " diff " + Math.abs(gradW2DisPhenix[i]-grad2WDis[i]),
                    gradW2DisPhenix[i],grad2WDis[i],NUMACC);
        }
        for(int i = 6; i < 9; i++){
            assertEquals("Failure in gradient element " + i + " reference element "
                    + gradW2DisPhenix[i] + " ours " + grad2WDis[i] + " diff " + Math.abs(gradW2DisPhenix[i]-grad2WDis[i]),
                    gradW2DisPhenix[i],grad2WDis[i],NUMACC);
        }
        
        
        // don't bother to put in individual elements...
        for(int i = 0; i < 120; i++){
            assertEquals("Failure in gradient element " + i + " reference element "
                    + 0.0 + " ours " + grad20W[i] + " diff " + Math.abs(grad20W[i]),
                    0.0,grad20W[i],NUMACC);
        }
    }
    
    @Test
    public void testNumericalGradientStepped() {
        
        System.out.println("numerical gradients (stepped)");

        // make sure that w2 and w20 have the right energy and gradient
        final double[] extCoord2Dis = rigidW2DisTIP4P.getActiveCoordinates(gWater2Distorted.clone());
        
        final double[] grad2WDis = new double[6*2];
        final double[] grad2WDisAna = new double[6*2];
        final double[] grad2WDisAna2 = new double[6*2];
                        
        rigidW2DisTIP4P.gradient(extCoord2Dis, grad2WDisAna, 0);
        // do it a second time, to be sure this works identically
        rigidW2DisTIP4P.gradient(extCoord2Dis, grad2WDisAna2, 0);
        
        for(int i = 0; i < 12; i++){
            assertEquals("(0) Failure in second invocation gradient element " + i + " analytical element (new) "
                    + grad2WDisAna2[i] + " analytical element (old) " + grad2WDisAna[i] + " diff " + Math.abs(grad2WDisAna2[i]-grad2WDisAna2[i]),
                    grad2WDisAna2[i],grad2WDisAna[i],NUMACC);
        }
        
        final int mols = 2;
        NumericalGradients.gradientForRigidBody(gWater2Distorted, rigidW2DisTIP4P.getCartesBackupCopy(), RigidBodyCoordinates.getCOMs(extCoord2Dis, mols), RigidBodyCoordinates.getEulers(extCoord2Dis, mols), grad2WDis, tip4p, 0);
        
        for(int i = 0; i < 12; i++){
            assertEquals("(I) Failure in numerical gradient element " + i + " analytical element "
                    + grad2WDisAna[i] + " ours " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna[i]-grad2WDis[i]),
                    grad2WDisAna[i],grad2WDis[i],NUMACC);
        }
        
        // do it a second time, to be sure this works identically
        rigidW2DisTIP4P.gradient(extCoord2Dis, grad2WDisAna2, 0);
        
        for(int i = 0; i < 12; i++){
            assertEquals("(II) Failure in second invocation gradient element " + i + " analytical element (new) "
                    + grad2WDisAna2[i] + " analytical element (old) " + grad2WDisAna[i] + " diff " + Math.abs(grad2WDisAna2[i]-grad2WDisAna2[i]),
                    grad2WDisAna2[i],grad2WDisAna[i],NUMACC);
        }
        
        NumericalGradients.gradientForRigidBody(gWater2Distorted, rigidW2DisTIP4P.getCartesBackupCopy(), RigidBodyCoordinates.getCOMs(extCoord2Dis, mols), RigidBodyCoordinates.getEulers(extCoord2Dis, mols), grad2WDis, tip4p, 0);
        
        for(int i = 0; i < 12; i++){
            assertEquals("(III) Failure in numerical gradient element " + i + " analytical element "
                    + grad2WDisAna2[i] + " ours " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna2[i]-grad2WDis[i]),
                    grad2WDisAna2[i],grad2WDis[i],NUMACC);
        }
        
        // step through the COM shifts of the second water, because we can!
        final double[] extCoord2DisBack = extCoord2Dis.clone();
        for(int ix = 0; ix < 20; ix++){
            extCoord2Dis[6] = ix*0.01+extCoord2DisBack[6];
            for(int iy = 0; iy < 20; iy++){
                extCoord2Dis[7] = iy*0.01+extCoord2DisBack[7];
                for(int iz = 0; iz < 20; iz++){
                    extCoord2Dis[8] = iz*0.01+extCoord2DisBack[8];
                    
                    final double eAna = rigidW2DisTIP4P.gradient(extCoord2Dis, grad2WDisAna, 0);
                    final double eNum = NumericalGradients.gradientForRigidBody(gWater2Distorted, rigidW2DisTIP4P.getCartesBackupCopy(), RigidBodyCoordinates.getCOMs(extCoord2Dis, mols), RigidBodyCoordinates.getEulers(extCoord2Dis, mols), grad2WDis, tip4p, 0);
        
                    // compare
                    boolean allFine = true;
                    if(Math.abs(eAna-eNum) > NUMACCLOOSE){
                            System.err.println("TRANSLATION: Failure in numerical vs analytical energy, analytical energy "
                            + eAna + " numerical " + eNum + " diff " + Math.abs(eAna-eNum));
                            allFine = false;
                    }
                    for(int i = 0; i < 12; i++){
                        if(Math.abs(grad2WDisAna[i]-grad2WDis[i]) > NUMACCLOOSE){
                            System.err.println("TRANSLATION: Failure in numerical vs analytical gradient element " + i + " analytical element "
                            + grad2WDisAna[i] + " numerical " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna[i]-grad2WDis[i]));
                            allFine = false;
                        }
                        /*assertEquals("Failure in numerical gradient element " + i + " analytical element "
                            + grad2WDisAna[i] + " ours " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna[i]-grad2WDis[i]),
                        grad2WDisAna[i],grad2WDis[i],NUMACCLOOSE);*/
                    }
                    assertTrue("TRANSLATION: Failure in one or more gradient elements, see above for details!",allFine);
                }
            }
        }
        
        
        // here, the Euler angles are all zero. Now, we want to see what happens with analytical vs numerical gradient when
        // we step through the allowed range of Euler angles one by one
        
        final double phiIncr = 2*Math.PI/20;
        final double omegaIncr = Math.PI/20;
        final double psiIncr = 2*Math.PI/20;
        
        final double phiStart = -Math.PI + 2*phiIncr;
        final double omegaStart = -0.5*Math.PI + 2*omegaIncr;
        final double psiStart = -Math.PI + 2*psiIncr;
        
        // as to not run into problems with the sanitizing (some gradient discontinuities may be expected if you go beyond the definition interval of the Eulers)
        // we start a litte ABOVE the lower bound and stop a little BELOW the upper one.
        
        double phi = phiStart;
        double omega = omegaStart;
        double psi = psiStart;
                
        final double phiEnd = Math.PI-2*phiIncr;
        final double omegaEnd = Math.PI/2-2*omegaIncr;
        final double psiEnd = Math.PI-2*psiIncr;
        
        while(phi <= phiEnd){
            while(omega <= omegaEnd){
                while(psi <= psiEnd){
                    
                    // put it in (first water)
                    extCoord2Dis[3] = phi;
                    extCoord2Dis[4] = omega;
                    extCoord2Dis[5] = psi;
                                        
                    // compute analytical and numerical
                    final double eAna = rigidW2DisTIP4P.gradient(extCoord2Dis, grad2WDisAna, 0);
                    final double eNum = NumericalGradients.gradientForRigidBody(gWater2Distorted, rigidW2DisTIP4P.getCartesBackupCopy(), RigidBodyCoordinates.getCOMs(extCoord2Dis, mols), RigidBodyCoordinates.getEulers(extCoord2Dis, mols), grad2WDis, tip4p, 0);
                    
                    // compare
                    boolean allFine = true;
                    if(Math.abs(eAna-eNum) > NUMACCLOOSE){
                            System.err.println("ROTATION: Failure in numerical vs analytical energy, analytical energy "
                            + eAna + " numerical " + eNum + " diff " + Math.abs(eAna-eNum) + " phi " + phi + " omega " + omega + " psi " + psi);
                            allFine = false;
                    }
                    for(int i = 0; i < 12; i++){
                        if(Math.abs(grad2WDisAna[i]-grad2WDis[i]) > NUMACCLOOSE){
                            System.err.println("ROTATION: Failure in numerical vs analytical gradient element " + i + " analytical element "
                            + grad2WDisAna[i] + " numerical " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna[i]-grad2WDis[i]) + " phi " + phi + " omega " + omega + " psi " + psi);
                            allFine = false;
                        }
                    }
                    assertTrue("ROTATION: Failure in one or more gradient elements, see above for details!",allFine);
                    
                    if(false){
                        System.out.println("One done!");
                    }
                    
                    psi += psiIncr;
                }
                psi = psiStart;
                
                omega += omegaIncr;
            }
            omega = omegaStart;
            
            phi += phiIncr;
        }
        
        // reset
        phi = phiStart;
        omega = omegaStart;
        psi = psiStart;
        
        while(phi <= phiEnd){
            while(omega <= omegaEnd){
                while(psi <= psiEnd){
                    
                    // put it in (second water)
                    extCoord2Dis[9] = phi;
                    extCoord2Dis[10] = omega;
                    extCoord2Dis[11] = psi;
                    
                    // compute analytical and numerical
                    final double eAna = rigidW2DisTIP4P.gradient(extCoord2Dis, grad2WDisAna, 0);
                    final double eNum = NumericalGradients.gradientForRigidBody(gWater2Distorted, rigidW2DisTIP4P.getCartesBackupCopy(), RigidBodyCoordinates.getCOMs(extCoord2Dis, mols), RigidBodyCoordinates.getEulers(extCoord2Dis, mols), grad2WDis, tip4p, 0);
                    
                    // compare
                    boolean allFine = true;
                    if(Math.abs(eAna-eNum) > NUMACCLOOSE){
                            System.err.println("ROTATION: Failure in numerical vs analytical energy, analytical energy "
                            + eAna + " numerical " + eNum + " diff " + Math.abs(eAna-eNum) + " phi " + phi + " omega " + omega + " psi " + psi);
                            allFine = false;
                    }
                    for(int i = 0; i < 12; i++){
                        if(Math.abs(grad2WDisAna[i]-grad2WDis[i]) > NUMACCLOOSE){
                            System.err.println("ROTATION: Failure in numerical vs analytical gradient element " + i + " analytical element "
                            + grad2WDisAna[i] + " numerical " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna[i]-grad2WDis[i]) + " phi " + phi + " omega " + omega + " psi " + psi);
                            allFine = false;
                        }
                    }
                    assertTrue("ROTATION: Failure in one or more gradient elements, see above for details!",allFine);
                    
                    if(false){
                        System.out.println("One done!");
                    }
                    
                    psi += psiIncr;
                }
                psi = psiStart;
                
                omega += omegaIncr;
            }
            omega = omegaStart;
            
            phi += phiIncr;
        }
        
    }
    
    
    @Test
    public void testNumericalGradientTIP3PStepped() {
        
        // a stupid copy of the TIP4P test above. it's ok though, we do not currently have anything better than this
        System.out.println("numerical gradients TIP3P (stepped)");

        // make sure that w2 and w20 have the right energy and gradient
        final double[] extCoord2Dis = rigidW2DisTIP3P.getActiveCoordinates(gWater2Distorted.clone());
        
        final double[] grad2WDis = new double[6*2];
        final double[] grad2WDisAna = new double[6*2];
        final double[] grad2WDisAna2 = new double[6*2];
                        
        rigidW2DisTIP3P.gradient(extCoord2Dis, grad2WDisAna, 0);
        // do it a second time, to be sure this works identically
        rigidW2DisTIP3P.gradient(extCoord2Dis, grad2WDisAna2, 0);
        
        for(int i = 0; i < 12; i++){
            assertEquals("(0) Failure in second invocation gradient element " + i + " analytical element (new) "
                    + grad2WDisAna2[i] + " analytical element (old) " + grad2WDisAna[i] + " diff " + Math.abs(grad2WDisAna2[i]-grad2WDisAna2[i]),
                    grad2WDisAna2[i],grad2WDisAna[i],NUMACC);
        }
        
        final int mols = 2;
        NumericalGradients.gradientForRigidBody(gWater2Distorted, rigidW2DisTIP3P.getCartesBackupCopy(), RigidBodyCoordinates.getCOMs(extCoord2Dis, mols), RigidBodyCoordinates.getEulers(extCoord2Dis, mols), grad2WDis, tip3p, 0);
        
        for(int i = 0; i < 12; i++){
            assertEquals("(I) Failure in numerical gradient element " + i + " analytical element "
                    + grad2WDisAna[i] + " ours " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna[i]-grad2WDis[i]),
                    grad2WDisAna[i],grad2WDis[i],NUMACC);
        }
        
        // do it a second time, to be sure this works identically
        rigidW2DisTIP3P.gradient(extCoord2Dis, grad2WDisAna2, 0);
        
        for(int i = 0; i < 12; i++){
            assertEquals("(II) Failure in second invocation gradient element " + i + " analytical element (new) "
                    + grad2WDisAna2[i] + " analytical element (old) " + grad2WDisAna[i] + " diff " + Math.abs(grad2WDisAna2[i]-grad2WDisAna2[i]),
                    grad2WDisAna2[i],grad2WDisAna[i],NUMACC);
        }
        
        NumericalGradients.gradientForRigidBody(gWater2Distorted, rigidW2DisTIP3P.getCartesBackupCopy(), RigidBodyCoordinates.getCOMs(extCoord2Dis, mols), RigidBodyCoordinates.getEulers(extCoord2Dis, mols), grad2WDis, tip3p, 0);
        
        for(int i = 0; i < 12; i++){
            assertEquals("(III) Failure in numerical gradient element " + i + " analytical element "
                    + grad2WDisAna2[i] + " ours " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna2[i]-grad2WDis[i]),
                    grad2WDisAna2[i],grad2WDis[i],NUMACC);
        }
        
        // step through the COM shifts of the second water, because we can!
        final double[] extCoord2DisBack = extCoord2Dis.clone();
        for(int ix = 0; ix < 20; ix++){
            extCoord2Dis[6] = ix*0.01+extCoord2DisBack[6];
            for(int iy = 0; iy < 20; iy++){
                extCoord2Dis[7] = iy*0.01+extCoord2DisBack[7];
                for(int iz = 0; iz < 20; iz++){
                    extCoord2Dis[8] = iz*0.01+extCoord2DisBack[8];
                    
                    final double eAna = rigidW2DisTIP3P.gradient(extCoord2Dis, grad2WDisAna, 0);
                    final double eNum = NumericalGradients.gradientForRigidBody(gWater2Distorted, rigidW2DisTIP3P.getCartesBackupCopy(), RigidBodyCoordinates.getCOMs(extCoord2Dis, mols), RigidBodyCoordinates.getEulers(extCoord2Dis, mols), grad2WDis, tip3p, 0);
        
                    // compare
                    boolean allFine = true;
                    if(Math.abs(eAna-eNum) > NUMACCLOOSE){
                            System.err.println("TRANSLATION: Failure in numerical vs analytical energy, analytical energy "
                            + eAna + " numerical " + eNum + " diff " + Math.abs(eAna-eNum));
                            allFine = false;
                    }
                    for(int i = 0; i < 12; i++){
                        if(Math.abs(grad2WDisAna[i]-grad2WDis[i]) > NUMACCLOOSE){
                            System.err.println("TRANSLATION: Failure in numerical vs analytical gradient element " + i + " analytical element "
                            + grad2WDisAna[i] + " numerical " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna[i]-grad2WDis[i]));
                            allFine = false;
                        }
                        /*assertEquals("Failure in numerical gradient element " + i + " analytical element "
                            + grad2WDisAna[i] + " ours " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna[i]-grad2WDis[i]),
                        grad2WDisAna[i],grad2WDis[i],NUMACCLOOSE);*/
                    }
                    assertTrue("TRANSLATION: Failure in one or more gradient elements, see above for details!",allFine);
                }
            }
        }
        
        
        // here, the Euler angles are all zero. Now, we want to see what happens with analytical vs numerical gradient when
        // we step through the allowed range of Euler angles one by one
        
        final double phiIncr = 2*Math.PI/20;
        final double omegaIncr = Math.PI/20;
        final double psiIncr = 2*Math.PI/20;
        
        final double phiStart = -Math.PI + 2*phiIncr;
        final double omegaStart = -0.5*Math.PI + 2*omegaIncr;
        final double psiStart = -Math.PI + 2*psiIncr;
        
        // as to not run into problems with the sanitizing (some gradient discontinuities may be expected if you go beyond the definition interval of the Eulers)
        // we start a litte ABOVE the lower bound and stop a little BELOW the upper one.
        
        double phi = phiStart;
        double omega = omegaStart;
        double psi = psiStart;
                
        final double phiEnd = Math.PI-2*phiIncr;
        final double omegaEnd = Math.PI/2-2*omegaIncr;
        final double psiEnd = Math.PI-2*psiIncr;
        
        while(phi <= phiEnd){
            while(omega <= omegaEnd){
                while(psi <= psiEnd){
                    
                    // put it in (first water)
                    extCoord2Dis[3] = phi;
                    extCoord2Dis[4] = omega;
                    extCoord2Dis[5] = psi;
                                        
                    // compute analytical and numerical
                    final double eAna = rigidW2DisTIP3P.gradient(extCoord2Dis, grad2WDisAna, 0);
                    final double eNum = NumericalGradients.gradientForRigidBody(gWater2Distorted, rigidW2DisTIP3P.getCartesBackupCopy(), RigidBodyCoordinates.getCOMs(extCoord2Dis, mols), RigidBodyCoordinates.getEulers(extCoord2Dis, mols), grad2WDis, tip3p, 0);
                    
                    // compare
                    boolean allFine = true;
                    if(Math.abs(eAna-eNum) > NUMACCLOOSE){
                            System.err.println("ROTATION: Failure in numerical vs analytical energy, analytical energy "
                            + eAna + " numerical " + eNum + " diff " + Math.abs(eAna-eNum) + " phi " + phi + " omega " + omega + " psi " + psi);
                            allFine = false;
                    }
                    for(int i = 0; i < 12; i++){
                        if(Math.abs(grad2WDisAna[i]-grad2WDis[i]) > NUMACCLOOSE){
                            System.err.println("ROTATION: Failure in numerical vs analytical gradient element " + i + " analytical element "
                            + grad2WDisAna[i] + " numerical " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna[i]-grad2WDis[i]) + " phi " + phi + " omega " + omega + " psi " + psi);
                            allFine = false;
                        }
                    }
                    assertTrue("ROTATION: Failure in one or more gradient elements, see above for details!",allFine);
                    
                    if(false){
                        System.out.println("One done!");
                    }
                    
                    psi += psiIncr;
                }
                psi = psiStart;
                
                omega += omegaIncr;
            }
            omega = omegaStart;
            
            phi += phiIncr;
        }
        
        // reset
        phi = phiStart;
        omega = omegaStart;
        psi = psiStart;
        
        while(phi <= phiEnd){
            while(omega <= omegaEnd){
                while(psi <= psiEnd){
                    
                    // put it in (second water)
                    extCoord2Dis[9] = phi;
                    extCoord2Dis[10] = omega;
                    extCoord2Dis[11] = psi;
                    
                    // compute analytical and numerical
                    final double eAna = rigidW2DisTIP3P.gradient(extCoord2Dis, grad2WDisAna, 0);
                    final double eNum = NumericalGradients.gradientForRigidBody(gWater2Distorted, rigidW2DisTIP3P.getCartesBackupCopy(), RigidBodyCoordinates.getCOMs(extCoord2Dis, mols), RigidBodyCoordinates.getEulers(extCoord2Dis, mols), grad2WDis, tip3p, 0);
                    
                    // compare
                    boolean allFine = true;
                    if(Math.abs(eAna-eNum) > NUMACCLOOSE){
                            System.err.println("ROTATION: Failure in numerical vs analytical energy, analytical energy "
                            + eAna + " numerical " + eNum + " diff " + Math.abs(eAna-eNum) + " phi " + phi + " omega " + omega + " psi " + psi);
                            allFine = false;
                    }
                    for(int i = 0; i < 12; i++){
                        if(Math.abs(grad2WDisAna[i]-grad2WDis[i]) > NUMACCLOOSE){
                            System.err.println("ROTATION: Failure in numerical vs analytical gradient element " + i + " analytical element "
                            + grad2WDisAna[i] + " numerical " + grad2WDis[i] + " diff " + Math.abs(grad2WDisAna[i]-grad2WDis[i]) + " phi " + phi + " omega " + omega + " psi " + psi);
                            allFine = false;
                        }
                    }
                    assertTrue("ROTATION: Failure in one or more gradient elements, see above for details!",allFine);
                    
                    if(false){
                        System.out.println("One done!");
                    }
                    
                    psi += psiIncr;
                }
                psi = psiStart;
                
                omega += omegaIncr;
            }
            omega = omegaStart;
            
            phi += phiIncr;
        }
        
    }
    
    private static Geometry getWater2TIP4P(){
        
        // plug in the known global minimum structure for water 2
        final CartesianCoordinates c = new CartesianCoordinates(6,2,new int[]{3,3});
        final double[][] xyz = c.getAllXYZCoord();
        xyz[0][0] = -0.1858140*ANGTOBOHR;
        xyz[1][0] = -1.1749469*ANGTOBOHR;
        xyz[2][0] = 0.7662596*ANGTOBOHR;
        xyz[0][1] = -0.1285513*ANGTOBOHR;
        xyz[1][1] = -0.8984365*ANGTOBOHR;
        xyz[2][1] = 1.6808606*ANGTOBOHR;
        xyz[0][2] = -0.0582782*ANGTOBOHR;
        xyz[1][2] = -0.3702550*ANGTOBOHR;
        xyz[2][2] = 0.2638279*ANGTOBOHR;
        xyz[0][3] = 0.1747051*ANGTOBOHR;
        xyz[1][3] = 1.1050002*ANGTOBOHR;
        xyz[2][3] = -0.7244430*ANGTOBOHR;
        xyz[0][4] = -0.5650842*ANGTOBOHR;
        xyz[1][4] = 1.3134964*ANGTOBOHR;
        xyz[2][4] = -1.2949455*ANGTOBOHR;
        xyz[0][5] = 0.9282185*ANGTOBOHR;
        xyz[1][5] = 1.0652990*ANGTOBOHR;
        xyz[2][5] = -1.313402*ANGTOBOHR;
        
        final String[] atoms = c.getAllAtomTypes();
        atoms[0] = "O";
        atoms[1] = "H";
        atoms[2] = "H";
        atoms[3] = "O";
        atoms[4] = "H";
        atoms[5] = "H";
        
        c.recalcAtomNumbers();
        
        final int noOfMols = 2;
        final int noOfAtoms = noOfMols*3;
        final int[] atsPerMol = new int[noOfMols];
        final String[] sids = new String[noOfMols];
        final boolean[] molFlexies = new boolean[noOfMols];
        final boolean[] molConstr = new boolean[noOfMols];
        final boolean[][] constrXYZ = new boolean[3][noOfAtoms];
        final BondInfo bonds = new SimpleBondInfo(noOfAtoms);
        final List<boolean[][]> allFlex = new ArrayList<>();
        for(int i = 0; i < noOfMols; i++){
            atsPerMol[i] = 3;
            sids[i] = "OHH";
            molFlexies[i] = false;
            molConstr[i] = false;
            final int o = i*3;
            final int h1 = i*3+1;
            final int h2 = i*3+2;
            bonds.setBond(o, h1, BondInfo.SINGLE);
            bonds.setBond(o, h2, BondInfo.SINGLE);
            allFlex.add(null);
        }
        for(int i = 0; i < noOfAtoms; i++){
            constrXYZ[0][i] = false;
            constrXYZ[1][i] = false;
            constrXYZ[2][i] = false;
        }
        
        final Geometry g = new Geometry(c, 0, noOfMols, atsPerMol,
            molFlexies, allFlex, molConstr, constrXYZ,
            sids, bonds);
            
        return g;
    }
    
    private static Geometry getWater2TIP4PDistorted(){
        
        final short[] spins = new short[6];
        final float[] charges = new float[6];
        final int[] atsPerMol = new int[]{3,3};
        
        final String[] w2_distorted = new String[]{
            "6",
            "Energy:              -6.5508",
            "O             -0.0053998              0.0223328              0.1261801",
            "H             -0.9200008             -0.0273239             -0.1517969",
            "H              0.4970319             -0.1584812             -0.6682384",
            "O              2.1931798              2.0700218             -2.1761271",
            "H              3.0481182              1.8086380             -1.8340878",
            "H              1.9936879              2.8892588             -1.7230422"
        };
        
        
        CartesianCoordinates w2 = null;
        try{
            w2 = Input.parseCartesFromFileData(w2_distorted, 2, atsPerMol, spins, charges);
        } catch(Exception e){
            throw new Error("Error to parse in known w2 distorted geometry.",e);
        }
        
        final int noOfMols = 2;
        final int noOfAtoms = noOfMols*3;
        final String[] sids = new String[noOfMols];
        final boolean[] molFlexies = new boolean[noOfMols];
        final boolean[] molConstr = new boolean[noOfMols];
        final boolean[][] constrXYZ = new boolean[3][noOfAtoms];
        final BondInfo bonds = new SimpleBondInfo(noOfAtoms);
        final List<boolean[][]> allFlex = new ArrayList<>();
        for(int i = 0; i < noOfMols; i++){
            atsPerMol[i] = 3;
            sids[i] = "OHH";
            molFlexies[i] = false;
            molConstr[i] = false;
            final int o = i*3;
            final int h1 = i*3+1;
            final int h2 = i*3+2;
            bonds.setBond(o, h1, BondInfo.SINGLE);
            bonds.setBond(o, h2, BondInfo.SINGLE);
            allFlex.add(null);
        }
        for(int i = 0; i < noOfAtoms; i++){
            constrXYZ[0][i] = false;
            constrXYZ[1][i] = false;
            constrXYZ[2][i] = false;
        }
        
        final Geometry g = new Geometry(w2, 0, noOfMols, atsPerMol,
            molFlexies, allFlex, molConstr, constrXYZ,
            sids, bonds);
            
        return g;
    }
    
    private static Geometry getWater20TIP4P(){
        final short[] spins = new short[60];
        final float[] charges = new float[60];
        final int[] atsPerMol = new int[]{3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
        
        final String[] w20_globMin = new String[]{
            "60",
            "25304   Energy =     -872.9887544952  kJ/mol",
            "O    0.3470450  -0.7178074  -0.3141646",
            "H   -0.4987771  -0.4023865  -0.6324831",
            "H    0.4310626  -1.5958882  -0.6858463",
            "O    1.6727305   1.3910313  -1.4036375",
            "H    1.3457584   0.5950317  -0.9844648",
            "H    2.0394098   1.0910128  -2.2353629",
            "O    2.3018465   0.3812069  -3.8845400",
            "H    3.1099899   0.3801893  -4.3975071",
            "H    2.0321756  -0.5368862  -3.8597484",
            "O    0.1025673  -3.1189675  -1.5785035",
            "H   -0.7300398  -2.8548031  -1.9699297",
            "H   -0.1423526  -3.7271640  -0.8811203",
            "O   -2.8168905   0.4280904   1.3812494",
            "H   -3.4076222  -0.3195473   1.4723865",
            "H   -1.9935236   0.1320979   1.7694340",
            "O   -4.0874498  -2.0205216   1.5782400",
            "H   -3.9284136  -2.5432513   0.7923054",
            "H   -4.9531663  -2.2967179   1.8790432",
            "O   -3.4229815  -3.5842142  -0.5206354",
            "H   -2.6411841  -4.0554473  -0.2326002",
            "H   -3.1633283  -3.1685169  -1.3428324",
            "O   -0.7506626   2.3661529  -2.1609725",
            "H   -0.9151074   3.1344968  -1.6143121",
            "H    0.1319230   2.0846668  -1.9200570",
            "O   -2.0420422   0.1630925  -1.2136159",
            "H   -1.7487716   1.0099069  -1.5499635",
            "H   -2.4304310   0.3644856  -0.3622482",
            "O   -0.2800640   1.1040261  -4.6117443",
            "H    0.6548920   0.9730038  -4.4538758",
            "H   -0.5698339   1.6605327  -3.8888571",
            "O    0.4395343   2.1254982   2.8973441",
            "H    0.8751504   2.4986854   2.1310530",
            "H    0.9424998   2.4534973   3.6427798",
            "O    1.4361964   3.2411346   0.6566335",
            "H    0.7004338   3.7381475   0.2990532",
            "H    1.6959489   2.6531787  -0.0526398",
            "O   -1.0265126  -4.5692310   0.4896893",
            "H   -0.8653801  -5.4371521   0.8597992",
            "H   -1.1597562  -4.0046194   1.2510643",
            "O    1.3071138  -2.1277620  -3.8956797",
            "H    0.4836151  -2.0150491  -4.3704194",
            "H    1.0632203  -2.6034941  -3.1016854",
            "O   -1.1385446  -1.4258139  -5.0114174",
            "H   -0.9184261  -0.5007828  -4.9014301",
            "H   -1.4524214  -1.4922285  -5.9132503",
            "O   -1.9905990   3.0299735   1.8988378",
            "H   -1.2100748   2.7681006   2.3871375",
            "H   -2.4725335   2.2144503   1.7613859",
            "O   -0.9695361   4.3238367  -0.2345309",
            "H   -1.2732431   5.2313581  -0.2145549",
            "H   -1.4117555   3.9064496   0.5046998",
            "O   -2.2412072  -2.1731324  -2.6014442",
            "H   -2.2859870  -1.3195776  -2.1705486",
            "H   -1.9941113  -1.9734685  -3.5043893",
            "O   -1.6276977  -2.9066441   2.5290618",
            "H   -2.5203412  -2.6143057   2.3447852",
            "H   -1.1228631  -2.0977796   2.6133951",
            "O   -0.3997161  -0.4285041   2.2859951",
            "H    0.0022440  -0.5697997   1.4288512",
            "H   -0.0190375   0.3927811   2.5971324"
        };
        
        CartesianCoordinates w20 = null;
        try{
            w20 = Input.parseCartesFromFileData(w20_globMin, 20, atsPerMol, spins, charges);
        } catch(Exception e){
            throw new Error("Error to parse in known w20 global minimum.",e);
        }
        
        final int noOfMols = 20;
        final int noOfAtoms = noOfMols*3;
        final String[] sids = new String[noOfMols];
        final boolean[] molFlexies = new boolean[noOfMols];
        final boolean[] molConstr = new boolean[noOfMols];
        final boolean[][] constrXYZ = new boolean[3][noOfAtoms];
        final BondInfo bonds = new SimpleBondInfo(noOfAtoms);
        final List<boolean[][]> allFlex = new ArrayList<>();
        for(int i = 0; i < noOfMols; i++){
            atsPerMol[i] = 3;
            sids[i] = "OHH";
            molFlexies[i] = false;
            molConstr[i] = false;
            final int o = i*3;
            final int h1 = i*3+1;
            final int h2 = i*3+2;
            bonds.setBond(o, h1, BondInfo.SINGLE);
            bonds.setBond(o, h2, BondInfo.SINGLE);
            allFlex.add(null);
        }
        for(int i = 0; i < noOfAtoms; i++){
            constrXYZ[0][i] = false;
            constrXYZ[1][i] = false;
            constrXYZ[2][i] = false;
        }
        
        final Geometry g = new Geometry(w20, 0, noOfMols, atsPerMol,
            molFlexies, allFlex, molConstr, constrXYZ,
            sids, bonds);
            
        return g;
    }
    
    private static Geometry getNoWaterTIP4P(){
        
        final short[] spins = new short[1];
        final float[] charges = new float[1];
        final int[] atsPerMol = new int[]{1};
        
        final String[] lj = new String[]{
            "1",
            "",
            "Ar 0.0 0.0 0.0"
        };
        
        CartesianCoordinates ar = null;
        try{
            ar = Input.parseCartesFromFileData(lj, 1, atsPerMol, spins, charges);
        } catch(Exception e){
            throw new Error("Error to parse in Ar single atom.",e);
        }
        
        final int noOfMols = 1;
        final int noOfAtoms = noOfMols*1;
        final String[] sids = new String[noOfMols];
        final boolean[] molFlexies = new boolean[noOfMols];
        final boolean[] molConstr = new boolean[noOfMols];
        final boolean[][] constrXYZ = new boolean[3][noOfAtoms];
        final BondInfo bonds = new SimpleBondInfo(noOfAtoms);
        final List<boolean[][]> allFlex = new ArrayList<>();
        for(int i = 0; i < noOfMols; i++){
            atsPerMol[i] = 1;
            sids[i] = "Ar";
            molFlexies[i] = false;
            molConstr[i] = false;
            allFlex.add(null);
        }
        for(int i = 0; i < noOfAtoms; i++){
            constrXYZ[0][i] = false;
            constrXYZ[1][i] = false;
            constrXYZ[2][i] = false;
        }
        
        final Geometry g = new Geometry(ar, 0, noOfMols, atsPerMol,
            molFlexies, allFlex, molConstr, constrXYZ,
            sids, bonds);
            
        return g;
    }
}
