/**
Copyright (c) 2020, J. M. Dieterich and B. Hartke
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

import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Johannes Dieterich
 * @version 2020-02-01
 */
public class LennardJonesFFTest {
    
    private static final String LJ38MIN = "38\n" +
"Energy:            -173.2453\n" +
"Ar              2.2635080             -1.3021673              0.2364510\n" +
"Ar             -4.0696589             -1.7833951             -1.0912553\n" +
"Ar              3.3465963              0.5473776              5.1173994\n" +
"Ar              4.6352036              2.9641120             -1.9064564\n" +
"Ar             -5.2396084              1.6412581             -2.1458985\n" +
"Ar              2.1975346              3.9699958              4.0662786\n" +
"Ar              2.9139336             -0.1048691             -3.2585750\n" +
"Ar             -0.1537578             -0.2936708              6.2036580\n" +
"Ar             -3.5310673              4.6998827             -0.7800296\n" +
"Ar             -1.8706342             -3.3672049              4.8312587\n" +
"Ar              3.9770026              1.7439961              1.5989158\n" +
"Ar              5.1470505             -1.6806710              2.6532521\n" +
"Ar              0.0611707              0.2542515             -5.6961461\n" +
"Ar              1.1218741              2.0986964             -0.8081952\n" +
"Ar             -5.8825641              0.4360893              1.3575178\n" +
"Ar              1.7779228              3.3278948             -4.3237295\n" +
"Ar              3.4384866             -4.7393483              1.2876312\n" +
"Ar              4.0814241             -3.5339516             -2.2157154\n" +
"Ar              5.7898690             -0.4751665             -0.8500914\n" +
"Ar             -2.9362859             -5.2206833             -0.0378087\n" +
"Ar             -1.2145583             -2.1381033              1.3158395\n" +
"Ar             -4.1741020              3.4946992              2.7233876\n" +
"Ar              2.8435532              5.1812951              0.5458881\n" +
"Ar             -0.5757456             -0.9405450             -2.1653462\n" +
"Ar              1.6313510             -2.5096021              3.7313159\n" +
"Ar             -1.3028293              3.1289101              5.1523874\n" +
"Ar             -3.0066453              0.0655822              3.7661092\n" +
"Ar             -0.0136789              5.5451888             -1.8714567\n" +
"Ar              0.5682738             -4.3585518             -1.1260230\n" +
"Ar             -0.6609501              4.3190723              1.6338157\n" +
"Ar             -2.2902497             -4.0096773             -3.5583525\n" +
"Ar              0.4831067              0.9010895              2.6729642\n" +
"Ar             -3.4390998             -0.5870227             -4.6098257\n" +
"Ar             -2.3561517              1.2627716              0.2711277\n" +
"Ar             -1.7239135              2.4701423             -3.2237536\n" +
"Ar             -4.7280805             -3.0031875              2.4141292\n" +
"Ar              1.2100750             -3.1683730             -4.6446476\n" +
"Ar             -0.0788494             -5.5846809              2.3792819";
    
    private static final String LJ55MIN = "55\n" +
"Energy:            -278.1517\n" +
"Ar              0.2637670              0.0750218             -0.0305859\n" +
"Ar              1.9400511             -5.9359266             -0.2182075\n" +
"Ar             -3.3887782              5.5878993             -3.0857245\n" +
"Ar              5.4645827              3.5286603             -0.0164234\n" +
"Ar             -0.0611506             -6.3446547             -3.4581804\n" +
"Ar             -4.7623515              2.0139627             -3.1856848\n" +
"Ar             -2.7311677              1.0755806              5.3550730\n" +
"Ar              5.1542978             -3.8014196             -0.2095618\n" +
"Ar             -3.0861314             -3.9969872             -3.3732330\n" +
"Ar              2.4421399              5.9244352              0.0941955\n" +
"Ar              0.1530409              0.2714680             -7.3117752\n" +
"Ar              3.2585466             -0.9253576             -5.4163715\n" +
"Ar              0.2094058              0.1714695             -3.6056473\n" +
"Ar              2.2225382              2.7895936             -5.3003453\n" +
"Ar             -1.6999051             -2.4743290              1.5308355\n" +
"Ar              6.5888663              1.7386441              3.1775426\n" +
"Ar              0.4233255              3.2271118              1.6523928\n" +
"Ar             -1.6951682             -2.6393958              5.2391891\n" +
"Ar             -2.8418301             -0.7418962             -1.6057981\n" +
"Ar             -2.7097004              1.1466854              1.6437936\n" +
"Ar              0.3743391             -0.1213402              7.2505633\n" +
"Ar              2.0571514             -2.6317879              1.4694758\n" +
"Ar              3.3693659              0.8918850              1.5446257\n" +
"Ar              3.9163409             -5.4378330              3.0245237\n" +
"Ar             -1.4124661              6.0859683              0.1570572\n" +
"Ar              0.4832106              3.2099543              5.3639202\n" +
"Ar             -5.9730344              0.3363606              0.0711651\n" +
"Ar             -4.6268008              3.9515551              0.1482164\n" +
"Ar             -2.6231431              4.4083720              3.4138865\n" +
"Ar             -2.9781957             -0.6641526             -5.3145255\n" +
"Ar              4.2630175              5.2673226             -3.2106167\n" +
"Ar             -1.9145285             -5.7743746             -0.1552363\n" +
"Ar             -4.9370769             -3.3786155             -0.0447196\n" +
"Ar             -6.0612815             -1.5887397             -3.2387501\n" +
"Ar             -1.5296300              2.7818574             -1.5306798\n" +
"Ar             -3.7355675             -5.1171473              3.1495165\n" +
"Ar              0.4384309              5.4677007             -3.1714606\n" +
"Ar              5.3289109              1.5912421             -3.3504794\n" +
"Ar              0.0890783             -5.3175663              3.1103721\n" +
"Ar             -4.8014991             -1.4410447              3.2891758\n" +
"Ar             -1.6319388              2.9510461             -5.2374343\n" +
"Ar              5.2899358             -1.8640178              3.1244180\n" +
"Ar             -5.7921455              2.2576283              3.3795587\n" +
"Ar              3.2371665             -0.9965492             -1.7050608\n" +
"Ar              0.0442854             -3.0598980             -5.4250516\n" +
"Ar              3.6137059              4.1469046              3.3121847\n" +
"Ar              0.3180569             -0.0213815              3.5444697\n" +
"Ar              2.1594096             -2.8009554              5.1762327\n" +
"Ar              0.5887912              6.4946966              3.3970661\n" +
"Ar              0.1042463             -3.0770498             -1.7135232\n" +
"Ar              6.3195618             -2.1073732             -3.4409865\n" +
"Ar              2.2274111              2.6244358             -1.5919829\n" +
"Ar              3.5057247              0.8140948              5.2533412\n" +
"Ar              3.1506760             -4.2582346             -3.4751765\n" +
"Ar              6.5005124             -0.1862658             -0.1324787";
    
    private static final double ENERGYLJ38 = -0.0659856619926505;
    private static final double ENERGYLJ55 = -0.10594240151944964;
    private static final double NUMACC = 1.0e-8;
    
    @Test
    public void testEnergyCalculationLJMins() {
        
        final int[] atsPerMol38 = new int[38];
        for(int i = 0; i < 38; i++){
            atsPerMol38[i] = 1;
        }
        
        final int[] atsPerMol55 = new int[55];
        for(int i = 0; i < 55; i++){
            atsPerMol55[i] = 1;
        }
        
        CartesianCoordinates lj38 = null;
        CartesianCoordinates lj55 = null;
        try {
            final String[] fileCont38 = LJ38MIN.split("\n");
            lj38 = Input.parseCartesFromFileData(fileCont38, 38, atsPerMol38, new short[38], new float[38]);
        
            final String[] fileCont55 = LJ55MIN.split("\n");
            lj55 = Input.parseCartesFromFileData(fileCont55, 55, atsPerMol55, new short[55], new float[55]);
        } catch(Exception e){
           fail(e.toString());
        }

        final BondInfo bonds38 = new SimpleBondInfo(38);
        final BondInfo bonds55 = new SimpleBondInfo(55);
        
        final LennardJonesFF ljFF = new LennardJonesFF(true);
       
        final double e38 = ljFF.energyCalculation(-1, 0, lj38.getAll1DCartes(),
                lj38.getAllAtomTypes(), lj38.getAllAtomNumbers(), atsPerMol38, new double[38],
                lj38.getNoOfAtoms(), lj38.getAllCharges(), lj38.getAllSpins(), bonds38);
        
        assertEquals(e38, ENERGYLJ38, NUMACC);
        
        final double e55 = ljFF.energyCalculation(-1, 0, lj55.getAll1DCartes(),
                lj55.getAllAtomTypes(), lj55.getAllAtomNumbers(), atsPerMol55, new double[55],
                lj55.getNoOfAtoms(), lj55.getAllCharges(), lj55.getAllSpins(), bonds55);
        
        assertEquals(e55, ENERGYLJ55, NUMACC);
    }
    
    @Test
    public void testGradientCalculationLJMins() {
        
        final int[] atsPerMol38 = new int[38];
        for(int i = 0; i < 38; i++){
            atsPerMol38[i] = 1;
        }
        
        final int[] atsPerMol55 = new int[55];
        for(int i = 0; i < 55; i++){
            atsPerMol55[i] = 1;
        }
        
        CartesianCoordinates lj38 = null;
        CartesianCoordinates lj55 = null;
        try {
            final String[] fileCont38 = LJ38MIN.split("\n");
            lj38 = Input.parseCartesFromFileData(fileCont38, 38, atsPerMol38, new short[38], new float[38]);
        
            final String[] fileCont55 = LJ55MIN.split("\n");
            lj55 = Input.parseCartesFromFileData(fileCont55, 55, atsPerMol55, new short[55], new float[55]);
        } catch(Exception e){
           fail(e.toString());
        }

        final BondInfo bonds38 = new SimpleBondInfo(38);
        final BondInfo bonds55 = new SimpleBondInfo(55);
        
        final LennardJonesFF ljFF = new LennardJonesFF(true);
       
        final Gradient grad38 = new Gradient(3, 38);
        ljFF.gradientCalculation(-1, 0, lj38.getAll1DCartes(),
                lj38.getAllAtomTypes(), lj38.getAllAtomNumbers(), atsPerMol38, new double[38],
                lj38.getNoOfAtoms(), lj38.getAllCharges(), lj38.getAllSpins(), bonds38,
                grad38);
        
        assertEquals(grad38.getTotalEnergy(), ENERGYLJ38, NUMACC);
        
        final double[][] gradient38 = grad38.getTotalGradient();
        
        double gradTot38 = 0.0;
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 38; j++) gradTot38 += gradient38[i][j];
        }
        
        assertEquals(gradTot38, 0.0, NUMACC);
        
        final Gradient grad55 = new Gradient(3, 55);
        ljFF.gradientCalculation(-1, 0, lj55.getAll1DCartes(),
                lj55.getAllAtomTypes(), lj55.getAllAtomNumbers(), atsPerMol55, new double[55],
                lj55.getNoOfAtoms(), lj55.getAllCharges(), lj55.getAllSpins(), bonds55,
                grad55);
        
        assertEquals(grad55.getTotalEnergy(), ENERGYLJ55, NUMACC);
        
        final double[][] gradient55 = grad55.getTotalGradient();
        
        double gradTot55 = 0.0;
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 55; j++) gradTot55 += gradient55[i][j];
        }
        
        assertEquals(gradTot55, 0.0, NUMACC);
    }
}
