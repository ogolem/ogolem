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
public class TIP4PForceFieldTest {
    
    private static final double NUMACC = 1e-6;

    /**
     * Test of energy method, of class TIP4PForceField.
     */
    @Test
    public void testEnergy() {
        System.out.println("energy");
        final Geometry ref = null; // can be null
        final List<CartesianCoordinates> cartes = new ArrayList<>();
        final double[][] coms = new double[2][];
        // plug in the known global minimum structure for water 2
        final CartesianCoordinates c1 = new CartesianCoordinates(3,1,new int[]{3});
        final double[][] xyz1 = c1.getAllXYZCoord();
        xyz1[0][0] = -0.1858140*ANGTOBOHR;
        xyz1[1][0] = -1.1749469*ANGTOBOHR;
        xyz1[2][0] = 0.7662596*ANGTOBOHR;
        xyz1[0][1] = -0.1285513*ANGTOBOHR;
        xyz1[1][1] = -0.8984365*ANGTOBOHR;
        xyz1[2][1] = 1.6808606*ANGTOBOHR;
        xyz1[0][2] = -0.0582782*ANGTOBOHR;
        xyz1[1][2] = -0.3702550*ANGTOBOHR;
        xyz1[2][2] = 0.2638279*ANGTOBOHR;
        final String[] atoms1 = c1.getAllAtomTypes();
        atoms1[0] = "O";
        atoms1[1] = "H";
        atoms1[2] = "H";
        c1.recalcAtomNumbers();
        
        final double[] com1 = c1.calculateTheCOM();
        c1.moveCoordsToCOM();
        coms[0] = com1.clone();
            
        final CartesianCoordinates c2 = new CartesianCoordinates(3,1,new int[]{3});
        final double[][] xyz2 = c2.getAllXYZCoord();
        xyz2[0][0] = 0.1747051*ANGTOBOHR;
        xyz2[1][0] = 1.1050002*ANGTOBOHR;
        xyz2[2][0] = -0.7244430*ANGTOBOHR;
        xyz2[0][1] = -0.5650842*ANGTOBOHR;
        xyz2[1][1] = 1.3134964*ANGTOBOHR;
        xyz2[2][1] = -1.2949455*ANGTOBOHR;
        xyz2[0][2] = 0.9282185*ANGTOBOHR;
        xyz2[1][2] = 1.0652990*ANGTOBOHR;
        xyz2[2][2] = -1.313402*ANGTOBOHR;
        final String[] atoms2 = c2.getAllAtomTypes();
        atoms2[0] = "O";
        atoms2[1] = "H";
        atoms2[2] = "H";
        c2.recalcAtomNumbers();
        
        final double[] com2 = c2.calculateTheCOM();
        c2.moveCoordsToCOM();
        coms[1] = com2.clone();
        
        final TIP4PForceField instance = new TIP4PForceField();
        final CartesianCoordinates c1Adj = instance.adjustCartesians(ref, c1, 0);
        final CartesianCoordinates c2Adj = instance.adjustCartesians(ref, c2, 0);
        
        cartes.add(c1Adj);
        cartes.add(c2Adj);
        
        final int counter = 0;
        
        final double expResult = -0.009936231756023471; // according to DJW data base
        final double result = instance.energy(ref, cartes, coms, counter);
        assertEquals(expResult, result, NUMACC);
    }
    
    @Test
    public void testEnergyDistorted(){
        System.out.println("energy (distorted)");
        // uses a distorted w2 geometry which proved cumbersome
        
        final short[] spins = new short[6];
        final float[] charges = new float[6];
        final int[] atsPerMol = new int[]{3,3};
        
        final String[] w2_distr = new String[]{
            "6",
            "Energy:             -29.3998(WRONG!)",
            "O              0.2561280             -0.7959166              1.0388088",
            "H             -0.4705858             -0.6439888              1.6513454",
            "H              0.8425446             -1.2880304              1.5995301",
            "O             -0.3597538              1.3336911             -0.6140416",
            "H             -1.0841008              0.7575674             -1.1219294",
            "H              0.0397342              0.5679421             -0.0935937"
        };
        
        CartesianCoordinates w2Dist = null;
        try{
            w2Dist = Input.parseCartesFromFileData(w2_distr, 2, atsPerMol, spins, charges);
        } catch(Exception e){
            throw new Error("Error to parse in known w2 distorted minimum.",e);
        }
       
        final TIP4PForceField instance = new TIP4PForceField();
        
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
        
        final Geometry ref = new Geometry(w2Dist, 0, noOfMols, atsPerMol,
            molFlexies, allFlex, molConstr, constrXYZ,
            sids, bonds);
        final TIPnPParameters.StandardTIP4PParameters tip4pParams = new TIPnPParameters.StandardTIP4PParameters();
        for(int i = 0; i < noOfMols; i++){
            // adjust
            instance.rigidify(ref.getMoleculeAtPosition(i));
            
            // check the OH ditance and H-O-H angle
            final double[][] xyz = ref.getMoleculeAtPosition(i).getReferenceCartesians();
            final double distOH1 = CoordTranslation.distance(xyz, 0, 1);
            final double distOH2 = CoordTranslation.distance(xyz, 0, 2);
            final double angHOH = CoordTranslation.calcAngle(xyz, 1, 0, 2);
            assertEquals(distOH1,tip4pParams.getOH(),NUMACC);
            assertEquals(distOH2,tip4pParams.getOH(),NUMACC);
            assertEquals(angHOH,Math.toRadians(tip4pParams.getHOH()),NUMACC);
        }
        
        final List<CartesianCoordinates> cartes = new ArrayList<>();
        final double[][] coms = new double[2][];
        final CartesianCoordinates rigidified = ref.getCartesians();
        
        for(int mol = 0; mol < 2; mol++){
            final CartesianCoordinates molC = rigidified.giveMolecularCartes(mol, false);
            final double[] com = molC.calculateTheCOM();
            coms[mol] = com;
            molC.moveCoordsToCOM();
            
            final CartesianCoordinates molCAdj = instance.adjustCartesians(null, molC, mol);
            
            cartes.add(molCAdj);
        }
        
        final int counter = 0;
        
        final double expResult = -23.168783044616152*Constants.KJTOHARTREE; // according to phenix
        final double result = instance.energy(ref, cartes, coms, counter);
        // we want to make sure, that our implementation AT LEAST is as "high energy" as phenix for this weird case
        assertTrue("Distorted TIP4P energy fails: " + expResult + " got " + result + " diff " + (expResult-result), (expResult-result) <= 0.0);
    }
    
    @Test
    public void testEnergyBig(){
        System.out.println("energy (big)");
        // get the w20 global minimum and calculate the energy. also stresses other subsystems.
        
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
       
        final TIP4PForceField instance = new TIP4PForceField();
        
        final List<CartesianCoordinates> cartes = new ArrayList<>();
        final double[][] coms = new double[20][];
        for(int mol = 0; mol < 20; mol++){
            final CartesianCoordinates molC = w20.giveMolecularCartes(mol, false);
            final double[] com = molC.calculateTheCOM();
            coms[mol] = com;
            molC.moveCoordsToCOM();
            
            final CartesianCoordinates molCAdj = instance.adjustCartesians(null, molC, mol);
            
            cartes.add(molCAdj);
        }
        
        final Geometry ref = null; // can be null
        final int counter = 0;
        
        final double expResult = -0.3325038695132564; // according to DJW data base
        final double result = instance.energy(ref, cartes, coms, counter);
        assertEquals("Big TIP4P energy fails: " + expResult + " got " + result + " diff " + Math.abs(expResult-result),expResult, result, NUMACC);
    }
}
