/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

import java.util.HashMap;

/**
 * This provides any other class with atomic properties for a given atom,
 * identified through its "name". The returned values are always in atomic
 * units.
 * Since these atomic features are accessed quite often, they are all stored in
 * HashMaps along side full arrays to provide a better performance.
 * @author Johannes Dieterich
 * @version 2015-07-27
 */
public class AtomicProperties {
    
    private static final boolean DEBUG = false;
    private static final double RADIUSDEFAULT = 5.0;
    private static final int MAXATOMNO = 113;
    
    private static final HashMap<String,Double> radiusTable = new HashMap<>();
    private static final double[] radiusArray = new double[MAXATOMNO];

    private static final HashMap<String,Double> weightTable = new HashMap<>();
    private static final double[] weightArray = new double[MAXATOMNO];

    private static final HashMap<String,Short> numberTable = new HashMap<>();

    private static final HashMap<String,Double> ljSigmaTable = new HashMap<>();

    private static final HashMap<String,Double> ljEpsilonTable = new HashMap<>();
    
    private static final String[] atomSymbolArray = new String[MAXATOMNO];

    // it's a singleton ;-)
    private static AtomicProperties props = new AtomicProperties();

    /**
     * A constructor populating both the HashMap containing radii and the one
     * containing weights. Private since we have a singleton design pattern.
     */
    private AtomicProperties(){
        
        // XXX as we miss some radii, use a safe "default" of 5.0 bohr for them
        for(int i = 0; i < 113; i++){
            radiusArray[i] = RADIUSDEFAULT;
        }
        radiusTable.put("Po", RADIUSDEFAULT);
        radiusTable.put("At", RADIUSDEFAULT);
        radiusTable.put("Rn", RADIUSDEFAULT);
        radiusTable.put("Fr", RADIUSDEFAULT);
        radiusTable.put("Ra", RADIUSDEFAULT);
        radiusTable.put("Ac", RADIUSDEFAULT);
        radiusTable.put("Th", RADIUSDEFAULT);
        radiusTable.put("Pa", RADIUSDEFAULT);
        radiusTable.put("Np", RADIUSDEFAULT);
        radiusTable.put("Pu", RADIUSDEFAULT);
        radiusTable.put("Am", RADIUSDEFAULT);
        radiusTable.put("Cm", RADIUSDEFAULT);
        radiusTable.put("Bk", RADIUSDEFAULT);
        radiusTable.put("Cf", RADIUSDEFAULT);
        radiusTable.put("Es", RADIUSDEFAULT);
        radiusTable.put("Fm", RADIUSDEFAULT);
        radiusTable.put("Md", RADIUSDEFAULT);
        radiusTable.put("No", RADIUSDEFAULT);
        radiusTable.put("Lr", RADIUSDEFAULT);
        radiusTable.put("Rf", RADIUSDEFAULT);
        radiusTable.put("Db", RADIUSDEFAULT);
        radiusTable.put("Sg", RADIUSDEFAULT);
        radiusTable.put("Bh", RADIUSDEFAULT);
        radiusTable.put("Hs", RADIUSDEFAULT);
        radiusTable.put("Mt", RADIUSDEFAULT);
        radiusTable.put("Ds", RADIUSDEFAULT);
        radiusTable.put("Rg", RADIUSDEFAULT);
        radiusTable.put("Cn", RADIUSDEFAULT);
        // defaults till here
        
        /*
         * first we populate the radius table
         * with [1] J. C. Slater, J. Chem. Phys., 41, 10, 1964
         * or [2] wikipedia.org: covalent radius
         */
        // ATTENTION: THIS IS IMPORTANT FOR BOND CREATION
        // a dummy of same radius as a H
        radiusTable.put("XX", 0.47); //XXX a bigger radius might be helpful to avoid problems in CD/DD
        // dummies for Tinker: TIP4P/5P etc
        radiusTable.put("LP", 1.0);
        radiusTable.put("M", 1.0);
        radiusTable.put("DM", 1.0);
        // [1]
        radiusTable.put("H", 0.47); //XXX a bigger radius might be helpful to avoid problems in CD/DD
        // [2]
        radiusTable.put("He", 0.53);
        // [2]
        radiusTable.put("Li", 2.42);
        // [2]
        radiusTable.put("Be", 1.81);
        // [2]
        radiusTable.put("B", 1.59);
        // [1]
        radiusTable.put("C", 1.32);
        // [1]
        radiusTable.put("N", 1.23);
        // [1]
        radiusTable.put("O", 1.13);
        // [1]
        radiusTable.put("F", 1.08);
        // [2]
        radiusTable.put("Ne", 1.10);
        // [2]
        radiusTable.put("Na", 3.15);
        // [2]
        radiusTable.put("Mg", 2.67);
        // [2]
        radiusTable.put("Al", 2.29);
        // [1]
        radiusTable.put("Si", 2.08);
        // [2]
        radiusTable.put("P", 2.02);
        // [1]
        radiusTable.put("S", 1.89);
        // [1]
        radiusTable.put("Cl", 1.89);
        // [2]
        radiusTable.put("Ar", 1.85);
        // [2]
        radiusTable.put("K", 3.19);
        // [2]
        radiusTable.put("Ca", 3.29);
        // [2]
        radiusTable.put("Sc", 3.21);
        // [2]
        radiusTable.put("Ti", 3.02);
        // [2]
        radiusTable.put("V", 2.36);
        // [2]
        radiusTable.put("Cr", 2.63);
        // [2] there is a big difference between high-spin and low-spin case,
        // we use the middle: 150 pm
        radiusTable.put("Mn", 2.83);
        // [2] and again a big low-spin/high-spin difference: 142 pm
        radiusTable.put("Fe", 2.68);
        // [2] use: 138 pm
        radiusTable.put("Co", 2.61);
        // [2]
        radiusTable.put("Ni", 2.34);
        // [1]
        radiusTable.put("Cu", 2.55);
        // [2]
        radiusTable.put("Zn", 2.31);
        // [2]
        radiusTable.put("Ga", 2.31);
        // [2]
        radiusTable.put("Ge", 2.31);
        // [2]
        radiusTable.put("As", 2.25);
        // [2]
        radiusTable.put("Se", 2.28);
        // [2]
        radiusTable.put("Br", 2.28);
        // [2]
        radiusTable.put("Kr", 2.19);
        // [2]
        radiusTable.put("Rb", 4.16);
        // [2]
        radiusTable.put("Sr", 3.68);
        // [2]
        radiusTable.put("Y", 3.59);
        // [2]
        radiusTable.put("Zr", 3.31);
        // [2]
        radiusTable.put("Nb", 3.10);
        // [2]
        radiusTable.put("Mo", 2.91);
        // [2]
        radiusTable.put("Tc", 2.78);
        // [2]
        radiusTable.put("Ru", 2.76);
        // [2]
        radiusTable.put("Rh", 2.68);
        // [2]
        radiusTable.put("Pd", 2.63);
        // [2]
        radiusTable.put("Ag", 2.74);
        // [2]
        radiusTable.put("Cd", 2.72);
        // [2]
        radiusTable.put("In", 2.68);
        // [2]
        radiusTable.put("Sn", 2.63);
        // [2]
        radiusTable.put("Sb", 2.61);
        // [2]
        radiusTable.put("Te", 2.61);
        // [2]
        radiusTable.put("I", 2.63);
        // [2]
        radiusTable.put("Xe", 2.65);
        // [2]
        radiusTable.put("Cs", 4.61);
        // [2]
        radiusTable.put("Ba", 4.06);
        // [2]
        radiusTable.put("La", 3.91);
        // [2]
        radiusTable.put("Ce", 3.86);
        // [2]
        radiusTable.put("Pr", 3.84);
        // [2]
        radiusTable.put("Nd", 3.80);
        // [2]
        radiusTable.put("Pm", 3.76);
        // [2]
        radiusTable.put("Sm", 3.74);
        // [2]
        radiusTable.put("Eu", 3.74);
        // [2]
        radiusTable.put("Gd", 3.71);
        // [2]
        radiusTable.put("Tb", 3.67);
        // [2]
        radiusTable.put("Dy", 3.63);
        // [2]
        radiusTable.put("Ho", 3.63);
        // [2]
        radiusTable.put("Er", 3.57);
        // [2]
        radiusTable.put("Tm", 3.59);
        // [2]
        radiusTable.put("Yb", 3.53);
        // [2]
        radiusTable.put("Lu", 3.53);
        // [2]
        radiusTable.put("Hf", 3.31);
        // [2]
        radiusTable.put("Ta", 3.21);
        // [1]
        radiusTable.put("W", 2.55);
        // [2]
        radiusTable.put("Re", 2.85);
        // [2]
        radiusTable.put("Os", 2.72);
        // [2]
        radiusTable.put("Ir", 2.66);
        // [2]
        radiusTable.put("Pt", 2.57);
        // [2]
        radiusTable.put("Au", 2.57);
        // [2]
        radiusTable.put("Hg", 2.49);
        // [2]
        radiusTable.put("Tl", 3.21);
        // [2]
        radiusTable.put("Pb", 2.76);
        // [2]
        radiusTable.put("Bi", 2.80);
        // [2]
        radiusTable.put("Po", 2.65);
        // [2]
        radiusTable.put("At", 2.83);
        // [2]
        radiusTable.put("Rn", 2.83);
        // [2]
        radiusTable.put("Fr", 4.91);
        // [2]
        radiusTable.put("Ra", 4.18);
        // [2]
        radiusTable.put("U", 3.70);
        
        /*
         * same as above, as an array
         */
        radiusArray[0] = 0.47; // all dummies are sized as an hydrogen
        radiusArray[ 1] = 0.47;
        radiusArray[ 2] = 0.53;
        radiusArray[ 3] = 2.42;
        radiusArray[ 4] = 1.81;
        radiusArray[ 5] = 1.59;
        radiusArray[ 6] = 1.32;
        radiusArray[ 7] = 1.23;
        radiusArray[ 8] = 1.13;
        radiusArray[ 9] = 1.08;
        radiusArray[10] = 1.10;
        radiusArray[11] = 3.15;
        radiusArray[12] = 2.67;
        radiusArray[13] = 2.29;
        radiusArray[14] = 2.08;
        radiusArray[15] = 2.02;
        radiusArray[16] = 1.89;
        radiusArray[17] = 1.89;
        radiusArray[18] = 1.85;
        radiusArray[19] = 3.19;
        radiusArray[20] = 3.29;
        radiusArray[21] = 3.21;
        radiusArray[22] = 3.02;
        radiusArray[23] = 2.36;
        radiusArray[24] = 2.63;
        radiusArray[25] = 2.83;
        radiusArray[26] = 2.68;
        radiusArray[27] = 2.61;
        radiusArray[28] = 2.34;
        radiusArray[29] = 2.55;
        radiusArray[30] = 2.31;
        radiusArray[31] = 2.31;
        radiusArray[32] = 2.31;
        radiusArray[33] = 2.25;
        radiusArray[34] = 2.28;
        radiusArray[35] = 2.28;
        radiusArray[36] = 2.19;
        radiusArray[37] = 4.16;
        radiusArray[38] = 3.68;
        radiusArray[39] = 3.59;
        radiusArray[40] = 3.31;
        radiusArray[41] = 3.10;
        radiusArray[42] = 2.91;
        radiusArray[43] = 2.78;
        radiusArray[44] = 2.76;
        radiusArray[45] = 2.68;
        radiusArray[46] = 2.63;
        radiusArray[47] = 2.74;
        radiusArray[48] = 2.72;
        radiusArray[49] = 2.68;
        radiusArray[50] = 2.63;
        radiusArray[51] = 2.61;
        radiusArray[52] = 2.61;
        radiusArray[53] = 2.63;
        radiusArray[54] = 2.65;
        radiusArray[55] = 4.61;
        radiusArray[56] = 4.06;
        radiusArray[57] = 3.91;
        radiusArray[58] = 3.86;
        radiusArray[59] = 3.84;
        radiusArray[60] = 3.80;
        radiusArray[61] = 3.76;
        radiusArray[62] = 3.74;
        radiusArray[63] = 3.74;
        radiusArray[64] = 3.71;
        radiusArray[65] = 3.67;
        radiusArray[66] = 3.63;
        radiusArray[67] = 3.63;
        radiusArray[68] = 3.57;
        radiusArray[69] = 3.59;
        radiusArray[70] = 3.53;
        radiusArray[71] = 3.53;
        radiusArray[72] = 3.31;
        radiusArray[73] = 3.21;
        radiusArray[74] = 2.55;
        radiusArray[75] = 2.85;
        radiusArray[76] = 2.72;
        radiusArray[77] = 2.66;
        radiusArray[78] = 2.57;
        radiusArray[79] = 2.57;
        radiusArray[80] = 2.49;
        radiusArray[81] = 3.21;
        radiusArray[82] = 2.76;
        radiusArray[83] = 2.80;
        radiusArray[84] = 2.65;
        radiusArray[85] = 2.83;
        radiusArray[86] = 2.83;
        radiusArray[87] = 4.91;
        radiusArray[88] = 4.18;
        // XXX missing
        radiusArray[92] = 3.70;
        

        /*
         * now we populate the weight table
         */
        weightTable.put("XX", 0.0);
        weightTable.put("LP", 0.0);
        weightTable.put("M", 0.0);
        weightTable.put("DM", 0.0);
        weightTable.put("H", 1.00794);
        weightTable.put("He", 4.0026);
        weightTable.put("Li", 6.941);
        weightTable.put("Be", 9.012182);
        weightTable.put("B", 10.811);
        weightTable.put("C", 12.010);
        weightTable.put("N", 14.006);
        weightTable.put("O", 15.999);
        weightTable.put("F", 18.9984032);
        weightTable.put("Ne", 20.1797);
        weightTable.put("Na", 22.98977);
        weightTable.put("Mg", 24.3050);
        weightTable.put("Al", 26.9815386);
        weightTable.put("Si", 28.0855);
        weightTable.put("P", 30.973762);
        weightTable.put("S", 32.065);
        weightTable.put("Cl", 35.453);
        weightTable.put("Ar", 39.948);
        weightTable.put("K", 39.0983);
        weightTable.put("Ca", 40.078);
        weightTable.put("Sc", 44.955912);
        weightTable.put("Ti", 47.867);
        weightTable.put("V", 50.9415);
        weightTable.put("Cr", 51.9961);
        weightTable.put("Mn", 54.938045);
        weightTable.put("Fe", 55.845);
        weightTable.put("Co", 58.933195);
        weightTable.put("Ni", 58.6934);
        weightTable.put("Cu", 63.546);
        weightTable.put("Zn", 65.38);
        weightTable.put("Ga", 69.723);
        weightTable.put("Ge", 72.64);
        weightTable.put("As", 74.92160);
        weightTable.put("Se", 78.96);
        weightTable.put("Br", 79.904);
        weightTable.put("Kr", 83.798);
        weightTable.put("Rb", 85.4678);
        weightTable.put("Sr", 87.62);
        weightTable.put("Y", 88.90585);
        weightTable.put("Zr", 91.224);
        weightTable.put("Nb", 92.90638);
        weightTable.put("Mo", 95.94);
        weightTable.put("Tc", 98.00);
        weightTable.put("Ru", 101.07);
        weightTable.put("Rh", 102.90550);
        weightTable.put("Pd", 106.42);
        weightTable.put("Ag", 107.8682);
        weightTable.put("Cd", 112.411);
        weightTable.put("In", 114.818);
        weightTable.put("Sn", 118.710);
        weightTable.put("Sb", 121.760);
        weightTable.put("Te", 127.60);
        weightTable.put("I", 126.90447);
        weightTable.put("Xe", 131.293);
        weightTable.put("Cs", 132.9054519);
        weightTable.put("Ba", 137.33);
        weightTable.put("La", 138.90547);
        weightTable.put("Ce", 140.116);
        weightTable.put("Pr", 140.90765);
        weightTable.put("Nd", 144.242);
        weightTable.put("Pm", 145.00);
        weightTable.put("Sm", 150.36);
        weightTable.put("Eu", 151.964);
        weightTable.put("Gd", 157.25);
        weightTable.put("Tb", 158.92535);
        weightTable.put("Dy", 162.500);
        weightTable.put("Ho", 164.93032);
        weightTable.put("Er", 167.259);
        weightTable.put("Tm", 168.93421);
        weightTable.put("Yb", 173.04);
        weightTable.put("Lu", 174.967);
        weightTable.put("Hf", 178.49);
        weightTable.put("Ta", 180.94788);
        weightTable.put("W", 183.84);
        weightTable.put("Re", 186.207);
        weightTable.put("Os", 190.23);
        weightTable.put("Ir", 192.217);
        weightTable.put("Pt", 195.084);
        weightTable.put("Au", 196.966569);
        weightTable.put("Hg", 200.59);
        weightTable.put("Tl", 204.3833);
        weightTable.put("Pb", 207.2);
        weightTable.put("Bi", 208.98040);
        weightTable.put("Po", 208.89);
        weightTable.put("At", 209.99);
        weightTable.put("Rn", 222.02);
        weightTable.put("Fr", 223.02);
        weightTable.put("Ra", 226.03);
        weightTable.put("Ac", 224.03);
        weightTable.put("Th", 232.04);
        weightTable.put("Pa", 231.04);
        weightTable.put("U", 238.02891);
        weightTable.put("Np", 237.05);
        weightTable.put("Pu", 244.06);
        weightTable.put("Am", 243.06);
        weightTable.put("Cm", 247.07);
        weightTable.put("Bk", 247.07);
        weightTable.put("Cf", 251.08);
        weightTable.put("Es", 252.08);
        weightTable.put("Fm", 257.10);
        weightTable.put("Md", 258.10);
        weightTable.put("No", 259.10);
        weightTable.put("Lr", 262.11);
        weightTable.put("Rf", 265.12);
        weightTable.put("Db", 268.13);
        weightTable.put("Sg", 271.13);
        weightTable.put("Bh", 270.);
        weightTable.put("Hs", 277.15);
        weightTable.put("Mt", 276.15);
        weightTable.put("Ds", 281.16);
        weightTable.put("Rg", 280.16);
        weightTable.put("Cn", 285.17);

        weightArray[ 0] = 0.0;
        weightArray[ 1] = 1.00794;
        weightArray[ 2] = 4.0026;
        weightArray[ 3] = 6.941;
        weightArray[ 4] = 9.012182;
        weightArray[ 5] = 10.811;
        weightArray[ 6] = 12.010;
        weightArray[ 7] = 14.006;
        weightArray[ 8] = 15.999;
        weightArray[ 9] = 18.9984032;
        weightArray[10] = 20.1797;
        weightArray[11] = 22.98977;
        weightArray[12] = 24.3050;
        weightArray[13] = 26.9815386;
        weightArray[14] = 28.0855;
        weightArray[15] = 30.973762;
        weightArray[16] = 32.065;
        weightArray[17] = 35.453;
        weightArray[18] = 39.948;
        weightArray[19] = 39.0983;
        weightArray[20] = 40.078;
        weightArray[21] = 44.955912;
        weightArray[22] = 47.867;
        weightArray[23] = 50.9415;
        weightArray[24] = 51.9961;
        weightArray[25] = 54.938045;
        weightArray[26] = 55.845;
        weightArray[27] = 58.933195;
        weightArray[28] = 58.6934;
        weightArray[29] = 63.546;
        weightArray[30] = 65.38;
        weightArray[31] = 69.723;
        weightArray[32] = 72.64;
        weightArray[33] = 74.92160;
        weightArray[34] = 78.96;
        weightArray[35] = 79.904;
        weightArray[36] = 83.798;
        weightArray[37] = 85.4678;
        weightArray[38] = 87.62;
        weightArray[39] = 88.90585;
        weightArray[40] = 91.224;
        weightArray[41] = 92.90638;
        weightArray[42] = 95.94;
        weightArray[43] = 98.00;
        weightArray[44] = 101.07;
        weightArray[45] = 102.90550;
        weightArray[46] = 106.42;
        weightArray[47] = 107.8682;
        weightArray[48] = 112.411;
        weightArray[49] = 114.818;
        weightArray[50] = 118.710;
        weightArray[51] = 121.760;
        weightArray[52] = 127.60;
        weightArray[53] = 126.90447;
        weightArray[54] = 131.293;
        weightArray[55] = 132.9054519;
        weightArray[56] = 137.33;
        weightArray[57] = 138.90547;
        weightArray[58] = 140.116;
        weightArray[59] = 140.90765;
        weightArray[60] = 144.242;
        weightArray[61] = 145.00;
        weightArray[62] = 150.36;
        weightArray[63] = 151.964;
        weightArray[64] = 157.25;
        weightArray[65] = 158.92535;
        weightArray[66] = 162.500;
        weightArray[67] = 164.93032;
        weightArray[68] = 167.259;
        weightArray[69] = 168.93421;
        weightArray[70] = 173.04;
        weightArray[71] = 174.967;
        weightArray[72] = 178.49;
        weightArray[73] = 180.94788;
        weightArray[74] = 183.84;
        weightArray[75] = 186.207;
        weightArray[76] = 190.23;
        weightArray[77] = 192.217;
        weightArray[78] = 195.084;
        weightArray[79] = 196.966569;
        weightArray[80] = 200.59;
        weightArray[81] = 204.3833;
        weightArray[82] = 207.2;
        weightArray[83] = 208.98040;
        weightArray[84] = 208.89;
        weightArray[85] = 209.99;
        weightArray[86] = 222.02;
        weightArray[87] = 223.02;
        weightArray[88] = 226.03;
        weightArray[89] = 224.03;
        weightArray[90] = 232.04;
        weightArray[91] = 231.04;
        weightArray[92] = 238.02891;
        weightArray[93] = 237.05;
        weightArray[94] = 244.06;
        weightArray[95] = 243.06;
        weightArray[96] = 247.07;
        weightArray[97] = 247.07;
        weightArray[98] = 251.08;
        weightArray[99] = 252.08;
        weightArray[100] = 257.10;
        weightArray[101] = 258.10;
        weightArray[102] = 259.10;
        weightArray[103] = 262.11;
        weightArray[104] = 265.12;
        weightArray[105] = 268.13;
        weightArray[106] = 271.13;
        weightArray[107] = 270.;
        weightArray[108] = 277.15;
        weightArray[109] = 276.15;
        weightArray[110] = 281.16;
        weightArray[111] = 280.16;
        weightArray[112] = 285.17;
        
        /*
         * then we populate the atomic number table
         */
        numberTable.put("XX", (short) 0);
        // dummies for Tinker: TIP4P/5P etc
        numberTable.put("LP", (short) 0);
        numberTable.put("M", (short) 0);
        numberTable.put("DM", (short) 0);
        numberTable.put("H", (short) 1);
        numberTable.put("He", (short) 2);
        numberTable.put("Li", (short) 3);
        numberTable.put("Be", (short) 4);
        numberTable.put("B", (short) 5);
        numberTable.put("C", (short) 6);
        numberTable.put("N", (short) 7);
        numberTable.put("O", (short) 8);
        numberTable.put("F", (short) 9);
        numberTable.put("Ne", (short) 10);
        numberTable.put("Na", (short) 11);
        numberTable.put("Mg", (short) 12);
        numberTable.put("Al", (short) 13);
        numberTable.put("Si", (short) 14);
        numberTable.put("P", (short) 15);
        numberTable.put("S", (short) 16);
        numberTable.put("Cl", (short) 17);
        numberTable.put("Ar", (short) 18);
        numberTable.put("K", (short) 19);
        numberTable.put("Ca", (short) 20);
        numberTable.put("Sc", (short) 21);
        numberTable.put("Ti", (short) 22);
        numberTable.put("V", (short) 23);
        numberTable.put("Cr", (short) 24);
        numberTable.put("Mn", (short) 25);
        numberTable.put("Fe", (short) 26);
        numberTable.put("Co", (short) 27);
        numberTable.put("Ni", (short) 28);
        numberTable.put("Cu", (short) 29);
        numberTable.put("Zn", (short) 30);
        numberTable.put("Ga", (short) 31);
        numberTable.put("Ge", (short) 32);
        numberTable.put("As", (short) 33);
        numberTable.put("Se", (short) 34);
        numberTable.put("Br", (short) 35);
        numberTable.put("Kr", (short) 36);
        numberTable.put("Rb", (short) 37);
        numberTable.put("Sr", (short) 38);
        numberTable.put("Y", (short) 39);
        numberTable.put("Zr", (short) 40);
        numberTable.put("Nb", (short) 41);
        numberTable.put("Mo", (short) 42);
        numberTable.put("Tc", (short) 43);
        numberTable.put("Ru", (short) 44);
        numberTable.put("Rh", (short) 45);
        numberTable.put("Pd", (short) 46);
        numberTable.put("Ag", (short) 47);
        numberTable.put("Cd", (short) 48);
        numberTable.put("In", (short) 49);
        numberTable.put("Sn", (short) 50);
        numberTable.put("Sb", (short) 51);
        numberTable.put("Te", (short) 52);
        numberTable.put("I", (short) 53);
        numberTable.put("Xe", (short) 54);
        numberTable.put("Cs", (short) 55);
        numberTable.put("Ba", (short) 56);
        numberTable.put("La", (short) 57);
        numberTable.put("Ce", (short) 58);
        numberTable.put("Pr", (short) 59);
        numberTable.put("Nd", (short) 60);
        numberTable.put("Pm", (short) 61);
        numberTable.put("Sm", (short) 62);
        numberTable.put("Eu", (short) 63);
        numberTable.put("Gd", (short) 64);
        numberTable.put("Tb", (short) 65);
        numberTable.put("Dy", (short) 66);
        numberTable.put("Ho", (short) 67);
        numberTable.put("Er", (short) 68);
        numberTable.put("Tm", (short) 69);
        numberTable.put("Yb", (short) 70);
        numberTable.put("Lu", (short) 71);
        numberTable.put("Hf", (short) 72);
        numberTable.put("Ta", (short) 73);
        numberTable.put("W", (short) 74);
        numberTable.put("Re", (short) 75);
        numberTable.put("Os", (short) 76);
        numberTable.put("Ir", (short) 77);
        numberTable.put("Pt", (short) 78);
        numberTable.put("Au", (short) 79);
        numberTable.put("Hg", (short) 80);
        numberTable.put("Tl", (short) 81);
        numberTable.put("Pb", (short) 82);
        numberTable.put("Bi", (short) 83);
        numberTable.put("Po", (short) 84);
        numberTable.put("At", (short) 85);
        numberTable.put("Rn", (short) 86);
        numberTable.put("Fr", (short) 87);
        numberTable.put("Ra", (short) 88);
        numberTable.put("Ac", (short) 89);
        numberTable.put("Th", (short) 90);
        numberTable.put("Pa", (short) 91);
        numberTable.put("U", (short) 92);
        numberTable.put("Np", (short) 93);
        numberTable.put("Pu", (short) 94);
        numberTable.put("Am", (short) 95);
        numberTable.put("Cm", (short) 96);
        numberTable.put("Bk", (short) 97);
        numberTable.put("Cf", (short) 98);
        numberTable.put("Es", (short) 99);
        numberTable.put("Fm", (short) 100);
        numberTable.put("Md", (short) 101);
        numberTable.put("No", (short) 102);
        numberTable.put("Lr", (short) 103);
        numberTable.put("Rf", (short) 104);
        numberTable.put("Db", (short) 105);
        numberTable.put("Sg", (short) 106);
        numberTable.put("Bh", (short) 107);
        numberTable.put("Hs", (short) 108);
        numberTable.put("Mt", (short) 109);
        numberTable.put("Ds", (short) 110);
        numberTable.put("Rg", (short) 111);
        numberTable.put("Cn", (short) 112);
        // that should be enough for the time being ;-)
        
        atomSymbolArray[0]   = "XX";
        atomSymbolArray[1]   = "H";
        atomSymbolArray[2]   = "He";
        atomSymbolArray[3]   = "Li";
        atomSymbolArray[4]   = "Be";
        atomSymbolArray[5]   = "B";
        atomSymbolArray[6]   = "C";
        atomSymbolArray[7]   = "N";
        atomSymbolArray[8]   = "O";
        atomSymbolArray[9]   = "F";
        atomSymbolArray[10]  = "Ne";
        atomSymbolArray[11]  = "Na";
        atomSymbolArray[12]  = "Mg";
        atomSymbolArray[13]  = "Al";
        atomSymbolArray[14]  = "Si";
        atomSymbolArray[15]  = "P";
        atomSymbolArray[16]  = "S";
        atomSymbolArray[17]  = "Cl";
        atomSymbolArray[18]  = "Ar";
        atomSymbolArray[19]  = "K";
        atomSymbolArray[20]  = "Ca";
        atomSymbolArray[21]  = "Sc";
        atomSymbolArray[22]  = "Ti";
        atomSymbolArray[23]  = "V";
        atomSymbolArray[24]  = "Cr";
        atomSymbolArray[25]  = "Mn";
        atomSymbolArray[26]  = "Fe";
        atomSymbolArray[27]  = "Co";
        atomSymbolArray[28]  = "Ni";
        atomSymbolArray[29]  = "Cu";
        atomSymbolArray[30]  = "Zn";
        atomSymbolArray[31]  = "Ga";
        atomSymbolArray[32]  = "Ge";
        atomSymbolArray[33]  = "As";
        atomSymbolArray[34]  = "Se";
        atomSymbolArray[35]  = "Br";
        atomSymbolArray[36]  = "Kr";
        atomSymbolArray[37]  = "Rb";
        atomSymbolArray[38]  = "Sr";
        atomSymbolArray[39]  = "Y";
        atomSymbolArray[40]  = "Zr";
        atomSymbolArray[41]  = "Nb";
        atomSymbolArray[42]  = "Mo";
        atomSymbolArray[43]  = "Tc";
        atomSymbolArray[44]  = "Ru";
        atomSymbolArray[45]  = "Rh";
        atomSymbolArray[46]  = "Pd";
        atomSymbolArray[47]  = "Ag";
        atomSymbolArray[48]  = "Cd";
        atomSymbolArray[49]  = "In";
        atomSymbolArray[50]  = "Sn";
        atomSymbolArray[51]  = "Sb";
        atomSymbolArray[52]  = "Te";
        atomSymbolArray[53]  = "I";
        atomSymbolArray[54]  = "Xe";
        atomSymbolArray[55]  = "Cs";
        atomSymbolArray[56]  = "Ba";
        atomSymbolArray[57]  = "La";
        atomSymbolArray[58]  = "Ce";
        atomSymbolArray[59]  = "Pr";
        atomSymbolArray[60]  = "Nd";
        atomSymbolArray[61]  = "Pm";
        atomSymbolArray[62]  = "Sm";
        atomSymbolArray[63]  = "Eu";
        atomSymbolArray[64]  = "Gd";
        atomSymbolArray[65]  = "Tb";
        atomSymbolArray[66]  = "Dy";
        atomSymbolArray[67]  = "Ho";
        atomSymbolArray[68]  = "Er";
        atomSymbolArray[69]  = "Tm";
        atomSymbolArray[70]  = "Yb";
        atomSymbolArray[71]  = "Lu";
        atomSymbolArray[72]  = "Hf";
        atomSymbolArray[73]  = "Ta";
        atomSymbolArray[74]  = "W";
        atomSymbolArray[75]  = "Re";
        atomSymbolArray[76]  = "Os";
        atomSymbolArray[77]  = "Ir";
        atomSymbolArray[78]  = "Pt";
        atomSymbolArray[79]  = "Au";
        atomSymbolArray[80]  = "Hg";
        atomSymbolArray[81]  = "Tl";
        atomSymbolArray[82]  = "Pb";
        atomSymbolArray[83]  = "Bi";
        atomSymbolArray[84]  = "Po";
        atomSymbolArray[85]  = "At";
        atomSymbolArray[86]  = "Rn";
        atomSymbolArray[87]  = "Fr";
        atomSymbolArray[88]  = "Ra";
        atomSymbolArray[89]  = "Ac";
        atomSymbolArray[90]  = "Th";
        atomSymbolArray[91]  = "Pa";
        atomSymbolArray[92]  = "U";
        atomSymbolArray[93]  = "Np";
        atomSymbolArray[94]  = "Pu";
        atomSymbolArray[95]  = "Am";
        atomSymbolArray[96]  = "Cm";
        atomSymbolArray[97]  = "Bk";
        atomSymbolArray[98]  = "Cf";
        atomSymbolArray[99]  = "Es";
        atomSymbolArray[100] = "Fm";
        atomSymbolArray[101] = "Md";
        atomSymbolArray[102] = "No";
        atomSymbolArray[103] = "Lr";
        atomSymbolArray[104] = "Rf";
        atomSymbolArray[105] = "Db";
        atomSymbolArray[106] = "Sg";
        atomSymbolArray[107] = "Bh";
        atomSymbolArray[108] = "Hs";
        atomSymbolArray[109] = "Mt";
        atomSymbolArray[110] = "Ds";
        atomSymbolArray[111] = "Rg";
        atomSymbolArray[112] = "Cn";
        
        /*
         * lennard-jones sigma table
         */
        ljSigmaTable.put("He", 4.830139627);
        ljSigmaTable.put("Ne", 5.194856743);
        ljSigmaTable.put("Ar", 6.434516991);
        ljSigmaTable.put("Kr", 6.803013559);
        ljSigmaTable.put("Xe", 7.747876553);
        // the following entries are just in here out of completeness.
        ljSigmaTable.put("H2", 5.533117695);
        ljSigmaTable.put("D2", 5.533117695);
        ljSigmaTable.put("O2", 6.765219039);
        ljSigmaTable.put("CO", 7.111038895);
        ljSigmaTable.put("N2", 6.988206706);
        ljSigmaTable.put("CH4", 7.213084098);

        /*
         * lennard-jones epsilon table
         */
        ljEpsilonTable.put("He", 0.000032365);
        ljEpsilonTable.put("Ne", 0.000112739);
        ljEpsilonTable.put("Ar", 0.000379384);
        ljEpsilonTable.put("Kr", 0.000541524);
        ljEpsilonTable.put("Xe", 0.000699857);
        // the following entries are just in here out of completeness.
        ljEpsilonTable.put("H2", 0.000117172);
        ljEpsilonTable.put("D2", 0.000117172);
        ljEpsilonTable.put("O2", 0.000372099);
        ljEpsilonTable.put("CO", 0.000317314);
        ljEpsilonTable.put("N2", 0.000301006);
        ljEpsilonTable.put("CH4", 0.000469321);
        
        // then we check for updates
        final String upFile = "atomic.properties";
        try{
            updateEntries(upFile);
        } catch(Exception | Error e){
            // this is empty because we do not necessarily have an update file
            if(DEBUG) e.printStackTrace(System.err);
        }
    }

    /**
     * Maps the atomic ID as a String to the atomic number.
     * @param atomID The atomic ID.
     * @return The atomic ID as defined by the PSE, plus zero for dummies and unknowns.
     */
    public static short giveAtomicNumber(final String atomID) {

        final Short number = numberTable.get(atomID);
        if (number != null) return number;
        else {
            System.err.println("WARNING: No number available for Atom " + atomID + " default set to 0.");
            return 0;
        }
    }

    /**
     * Provides the weight for an atom and checks for errors.
     * @param atomID The atomic ID.
     * @return the weight of this atom in u, not a.u.
     */
    public static double giveWeight(final String atomID) {
    
        final Double weight = weightTable.get(atomID);
        if (weight != null) {
            return weight;
        } else {
            System.err.println("WARNING: No mass available for Atom " + atomID + " default set to 1.0 u.");
            return 1.0;
        }
    }
    
    /**
     * Faster weight version but errors will NOT be checked or caught!
     * @param which Which Atom
     * @return the weight of this atom in u, not a.u. 
     */
    public static double giveWeight(final short which) {
        return weightArray[which];
    }

    /**
     * This returns empirical atomic radii as found in: J. C. Slater, J. Chem. Phys.,
     * 41, 10, 1964. Or, for the majority: wikipedia.org, covalent radii.
     * @param atomID The atomic ID.
     * @return radius The radius in atomic units (meaning: bohr).
     */
    public static double giveRadius(final String atomID) {
       
        final Double radius = radiusTable.get(atomID);
        if (radius != null) {
            // return it
            return radius;
        } else {
            System.err.println("WARNING: No radius available for Atom " + atomID + " default set to 1.0 bohr.");
            return RADIUSDEFAULT;
        }
    }
    
    /**
     * This returns empirical atomic radii as found in: J. C. Slater, J. Chem. Phys.,
     * 41, 10, 1964. Or, for the majority: wikipedia.org, covalent radii.
     * It is a faster access but errros will NOT be checked or caught!
     * @param which The atom number
     * @return The radius in atomic units (meaning: bohr).
     */
    public static double giveRadius(final short which) {
        return radiusArray[which];
    }

    /**
     * The returned sigma value is in a. u.'s (here: bohr). Values were obtained from
     * Berry, Rice, Ross: Physical Chemistry, 1980, Wiley.
     * @param atomID The atomic ID.
     * @return sigma The sigma parameter of a Lennard-Jones force field.
     */
    public static double giveLennardJonesSigma(final String atomID) {

        final Double sigma = ljSigmaTable.get(atomID);
        if (sigma != null) {
            // return it
            return sigma;
        } else {
            System.err.println("WARNING: No Lennard-Jones sigma available for Atom " + atomID + " default set to 1.0.");
            return 1.0;
        }
    }

    /**
     * The returned epsilon value is in a. u.'s (here: hartree). Values were obtained from
     * Berry, Rice, Ross: Physical Chemistry, 1980, Wiley.
     * @param atomID The atomic ID.
     * @return epsilon The epsilon parameter of a Lennard-Jones force field.
     */
    public static double giveLennardJonesEpsilon(final String atomID) {

        final Double epsilon = ljEpsilonTable.get(atomID);
        if (epsilon != null) {
            // return it
            return epsilon;
        } else {
            System.err.println("WARNING: No Lennard-Jones epsilon available for Atom " + atomID + " default set to 0.0001.");
            return 0.0001;
        }
    }
    
    /**
     * Give the atom symbol for a given atom number.
     * @param atomNo the atom number. must be between 0 and 112
     * @return the atom symbol, may throw a runtime exception for illegal atom numbers. 
     */
    public static String giveAtomSymbol(final short atomNo){
        return atomSymbolArray[atomNo];
    }
    
    private static void updateEntries(final String upFile) throws Exception {
        final String[] sa = org.ogolem.io.InputPrimitives.readFileIn(upFile);
        for (final String s : sa) {
            final String[] sa2 = s.trim().split("\\s+");
            final double d = Double.parseDouble(sa2[2]);
            if (sa2[0].equalsIgnoreCase("WEIGHT")) {
                if(weightTable.containsKey(sa2[1])){
                    weightTable.remove(sa2[1]);
                    final int no = numberTable.get(sa2[1]);
                    /*if(no >= weightArray.length){
                        weightArray = new double[no+1];
                    }*/ //XXX we would like to refresh perhaps then...
                    weightArray[no] = d;
                }
                weightTable.put(sa2[1], d);
            } else if (sa2[0].equalsIgnoreCase("RADIUS")) {
                if(radiusTable.containsKey(sa2[1])){
                    radiusTable.remove(sa2[1]);
                    final int no = numberTable.get(sa2[1]);
                    //XXX same here as above
                    radiusArray[no] = d;
                }
                radiusTable.put(sa2[1], d);
            } else if (sa2[0].equalsIgnoreCase("NUMBER")) {
                final short i = Short.parseShort(sa2[2]);
                if(numberTable.containsKey(sa2[1]))numberTable.remove(sa2[1]);
                numberTable.put(sa2[1], i);
            } else if (sa2[0].equalsIgnoreCase("LJSIGMA")) {
                if(ljSigmaTable.containsKey(sa2[1]))ljSigmaTable.remove(sa2[1]);
                ljSigmaTable.put(sa2[1], d);
            } else if (sa2[0].equalsIgnoreCase("LJEPSILON")) {
                if(ljEpsilonTable.containsKey(sa2[1]))ljEpsilonTable.remove(sa2[1]);
                ljEpsilonTable.put(sa2[1], d);
            } else {
                System.err.println("WARNING: Couldn't parse entry " + s + " from atomic.properties file. Continue.");
            }
        }
    }
}
