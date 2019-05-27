/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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
package org.ogolem.switches;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CollisionDetection;
import org.ogolem.core.CollisionInfo;
import org.ogolem.core.CollisionInfo.Collision;
import org.ogolem.core.SimpleBondInfo;
import org.ogolem.core.ZMatrix;

/**
 * Magically glues together backbone and sidechains, hopefully resulting in somewhat of a reasonable
 * structure.
 * @author Johannes Dieterich
 * @version 2015-07-20
 */
final class MagicGlue {
    //TODO in general, we should perhaps consider whether we ask for a random sign and increment or decrement
    // the dihedral then. otherwise we potentially(!) could run into problems with bigger sidechains
    // but then we for sure should apply the same randomized values to cis and trans
    // take into consideration that we do that individually for the sidechains -> ALs of these infos

    private static final boolean bDebug = false;

    /**
     * Glues (or at least tries to glue) the backbone and the sidechains together.
     * @param back
     * @param alSidechains
     * @return
     */
    static ArrayList<CartesianCoordinates> glueTogether(final Backbone back,
            final ArrayList<Sidechain> alSidechains, final double dBlowBondDetect){

        final ArrayList<CartesianCoordinates> alCartes = new ArrayList<>(2);

        /*
         * we NEED a collision detection, we take a standard one
         */
        final CollisionDetection collDetect = new CollisionDetection(CollisionDetection.CDTYPE.SIMPLEPAIRWISE);

        final int iNoOfSidechains = alSidechains.size();

        final int[] iaConnectsCis = back.getConnectsCopy(true);

        final BondInfo baBondsBackCis = back.getBondsCopy(true);

        /*
         * we handle one after the other, first cis then trans
         */
        final CartesianCoordinates caBackCis = back.getCartesCopy(true);
        CartesianCoordinates caPrev = new CartesianCoordinates(caBackCis);
        CartesianCoordinates caNext;

        BondInfo baBondsPrev = baBondsBackCis.clone();

        final double dDihedralIncrement = 10.0*Math.PI/180.0;

        SidechainLoop:
        for(int iSide = 0; iSide < iNoOfSidechains; iSide++){

            final ZMatrix zmatCurrent = alSidechains.get(iSide).returnZMatrixCopy();
            final boolean[][] baCurrBonds = alSidechains.get(iSide).returnBondingCopy();
            //reduced by two since the two dummies (one in cartes, one in zmat) will be gone
            final int iCurrNoOfAts = caPrev.getNoOfAtoms() + zmatCurrent.getNoOfAtoms()-2;

            final int[] iaCurrAtsPerMol = {iCurrNoOfAts};

            // the dummy
            final int iCurrConnect = iaConnectsCis[iSide];

            // which atom the dummy is connected to
            int iHostAtom = -1;
            for(int i = 0; i < baBondsPrev.getNoOfAtoms(); i++){
                if(baBondsPrev.hasBond(iCurrConnect, i) && i != iCurrConnect){
                    iHostAtom = i;
                    break;
                }
            }

            if(iHostAtom == -1){
                // this is not good
                System.err.println("ERROR: No connected atom found to dummy. This is a problem.");
                alCartes.add(null);
                alCartes.add(null);

                return alCartes;
            }

            int iDihedralAt = -1;
            for(int i = 0; i < baBondsPrev.getNoOfAtoms(); i++){
                if(baBondsPrev.hasBond(iHostAtom, i) && i != iHostAtom && i != iCurrConnect){
                    iDihedralAt = i;
                    break;
                }
            }

            if(iDihedralAt == -1){
                // this is also not good
                System.err.println("ERROR: No connected atom found for creating dihedral. This is a problem.");
                alCartes.add(null);
                alCartes.add(null);

                return alCartes;
            }
            
            // dynamically adjust the bonding information
            final BondInfo baBonds = createNewBonding(baBondsPrev, baCurrBonds,
                    iCurrNoOfAts, iCurrConnect);

            if(bDebug){
                System.out.println("DEBUG: Dynamically created bonding info coming.");
                for(int i = 0; i < iCurrNoOfAts; i++){
                    String sTemp = " ";
                    for(int j = 0; j < iCurrNoOfAts; j++){
                        sTemp += baBonds.bondType(i, j) + " ";
                    }
                    System.out.println("DEBUG: " + sTemp);
                }
            }

            caNext = new CartesianCoordinates(iCurrNoOfAts, 1, iaCurrAtsPerMol);

            boolean bSuccess = false;
            double dDihedral = 0.0;

            final double[][] daNextXYZ = caNext.getAllXYZCoordsCopy();
            final double[][] daPrevXYZ = caPrev.getAllXYZCoordsCopy();
            for(int i = 0; i < 3; i++){
                System.arraycopy(daPrevXYZ[i], 0, daNextXYZ[i], 0, daPrevXYZ[0].length);
            }

            // put it back in
            caNext.setAllXYZ(daNextXYZ);

            final String[] saNextAts = caNext.getAllAtomTypes();
            final String[] saPrevAts = caPrev.getAllAtomTypes();
            System.arraycopy(saPrevAts, 0, saNextAts, 0, saPrevAts.length);

            // now put the second atom of the sidechain into the correct spot
            final String[] saCurrSide = zmatCurrent.getAllAtomNames();
            saNextAts[iCurrConnect] = saCurrSide[1];

            // now the rest
            for(int i = 2; i < saCurrSide.length; i++){
                saNextAts[i + saPrevAts.length - 2] = saCurrSide[i];
            }

            final double[] daCoordsDummy = caPrev.getXYZCoordinatesOfAtom(iCurrConnect);
            final double[] daCoordsHost = caPrev.getXYZCoordinatesOfAtom(iHostAtom);
            final double[] daCoordsDihedral = caPrev.getXYZCoordinatesOfAtom(iDihedralAt);

            final int[] iaBondConnects = zmatCurrent.getAllBondConnects();
            final int[] iaAngleConnects = zmatCurrent.getAllAnglesConnects();
            final int[] iaDihedralConnects = zmatCurrent.getAllDihedralConnects();

            final double[] daBondLengths = zmatCurrent.getAllBondLengths();
            final double[] daBondAngles = zmatCurrent.getAllBondAngles();
            final double[] daDihedrals = zmatCurrent.getAllDihedrals();

            final int iOffset = caPrev.getNoOfAtoms();

            // of course, this just applies if we have more than one atom
            if(daBondAngles.length > 2){
                DegreesLoop:
                while (!bSuccess && dDihedral < (360.0 * Math.PI / 180.0)) {
                    /*
                     * populate the cartesian using the old one and the z-Matrix
                     * 1) the third atom and add the dihedral there
                     */
                    final double[] daFirstReal = giveCoordsWithDihedral(daCoordsDummy,
                            daCoordsHost, daCoordsDihedral, daBondLengths[2],
                            daBondAngles[2], dDihedral);

                    if(bDebug){
                        System.out.println("DEBUG: Tried to rotate the atom with " + dDihedral + " at atoms "
                                + iCurrConnect + " " + iHostAtom + " " + iDihedralAt);
                    }

                    // put it into the right spot
                    daNextXYZ[0][iOffset] = daFirstReal[0];
                    daNextXYZ[1][iOffset] = daFirstReal[1];
                    daNextXYZ[2][iOffset] = daFirstReal[2];

                    // loop over the rest
                    for (int i = 3; i < daBondLengths.length; i++) {

                        int iTmpBondConnect;
                        if(iaBondConnects[i] == 0){
                            iTmpBondConnect = iHostAtom;
                        } else if(iaBondConnects[i] == 1){
                            iTmpBondConnect = iCurrConnect;
                        } else {
                            iTmpBondConnect = iaBondConnects[i] + iOffset - 2;
                        }

                        int iTmpAngleConnect;
                        if(iaAngleConnects[i] == 0){
                            iTmpAngleConnect = iHostAtom;
                        } else if(iaAngleConnects[i] == 1){
                            iTmpAngleConnect = iCurrConnect;
                        } else {
                            iTmpAngleConnect = iaAngleConnects[i] + iOffset - 2;
                        }

                        int iTmpDihedralConnect;
                        if(iaDihedralConnects[i] == 0){
                            iTmpDihedralConnect = iHostAtom;
                        } else if(iaDihedralConnects[i] == 1){
                            iTmpDihedralConnect = iCurrConnect;
                        } else {
                            iTmpDihedralConnect = iaDihedralConnects[i] + iOffset - 2;
                        }


                        final double[] daXYZBond = new double[3];
                        final double[] daXYZAngle = new double[3];
                        final double[] daXYZDihedral = new double[3];

                        daXYZBond[0] = daNextXYZ[0][iTmpBondConnect];
                        daXYZBond[1] = daNextXYZ[1][iTmpBondConnect];
                        daXYZBond[2] = daNextXYZ[2][iTmpBondConnect];

                        daXYZAngle[0] = daNextXYZ[0][iTmpAngleConnect];
                        daXYZAngle[1] = daNextXYZ[1][iTmpAngleConnect];
                        daXYZAngle[2] = daNextXYZ[2][iTmpAngleConnect];

                        daXYZDihedral[0] = daNextXYZ[0][iTmpDihedralConnect];
                        daXYZDihedral[1] = daNextXYZ[1][iTmpDihedralConnect];
                        daXYZDihedral[2] = daNextXYZ[2][iTmpDihedralConnect];

                        final double[] daTempXYZ = giveCoordsWithDihedral(daXYZBond,
                                daXYZAngle, daXYZDihedral, daBondLengths[i],
                                daBondAngles[i], daDihedrals[i]);

                        // put it into the right spot
                        daNextXYZ[0][iOffset + i - 2] = daTempXYZ[0];
                        daNextXYZ[1][iOffset + i - 2] = daTempXYZ[1];
                        daNextXYZ[2][iOffset + i - 2] = daTempXYZ[2];
                    }

                    // check for collisions
                    caNext.setAllXYZ(daNextXYZ);

                    if(bDebug){
                        System.out.println("DEBUG: We are at " + dDihedral + " degrees. Cartesian coming.");
                        final String[] saCartes = caNext.createPrintableCartesians();
                        for(final String s : saCartes){
                            System.out.println(s);
                        }
                    }

                    final CollisionInfo  collInfo = collDetect.checkForCollision(caNext, dBlowBondDetect, baBonds);
                    final boolean bColl = collInfo.hasCollision();

                    if (!bColl) {
                        bSuccess = true;
                        break DegreesLoop;
                    } else if (bColl && bDebug){
                        final List<Collision> colls = collInfo.getCollisions();
                        String st = "";
                        for(final Collision coll : colls){
                            st += "  " + coll.getAtomOne() + "/" + coll.getAtomTwo();
                        }
                        System.out.println("DEBUG: Collision found at " + st);
                    }
                    dDihedral += dDihedralIncrement;
                }
                if (!bSuccess) {
                    // react on this... meaning: return null'd cartesians
                    if(bDebug){
                        System.out.println("DEBUG: Didn't manage to glue. My last attempt comes");
                        final String[] saDebugStuff = caNext.createPrintableCartesians();
                        for(final String sTemp: saDebugStuff){
                            System.out.println(sTemp);
                        }
                    }
                    alCartes.add(null);
                    alCartes.add(null);

                    return alCartes;
                }
            }

            caPrev = caNext;
            baBondsPrev = baBonds;
        }


        // if we reach here, then we must have one valid cartesian set
        final CartesianCoordinates caFinalCis = new CartesianCoordinates(caPrev);
        alCartes.add(caFinalCis);

        /*
         * now we handle the trans part
         */
        final CartesianCoordinates caBackTrans = back.getCartesCopy(false);
        final int[] iaConnectsTrans = back.getConnectsCopy(false);
        final BondInfo baBondsBackTrans = back.getBondsCopy(false);

        caPrev = new CartesianCoordinates(caBackTrans);

        baBondsPrev = baBondsBackTrans.clone();

        SidechainLoop:
        for(int iSide = 0; iSide < iNoOfSidechains; iSide++){

            final ZMatrix zmatCurrent = alSidechains.get(iSide).returnZMatrixCopy();
            final boolean[][] baCurrBonds = alSidechains.get(iSide).returnBondingCopy();
            //reduced by two since the two dummies (one in cartes, one in zmat) will be gone
            final int iCurrNoOfAts = caPrev.getNoOfAtoms() + zmatCurrent.getNoOfAtoms()-2;

            final int[] iaCurrAtsPerMol = {iCurrNoOfAts};

            // the dummy
            final int iCurrConnect = iaConnectsTrans[iSide];

            // which atom the dummy is connected to
            int iHostAtom = -1;
            for(int i = 0; i < baBondsPrev.getNoOfAtoms(); i++){
                if(baBondsPrev.hasBond(iCurrConnect, i) && i != iCurrConnect){
                    iHostAtom = i;
                    break;
                }
            }

            if(iHostAtom == -1){
                // this is not good
                System.err.println("ERROR: No connected atom found to dummy. This is a problem.");
                alCartes.add(null);
                alCartes.add(null);

                return alCartes;
            }

            int iDihedralAt = -1;
            for(int i = 0; i < baBondsPrev.getNoOfAtoms(); i++){
                if(baBondsPrev.hasBond(iHostAtom, i) && i != iHostAtom && i != iCurrConnect){
                    iDihedralAt = i;
                    break;
                }
            }

            if(iDihedralAt == -1){
                // this is also not good
                System.err.println("ERROR: No connected atom found for creating dihedral. This is a problem.");
                alCartes.add(null);
                alCartes.add(null);

                return alCartes;
            }

            // dynamically adjust the bonding information
            final BondInfo bonds = createNewBonding(baBondsPrev, baCurrBonds,
                    iCurrNoOfAts, iCurrConnect);

            caNext = new CartesianCoordinates(iCurrNoOfAts, 1, iaCurrAtsPerMol);

            boolean bSuccess = false;
            double dDihedral = 0.0;

            final double[][] daNextXYZ = caNext.getAllXYZCoordsCopy();
            final double[][] daPrevXYZ = caPrev.getAllXYZCoordsCopy();
            for(int i = 0; i < 3; i++){
                System.arraycopy(daPrevXYZ[i], 0, daNextXYZ[i], 0, daPrevXYZ[0].length);
            }

            // put it back in
            caNext.setAllXYZ(daNextXYZ);

            final String[] saNextAts = caNext.getAllAtomTypes();
            final String[] saPrevAts = caPrev.getAllAtomTypes();
            System.arraycopy(saPrevAts, 0, saNextAts, 0, saPrevAts.length);

            // now put the second atom of the sidechain into the correct spot
            final String[] saCurrSide = zmatCurrent.getAllAtomNames();
            saNextAts[iCurrConnect] = saCurrSide[1];

            // now the rest
            for(int i = 2; i < saCurrSide.length; i++){
                saNextAts[i + saPrevAts.length - 2] = saCurrSide[i];
            }

            final double[] daCoordsDummy = caPrev.getXYZCoordinatesOfAtom(iCurrConnect);
            final double[] daCoordsHost = caPrev.getXYZCoordinatesOfAtom(iHostAtom);
            final double[] daCoordsDihedral = caPrev.getXYZCoordinatesOfAtom(iDihedralAt);

            final int[] iaBondConnects = zmatCurrent.getAllBondConnects();
            final int[] iaAngleConnects = zmatCurrent.getAllAnglesConnects();
            final int[] iaDihedralConnects = zmatCurrent.getAllDihedralConnects();

            final double[] daBondLengths = zmatCurrent.getAllBondLengths();
            final double[] daBondAngles = zmatCurrent.getAllBondAngles();
            final double[] daDihedrals = zmatCurrent.getAllDihedrals();

            final int iOffset = caPrev.getNoOfAtoms();

            // of course, this just applies if we have more than one atom
            if(daBondAngles.length > 2){
                DegreesLoop:
                while (!bSuccess && dDihedral < (360.0 * Math.PI / 180.0)) {
                    /*
                     * populate the cartesian using the old one and the z-Matrix
                     * 1) the third atom and add the dihedral there
                     */
                    final double[] daFirstReal = giveCoordsWithDihedral(daCoordsDummy,
                            daCoordsHost, daCoordsDihedral, daBondLengths[2],
                            daBondAngles[2], dDihedral);

                    // put it into the right spot
                    daNextXYZ[0][iOffset] = daFirstReal[0];
                    daNextXYZ[1][iOffset] = daFirstReal[1];
                    daNextXYZ[2][iOffset] = daFirstReal[2];

                    // loop over the rest
                    for (int i = 3; i < daBondLengths.length; i++) {

                        int iTmpBondConnect;
                        if(iaBondConnects[i] == 0){
                            iTmpBondConnect = iHostAtom;
                        } else if(iaBondConnects[i] == 1){
                            iTmpBondConnect = iCurrConnect;
                        } else {
                            iTmpBondConnect = iaBondConnects[i] + iOffset - 2;
                        }

                        int iTmpAngleConnect;
                        if(iaAngleConnects[i] == 0){
                            iTmpAngleConnect = iHostAtom;
                        } else if(iaAngleConnects[i] == 1){
                            iTmpAngleConnect = iCurrConnect;
                        } else {
                            iTmpAngleConnect = iaAngleConnects[i] + iOffset - 2;
                        }

                        int iTmpDihedralConnect;
                        if(iaDihedralConnects[i] == 0){
                            iTmpDihedralConnect = iHostAtom;
                        } else if(iaDihedralConnects[i] == 1){
                            iTmpDihedralConnect = iCurrConnect;
                        } else {
                            iTmpDihedralConnect = iaDihedralConnects[i] + iOffset - 2;
                        }


                        final double[] daXYZBond = new double[3];
                        final double[] daXYZAngle = new double[3];
                        final double[] daXYZDihedral = new double[3];

                        daXYZBond[0] = daNextXYZ[0][iTmpBondConnect];
                        daXYZBond[1] = daNextXYZ[1][iTmpBondConnect];
                        daXYZBond[2] = daNextXYZ[2][iTmpBondConnect];

                        daXYZAngle[0] = daNextXYZ[0][iTmpAngleConnect];
                        daXYZAngle[1] = daNextXYZ[1][iTmpAngleConnect];
                        daXYZAngle[2] = daNextXYZ[2][iTmpAngleConnect];

                        daXYZDihedral[0] = daNextXYZ[0][iTmpDihedralConnect];
                        daXYZDihedral[1] = daNextXYZ[1][iTmpDihedralConnect];
                        daXYZDihedral[2] = daNextXYZ[2][iTmpDihedralConnect];

                        final double[] daTempXYZ = giveCoordsWithDihedral(daXYZBond,
                                daXYZAngle, daXYZDihedral, daBondLengths[i],
                                daBondAngles[i], daDihedrals[i]);

                        // put it into the right spot
                        daNextXYZ[0][iOffset + i - 2] = daTempXYZ[0];
                        daNextXYZ[1][iOffset + i - 2] = daTempXYZ[1];
                        daNextXYZ[2][iOffset + i - 2] = daTempXYZ[2];
                    }

                    // check for collisions
                    caNext.setAllXYZ(daNextXYZ);
                    final boolean bColl = collDetect.checkForCollision(caNext, dBlowBondDetect, bonds).hasCollision();

                    if (!bColl) {
                        bSuccess = true;
                        break DegreesLoop;
                    }
                    dDihedral += dDihedralIncrement;
                }
                if (!bSuccess) {
                    // react on this... meaning: return null'd cartesians
                    if(bDebug){
                        System.out.println("DEBUG: Didn't manage to glue. My last attempt comes. Means one went through fine.");
                        final String[] saDebugStuff = caNext.createPrintableCartesians();
                        for(final String sTemp: saDebugStuff){
                            System.out.println(sTemp);
                        }
                    }
                    alCartes.add(null);

                    return alCartes;
                }
            }
            caPrev = caNext;
            baBondsPrev = bonds;
        }


        // if we reach here, then we must have a valid cartesian set
        final CartesianCoordinates caFinalTrans = new CartesianCoordinates(caPrev);
        alCartes.add(caFinalTrans);


        return alCartes;
    }
    
    private static BondInfo createNewBonding(final BondInfo original,
            final boolean[][] sidechain, final int totNoOfAts,
            final int whichDummy){
        
        if(false){
            System.out.println("DEBUG: Creating new bonding with " + original.getNoOfAtoms() + "  " + sidechain.length
                    + "  " + totNoOfAts + "  " + whichDummy);
        }

        if(totNoOfAts != (original.getNoOfAtoms() + sidechain.length -2)){
            System.err.println("ERROR: Non-equal number of atoms: " + totNoOfAts
                    + " and " + (original.getNoOfAtoms() + sidechain.length -2)
                    + ". Returning null. Contact the author(s).");
            return null;
        }

        final BondInfo newBonds = new SimpleBondInfo(totNoOfAts);

        // put original in
        for(int i = 0; i < original.getNoOfAtoms(); i++){
            for(int j = 0; j < original.getNoOfAtoms(); j++){
                newBonds.setBond(i, j, original.bondType(i, j));
            }
        }

        // the offset for the sidechain
        final int offset = original.getNoOfAtoms();

        // copy the bonding of the sidechain except for the dummy of the sidechain and the
        // replacement for the dummy at the backbone, therefore it starts always with 
        // number two: 0) XX, 1) dummy replacement
        for(int i = 2; i < sidechain.length; i++){
            for(int j = 2; j < sidechain.length; j++){
                if(false) {System.out.println("DEBUG: working on " + i + "  " + j);}
                final short bond = (sidechain[i][j]) ? BondInfo.UNCERTAIN : BondInfo.NOBOND;
                newBonds.setBond(i+offset-2, j+offset-2, bond);
            }
        }

        // adjust the former dummy and the bonds towards it
        for(int i = offset; i < newBonds.getNoOfAtoms(); i++){
            final short bond = (sidechain[1][i-offset+2]) ? BondInfo.UNCERTAIN : BondInfo.NOBOND;
            newBonds.setBond(i, whichDummy, bond);
        }

        return newBonds;
    }

    private static double[] giveCoordsWithDihedral(final double[] daXYZFirst,
            final double[] daXYZSecond, final double[] daXYZThird,
            final double dBondLength, final double dAngle, final double dDihedral){

        final double[] daAVec = daXYZFirst;
        final double[] daBVec = daXYZSecond;
        final double[] daCVec = daXYZThird;

        final double[] daVector1 = new double[3];
        final double[] daVector2 = new double[3];

        for (int j = 0; j < 3; j++) {
            daVector1[j] = daAVec[j] - daBVec[j];
            daVector2[j] = daAVec[j] - daCVec[j];
        }

        final double[] daNVec = new double[3];
        final double[] daNNVec = new double[3];

        daNVec[0] = daVector1[1] * daVector2[2] - daVector1[2] * daVector2[1];
        daNVec[1] = daVector1[2] * daVector2[0] - daVector1[0] * daVector2[2];
        daNVec[2] = daVector1[0] * daVector2[1] - daVector1[1] * daVector2[0];

        daNNVec[0] = daVector1[1] * daNVec[2] - daVector1[2] * daNVec[1];
        daNNVec[1] = daVector1[2] * daNVec[0] - daVector1[0] * daNVec[2];
        daNNVec[2] = daVector1[0] * daNVec[1] - daVector1[1] * daNVec[0];

        // normalize the new vectors
        final double dNormalN = 1 / Math.sqrt(Math.pow(daNVec[0], 2) + Math.pow(daNVec[1], 2) + Math.pow(daNVec[2], 2));
        final double dNormalNN = 1 / Math.sqrt(Math.pow(daNNVec[0], 2) + Math.pow(daNNVec[1], 2) + Math.pow(daNNVec[2], 2));
        for (int j = 0; j < 3; j++) {
            daNVec[j] *= dNormalN;
            daNNVec[j] *= dNormalNN;
        }

        // manipulate them with the dihedral
        final double dTempN = -Math.sin(dDihedral);
        final double dTempNN = Math.cos(dDihedral);

        for (int j = 0; j < 3; j++) {
            daNVec[j] *= dTempN;
            daNNVec[j] *= dTempNN;
        }

        //create the next vector
        final double[] daVector3 = new double[3];
        for (int j = 0; j < 3; j++) {
            daVector3[j] = daNVec[j] + daNNVec[j];
        }

        // normalize that one as well
        final double dTemp3 = 1 / Math.sqrt(Math.pow(daVector3[0], 2) + Math.pow(daVector3[1], 2) + Math.pow(daVector3[2], 2));
        for (int j = 0; j < 3; j++) {
            daVector3[j] *= dTemp3;
        }

        // manipulate it with the bond length and angle
        for (int j = 0; j < 3; j++) {
            daVector3[j] *= dBondLength * Math.sin(dAngle);
        }

        // normalize the first vector
        final double dTemp1 = 1 / Math.sqrt(Math.pow(daVector1[0], 2) + Math.pow(daVector1[1], 2) + Math.pow(daVector1[2], 2));
        for (int j = 0; j < 3; j++) {
            daVector1[j] *= dTemp1;
        }

        // manipulate it with angle and length
        for (int j = 0; j < 3; j++) {
            daVector1[j] *= dBondLength * Math.cos(dAngle);
        }

        // calculate the actual postion
        final double[] daXYZ = new double[3];
        for (int j = 0; j < 3; j++) {
            daXYZ[j] = daAVec[j] + daVector3[j] - daVector1[j];
        }
        
        return daXYZ;
    }
}
