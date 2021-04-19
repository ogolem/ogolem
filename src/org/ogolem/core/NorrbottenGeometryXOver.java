/*
Copyright (c) 2014, J. M. Dieterich
              2016, J. M. Dieterich and B. Hartke
              2018-2020, J. M. Dieterich and B. Hartke
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

import static org.ogolem.core.GlobOptAtomics.assignMolecularTypes;
import static org.ogolem.core.GlobOptAtomics.findOptimalXPlaneHeight;
import static org.ogolem.core.GlobOptAtomics.findOptimalYPlaneHeight;
import static org.ogolem.core.GlobOptAtomics.findOptimalZPlaneHeight;
import static org.ogolem.core.GlobOptAtomics.randomCuttingXPlane;
import static org.ogolem.core.GlobOptAtomics.randomCuttingYPlane;
import static org.ogolem.core.GlobOptAtomics.randomCuttingZPlane;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.Lottery;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Our basic N-type phenotype X-Over in a user-specified plane. Rotates the input structures.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class NorrbottenGeometryXOver implements GenericCrossover<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20160403;
  private static final boolean DEBUG = false;
  private static final Logger log = LoggerFactory.getLogger(NorrbottenGeometryXOver.class);

  public static enum ROTATIONPLANE {
    XY,
    XZ,
    YZ
  };

  private final Lottery random = Lottery.getInstance();
  private final GlobOptAtomics.CUTTINGMODE whichGlobOpt;

  private final ROTATIONPLANE plane;

  NorrbottenGeometryXOver(
      final GlobOptAtomics.CUTTINGMODE whichGlobOpt, final ROTATIONPLANE plane) {
    assert (whichGlobOpt != null);
    assert (plane != null);

    this.whichGlobOpt = whichGlobOpt;
    this.plane = plane;
  }

  NorrbottenGeometryXOver(final NorrbottenGeometryXOver orig) {
    this.whichGlobOpt = orig.whichGlobOpt;
    this.plane = orig.plane;
  }

  @Override
  public NorrbottenGeometryXOver copy() {
    return new NorrbottenGeometryXOver(this);
  }

  @Override
  public String getMyID() {
    return "NORRBOTTEN\nphenotype X-Over\n\tglobopt style: "
        + whichGlobOpt.toString()
        + "\n\trotation plane: "
        + plane.name();
  }

  @Override
  public Tuple<Geometry, Geometry> crossover(
      final Geometry mother, final Geometry father, final long futureID) {

    final Geometry gChild1 = new Geometry(father);
    gChild1.setID(futureID);
    final Geometry gChild2 = new Geometry(mother);
    gChild2.setID(futureID);

    final CartesianCoordinates cartesMother = gChild2.getCartesians();
    final CartesianCoordinates cartesFather = gChild1.getCartesians();

    cartesMother.moveCoordsToCOM();
    cartesFather.moveCoordsToCOM();

    // rotate the two
    final double rot1 = 2 * Math.PI * random.nextDouble();
    final double rot2 = 2 * Math.PI * random.nextDouble();

    switch (plane) {
      case XY:
        CoordTranslation.rotateXYZAroundZ(
            cartesMother.getAllXYZCoordsCopy(), rot1, cartesMother.getAllXYZCoord());
        CoordTranslation.rotateXYZAroundZ(
            cartesFather.getAllXYZCoordsCopy(), rot2, cartesFather.getAllXYZCoord());
        break;
      case XZ:
        CoordTranslation.rotateXYZAroundY(
            cartesMother.getAllXYZCoordsCopy(), rot1, cartesMother.getAllXYZCoord());
        CoordTranslation.rotateXYZAroundY(
            cartesFather.getAllXYZCoordsCopy(), rot2, cartesFather.getAllXYZCoord());
        break;
      case YZ:
        CoordTranslation.rotateXYZAroundX(
            cartesMother.getAllXYZCoordsCopy(), rot1, cartesMother.getAllXYZCoord());
        CoordTranslation.rotateXYZAroundX(
            cartesFather.getAllXYZCoordsCopy(), rot2, cartesFather.getAllXYZCoord());
        break;
    }

    if (DEBUG) {
      final String[] printCart1 = cartesMother.createPrintableCartesians();
      final String[] printCart2 = cartesFather.createPrintableCartesians();
      System.out.println("DEBUG: Cartesians after rotation: mother");
      for (final String s : printCart1) {
        System.out.println(s);
      }
      System.out.println("DEBUG: Cartesians after rotation: father");
      for (final String s : printCart2) {
        System.out.println(s);
      }
    }

    // translate it to geometries again
    gChild1.updateAllCoordinates(cartesFather.getAllXYZCoord());
    gChild2.updateAllCoordinates(cartesMother.getAllXYZCoord());

    // set the IDs of mother and father
    gChild1.setMother(mother.getID());
    gChild1.setFather(father.getID());
    gChild2.setMother(mother.getID());
    gChild2.setFather(father.getID());

    // in case of environments, add it in at this point
    if (mother.containsEnvironment()) {
      // the environment
      final List<Environment> envChilds =
          mother.getEnvironment().createOffspring(father.getEnvironment());
      gChild1.setEnvironment(envChilds.get(0));
      gChild2.setEnvironment(envChilds.get(1));
    }

    // first get all COMs and put them into matrices/double arrays
    final double[][] daMotherCOMs = new double[3][gChild1.getNumberOfIndieParticles()];
    final double[][] daFatherCOMs = new double[3][gChild1.getNumberOfIndieParticles()];

    for (int i = 0; i < gChild1.getNumberOfIndieParticles(); i++) {
      final double[] daTemp1 = gChild2.getMoleculeAtPosition(i).getExternalCenterOfMass();
      final double[] daTemp2 = gChild1.getMoleculeAtPosition(i).getExternalCenterOfMass();
      for (int j = 0; j < 3; j++) {
        daMotherCOMs[j][i] = daTemp1[j];
        daFatherCOMs[j][i] = daTemp2[j];
      }
    }

    /*
     * now we want to take a plane cutting through both
     * 1) which height on the z axis?
     */
    final int noMols = gChild1.getNumberOfIndieParticles();
    final List<Integer> fatherAbove = new ArrayList<>(noMols);
    final List<Integer> fatherUnder = new ArrayList<>(noMols);

    double dPlaneHeightFather = 0.0;
    short lookAtCoord = 0;
    int ec = 0;
    boolean cont;
    int noFatherUnder = -1;
    do {
      switch (plane) {
        case XY:
          dPlaneHeightFather = randomCuttingXPlane(whichGlobOpt, daFatherCOMs, random);
          lookAtCoord = 0;
          break;
        case XZ:
          dPlaneHeightFather = randomCuttingZPlane(whichGlobOpt, daFatherCOMs, random);
          lookAtCoord = 2;
          break;
        case YZ:
          dPlaneHeightFather = randomCuttingYPlane(whichGlobOpt, daFatherCOMs, random);
          lookAtCoord = 1;
          break;
      }
      // check which COMs are above and underneath the plane (for the father)
      for (int i = 0; i < noMols; i++) {
        if (daFatherCOMs[lookAtCoord][i] <= dPlaneHeightFather) {
          // "underneath"
          fatherUnder.add(i);
        } else {
          // "above"
          fatherAbove.add(i);
        }
      }

      noFatherUnder = fatherUnder.size();
      if (noFatherUnder == 0) {
        // then there was none underneath the plane, this is NOT ok.
        cont = true;
        fatherUnder.clear();
        fatherAbove.clear();
      } else if (noFatherUnder >= noMols) {
        // all are above the plane, also not OK.
        cont = true;
        fatherUnder.clear();
        fatherAbove.clear();
      } else {
        cont = false;
      }

      ec++;
    } while (cont && ec < FixedValues.MAXTOEMERGENCY);

    if (cont) {
      System.err.println("Too many attemps in finding a first plane height!");
      return new Tuple<>(null, null);
    }

    // now move the plane in the mother geometry till enough COMs are above/underneath
    double dPlaneHeightMother = 0.0;
    switch (plane) {
      case XY:
        dPlaneHeightMother = findOptimalXPlaneHeight(daMotherCOMs, noFatherUnder);
        break;
      case XZ:
        dPlaneHeightMother = findOptimalZPlaneHeight(daMotherCOMs, noFatherUnder);
        break;
      case YZ:
        dPlaneHeightMother = findOptimalYPlaneHeight(daMotherCOMs, noFatherUnder);
        break;
    }

    if (Double.isNaN(dPlaneHeightMother)) {
      // two COMs were too close together, return null'd geometries
      System.err.println("WARNING: Couldn't find cutting plane for mother, two COMs too close.");
      return new Tuple<>(null, null);
    }

    if (DEBUG) {
      System.out.println(
          "DEBUG: Plane height father: "
              + dPlaneHeightFather
              + " and mother "
              + dPlaneHeightMother);
      System.out.println("DEBUG Father coming .... ");
      for (final String s : gChild1.makePrintableAbsoluteCoord(true)) {
        System.out.println(s);
      }
      System.out.println("DEBUG Mother coming .... ");
      for (final String s : gChild2.makePrintableAbsoluteCoord(true)) {
        System.out.println(s);
      }
    }

    final List<Integer> motherUnder = new ArrayList<>(noMols);
    final List<Integer> motherAbove = new ArrayList<>(noMols);

    for (int ii = 0; ii < noMols; ii++) {
      if (daMotherCOMs[lookAtCoord][ii] <= dPlaneHeightMother) {
        // "underneath"
        motherUnder.add(ii);
      } else {
        // "above"
        motherAbove.add(ii);
      }
    }

    if (noFatherUnder != motherUnder.size()
        || motherAbove.size() != fatherAbove.size()
        || fatherUnder.size() != motherUnder.size()) {
      System.err.println(
          "DEBUG: Plane height is bad. Returning. Mismatch: "
              + noFatherUnder
              + " vs "
              + motherUnder.size()
              + " or "
              + motherAbove.size()
              + " vs "
              + fatherAbove.size());
      return new Tuple<>(null, null);
    }

    /*
     * now we want to see which molecular types exist above and underneath for both geometries
     * first: the father
     */
    final ArrayList<Molecule> alFatherMolecules = new ArrayList<>(noMols);

    for (int i = 0; i < noMols; i++) {
      alFatherMolecules.add(gChild2.getMoleculeAtPosition(i));
    }

    final int[] iaMolTypesFather = assignMolecularTypes(alFatherMolecules);

    // now we want to figure out how many of each are above and underneath

    // find the "highest" molecular type
    int iHighest = 0;
    for (int i = 0; i < noMols; i++) {
      iHighest = Math.max(iHighest, iaMolTypesFather[i]);
    }

    final int[] iaTotTypeCounter = new int[iHighest + 1];
    final int[] iaFatherAboveCounter = new int[iHighest + 1];
    final int[] iaMotherAboveCounter = new int[iHighest + 1];

    for (int i = 0; i < noMols; i++) {
      // total
      iaTotTypeCounter[iaMolTypesFather[i]]++;
    }

    for (int i = 0; i < fatherAbove.size(); i++) {
      // father
      iaFatherAboveCounter[iaMolTypesFather[fatherAbove.get(i)]]++;

      // mother, has p.d. the same overall molecular types as father
      iaMotherAboveCounter[iaMolTypesFather[motherAbove.get(i)]]++;
    }

    final int[] iaDiscrepancy = new int[iaFatherAboveCounter.length];
    for (int i = 0; i < iaDiscrepancy.length; i++) {
      /*
       * this means, that if both are same the result is 0, if mother has more it is negative
       * and if father has more it is positive
       */
      iaDiscrepancy[i] = iaFatherAboveCounter[i] - iaMotherAboveCounter[i];
      if (DEBUG) {
        System.out.println("DEBUG: Discrepancy " + i + "\t" + iaDiscrepancy[i]);
      }
    }

    for (int i = 0; i < iaDiscrepancy.length; i++) {
      if (iaDiscrepancy[i] == 0) {
        continue;
      } else if (iaDiscrepancy[i] > 0) {
        // from the mother needs to come up
        while (iaDiscrepancy[i] > 0) {
          // which one of the under goes up?
          final int iWhichUp = random.nextInt(iaTotTypeCounter[i] - iaMotherAboveCounter[i]);

          // which type to change with?
          int iNoOfPoss = 0;
          for (int iDis : iaDiscrepancy) {
            if (iDis < 0) {
              iNoOfPoss++;
            }
          }

          final int iWhichKindDown = random.nextInt(iNoOfPoss);

          // map the random to the type
          int iCounter = -1;
          int iType = -1;
          for (int j = 0; j < iaDiscrepancy.length; j++) {
            if (iaDiscrepancy[j] < 0) {
              iCounter++;
            }

            if (iCounter == iWhichKindDown) {
              iType = j;
              break;
            }
          }

          // which molecule of type iType?
          final int iWhichMolOfType = random.nextInt(iaMotherAboveCounter[iType]);

          // now we need to have the ID for the one coming up
          int iMolIDUp = -1;
          int iMolCounter = -1;
          for (int iMol = 0; iMol < iaMolTypesFather.length; iMol++) {
            if (iaMolTypesFather[iMol] == i && !motherAbove.contains(iMol)) {
              iMolCounter++;
            }

            if (iMolCounter == iWhichUp) {
              iMolIDUp = iMol;
              break;
            }
          }

          // and for the one going down
          int iMolIDDown = -1;
          iMolCounter = -1;
          for (int iMol = 0; iMol < iaMolTypesFather.length; iMol++) {
            if (iaMolTypesFather[iMol] == iType && motherAbove.contains(iMol)) {
              iMolCounter++;
            }

            if (iMolCounter == iWhichMolOfType) {
              iMolIDDown = iMol;
              break;
            }
          }

          // now exchange the coordinates and angles of the two
          final double[] daTempCOM =
              gChild2.getMoleculeAtPosition(iMolIDUp).getExternalCenterOfMass().clone();
          final double[] daTempOrient =
              gChild2.getMoleculeAtPosition(iMolIDUp).getOrientation().clone();

          final double[] daTempCOM2 =
              gChild2.getMoleculeAtPosition(iMolIDDown).getExternalCenterOfMass().clone();
          final double[] daTempOrient2 =
              gChild2.getMoleculeAtPosition(iMolIDDown).getOrientation().clone();

          gChild2.getMoleculeAtPosition(iMolIDUp).setExternalCenterOfMass(daTempCOM2);
          gChild2.getMoleculeAtPosition(iMolIDUp).setOrientation(daTempOrient2);

          gChild2.getMoleculeAtPosition(iMolIDDown).setExternalCenterOfMass(daTempCOM);
          gChild2.getMoleculeAtPosition(iMolIDDown).setOrientation(daTempOrient);

          // adjust things a little
          iaDiscrepancy[i]--;
          iaDiscrepancy[iType]++;
          iaMotherAboveCounter[i]++;
          iaMotherAboveCounter[iType]--;

          // and also the arraylists
          motherAbove.remove(Integer.valueOf(iMolIDDown));
          motherAbove.add(iMolIDUp);
          motherUnder.remove(Integer.valueOf(iMolIDUp));
          motherUnder.add(iMolIDDown);
        }
      }

      /*
       * we do not need to explicitly handle the iaDiscrepancy < 0 case, the >0 does that for us :-)
       */
    }

    /*
     * a quick check whether everything is ok
     */
    for (int i : iaDiscrepancy) {
      if (i != 0) {
        System.out.println("ERROR: Discrepancy is still existing: " + i + ".");
      }
    }

    for (int i = 0; i < iaMotherAboveCounter.length; i++) {
      if (iaMotherAboveCounter[i] != iaFatherAboveCounter[i]) {
        System.out.println(
            "DEBUG: Mismatch father mother: "
                + iaMotherAboveCounter[i]
                + "   "
                + iaFatherAboveCounter[i]);
      }
    }

    /*
     * now do the real exchanging
     */

    /* still we need to exchange parts though
     * remember that we need to adjust with the cutting height (so
     * substracting it from the z component of the mother COMs)
     */
    // let's adjust the mother COMs again, since we swapped them round
    for (int i = 0; i < noMols; i++) {
      final double[] daTemp1 = gChild2.getMoleculeAtPosition(i).getExternalCenterOfMass();
      for (int j = 0; j < 3; j++) {
        daMotherCOMs[j][i] = daTemp1[j];
      }
    }

    // first always null the ArrayLists
    fatherUnder.clear();
    fatherAbove.clear();
    motherUnder.clear();
    motherAbove.clear();
    for (int kk = 0; kk < noMols; kk++) {
      if (daMotherCOMs[lookAtCoord][kk] <= dPlaneHeightMother) {
        // underneath
        motherUnder.add(kk);
      } else {
        // above
        motherAbove.add(kk);
      }
    }

    for (int kk = 0; kk < noMols; kk++) {
      if (daFatherCOMs[lookAtCoord][kk] <= dPlaneHeightFather) {
        // underneath
        fatherUnder.add(kk);
      } else {
        // above
        fatherAbove.add(kk);
      }
    }

    final int iMotherCOMs = motherUnder.size();
    if (iMotherCOMs != fatherUnder.size()) {
      // we f***** it up
      System.err.println(
          "ERROR: We f***** it up in Norrbotten. Mother: "
              + iMotherCOMs
              + " father: "
              + fatherUnder.size());
      return new Tuple<>(null, null);
    }

    // and again we need to figure out the molecular types above
    final int[] iaFatherAboveCounter2 = new int[iHighest + 1];
    final int[] iaMotherAboveCounter2 = new int[iHighest + 1];
    for (int i = 0; i < motherAbove.size(); i++) {
      iaFatherAboveCounter2[iaMolTypesFather[fatherAbove.get(i)]]++;
      // mother has p.d. the same overall molecular types as the father
      iaMotherAboveCounter2[iaMolTypesFather[motherAbove.get(i)]]++;
    }

    // check whether the exchange worked
    for (int i = 0; i < iaMotherAboveCounter2.length; i++) {
      if (iaMotherAboveCounter2[i] != iaFatherAboveCounter2[i]) {
        System.err.println(
            "ERROR: Exchange didn't work out. At position "
                + i
                + " mismatch "
                + iaMotherAboveCounter2[i]
                + " to "
                + iaFatherAboveCounter2[i]
                + ".");
        return new Tuple<>(null, null);
      }
    }

    // check whether we messed something up badly
    final int[] iaTypesChild1 = assignMolecularTypes(gChild1.getMolecules());
    final int[] iaTypesChild2 = assignMolecularTypes(gChild2.getMolecules());
    for (int i = 0; i < iaTypesChild1.length; i++) {
      if (iaTypesChild1[i] != iaMolTypesFather[i] || iaTypesChild2[i] != iaMolTypesFather[i]) {
        // darn, didn't work
        System.err.println(
            "ERROR: Problem in Sweden after internal exchange "
                + "at position "
                + i
                + " father "
                + iaMolTypesFather[i]
                + " child 1 "
                + iaTypesChild1[i]
                + " child 2 "
                + iaTypesChild2[i]);
        return new Tuple<>(null, null);
      }
    }

    if (DEBUG) {
      final GeometryConfig gFaUnder = new GeometryConfig();
      gFaUnder.noOfParticles = fatherUnder.size();
      gFaUnder.geomMCs = new ArrayList<>();
      for (final int under : fatherUnder) {
        final Molecule m = new Molecule(gChild1.getMoleculeAtPosition(under));
        gFaUnder.geomMCs.add(m.returnMyConfig());
      }

      final GeometryConfig gFaAbove = new GeometryConfig();
      gFaAbove.noOfParticles = fatherAbove.size();
      gFaAbove.geomMCs = new ArrayList<>();
      for (final int above : fatherAbove) {
        final Molecule m = new Molecule(gChild1.getMoleculeAtPosition(above));
        gFaAbove.geomMCs.add(m.returnMyConfig());
      }

      final GeometryConfig gMaUnder = new GeometryConfig();
      gMaUnder.noOfParticles = motherUnder.size();
      gMaUnder.geomMCs = new ArrayList<>();
      for (final int under : motherUnder) {
        final Molecule m = new Molecule(gChild2.getMoleculeAtPosition(under));
        gMaUnder.geomMCs.add(m.returnMyConfig());
      }

      final GeometryConfig gMaAbove = new GeometryConfig();
      gMaAbove.noOfParticles = motherAbove.size();
      gMaAbove.geomMCs = new ArrayList<>();
      for (final int above : motherAbove) {
        final Molecule m = new Molecule(gChild2.getMoleculeAtPosition(above));
        gMaAbove.geomMCs.add(m.returnMyConfig());
      }

      final Geometry g1 = new Geometry(gFaAbove);
      final Geometry g2 = new Geometry(gFaUnder);
      final Geometry g3 = new Geometry(gMaAbove);
      final Geometry g4 = new Geometry(gMaUnder);

      System.out.println("DEBUG: Father above");
      for (final String s : g1.makePrintableAbsoluteCoord(true)) {
        System.out.println(s);
      }
      System.out.println("DEBUG: Father under");
      for (final String s : g2.makePrintableAbsoluteCoord(true)) {
        System.out.println(s);
      }
      System.out.println("DEBUG: Mother above");
      for (final String s : g3.makePrintableAbsoluteCoord(true)) {
        System.out.println(s);
      }
      System.out.println("DEBUG: Mother under");
      for (final String s : g4.makePrintableAbsoluteCoord(true)) {
        System.out.println(s);
      }
    }

    // exchange the molecules but keep in mind that the order needs to
    // be conserved!
    final int[] iaMolTypesMother = iaMolTypesFather.clone();
    final List<Molecule> allMols1 = gChild1.getMolecules();
    final List<Molecule> allMols2 = gChild2.getMolecules();
    for (final int whichMol : fatherAbove) {
      final int iWhichType = iaMolTypesFather[whichMol];
      for (int j = 0; j < motherAbove.size(); j++) {
        final int iTempMol = motherAbove.get(j);
        final int iTempType = iaMolTypesMother[iTempMol];
        if (iTempType == iWhichType) {
          // exchange
          final Molecule m1 = allMols1.set(whichMol, null); // replace w/ null temporarily
          final Molecule m2 = allMols2.set(iTempMol, null);

          // adjust the COM heights
          final double[] daTempCOM = m1.getExternalCenterOfMass();
          final double[] daTempCOM2 = m2.getExternalCenterOfMass();

          daTempCOM[lookAtCoord] =
              daTempCOM[lookAtCoord] - dPlaneHeightMother + dPlaneHeightFather; // + heightIncr;
          daTempCOM2[lookAtCoord] =
              daTempCOM2[lookAtCoord] + dPlaneHeightMother - dPlaneHeightFather; // + heightIncr;

          gChild1.setMoleculeAtPosition(whichMol, m2);
          gChild2.setMoleculeAtPosition(iTempMol, m1);
          iaMolTypesMother[iTempMol] = -10000;

          break;
        }
        if (j == motherAbove.size() - 1) {
          // if we end here, we messed it up
          System.err.println("ERROR: We messed it up in Sweden.");
          return new Tuple<>(null, null);
        }
      }
    }

    // again check whether we messed something up badly
    final int[] iaTypesChild12 = assignMolecularTypes(gChild1.getMolecules());
    final int[] iaTypesChild22 = assignMolecularTypes(gChild2.getMolecules());
    for (int i = 0; i < iaTypesChild12.length; i++) {
      if (iaTypesChild12[i] != iaMolTypesFather[i] || iaTypesChild22[i] != iaMolTypesFather[i]) {
        // darn, didn't work
        System.err.println(
            "ERROR: Problem in Sweden after internal exchange "
                + "at position "
                + i
                + " father "
                + iaMolTypesFather[i]
                + " child 1 "
                + iaTypesChild12[i]
                + " child 2 "
                + iaTypesChild22[i]);
        return new Tuple<>(null, null);
      }
    }

    // ensure that the molecular IDs are again correct
    for (int i = 0; i < noMols; i++) {
      gChild1.getMoleculeAtPosition(i).setID(i);
      gChild2.getMoleculeAtPosition(i).setID(i);
    }

    // adding huge fitness values, just in case something goes really wrong
    gChild1.setFitness(Double.MAX_VALUE);
    gChild2.setFitness(Double.MAX_VALUE);

    if (DEBUG) {
      System.out.println("DEBUG: First child coming:");
      for (final String s : gChild1.makePrintableAbsoluteCoord(true)) {
        System.out.println(s);
      }

      System.out.println("DEBUG: Second child coming:");
      for (final String s : gChild2.makePrintableAbsoluteCoord(true)) {
        System.out.println(s);
      }
    }

    final boolean child1Fine = GlobOptAtomics.areGeometriesStillCorrect(mother, gChild1);
    final boolean child2Fine = GlobOptAtomics.areGeometriesStillCorrect(mother, gChild2);
    if (!child1Fine || !child2Fine) {
      System.err.println(
          "ERROR: After crossover, child1 is " + child1Fine + " and child2 is " + child2Fine);
      return new Tuple<>(null, null);
    }

    assert (gChild1.containsEnvironment() == mother.containsEnvironment());
    assert (gChild2.containsEnvironment() == father.containsEnvironment());
    assert (gChild1.containsEnvironment() == gChild2.containsEnvironment());

    return new Tuple<>(gChild1, gChild2);
  }

  @Override
  public short hasPriority() {
    return -1;
  }
}
