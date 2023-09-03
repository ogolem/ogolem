/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive;

import static org.ogolem.core.Constants.ANGTOBOHR;
import static org.ogolem.core.FixedValues.NONCONVERGEDENERGY;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import org.ogolem.adaptive.genericfitness.PropertyCalculator;
import org.ogolem.adaptive.genericfitness.ReferenceInputData;
import org.ogolem.core.*;
import org.ogolem.properties.Property;

/**
 * Calls orca and allows the manipulation of semiempirical parameters.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public final class AdaptiveOrcaCaller extends AbstractAdaptivable implements Newton {

  // the ID
  private static final long serialVersionUID = (long) 20160121;

  private static long noOfGeomLocalOpts = (long) 0;
  private static long noOfMolLocalOpts = (long) 0;

  private AdaptiveParameters parameters = null;

  private final double dBlowBonds;

  private final double dBlowBondsEnv;

  private final boolean bDoSanityCheck;

  public AdaptiveOrcaCaller(
      final double blowFacBondDetect,
      final double blowFacClusterEnvClashDetect,
      final boolean doSanityCheck) {
    // TODO fix this!
    this.dBlowBonds = blowFacBondDetect;
    this.dBlowBondsEnv = blowFacClusterEnvClashDetect;
    this.bDoSanityCheck = doSanityCheck;
  }

  private AdaptiveOrcaCaller(AdaptiveOrcaCaller orig) {
    this.bDoSanityCheck = orig.bDoSanityCheck;
    this.dBlowBonds = orig.dBlowBonds;
    this.dBlowBondsEnv = orig.dBlowBondsEnv;
  }

  @Override
  public AdaptiveOrcaCaller copy() {
    return new AdaptiveOrcaCaller(this);
  }

  @Override
  public String myIDandMethod() {
    return "adaptive Orca";
  }

  @Override
  public Geometry localOptimization(Geometry gStartGeom) {

    // first read the parameters in if needed, from aux file
    if (parameters == null) {
      String sFile = "adaptive-orca.param";
      String[] saData = null;
      try {
        saData = Input.ReadFile(sFile);
      } catch (Exception e) {
        System.err.println("ERROR: Couldn't read parameter file in: " + e.toString());
        gStartGeom.setFitness(NONCONVERGEDENERGY);
        return gStartGeom;
      }
      int iID = 0;
      parameters = new AdaptiveParameters(saData, iID);
    }

    CartesianCoordinates cartes = gStartGeom.getCartesians();

    final ZMatrix[] zmatVec = cartes.getZMatrices();
    try {
      cartes = CartesToCartes(gStartGeom.getID(), cartes, parameters);
    } catch (Exception e) {
      System.err.println(
          "Error in adaptive orca locopt. Returning non-optimized geometry." + e.toString());
      gStartGeom.setFitness(NONCONVERGEDENERGY);
      return gStartGeom;
    }

    // we anyway set the old zmatrices first in, so that we have something in there
    cartes.setZMatrices(zmatVec);

    for (int i = 0; i < zmatVec.length; i++) {
      if (zmatVec[i] == null) {
        // nothing needs to happen, unflexible molecule
      } else {
        CartesianCoordinates cartesTemp = cartes.giveMolecularCartes(i, true);
        zmatVec[i] = cartesTemp.calculateZMatrix();
      }
    }

    cartes.setZMatrices(zmatVec);

    if (bDoSanityCheck) {
      // check the cartesian for sanity
      final boolean bSanity =
          GeometrySanityCheck.checkSanity(
              cartes, gStartGeom.getBondInfo(), dBlowBonds, dBlowBondsEnv);

      if (!bSanity) {
        // something's wrong
        gStartGeom.setLocalOptimized(true);
        gStartGeom.setFitness(FixedValues.NONCONVERGEDENERGY);

        return gStartGeom;
      }
    }

    Geometry gEndGeom =
        new Geometry(
            cartes,
            gStartGeom.getID(),
            gStartGeom.getNumberOfIndieParticles(),
            cartes.getAllAtomsPerMol(),
            gStartGeom.getAllFlexies(),
            gStartGeom.getExplicitDoFs(),
            gStartGeom.getAllConstraints(),
            gStartGeom.getAllConstraintsXYZ(),
            gStartGeom.getSIDs(),
            gStartGeom.getBondInfo());
    gEndGeom.setFitness(cartes.getEnergy());
    gEndGeom.setFather(gStartGeom.getFatherID());
    gEndGeom.setMother(gStartGeom.getMotherID());
    gEndGeom.setLocalOptimized(true);

    // return it
    return gEndGeom;
  }

  @Override
  public Molecule localOptimization(Molecule mStartMolecule) {

    // first read the parameters in if needed, from aux file
    if (parameters == null) {
      String sFile = "adaptive-orca.param";
      String[] saData = null;
      try {
        saData = Input.ReadFile(sFile);
      } catch (Exception e) {
        System.err.println("ERROR: Couldn't read parameter file in: " + e.toString());
        mStartMolecule.setEnergy(NONCONVERGEDENERGY);
        return mStartMolecule;
      }
      int iID = 0;
      parameters = new AdaptiveParameters(saData, iID);
    }

    CartesianCoordinates cartes = mStartMolecule.getCartesians();

    final ZMatrix[] zmatVec = cartes.getZMatrices();
    try {
      // since we can't be certain that there is something really "unique" in the molecules ID, we
      // hash it :-)
      int iHashCode = mStartMolecule.hashCode();
      cartes = CartesToCartes((long) iHashCode, cartes, parameters);
    } catch (Exception e) {
      System.err.println("Error in orca locopt. Returning non-optimized molecule." + e.toString());
      mStartMolecule.setEnergy(NONCONVERGEDENERGY);
      return mStartMolecule;
    }

    // we anyway set the old zmatrix first in, so that we have something in there
    cartes.setZMatrices(zmatVec);

    for (int i = 0; i < zmatVec.length; i++) {
      if (zmatVec[i] == null) {
        // nothing needs to happen, unflexible molecule
      } else {
        CartesianCoordinates cartesTemp = cartes.giveMolecularCartes(i, true);
        zmatVec[i] = cartesTemp.calculateZMatrix();
      }
    }

    cartes.setZMatrices(zmatVec);

    Molecule mEndMolecule =
        new Molecule(
            cartes,
            mStartMolecule.getMolPosition(),
            mStartMolecule.getSID(),
            mStartMolecule.getFlexy(),
            mStartMolecule.getDegreesOfFreedom(),
            mStartMolecule.isConstricted(),
            mStartMolecule.getConstraints());
    mEndMolecule.setID(mStartMolecule.getMolPosition());
    mEndMolecule.setSID(mStartMolecule.getSID());

    return mEndMolecule;
  }

  @Override
  public double energyOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds) {

    long lTemp = params.getID() + System.currentTimeMillis() + cartes.hashCode();

    String sOrcaInput = "orca_energy" + lTemp + ".inp";
    String sOrcaBasis = "orca_energy" + lTemp;

    /*
     * write the output aka input for the to be called program
     */
    try {
      Output.writeAdaptOrcaInput(sOrcaInput, cartes, false, false, true, params);
    } catch (IOException e) {
      System.err.println(
          "Problem in writing geometry for orca input (energy calculation). Returning NONCONVERGEDENERGY"
              + e.toString());
      return NONCONVERGEDENERGY;
    }

    /*
     * call orca
     */
    String[] saOutput;
    try {
      Runtime rt = Runtime.getRuntime();
      String[] saCmd = new String[] {"orca", sOrcaInput};
      Process proc = rt.exec(saCmd);

      // any error message?
      StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      int iExitValue = proc.waitFor();
      if (iExitValue != 0) {
        throw new ConvergenceException("Orca returns non-zero return value (energy calculation).");
      } else {
        // orca should(!) have completed normally...
        saOutput = outputGobbler.getData();
      }

    } catch (Exception e) {
      System.err.println("Orca has a problem (energy calculation). " + e.toString());
      return NONCONVERGEDENERGY;
    }

    /*
     * read orcas output
     */

    double dEnergy = getEnergyFromOutput(saOutput);

    /*
     * clean up
     */
    try {
      org.ogolem.core.Input.RemoveOrcaFiles(sOrcaBasis);
    } catch (Exception e) {
      System.err.println("Problem cleaning orca files up." + e.toString());
    }

    return dEnergy;
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds,
      final double[] grad) {

    long lTemp = params.getID() + System.currentTimeMillis() + cartes.hashCode();

    String sOrcaInput = "orca_gradient" + lTemp + ".inp";
    String sOrcaBasis = "orca_gradient" + lTemp;

    /*
     * write the output aka input for the to be called program
     */
    try {
      Output.writeAdaptOrcaInput(sOrcaInput, cartes, false, true, false, params);
    } catch (IOException e) {
      System.err.println(
          "Problem in writing geometry for orca input (gradient calculation). Returning empty gradient"
              + e.toString());
      return 1000.0;
    }

    // System.OUT.println("DEBUG: Successfully written input.");

    /*
     * call orca
     */
    String[] saOutput;
    try {
      Runtime rt = Runtime.getRuntime();
      String[] saCmd = new String[] {"orca", sOrcaInput};
      Process proc = rt.exec(saCmd);

      // any error message?
      StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      int iExitValue = proc.waitFor();
      if (iExitValue != 0) {
        throw new ConvergenceException(
            "Orca returns non-zero return value (gradient calculation).");
      } else {
        // orca should(!) have completed normally...
        saOutput = outputGobbler.getData();
      }

    } catch (Exception e) {
      System.err.println("Orca has a problem (gradient calculation). " + e.toString());
      return 1000.0;
    }

    /*
     * read orcas output
     */

    final double en = getGradientFromOutput(saOutput, params.getNumberOfParamters(), grad);

    /*
     * clean up
     */
    try {
      org.ogolem.core.Input.RemoveOrcaFiles(sOrcaBasis);
    } catch (Exception e) {
      System.err.println("Problem cleaning orca files up." + e.toString());
    }

    return en;
  }

  private int numberParamsRequiredForAtom(final String sAtomID, final String sMethod) {
    if (sMethod.equalsIgnoreCase("AM1")) {
      int iNumber = 0;
      if (sAtomID.equalsIgnoreCase("H")) {
        iNumber = 2;
      } else {
        System.err.print(
            "ERROR: No information about the number of parameters "
                + "for atom "
                + sAtomID
                + ". Returning zero.");
        iNumber = 0;
      }
      return iNumber;
    } else {
      System.err.print(
          "ERROR: No information about the number of parameters "
              + "for atom "
              + sAtomID
              + ". Returning zero.");
      return 0;
    }
    // TODO more on parameters
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      final ArrayList<CartesianCoordinates> refCartes, final String sMethod) {

    // loop through all reference geometries and figure a list of non-redundant atom tyes out
    final LinkedList<String> llAtoms = new LinkedList<>();

    final Iterator<CartesianCoordinates> itRefGeoms = refCartes.iterator();

    while (itRefGeoms.hasNext()) {
      final CartesianCoordinates cartesTemp = itRefGeoms.next();
      final String[] saAtoms = cartesTemp.getAllAtomTypes();

      for (final String saAtom : saAtoms) {
        if (!llAtoms.contains(saAtom)) {
          // add it to the list
          llAtoms.add(saAtom);
        }
      }
    }

    int iParamSum = 0;

    String[] saAtoms = new String[llAtoms.size()];
    int[] iaParamsPerAt = new int[llAtoms.size()];

    // for each atom we need a set of parameters
    for (int i = 0; i < saAtoms.length; i++) {
      String sTempAtom = llAtoms.get(i);
      iaParamsPerAt[i] = numberParamsRequiredForAtom(sTempAtom, sMethod);
      iParamSum += iaParamsPerAt[i];
      saAtoms[i] = sTempAtom;
    }

    final AdaptiveParameters paramStub =
        new AdaptiveParameters(iParamSum, -1, saAtoms, iaParamsPerAt, sMethod);

    return paramStub;
  }

  @Override
  public double[][] minMaxBordersForParams(final AdaptiveParameters params) {

    final double[][] daBorders = new double[2][params.getNumberOfParamters()];

    final String[] saAtoms = params.getForWhichAtoms();

    int iCounter = 0;

    for (final String saAtom : saAtoms) {
      if (saAtom.equalsIgnoreCase("H")) {
        daBorders[0][iCounter] = -100;
        daBorders[1][iCounter] = +100;
        iCounter++;

        daBorders[0][iCounter] = -100;
        daBorders[1][iCounter] = +100;
        iCounter++;
      } else {
        System.err.println(
            "WARNING: No informations on borders for " + saAtom + " found. Using big ones.");
        for (int j = 0; j < params.getAmountOfParametersForKey(saAtom); j++) {
          daBorders[0][iCounter] = Double.MIN_VALUE;
          daBorders[0][iCounter] = Double.MAX_VALUE;
          iCounter++;
        }
        // TODO atom specific
      }
    }

    return daBorders;
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      long id,
      CartesianCoordinates cartes,
      boolean[][] constraints,
      boolean isConstricted,
      final BondInfo bonds)
      throws Exception {
    throw new Exception("cartes2cartes not implemented in adaptive orca.");
  }

  private static CartesianCoordinates CartesToCartes(
      long lID, CartesianCoordinates cartes, AdaptiveParameters params)
      throws ConvergenceException, InitIOException, CastException {

    String sOrcaInput = "orca" + lID + ".inp";
    String sOrcaBasis = "orca" + lID;
    float[] iaCharges = cartes.getAllCharges();
    short[] iaSpins = cartes.getAllSpins();
    /*
     * write the output aka input for the to be called program
     */
    try {
      Output.writeAdaptOrcaInput(sOrcaInput, cartes, true, false, false, params);
    } catch (IOException e) {
      throw new ConvergenceException(
          "Problem in writing geometry for orca input (local optimization).", e);
    }

    // System.OUT.println("DEBUG: Successfully written input.");

    /*
     * call orca
     */
    String[] saOutput;
    try {
      Runtime rt = Runtime.getRuntime();
      String[] saCmd = new String[] {"orca", sOrcaInput};
      Process proc = rt.exec(saCmd);

      // any error message?
      StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      int iExitValue = proc.waitFor();
      if (iExitValue != 0) {
        throw new ConvergenceException("Orca returns non-zero return value (local optimization).");
      } else {
        // orca should(!) have completed normally...
        saOutput = outputGobbler.getData();
      }

    } catch (Exception e) {
      throw new ConvergenceException("Orca has a problem (local optimization).", e);
    }

    /*
     * read orcas output
     */

    cartes =
        createCartesFromOutput(
            saOutput, cartes.getNoOfAtoms(), cartes.getNoOfMolecules(), cartes.getAllAtomsPerMol());

    /*
     * clean up
     */
    try {
      org.ogolem.core.Input.RemoveOrcaFiles(sOrcaBasis);
    } catch (Exception e) {
      System.err.println("Problem cleaning orca files up." + e.toString());
    }

    /*
     * return it
     */
    cartes.setAllCharges(iaCharges);
    cartes.setAllSpins(iaSpins);

    return cartes;
  }

  private static CartesianCoordinates createCartesFromOutput(
      String[] saOutput, int iNoOfAtoms, int iNoOfMoles, int[] iaNoOfAtsPerMol)
      throws InitIOException, CastException {
    CartesianCoordinates cartes = new CartesianCoordinates(iNoOfAtoms, iNoOfMoles, iaNoOfAtsPerMol);

    /*
     * check where the coordinates start from:
     * Last occurence of "CARTESIAN COORDINATES (ANGSTROEM)"
     * check where the energy is:
     * Last occurence of "Total Energy"
     */
    int iCartesStart = 0;
    int iEnergyIndex = 0;
    for (int i = 0; i < saOutput.length; i++) {
      if (saOutput[i].contains("CARTESIAN COORDINATES (ANGSTROEM)")) {
        iCartesStart = i + 2;
      } else if (saOutput[i].contains("Total Energy")) {
        iEnergyIndex = i;
      } else {
        // really nothing
      }
    }
    // get all coordinates
    String sTemp2;
    int iIndex;
    for (int i = 0; i < iNoOfAtoms; i++) {
      String sTemp = saOutput[i + iCartesStart];
      sTemp = sTemp.trim();
      String sAtom = sTemp.substring(0, sTemp.indexOf(" "));
      sTemp = sTemp.substring(sTemp.indexOf(" "));
      sTemp = sTemp.trim();

      // get the xyz coordinates
      double[] daCoords = new double[3];

      iIndex = sTemp.indexOf(" ");
      sTemp2 = sTemp.substring(0, iIndex);
      sTemp = sTemp.substring(iIndex);
      sTemp = sTemp.trim();
      try {
        daCoords[0] = Double.parseDouble(sTemp2) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException(e);
      }

      iIndex = sTemp.indexOf(" ");
      sTemp2 = sTemp.substring(0, iIndex);
      sTemp = sTemp.substring(iIndex);
      sTemp = sTemp.trim();
      try {
        daCoords[1] = Double.parseDouble(sTemp2) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException(e);
      }

      try {
        daCoords[2] = Double.parseDouble(sTemp) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException(e);
      }

      // set the coordinates and the atom to the cartes
      cartes.setAtom(sAtom, i);
      cartes.setXYZCoordinatesOfAtom(daCoords, i);
    }

    // we know where the energy is from before
    String sTemp = saOutput[iEnergyIndex];
    sTemp = sTemp.trim();
    sTemp = sTemp.substring(sTemp.indexOf(":") + 1);
    sTemp = sTemp.trim();
    sTemp = sTemp.substring(0, sTemp.indexOf(" "));
    double dEnergy;
    try {
      dEnergy = Double.parseDouble(sTemp);
    } catch (Exception e) {
      throw new CastException(e);
    }

    cartes.setEnergy(dEnergy);
    return cartes;
  }

  private static double getEnergyFromOutput(final String[] saOutput) {

    double dEnergy = NONCONVERGEDENERGY;

    // start from the end, that way we'll reach the energy quicker
    for (int i = saOutput.length - 1; i >= 0; i--) {
      if (saOutput[i].startsWith("FINAL SINGLE POINT ENERGY")) {
        String sTemp = saOutput[i].substring(25).trim();
        try {
          dEnergy = Double.parseDouble(sTemp);
        } catch (Exception e) {
          System.err.println(
              "ERROR: Failed to parse the energy from the orca output. "
                  + "Returning high energy. "
                  + e.toString());
          return NONCONVERGEDENERGY;
        }
        break;
      }
    }

    return dEnergy;
  }

  private static double getGradientFromOutput(
      final String[] saOutput, final int iNoOfParams, final double[] grad) {

    int iGradBegin = -1;
    for (int i = 0; i < saOutput.length; i++) {
      if (saOutput[i].contains("   GRADGRAD   ")) {
        iGradBegin = i;
        break;
      }
    }

    if (iGradBegin == -1) {
      System.err.println("ERROR: No gradient found in the orca output.");
      return NONCONVERGEDENERGY;
    }

    // TODO once the orca implementation is complete

    final double dEnergy = getEnergyFromOutput(saOutput);

    return dEnergy;
  }

  @Override
  public CartesianFullBackend getBackend() {
    return null;
  }

  private static synchronized void incrGeom() {
    noOfGeomLocalOpts++;
  }

  private static synchronized void incrMol() {
    noOfMolLocalOpts++;
  }

  @Override
  public long getNumberOfGeomLocalOpts() {
    return noOfGeomLocalOpts;
  }

  @Override
  public long getNumberOfMolLocalOpts() {
    return noOfMolLocalOpts;
  }

  @Override
  public <T extends Property, V extends ReferenceInputData<T>>
      PropertyCalculator<T, V> getCalculatorForProperty(final T property, final V data) {
    return null;
  }
}
