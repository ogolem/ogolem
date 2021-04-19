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

import static org.ogolem.core.FixedValues.NONCONVERGEDENERGY;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import org.ogolem.adaptive.genericfitness.PropertyCalculator;
import org.ogolem.adaptive.genericfitness.ReferenceInputData;
import org.ogolem.core.*;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.properties.Property;

/**
 * Provides global optimization of MOPAC parameters.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public final class AdaptiveMopacCaller extends AbstractAdaptivable implements Newton {

  // the ID
  private static final long serialVersionUID = (long) 20200729;

  private static long noOfGeomLocalOpts = (long) 0;
  private static long noOfMolLocalOpts = (long) 0;

  private final AdaptiveParameters parameters;

  private final boolean azoBenzene;

  private final boolean doSanityCheck;

  private final double blowBonds;

  private final double blowBondsEnv;

  public AdaptiveMopacCaller(
      final boolean doSanityCheck,
      final boolean bAzobenzeneStyle,
      final double blowFacBondDetect,
      final double blowFacClusterEnvClashDetect,
      final boolean inAdaptive) {
    this.azoBenzene = bAzobenzeneStyle;
    this.blowBonds = blowFacBondDetect;
    this.blowBondsEnv = blowFacClusterEnvClashDetect;
    this.doSanityCheck = doSanityCheck;
    if (!inAdaptive) {
      final String sFile = "adaptive-mopac.param";
      String[] saData;
      try {
        saData = Input.ReadFile(sFile);
      } catch (Exception e) {
        throw new RuntimeException("ERROR: Couldn't read parameter file in: " + e.toString());
      }
      parameters = new AdaptiveParameters(saData, -1);
    } else {
      parameters = null;
    }
  }

  private AdaptiveMopacCaller(final AdaptiveMopacCaller orig) {
    this.azoBenzene = orig.azoBenzene;
    this.doSanityCheck = orig.doSanityCheck;
    this.blowBonds = orig.blowBonds;
    this.blowBondsEnv = orig.blowBondsEnv;
    this.parameters = (orig.parameters == null) ? null : orig.parameters.copy();
  }

  @Override
  public AdaptiveMopacCaller copy() {
    return new AdaptiveMopacCaller(this);
  }

  @Override
  public String myIDandMethod() {
    return "adaptive MOPAC";
  }

  @Override
  public Geometry localOptimization(final Geometry gStartGeom) {

    incrGeom();

    if (parameters == null) {
      throw new RuntimeException(
          "Parameters in AdaptiveMopacCaller null. No parameters read in in consrtructor?");
    }

    CartesianCoordinates cartes = gStartGeom.getCartesians();

    final ZMatrix[] zmatVec = cartes.getZMatrices();
    try {
      cartes = CartesToCartes(gStartGeom.getID(), cartes, parameters);
    } catch (Exception e) {
      System.err.println(
          "Error in adaptive mopac locopt. Returning non-optimized geometry." + e.toString());
      gStartGeom.setFitness(NONCONVERGEDENERGY);
      return gStartGeom;
    }

    // we anyway set the old zmatrices first in, so that we have something in there
    cartes.setZMatrices(zmatVec);

    for (int i = 0; i < zmatVec.length; i++) {
      if (zmatVec[i] != null) {
        final CartesianCoordinates cartesTemp = cartes.giveMolecularCartes(i, true);
        zmatVec[i] = cartesTemp.calculateZMatrix();
      } // else: nothing needs to happen, molecule not flexible
    }

    cartes.setZMatrices(zmatVec);

    if (doSanityCheck) {
      // check the cartesian for sanity
      final boolean bSanity =
          GeometrySanityCheck.checkSanity(
              cartes, gStartGeom.getBondInfo(), blowBonds, blowBondsEnv);

      if (!bSanity) {
        // something's wrong
        gStartGeom.setLocalOptimized(true);
        gStartGeom.setFitness(FixedValues.NONCONVERGEDENERGY);

        return gStartGeom;
      }
    }

    final Geometry gEndGeom =
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
            gStartGeom.getBondInfo().copy());
    gEndGeom.setFitness(cartes.getEnergy());
    gEndGeom.setFather(gStartGeom.getFatherID());
    gEndGeom.setMother(gStartGeom.getMotherID());
    gEndGeom.setLocalOptimized(true);

    // return it
    return gEndGeom;
  }

  @Override
  public Molecule localOptimization(Molecule mStartMolecule) {

    incrMol();

    // first read the parameters in if needed, from aux file
    if (parameters == null) {
      throw new RuntimeException(
          "Parameters in AdaptiveMopacCaller null. Perhaps not properly constructed?");
    }

    CartesianCoordinates cartes = mStartMolecule.getCartesians();

    final ZMatrix[] zmatVec = cartes.getZMatrices();
    try {
      // since we can't be certain that there is something really "unique" in the molecules ID, we
      // hash it :-)
      int iHashCode = mStartMolecule.hashCode();
      cartes = CartesToCartes((long) iHashCode, cartes, parameters);
    } catch (Exception e) {
      System.err.println("Error in mopac locopt. Returning non-optimized molecule." + e.toString());
      mStartMolecule.setEnergy(NONCONVERGEDENERGY);
      return mStartMolecule;
    }

    // we anyway set the old zmatrices first in, so that we have something in there
    cartes.setZMatrices(zmatVec);

    for (int i = 0; i < zmatVec.length; i++) {
      if (zmatVec[i] != null) {
        CartesianCoordinates cartesTemp = cartes.giveMolecularCartes(i, true);
        zmatVec[i] = cartesTemp.calculateZMatrix();
      } // else: nothing needs to happen, unflexible molecule
    }

    cartes.setZMatrices(zmatVec);

    final Molecule mEndMolecule =
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

    final long id = params.getID();

    final String sMopacInput = "mopac_e" + id + ".dat";
    final String sMopacBasis = "mopac_e" + id;
    final String sParamFile = "mopac_e" + id + ".prm";
    final String sMopacFolder = "mopac_e" + id;

    /*
     * write the output aka input for the to be called program
     */
    try {
      Output.writeAdaptMopacInput(
          sMopacFolder, sMopacInput, cartes, false, true, params, azoBenzene, sParamFile);
    } catch (IOException e) {
      System.err.println(
          "Problem in writing geometry for mopac input (energy calculation). Returning NONCONVERGEDENERGY"
              + e.toString());
      return NONCONVERGEDENERGY;
    }

    /*
     * call mopac
     */
    Process proc;
    try {
      final Runtime rt = Runtime.getRuntime();
      final String sCommand = "mopac " + sMopacBasis;
      final File dir = new File(sMopacFolder);
      proc = rt.exec(sCommand, null, dir);

      // any error message?
      final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      if (proc.waitFor() != 0) {
        throw new ConvergenceException("Mopac returns non-zero return value (energy calculation).");
      } // else: mopac should(!) have completed normally...

    } catch (Exception e) {
      System.err.println("Mopac has a problem (energy calculation). ");
      e.printStackTrace(System.err);
      try {
        removeAllMopacFiles(sMopacFolder);
      } catch (Exception e1) {
        System.err.println("WARNING: Failed to clean up Mopac files. ");
        e1.printStackTrace(System.err);
      }
      return NONCONVERGEDENERGY;
    }

    /*
     * read mopac's output
     */

    double energy;

    final String sMopacOutput =
        sMopacFolder + System.getProperty("file.separator") + sMopacBasis + ".out";
    try {
      energy = Input.ReadEnergyMopacOutput(sMopacOutput, azoBenzene);
    } catch (Exception e) {
      System.err.println(
          "WARNING: Mopac has a problem to read the energy in (energy calculation). ");
      e.printStackTrace(System.err);
      try {
        removeAllMopacFiles(sMopacFolder);
      } catch (Exception e1) {
        System.err.println("WARNING: Failed to clean up Mopac files. ");
        e1.printStackTrace(System.err);
      }
      return NONCONVERGEDENERGY;
    }

    /*
     * clean up
     */
    try {
      removeAllMopacFiles(sMopacFolder);
    } catch (Exception e1) {
      System.err.println("WARNING: Failed to clean up Mopac files. ");
      e1.printStackTrace(System.err);
    }

    return energy;
  }

  @Override
  public double gradientOfStructWithParams(
      final CartesianCoordinates cartes,
      final AdaptiveParameters params,
      final int geomID,
      final BondInfo bonds,
      final double[] grad) {

    // XXX we just have a numerical gradient ATM
    return NumericalGradients.calculateParamGrad(cartes, params, this, geomID, bonds, grad);
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      final ArrayList<CartesianCoordinates> refCartes, final String sMethod) {

    // loop through all reference geometries and figure a list of non-redundant atom tyes out
    final List<String> llAtoms = new LinkedList<>();

    refCartes.stream()
        .map((cartesTemp) -> cartesTemp.getAllAtomTypes())
        .forEachOrdered(
            (saAtoms) -> {
              for (final String saAtom : saAtoms) {
                if (!llAtoms.contains(saAtom)) {
                  // add it to the list
                  llAtoms.add(saAtom);
                }
              }
            });

    int iParamSum = 0;

    final String[] saAtoms = new String[llAtoms.size()];
    final int[] iaParamsPerAt = new int[llAtoms.size()];

    // for each atom we need a set of parameters
    for (int i = 0; i < saAtoms.length; i++) {
      final String sTempAtom = llAtoms.get(i);
      iaParamsPerAt[i] = numberParamsRequiredForAtom(sTempAtom, sMethod);
      iParamSum += iaParamsPerAt[i];
      saAtoms[i] = sTempAtom;
    }

    final AdaptiveParameters paramStub =
        new AdaptiveParameters(iParamSum, -1, saAtoms, iaParamsPerAt, sMethod);

    return paramStub;
  }

  private int numberParamsRequiredForAtom(final String sAtomID, final String sMethod) {

    if (azoBenzene) {
      int iNumber;
      if (sAtomID.equalsIgnoreCase("H")) {
        iNumber = 0;
      } else if (sAtomID.equalsIgnoreCase("C")) {
        iNumber = 14;
      } else if (sAtomID.equalsIgnoreCase("N")) {
        iNumber = 21;
      } else {
        System.err.print(
            "ERROR: No information about the number of parameters "
                + "for atom "
                + sAtomID
                + " in azobenzene context. Returning zero.");
        iNumber = 0;
      }
      return iNumber;
    } else {
      if (sMethod.equalsIgnoreCase("AM1")) {
        int iNumber;
        if (sAtomID.equalsIgnoreCase("H")) {
          iNumber = 15;
        } else {
          // all non-hydrogens are treated the same
          iNumber = 21;
        }
        return iNumber;
      } else {
        System.err.print(
            "ERROR: No information about the number of parameters "
                + "for atom "
                + sAtomID
                + " in context of method "
                + sMethod
                + ". Returning zero.");
        return 0;
      }
    }
  }

  @Override
  public double[][] minMaxBordersForParams(final AdaptiveParameters params) {

    final double[][] daBorders = new double[2][params.getNumberOfParamters()];

    final String[] saAtoms = params.getForWhichAtoms();

    int iCounter = 0;

    if (azoBenzene) {
      for (final String saAtom : saAtoms) {
        if (saAtom.equalsIgnoreCase("H")) {
          // nothing, since we do NOT reoptimize hydrogen parameters p.d.
        } else if (saAtom.equalsIgnoreCase("C")) {
          // 1: USS
          daBorders[0][iCounter] = -163.32;
          daBorders[1][iCounter] = -9.11;
          iCounter++;

          // 2: UPP
          daBorders[0][iCounter] = -125.87;
          daBorders[1][iCounter] = 0.00;
          iCounter++;

          // 3: BETST this speculates a bit
          daBorders[0][iCounter] = -100;
          daBorders[1][iCounter] = 0;
          iCounter++;

          // 4: BETPT this speculates a bit
          daBorders[0][iCounter] = -100;
          daBorders[1][iCounter] = 0;
          iCounter++;

          // 5: BETAS
          daBorders[0][iCounter] = -83.51;
          daBorders[1][iCounter] = -0.71;
          iCounter++;

          // 6: BETAP
          daBorders[0][iCounter] = -35.12;
          daBorders[1][iCounter] = 0.00;
          iCounter++;

          // 7: ZS
          daBorders[0][iCounter] = 0.95;
          daBorders[1][iCounter] = 4.52;
          iCounter++;

          // 8: ZP
          daBorders[0][iCounter] = 0.80;
          daBorders[1][iCounter] = 3.02;
          iCounter++;

          // 9: ALP
          daBorders[0][iCounter] = 1.18;
          daBorders[1][iCounter] = 6.61;
          iCounter++;

          // 10: GSS
          daBorders[0][iCounter] = 6.47;
          daBorders[1][iCounter] = 20.30;
          iCounter++;

          // 11: GSP
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 20.70;
          iCounter++;

          // 12: GPP
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 20.05;
          iCounter++;

          // 13: GP2
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 17.89;
          iCounter++;

          // 14: HSP
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 5.80;
          iCounter++;
        } else if (saAtom.equalsIgnoreCase("N")) {
          // 1: USS
          daBorders[0][iCounter] = -163.32;
          daBorders[1][iCounter] = -9.11;
          iCounter++;

          // 2: UPP
          daBorders[0][iCounter] = -125.87;
          daBorders[1][iCounter] = 0.00;
          iCounter++;

          // 3: BETAS
          daBorders[0][iCounter] = -83.51;
          daBorders[1][iCounter] = -0.71;
          iCounter++;

          // 4: BETAP
          daBorders[0][iCounter] = -35.12;
          daBorders[1][iCounter] = 0.00;
          iCounter++;

          // 5: GSS
          daBorders[0][iCounter] = 6.47;
          daBorders[1][iCounter] = 20.30;
          iCounter++;

          // 6: GSP
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 20.70;
          iCounter++;

          // 7: GPP
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 20.05;
          iCounter++;

          // 8: GP2
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 17.89;
          iCounter++;

          // 9: HSP
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 5.80;
          iCounter++;

          // 10: ZS
          daBorders[0][iCounter] = 0.95;
          daBorders[1][iCounter] = 4.52;
          iCounter++;

          // 11: ZP
          daBorders[0][iCounter] = 0.80;
          daBorders[1][iCounter] = 3.02;
          iCounter++;

          // 12: ALP
          daBorders[0][iCounter] = 1.18;
          daBorders[1][iCounter] = 6.61;
          iCounter++;

          // 13-21: FN11-FN33
          for (int k = 0; k < 9; k++) {
            daBorders[0][iCounter] = -9.9;
            daBorders[1][iCounter] = 9.9;
            iCounter++;
          }
        } else {
          System.err.println(
              "WARNING: No informations on borders for " + saAtom + " found. Using big ones.");
          for (int j = 0; j < params.getAmountOfParametersForKey(saAtom); j++) {
            daBorders[0][iCounter] = -500;
            daBorders[1][iCounter] = 100;
            iCounter++;
          }
        }
      }
    } else {

      for (final String saAtom : saAtoms) {
        if (saAtom.equalsIgnoreCase("H")) {
          // 1. USS
          daBorders[0][iCounter] = -13.00;
          daBorders[1][iCounter] = -9.00;
          iCounter++;

          // 2. BETAS
          daBorders[0][iCounter] = -8.00;
          daBorders[1][iCounter] = -4.00;
          iCounter++;

          // 3. ZS
          daBorders[0][iCounter] = 0.70;
          daBorders[1][iCounter] = 1.50;
          iCounter++;

          // 4. ALP
          daBorders[0][iCounter] = 2.00;
          daBorders[1][iCounter] = 4.00;
          iCounter++;

          // 5. GSS
          daBorders[0][iCounter] = 10.00;
          daBorders[1][iCounter] = 15.00;
          iCounter++;

          // 6.-15. FN11-FN33
          for (int k = 0; k < 9; k++) {
            daBorders[0][iCounter] = -9.9;
            daBorders[1][iCounter] = 9.9;
            iCounter++;
          }
        } else {
          // 1: USS
          daBorders[0][iCounter] = -163.32;
          daBorders[1][iCounter] = -9.11;
          iCounter++;

          // 2: UPP
          daBorders[0][iCounter] = -125.87;
          daBorders[1][iCounter] = 0.00;
          iCounter++;

          // 3: BETAS
          daBorders[0][iCounter] = -83.51;
          daBorders[1][iCounter] = -0.71;
          iCounter++;

          // 4: BETAP
          daBorders[0][iCounter] = -35.12;
          daBorders[1][iCounter] = 0.00;
          iCounter++;

          // 5: GSS
          daBorders[0][iCounter] = 6.47;
          daBorders[1][iCounter] = 20.30;
          iCounter++;

          // 6: GSP
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 20.70;
          iCounter++;

          // 7: GPP
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 20.05;
          iCounter++;

          // 8: GP2
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 17.89;
          iCounter++;

          // 9: HSP
          daBorders[0][iCounter] = 0.00;
          daBorders[1][iCounter] = 5.80;
          iCounter++;

          // 10: ZS
          daBorders[0][iCounter] = 0.95;
          daBorders[1][iCounter] = 4.52;
          iCounter++;

          // 11: ZP
          daBorders[0][iCounter] = 0.80;
          daBorders[1][iCounter] = 3.02;
          iCounter++;

          // 12: ALP
          daBorders[0][iCounter] = 1.18;
          daBorders[1][iCounter] = 6.61;
          iCounter++;

          // 13-21: FN11-FN33
          for (int k = 0; k < 9; k++) {
            daBorders[0][iCounter] = -9.9;
            daBorders[1][iCounter] = 9.9;
            iCounter++;
          }
        }
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
    throw new Exception("cartes2cartes not implemented in adaptive mopac.");
  }

  private static CartesianCoordinates CartesToCartes(
      final long lID, final CartesianCoordinates cartes, final AdaptiveParameters params)
      throws ConvergenceException, InitIOException, CastException {

    // TODO just for the general case needed
    return null;
  }

  @Override
  public CartesianFullBackend getBackend() {
    return null;
  }

  private static void removeAllMopacFiles(final String sMopacFolder) throws IOException {
    // just remove the complete folder
    ManipulationPrimitives.remove(sMopacFolder);
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
