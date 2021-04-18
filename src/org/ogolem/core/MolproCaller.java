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
package org.ogolem.core;

import org.ogolem.io.OutputPrimitives;

/**
 * This calls MOLPRO for energy and/or gradient calculation as well as local optimizations using the
 * famous MOLPRO optimizer.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
final class MolproCaller implements CartesianFullBackend, Newton {

  // the ID
  private static final long serialVersionUID = (long) 20200729;

  static enum METHOD {
    AM1,
    HFVDZ,
    B86VDZ,
    MP2AVTZ,
    CUSTOM,
    CUSTOMNOLOCOPT
  };

  private static long noOfGeomLocalOpts = (long) 0;
  private static long noOfMolLocalOpts = (long) 0;

  /*
   * The choice of which method should provide energy and gradient
   * -2: use custom file and do not read cartes in fresh
   * -1: use custom file
   * 0: am1
   * 1: hf/vdz
   * 2: b86/vdz
   * 3: mp2/avtz
   */
  private final METHOD whichMethod;

  private final String sCustomInp;

  private double dBlowBonds;

  private double dBlowBondsEnv;

  private boolean bDoSanityCheck;

  /**
   * The constructor for usage as a backend.
   *
   * @param whichMethod The choice of method.
   */
  MolproCaller(final METHOD whichMethod) {
    // we do not allow custom input here. Since this "Interface" makes little sense anyways
    this.sCustomInp = null;
    this.whichMethod = whichMethod;
  }

  /**
   * The constructor for usage as a local optimizing engine.
   *
   * @param globconf A complete configuration set.
   * @param whichMethod The choice of method.
   */
  MolproCaller(GlobalConfig globconf, METHOD whichMethod) {

    if (whichMethod == METHOD.CUSTOM || whichMethod == METHOD.CUSTOMNOLOCOPT) {
      final String sWhichAuxFile = globconf.outputFolder + "-molpro.aux";
      this.sCustomInp = readCustomInp(sWhichAuxFile);
    } else {
      this.sCustomInp = null;
    }

    this.whichMethod = whichMethod;
    this.dBlowBonds = globconf.blowFacBondDetect;
    this.dBlowBondsEnv = globconf.blowFacEnvClusterClashes;
    this.bDoSanityCheck = globconf.doPostSanityCheck;
  }

  private MolproCaller(final MolproCaller orig) {
    this.whichMethod = orig.whichMethod;
    this.dBlowBonds = orig.dBlowBonds;
    this.dBlowBondsEnv = orig.dBlowBondsEnv;
    this.bDoSanityCheck = orig.bDoSanityCheck;
    if (orig.sCustomInp != null) {
      this.sCustomInp = orig.sCustomInp;
    } else {
      this.sCustomInp = null;
    }
  }

  @Override
  public MolproCaller copy() {
    return new MolproCaller(this);
  }

  @Override
  public String getMethodID() {
    String method = "";
    switch (whichMethod) {
      case CUSTOMNOLOCOPT:
        method = "custom input, no locopt";
        break;
      case CUSTOM:
        method = "custom input";
        break;
      case AM1:
        method = "am1";
        break;
      case HFVDZ:
        method = "hf/vdz";
        break;
      case B86VDZ:
        method = "b86/vdz";
        break;
      case MP2AVTZ:
        method = "mp2/avtz";
        break;
      default:
        break;
    }

    return "MOLPRO: " + method;
  }

  @Override
  public String myIDandMethod() {

    String method = "";
    switch (whichMethod) {
      case CUSTOMNOLOCOPT:
        method = "custom input, no locopt";
        break;
      case CUSTOM:
        method = "custom input";
        break;
      case AM1:
        method = "am1";
        break;
      case HFVDZ:
        method = "hf/vdz";
        break;
      case B86VDZ:
        method = "b86/vdz";
        break;
      case MP2AVTZ:
        method = "mp2/avtz";
        break;
      default:
        break;
    }

    return "MOLPRO: " + method;
  }

  @Override
  public double energyCalculation(
      long lID,
      int iIteration,
      double[] da1DXYZ,
      String[] saAtomType,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      int iNoOfAtoms,
      float[] faCharges,
      short[] iaSpins,
      final BondInfo bonds,
      final boolean hasRigidEnv) {

    String sMolproInput = "molpro" + lID + "i" + iIteration + "energy.inp";
    float fTotalCharge = 0;
    for (int i = 0; i < faCharges.length; i++) {
      fTotalCharge += faCharges[i];
    }
    /*
     * write the input file
     */
    try {
      WriteOutput(sMolproInput, false, false, da1DXYZ, saAtomType, Math.round(fTotalCharge));
    } catch (InitIOException e) {
      System.err.println(
          "Problem in writing geometry for molpro input (energy calculation)." + e.toString());
      return FixedValues.NONCONVERGEDENERGY;
    }

    /*
     * call molpro
     */
    try {
      Runtime rt = Runtime.getRuntime();
      String sMolproCmd = System.getenv("OGO_MOLPROCMD");
      if (sMolproCmd == null) {
        // default
        sMolproCmd = "molpro";
      }
      Process proc = rt.exec(sMolproCmd + " " + sMolproInput);

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
        return 1001.1;
      } else {
        // molpro should(!) have completed normally...
      }

    } catch (Exception e) {
      System.err.println("Problem calling molpro!" + e.toString());
      return FixedValues.NONCONVERGEDENERGY;
    }

    /*
     * read molpros output
     */
    double dEnergy;
    String sMolproOutput = sMolproInput.substring(0, sMolproInput.indexOf(".")) + ".out";
    try {
      dEnergy = Input.ReadEnergyMolproOutput(sMolproOutput);
    } catch (Exception e) {
      System.err.println("Problem reading molpros output!" + e.toString());
      return FixedValues.NONCONVERGEDENERGY;
    }

    /*
     * clean up
     */
    try {
      Input.RemoveMolproFiles(sMolproInput, false);
    } catch (Exception e) {
      System.err.println("Problem cleaning molpro files up." + e.toString());
    }

    /*
     * return it
     */
    return dEnergy;
  }

  @Override
  public void gradientCalculation(
      long lID,
      int iIteration,
      double[] da1DXYZ,
      String[] saAtomType,
      short[] atomNos,
      int[] atsPerMol,
      double[] energyparts,
      int iNoOfAtoms,
      float[] faCharges,
      short[] iaSpins,
      final BondInfo bonds,
      final Gradient gradient,
      final boolean hasRigidEnv) {

    String sMolproInput = "molpro" + lID + "i" + iIteration + "grad.inp";
    float fTotalCharge = 0;
    for (int i = 0; i < faCharges.length; i++) {
      fTotalCharge += faCharges[i];
    }
    /*
     * write the input file
     */
    try {
      WriteOutput(sMolproInput, false, true, da1DXYZ, saAtomType, Math.round(fTotalCharge));
    } catch (InitIOException e) {
      System.err.println(
          "Problem in writing geometry for molpro input (energy calculation)." + e.toString());
      gradient.markProblem();
      return;
    }

    /*
     * call molpro
     */
    try {
      Runtime rt = Runtime.getRuntime();
      Process proc = rt.exec("molpro " + sMolproInput);

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
        gradient.markProblem();
        return;
      } else {
        // molpro should(!) have completed normally...
      }

    } catch (Exception e) {
      System.err.println("Problem calling molpro!" + e.toString());
      gradient.markProblem();
      return;
    }

    /*
     * read molpros output
     */
    String sMolproOutput = sMolproInput.substring(0, sMolproInput.indexOf(".")) + ".out";
    try {
      final Gradient gradientTmp = Input.ReadGradientMolproOutput(sMolproOutput, iNoOfAtoms);
      gradient.copyDataIn(gradientTmp);
    } catch (Exception e) {
      System.err.println("Problem reading molpros output for gradient!" + e.toString());
      gradient.markProblem();
      return;
    }

    /*
     * clean up
     */
    try {
      Input.RemoveMolproFiles(sMolproInput, false);
    } catch (Exception e) {
      System.err.println("Problem cleaning molpro files up." + e.toString());
    }
  }

  @Override
  public Geometry localOptimization(Geometry gStartGeom) {

    incrGeom();

    CartesianCoordinates cartes = gStartGeom.getCartesians();
    ZMatrix[] zmatVec = cartes.getZMatrices();
    try {
      cartes =
          cartesToCartes(
              gStartGeom.getID(),
              cartes,
              gStartGeom.getAllConstraintsXYZ(),
              gStartGeom.isThereAConstraint(),
              new SimpleBondInfo(1)); // doesn't matter really
    } catch (Exception e) {
      System.err.println(
          "Error in molpro locopt. Returning non-optimized geometry." + e.toString());
      gStartGeom.setFitness(FixedValues.NONCONVERGEDENERGY);
      return gStartGeom;
    }

    // we anyway set the old zmatrices first in, so that we have something in there
    cartes.setZMatrices(zmatVec);

    for (int i = 0; i < zmatVec.length; i++) {
      if (zmatVec[i] != null) {
        final CartesianCoordinates cartesTemp = cartes.giveMolecularCartes(i, true);
        zmatVec[i] = cartesTemp.calculateZMatrix();
      } // else: nothing needs to happen, unflexible molecule
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
            gStartGeom.getBondInfo().copy());
    gEndGeom.setFitness(cartes.getEnergy());
    gEndGeom.setFather(gStartGeom.getFatherID());
    gEndGeom.setMother(gStartGeom.getMotherID());
    gEndGeom.setLocalOptimized(true);
    return gEndGeom;
  }

  @Override
  public Molecule localOptimization(Molecule mStartMolecule) {

    incrMol();

    CartesianCoordinates cartes = mStartMolecule.getCartesians();
    final ZMatrix[] zmatVec = cartes.getZMatrices();
    try {
      // since we can't be certain that there is something really "unique" in the molecules ID, we
      // hash it :-)
      int iHashCode = mStartMolecule.hashCode();
      cartes =
          cartesToCartes(
              (long) iHashCode,
              cartes,
              mStartMolecule.getConstraints(),
              mStartMolecule.isConstricted(),
              new SimpleBondInfo(1)); // doesn't matter really
    } catch (Exception e) {
      System.err.println(
          "Error in molpro locopt. Returning non-optimized molecule." + e.toString());
      mStartMolecule.setEnergy(FixedValues.NONCONVERGEDENERGY);
      return mStartMolecule;
    }

    // we anyway set the old zmatrix first in, so that we have something in there
    cartes.setZMatrices(zmatVec);

    if (zmatVec[0] != null) {
      final CartesianCoordinates cartesTemp = cartes.giveMolecularCartes(0, true);
      zmatVec[0] = cartesTemp.calculateZMatrix();
    } // else: nothing needs to happen, unflexible molecule

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

  private void WriteOutput(
      String sMolproInput,
      boolean bLocOpt,
      boolean bGradient,
      double[] daXYZ1D,
      String[] saAtoms,
      int iTotalCharge)
      throws InitIOException {
    int iAtoms = daXYZ1D.length / 3;
    double[][] daXYZ = new double[3][iAtoms];
    System.arraycopy(daXYZ1D, 0, daXYZ[0], 0, iAtoms);
    System.arraycopy(daXYZ1D, iAtoms, daXYZ[1], 0, iAtoms);
    System.arraycopy(daXYZ1D, iAtoms * 2, daXYZ[2], 0, iAtoms);
    switch (whichMethod) {
      case CUSTOMNOLOCOPT:
        // also custom, same behaviour
      case CUSTOM:
        // custom
        final String[] geom = new String[iAtoms + 2];
        geom[0] = "" + iAtoms;
        geom[1] = "";
        for (int i = 0; i < iAtoms; i++) {
          geom[i + 2] =
              saAtoms[i]
                  + "    "
                  + daXYZ[0][i] * Constants.BOHRTOANG
                  + "    "
                  + daXYZ[1][i] * Constants.BOHRTOANG
                  + "    "
                  + daXYZ[2][i] * Constants.BOHRTOANG;
        }
        final String[] out = glue(sCustomInp, geom);
        try {
          OutputPrimitives.writeOut(sMolproInput, out, false);
        } catch (Exception e) {
          throw new InitIOException(e);
        }
        break;
      case AM1:
        try {
          Output.WriteMolproInput(
              sMolproInput, bLocOpt, bGradient, "semi,am1", "vdz", daXYZ, saAtoms, iTotalCharge);
        } catch (InitIOException e) {
          throw e;
        }
        ;
        break;
      case HFVDZ:
        try {
          Output.WriteMolproInput(
              sMolproInput, bLocOpt, bGradient, "rhf", "vdz", daXYZ, saAtoms, iTotalCharge);
        } catch (InitIOException e) {
          throw e;
        }
        ;
        break;
      case B86VDZ:
        try {
          Output.WriteMolproInput(
              sMolproInput, bLocOpt, bGradient, "ks,b86", "vdz", daXYZ, saAtoms, iTotalCharge);
        } catch (InitIOException e) {
          throw e;
        }
        ;
        break;
      case MP2AVTZ:
        try {
          Output.WriteMolproInput(
              sMolproInput, bLocOpt, bGradient, "mp2", "aug-cc-pVTZ", daXYZ, saAtoms, iTotalCharge);
        } catch (InitIOException e) {
          throw e;
        }
        ;
        break;
      default:
        try {
          Output.WriteMolproInput(
              sMolproInput, bLocOpt, bGradient, "semi,am1", "vdz", daXYZ, saAtoms, iTotalCharge);
        } catch (InitIOException e) {
          throw e;
        }
        ;
        break;
    }
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      long lID,
      CartesianCoordinates cartes,
      boolean[][] constraints,
      boolean isConstricted,
      final BondInfo bonds)
      throws ConvergenceException {
    String sMolproInput = "molpro" + lID + ".inp";

    /*
     * write the input file
     */
    try {
      WriteOutput(
          sMolproInput,
          true,
          false,
          cartes.getAll1DCartes(),
          cartes.getAllAtomTypes(),
          cartes.getTotalCharge());
    } catch (InitIOException e) {
      throw new ConvergenceException(
          "Problem in writing geometry for molpro input (local optimisation).", e);
    }

    float[] faCharges = cartes.getAllCharges();
    short[] iaSpins = cartes.getAllSpins();

    /*
     * call molpro
     */
    try {
      Runtime rt = Runtime.getRuntime();
      String sMolproCmd = System.getenv("OGO_MOLPROCMD");
      if (sMolproCmd == null) {
        sMolproCmd = "molpro";
      }
      Process proc = rt.exec(sMolproCmd + " " + sMolproInput);

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
            "Molpro returns non-zero return value (local optimisation).");
      } else {
        // molpro should(!) have completed normally...
      }

    } catch (Exception e) {
      throw new ConvergenceException("Molpro has a problem (local optimisation).", e);
    }

    /*
     * read molpros output
     */
    if (whichMethod != METHOD.CUSTOMNOLOCOPT) {
      String sMolproOutput = sMolproInput.substring(0, sMolproInput.indexOf(".")) + ".log";
      try {
        cartes =
            Input.ReadXYZMolproOutput(
                sMolproOutput,
                cartes.getNoOfAtoms(),
                cartes.getNoOfMolecules(),
                cartes.getAllAtomsPerMol());
      } catch (Exception e) {
        throw new ConvergenceException("Problem in reading the xyz output of molpro.", e);
      }
    } else {
      String sMolproOutput = sMolproInput.substring(0, sMolproInput.indexOf(".")) + ".out";
      double energy = FixedValues.NONCONVERGEDENERGY;
      try {
        energy = Input.ReadEnergyMolproOutput(sMolproOutput);
      } catch (Exception e) {
        System.err.println("WARNING: Problem reading molpros energy!" + e.toString());
      }
      cartes.setEnergy(energy);
    }

    /*
     * clean up
     */
    try {
      Input.RemoveMolproFiles(sMolproInput, true);
    } catch (Exception e) {
      System.err.println("Problem cleaning molpro files up." + e.toString());
    }

    /*
     * return it
     */
    cartes.setAllCharges(faCharges);
    cartes.setAllSpins(iaSpins);

    return cartes;
  }

  private static String readCustomInp(final String file) {

    String[] sa = null;
    try {
      sa = Input.ReadFile(file);
    } catch (Exception e) {
      System.err.println(
          "ERROR: Couldn't read in reference input. This will fail! " + e.getMessage());
      return "";
    }

    String stot = "";
    for (final String s : sa) {
      stot += s + System.getProperty("line.separator");
    }

    return stot;
  }

  private static String[] glue(final String ref, final String[] geom) {

    // chop the ref into two
    final String[] refs = ref.split("OGOCOORDS");

    if (refs.length != 2) {
      System.err.println("ERROR: Reference input wrong. Returning empty input.");
      return new String[0];
    }

    // put the geom in between
    final String[] out = new String[2 + geom.length];
    out[0] = refs[0];
    int count = 1;
    for (final String s : geom) {
      out[count] = s;
      count++;
    }
    out[count] = refs[1];

    return out;
  }

  @Override
  public CartesianFullBackend getBackend() {
    return this;
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
}
