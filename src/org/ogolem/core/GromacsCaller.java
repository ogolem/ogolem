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

import java.io.File;
import java.util.Locale;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * This calls Gromacs for local optimization. Supports QM/MM hybrid schemes. Configuration mainly
 * done externally.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
class GromacsCaller extends AbstractLocOpt {

  // the ID
  private static final long serialVersionUID = (long) 20150421;

  static enum QMENGINE {
    NONE,
    MOPAC,
    Orca,
    Gaussian,
    GamessUK
  };

  private final String[] groFirstHalf;
  private final String topFile;
  private final String paramFile;
  private final String ndxFile;

  private final QMENGINE whichQM;
  private final double[] box;

  GromacsCaller(GlobalConfig globconf, QMENGINE whichQMEngine) {
    super(globconf);

    final String sAuxFile = globconf.outputFolder + "-gromacs.aux";
    String[] sa = null;
    double[] da = null;
    try {
      sa = Input.ReadFile(sAuxFile);
      if (!sa[0].trim().equalsIgnoreCase("###OGOLEMAUX###")) {
        throw new Exception("Wrong aux file in Gromacs.");
      }

      final String[] sa2 = sa[1].split("\\s");
      da = new double[3];
      try {
        for (int i = 0; i < 3; i++) {
          da[i] = Double.parseDouble(sa2[i]);
        }
      } catch (Exception e) {
        System.err.println("INFO: Box usage disabled!");
        da = null;
      }
    } catch (Exception e) {
      System.err.println(
          "ERROR: Not capable of reading gromacs aux file " + sAuxFile + " in. " + e.toString());
      sa = new String[0];
    }
    this.box = da;
    this.groFirstHalf = sa;
    this.topFile = "gromacs-topology.top";
    this.paramFile = "gromacs-parameters.mdp";
    this.ndxFile = "gromacs-index.ndx";
    this.whichQM = whichQMEngine;
  }

  private GromacsCaller(final GromacsCaller orig) {
    super(orig);
    this.groFirstHalf = orig.groFirstHalf.clone();
    this.topFile = orig.topFile;
    this.paramFile = orig.paramFile;
    this.ndxFile = orig.ndxFile;
    if (orig.box != null) {
      this.box = orig.box.clone();
    } else {
      this.box = null;
    }
    this.whichQM = orig.whichQM;
  }

  @Override
  public GromacsCaller copy() {
    return new GromacsCaller(this);
  }

  @Override
  public String myID() {
    return "Gromacs";
  }

  @Override
  public String myIDandMethod() {
    switch (whichQM) {
      case NONE:
        return "Gromacs";
      case MOPAC:
        return "Gromacs with MOPAC";
      case Orca:
        return "Gromacs with Orca";
      case Gaussian:
        return "Gromacs with Gaussian";
      case GamessUK:
        return "Gromacs with Gamess-UK";
      default:
        return "Gromacs unknown QM";
    }
  }

  @Override
  public Molecule localOptimization(Molecule mStartMolecule) {

    System.out.println(
        "WARNING: Gromacs caller does not support molecule optimization. Returning non-optimized molecule.");
    mStartMolecule.setEnergy(FixedValues.NONCONVERGEDENERGY);

    return mStartMolecule;
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      final long lID,
      final CartesianCoordinates cartes,
      boolean[][] constraints,
      boolean isConstricted,
      final BondInfo bonds)
      throws ConvergenceException {

    final String sGromacsFolder = "gromacs" + lID;

    final float[] faCharges = cartes.getAllCharges();
    final short[] iaSpins = cartes.getAllSpins();

    // create the folder
    try {
      OutputPrimitives.createAFolder(sGromacsFolder);
    } catch (Exception e) {
      throw new ConvergenceException(e);
    }

    // copy the topology and parameters
    try {
      OutputPrimitives.createLink(
          topFile, sGromacsFolder + System.getProperty("file.separator") + topFile);
      OutputPrimitives.createLink(
          paramFile, sGromacsFolder + System.getProperty("file.separator") + paramFile);
      OutputPrimitives.createLink(
          ndxFile, sGromacsFolder + System.getProperty("file.separator") + ndxFile);
    } catch (Exception e) {
      throw new ConvergenceException(e);
    }

    // depending upon what QM part we want to use, switch statement
    switch (whichQM) {
        // only orca needs special attention
      case Orca:
        setupForOrca(sGromacsFolder);
        break;
      default:
        break;
    }

    // set the input coordinates up
    cartes.moveCoordsToCOM();
    final String[] coords = createG96File(cartes.getAllXYZCoord(), groFirstHalf, box);
    final String coordFile = sGromacsFolder + System.getProperty("file.separator") + "coords.g96";
    try {
      OutputPrimitives.writeOut(coordFile, coords, false);
    } catch (Exception e) {
      throw new ConvergenceException(e);
    }

    // minimize
    doTheMinimization(sGromacsFolder, "coords.g96");

    // read the output
    final CartesianCoordinates localOptCartes = readXYZGromacsOutput(sGromacsFolder, cartes);

    // remove all traces
    try {
      ManipulationPrimitives.remove(sGromacsFolder);
    } catch (Exception e) {
      System.err.println("Failed to delete gromacs folder." + e.toString());
    }

    localOptCartes.setAllAtomTypes(cartes.getAllAtomTypes());
    localOptCartes.setAllCharges(faCharges);
    localOptCartes.setAllSpins(iaSpins);

    final String[] sa = localOptCartes.createPrintableCartesians();
    for (String s : sa) {
      System.out.println(s);
    }

    return localOptCartes;
  }

  /*    private static void initialSetUp(String sGromacsFilePrefix) throws ConvergenceException {
      // spawn the thread running pdb2gmx to create the topology
      try{
          Runtime rt = Runtime.getRuntime();
          String[] saEnvp = new String[0];
          File dir = new File(sGromacsFilePrefix);
          String sCommand = "pdb2gmx -f geometry.pdb -o conf.gro -p topol.top";
          Process proc = rt.exec(sCommand,saEnvp,dir);

          // any error message?
          StreamGobbler errorGobbler = new
              StreamGobbler(proc.getErrorStream(), "ERROR");

          // any output?
          StreamGobbler outputGobbler = new
              StreamGobbler(proc.getInputStream(), "OUTPUT");

          // kick them off
          errorGobbler.start();
          outputGobbler.start();

          // any error???
          int iExitValue = proc.waitFor();
          if(iExitValue != 0){
              throw new ConvergenceException("pdb2gmx returns non-zero return value (initial input creation).");
          } else{
              // pdb2gmx should(!) have completed normally...
          }

      } catch(Exception e){
          throw new ConvergenceException("pdb2gmx has a problem (initial input creation).",e);
      }
  }

  private static void solvateWithWater() throws ConvergenceException{
      // add solvation using gromacs genbox
          // editconf and genbox, I guess
          /*
           * editconf -f speptide -o -d 0.5
           * genbox -cp out -cs -p speptide -o b4em
           */
  // }

  private void doTheMinimization(final String folder, final String g96File)
      throws ConvergenceException {

    // spawn the thread running first grompp and then mdrun
    try {
      final Runtime rt = Runtime.getRuntime();
      final String[] saEnvp = null;
      final File dir = new File(folder);

      String grompp = System.getenv("OGO_GROMACSMPP");
      if (grompp == null) grompp = "grompp";

      final String[] cmd1 =
          new String[] {grompp, "-f", paramFile, "-c", g96File, "-p", topFile, "-n", ndxFile};

      Process proc = rt.exec(cmd1, saEnvp, dir);

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
            "grompp returns non-zero return value (initial input creation).");
      }

      String mdrun = System.getenv("OGO_GROMACSCMD");
      if (mdrun == null) mdrun = "mdrun";

      final String[] cmd2 = new String[] {mdrun, "-c", "coords.g96", "-nt", "1"};

      proc = rt.exec(cmd2, saEnvp, dir);

      // any error message?
      errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      // any output?
      outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      // kick them off
      errorGobbler.start();
      outputGobbler.start();

      // any error???
      iExitValue = proc.waitFor();
      if (iExitValue != 0) {
        throw new ConvergenceException("mdrun returns non-zero return value (minimization).");
      }

    } catch (Exception e) {
      throw new ConvergenceException("gromacs has a problem.", e);
    }
  }

  private static String[] createG96File(
      final double[][] xyz, final String[] firstHalf, final double[] box) {

    int length;
    if (box == null) {
      length = firstHalf.length + 3;
    } else {
      length = firstHalf.length + 6;
    }
    final String[] g96 = new String[length];
    g96[0] = "TITLE";
    g96[1] = "OGOLEM generated coordinates";
    g96[2] = "END";
    g96[3] = "POSITION";
    for (int i = 0; i < firstHalf.length - 2; i++) {
      // coordinates are in nm, not angstrom, file is fixed-format *grrr*
      g96[i + 4] =
          firstHalf[i + 2]
              + String.format(Locale.US, "%9d", (i + 1))
              + String.format(Locale.US, "%15.9f", xyz[0][i] * Constants.BOHRTOANG * 0.1)
              + String.format(Locale.US, "%15.9f", xyz[1][i] * Constants.BOHRTOANG * 0.1)
              + String.format(Locale.US, "%15.9f", xyz[2][i] * Constants.BOHRTOANG * 0.1);
    }
    g96[firstHalf.length + 2] = "END";
    if (box != null) {
      g96[g96.length - 3] = "BOX";
      g96[g96.length - 2] =
          String.format(Locale.US, "%15.9f", box[0])
              + String.format(Locale.US, "%15.9f", box[1])
              + String.format(Locale.US, "%15.9f", box[2]);
      g96[g96.length - 1] = "END";
    }

    return g96;
  }

  private static CartesianCoordinates readXYZGromacsOutput(
      final String folder, final CartesianCoordinates refCartes) throws ConvergenceException {

    // we can't just copy here since the refCartes does not contain any ZMatrix info e.g.
    final CartesianCoordinates optCartes =
        new CartesianCoordinates(
            refCartes.getNoOfAtoms(),
            refCartes.getNoOfMolecules(),
            refCartes.getAllAtomsPerMol().clone());

    String[] saCoords = null;
    String[] saLogFile = null;
    try {
      saCoords = Input.ReadFile(folder + System.getProperty("file.separator") + "coords.g96");
      saLogFile = Input.ReadFile(folder + System.getProperty("file.separator") + "md.log");
    } catch (Exception e) {
      throw new ConvergenceException(e);
    }

    // coordinates
    final double[][] xyz = optCartes.getAllXYZCoord();
    for (int i = 4; i < xyz[0].length + 4; i++) {
      final String s = saCoords[i];
      // fixed-formatting once more...
      final String s1 = s.substring(25, 39).trim();
      final String s2 = s.substring(40, 54).trim();
      final String s3 = s.substring(55).trim();
      try {
        // remember: nanometer!
        xyz[0][i - 4] = Double.parseDouble(s1) * 10.0 * Constants.ANGTOBOHR;
        xyz[1][i - 4] = Double.parseDouble(s2) * 10.0 * Constants.ANGTOBOHR;
        xyz[2][i - 4] = Double.parseDouble(s3) * 10.0 * Constants.ANGTOBOHR;
      } catch (Exception e) {
        throw new ConvergenceException(e);
      }
    }

    // energy
    double energy = FixedValues.NONCONVERGEDENERGY;
    for (int i = saLogFile.length - 1; i >= 0; i--) {
      if (saLogFile[i].startsWith("Potential Energy  =")) {
        final String sEnergy = saLogFile[i].substring(19).trim();
        try {
          energy = Double.parseDouble(sEnergy) * Constants.KCALTOHARTREE;
        } catch (Exception e) {
          throw new ConvergenceException(e);
        }
        break;
      }
    }

    if (energy == FixedValues.NONCONVERGEDENERGY)
      throw new ConvergenceException("No or too high energy found in Gromacs output.");
    optCartes.setEnergy(energy);

    return optCartes;
  }

  private void setupForOrca(final String folder) throws ConvergenceException {

    final String orcaRef = "gromacs-orca.aux";
    final String target = folder + System.getProperty("file.separator") + "topol.ORCAINFO";
    try {
      OutputPrimitives.createLink(orcaRef, target);
    } catch (Exception e) {
      throw new ConvergenceException(e);
    }
  }
}
