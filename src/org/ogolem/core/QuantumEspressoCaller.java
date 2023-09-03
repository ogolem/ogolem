/*
Copyright (c) 2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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

import java.util.Iterator;
import java.util.LinkedList;
import org.ogolem.io.OutputPrimitives;

/**
 * Calls the Quantum Espresso program suite for local optimizations.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
final class QuantumEspressoCaller extends AbstractLocOpt {

  // the ID
  private static final long serialVersionUID = (long) 20101108;

  QuantumEspressoCaller(final GlobalConfig globconf) {
    super(globconf);
  }

  // TODO everything

  private QuantumEspressoCaller(final QuantumEspressoCaller orig) {
    super(orig);
  }

  @Override
  public QuantumEspressoCaller copy() {
    return new QuantumEspressoCaller(this);
  }

  @Override
  public String myIDandMethod() {
    return "Quantum Espresso";
  }

  @Override
  public String myID() {
    return "Quantum Espresso";
  }

  @Override
  public CartesianCoordinates cartesToCartes(
      final long id,
      final CartesianCoordinates cartes,
      boolean[][] constraints,
      boolean isConstricted,
      final BondInfo bonds)
      throws Exception {

    // file/folder setup (might throw exception)
    final String sFolder = "quantumespresso" + id;
    OutputPrimitives.createAFolder(sFolder);

    // create program input
    final String[] saOutput;

    // call program

    // parse output

    return null;
  }

  private static String[] createOutput(
      final String[] atoms,
      final double[][] xyz,
      final double convThresh,
      final float[] cell,
      final String[] potentials,
      final String outdir,
      final double beta,
      final double cutWave,
      final double cutRho,
      final int charge) {

    final LinkedList<String> llOutput = new LinkedList<>();
    final int iNoOfAtoms = atoms.length;

    final LinkedList<String> llTypes = new LinkedList<>();
    for (String atom : atoms) {
      if (!llTypes.contains(atom)) {
        llTypes.add(atom);
      }
    }

    final int iNoOfTypes = llTypes.size();
    final String sWorkDir = System.getProperty("user.dir");
    final String sOutDir = sWorkDir + System.getProperty("file.separator") + outdir;
    final String sPseudoDir = sWorkDir + System.getProperty("file.separator") + "pseudo";

    // general system setup
    llOutput.add("&CONTROL");
    llOutput.add("calculation = \"relax\" ");
    llOutput.add("prefix = espresso");
    llOutput.add("pseudo_dir = " + sPseudoDir);
    llOutput.add("outdir = " + sOutDir);
    llOutput.add("/");

    // actual system setup
    llOutput.add("&SYSTEM");
    llOutput.add(" ibrav = 0"); // lattice index "free"
    llOutput.add(" nat = " + iNoOfAtoms);
    llOutput.add(" ntyp = " + iNoOfTypes);
    llOutput.add(" ecutwfc = " + cutWave);
    llOutput.add(" ecutrho = " + cutRho);
    llOutput.add(" tot_charge = " + charge);
    llOutput.add("/");

    // electrons stuff
    llOutput.add("&ELECTRONS");
    llOutput.add("  conv_thr = " + convThresh);
    llOutput.add("  mixing_beta = " + beta);
    llOutput.add("/");

    // ions
    llOutput.add("&IONS");
    llOutput.add("/");

    // cell (assuming cubic)
    llOutput.add("CELL_PARAMETERS cubic");
    llOutput.add(cell[0] + " 0.0  0.0");
    llOutput.add("0.0  " + cell[1] + " 0.0");
    llOutput.add("0.0  0.0  " + cell[2]);

    // pseudopotentials
    llOutput.add("ATOMIC_SPECIES");
    // TODO

    // coordinates
    llOutput.add("ATOMIC_POSITIONS {bohr}");
    for (int i = 0; i < iNoOfAtoms; i++) {
      llOutput.add("" + atoms[i] + "  " + xyz[0][i] + "  " + xyz[1][i] + "  " + xyz[2][i]);
    }

    // k-point
    llOutput.add("K_POINTS {gamma}");

    // copy to String[]
    final String[] saOutput = new String[llOutput.size()];
    final Iterator<String> it = llOutput.iterator();
    int iCount = 0;
    while (it.hasNext()) {
      saOutput[iCount] = it.next();
      iCount++;
    }

    return saOutput;
  }
}
