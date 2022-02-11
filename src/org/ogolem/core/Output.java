/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2014-2020, J. M. Dieterich and B. Hartke
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

import java.io.IOException;
import org.ogolem.io.OutputPrimitives;

/**
 * All kind of output related features using the primitive I/O classes.
 *
 * @author Johannes Dieterich
 * @version 2022-02-05
 */
public final class Output {

  public static String[] getHeader() {

    final String[] header =
        new String[] {
          "",
          "",
          "THIS IS OGOLEM",
          "A framework for genetic algorithms-based global optimization.",
          "",
          "",
          "     OOO OOO    GGG GGGG   OOO OOO   LLLL       EEEE EEEE  MMM   MMM",
          "    OOOO OOOO  GGGG  GGG  OOOO OOOO  LLLL       EEEE  EEE  MMMM MMMM",
          "    OOOO OOOO  GGGG       OOOO OOOO  LLLL       EEEE E     MMMM MMMM",
          "    OOOO OOOO  GGGG GGGG  OOOO OOOO  LLLL  LLL  EEEE  EEE  MMMM MMMM",
          "     OOO OOO    GGG GGGG   OOO OOO   LLLL LLLL  EEEE EEEE  MMMM MMMM",
          "",
          "THE OGOLEM.ORG PROJECT",
          "Developed at the Institute for Physical Chemistry, University of Kiel,",
          "          at the Institute for Physical Chemistry, University of Goettingen and",
          "          at the Department of Mechanical&Aerospace Engineering, Princeton University",
          "",
          "",
          "",
          "Copyright (c) 2009-2022, J. M. Dieterich and B. Hartke",
          "All rights reserved.",
          "",
          "Redistribution and use in source and binary forms, with or without",
          "modification, are permitted provided that the following conditions are met:",
          "",
          "    * Redistributions of source code must retain the above copyright",
          "      notice, this list of conditions and the following disclaimer.",
          "",
          "    * Redistributions in binary form must reproduce the above copyright",
          "      notice, this list of conditions and the following disclaimer in the",
          "      documentation and/or other materials provided with the distribution.",
          "",
          "    * All advertising materials mentioning features or use of this software",
          "      must display the following acknowledgement:",
          "",
          "      This product includes software of the ogolem.org project developed by",
          "      J. M. Dieterich and B. Hartke (Christian-Albrechts-University Kiel, Germany)",
          "      and contributors.",
          "",
          "    * Neither the name of the ogolem.org project, the University of Kiel",
          "      nor the names of its contributors may be used to endorse or promote products",
          "      derived from this software without specific prior written permission.",
          "",
          "THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ''AS IS'' AND ANY",
          "EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED",
          "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE",
          "DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY",
          "DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES",
          "(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;",
          "LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND",
          "ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT",
          "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS",
          "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.",
          "",
          "",
          "With contributions from:",
          "* Dominik Behrens",
          "* Mark Dittner",
          "",
          "OGOLEM uses (and wants to thank):",
          "* Atomdroid (J. Feldt, R. A. Mata and J. M. Dieterich)",
          "  license: proprietary",
          "* Lionbench (M. Kammler and J. M. Dieterich)",
          "  license: 4-clause BSD",
          "* L-BFGS from RISO (J. Nocedal/R. Dodier)",
          "  license: GPL2 / public domain",
          "* BOBYQA (Apache Software Foundation)",
          "  license: Apache",
          "* NEWUOA (code by J. Feldt and J. M. Dieterich, original M. J. D. Powell)",
          "  license: public domain",
          "* JAMA",
          "* Apache Commons Math",
          "* Scala",
          "* Numal",
          "* Arpack/Netlib/MTJ/EJML",
          "* SLF4J",
          "",
          "OGOLEM was partly inspired by:",
          "* phenix v0-42",
          "* Tinker",
          ""
        };

    return header;
  }

  public static void printXYZSingleGeometry(
      final Geometry geometry, final boolean withEnv, final String folder) throws IOException {

    final boolean localOpt = geometry.isLocalOptimized();

    String out;
    if (localOpt) {
      out =
          folder
              + System.getProperty("file.separator")
              + "geometry"
              + geometry.getID()
              + "locopt.xyz";
    } else {
      out = folder + System.getProperty("file.separator") + "geometry" + geometry.getID() + ".xyz";
    }

    final String[] cart = geometry.makePrintableAbsoluteCoord(withEnv);

    OutputPrimitives.writeOut(out, cart, true);
  }

  static void WriteMolproInput(
      final String sMolproInput,
      final boolean bLocalOpt,
      final boolean bGradient,
      final String sMethod,
      final String sBasis,
      final double[][] daXYZ,
      final String[] saAtoms,
      final int iTotalCharge)
      throws InitIOException {
    final int iNoAtoms = saAtoms.length;
    final int iLinesOfInput = iNoAtoms + 13;
    final String[] saInput = new String[iLinesOfInput];
    // the first line...
    saInput[0] = "***,OGOLEM generated input";

    /*
     * calculate the amount of memory needed
     */
    int iHowManyWords = iNoAtoms;
    if (sMethod.equalsIgnoreCase("mp2")) {
      iHowManyWords *= 3;
    } else if (sMethod.equalsIgnoreCase("hf")) {
      iHowManyWords *= 2;
    } else if (sMethod.equalsIgnoreCase("ks,bp96")) {
      iHowManyWords *= 2;
    } else {
      // ok...
    }
    if (sBasis.equalsIgnoreCase("vdz")) {
      iHowManyWords *= 3;
    } else if (sBasis.equalsIgnoreCase("aug-cc-pVTZ")) {
      iHowManyWords *= 2;
    } else {
      // ok....
    }
    iHowManyWords /= 8;
    if (iHowManyWords < 128) {
      iHowManyWords = 128;
    } else {
      iHowManyWords += 20;
    }
    saInput[1] = "memory," + iHowManyWords + ",m";

    /*
     * geometry
     */
    saInput[2] = "geomtyp=xyz";
    saInput[3] = "geometry={";
    saInput[4] = "nosym,bohr";
    saInput[5] = " " + iNoAtoms;
    saInput[6] = "This was created by OGOLEM.";

    for (int i = 7; i < iNoAtoms + 7; i++) {
      saInput[i] =
          saAtoms[i - 7] + "," + daXYZ[0][i - 7] + "," + daXYZ[1][i - 7] + "," + daXYZ[2][i - 7];
    }
    saInput[iNoAtoms + 7] = "}";
    saInput[iNoAtoms + 8] = "basis=" + sBasis;

    /*
     * charge
     */
    if (iTotalCharge != 0) {
      saInput[iNoAtoms + 9] = "set,charge=" + iTotalCharge;
    } else {
      saInput[iNoAtoms + 9] = "";
    }

    /*
     * methods
     */
    if (sMethod.equalsIgnoreCase("mp2")) {
      saInput[iNoAtoms + 10] = "hf";
      saInput[iNoAtoms + 11] = "mp2";
    } else {
      saInput[iNoAtoms + 10] = sMethod;
      saInput[iNoAtoms + 11] = "";
    }

    /*
     * gradient? locopt? anyone?
     */
    if (bGradient) {
      saInput[iNoAtoms + 12] = "forces";
    } else if (bLocalOpt) {
      saInput[iNoAtoms + 12] = "optg";
    } else {
      saInput[iNoAtoms + 12] = "";
    }

    // write it actually
    try {
      OutputPrimitives.writeOut(sMolproInput, saInput, false);
    } catch (IOException e) {
      throw new InitIOException(e);
    }
  }
}
