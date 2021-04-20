/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
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
package org.ogolem.beswitched;

import java.io.File;
import java.io.IOException;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.switches.Color;
import org.ogolem.switches.Switch;
import org.ogolem.switches.SwitchesInput;

/**
 * Allows for intermediate and final analyzation of an ogolem.switches run.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class MainBeswitched {

  /**
   * The actual entry point to analyzation. -i inputfile (default: switchespool.bin) -collblow blow
   * factor for the dissociation detection (default: 1.4) -getswitches reads all switches out of a
   * pool (default: off)
   *
   * @param args
   */
  public static void main(String[] args) {

    String sPoolBinFile = "switchespool.bin";
    String sSwitchesFolder = "switches";

    boolean bGetSwitches = false;

    double dBlowCollDetect = 1.4;

    for (int i = 0; i < args.length; i++) {
      if (args[i].equalsIgnoreCase("-i")) {
        // the next entry in the array contains the inputfile
        if (args[i + 1].startsWith("-")) {
          System.err.println(
              "ERROR: You need to specify an input file after the -i switch. Continuing.");
        } else {
          sPoolBinFile = args[i + 1];
        }
        i++;
        continue;
      } else if (args[i].equalsIgnoreCase("-collblow")) {
        // the next entry in the array contains the blow factor
        if (args[i + 1].startsWith("-")) {
          System.err.println(
              "ERROR: You need to specify a blow factor after "
                  + "the -collblow switch. Continuing.");
        } else {
          try {
            dBlowCollDetect = Double.parseDouble(args[i + 1]);
          } catch (Exception e) {
            System.err.println(
                "ERROR: Failure to parse the blow factor for "
                    + "the collision detection. Continuing with default. "
                    + e.toString());
          }
        }
        i++;
        continue;
      } else if (args[i].equalsIgnoreCase("-getswitches")) {
        bGetSwitches = true;
        continue;
      } else {
        // if we end up here, there was a non valid option
        System.err.println("ERROR: Unrecognized option " + args[i] + " Continuing.");
      }
    }

    // read the pool in
    GenericPool<Color, Switch> pool = null;
    try {
      pool = SwitchesInput.readSwitchPool(sPoolBinFile);
    } catch (Exception e) {
      System.err.println("ERROR: Can't read SwitchPopulation in. " + e.toString());
    }

    if (pool == null) {
      System.exit(42);
    }

    if (bGetSwitches) {

      // create the switches folder
      try {
        File folder = new File("switches");
        if (!folder.exists()) {
          if (!folder.mkdir()) {
            throw new IOException("Something went wrong.");
          }
        }
      } catch (Exception e) {
        System.err.println(
            "ERROR Not capable of creating folder switches. Exiting. " + e.toString());
        System.exit(1);
      }

      // get the switches and print them out
      for (int i = 0; i < pool.getPoolSize(); i++) {
        final Switch sw = pool.getIndividualAtPosition(i);

        final String[][] saCisTrans = sw.createPrintableCisTrans(dBlowCollDetect);
        final String[] saColorFile = sw.createPrintableColors();

        final String sFileBase =
            sSwitchesFolder
                + System.getProperty("file.separator")
                + "rank"
                + i
                + "switch"
                + sw.getID();

        final String sColorFile = sFileBase + "-colors.col";
        final String sCisFile = sFileBase + "-cis.xyz";
        final String sTransFile = sFileBase + "-trans.xyz";

        // actually write it out
        try {
          WriteNow(sColorFile, saColorFile);
          WriteNow(sCisFile, saCisTrans[0]);
          WriteNow(sTransFile, saCisTrans[1]);
        } catch (Exception e) {
          System.err.println(
              "ERROR: Failed to write results for rank "
                  + i
                  + " continuing anyway. "
                  + e.toString());
        }
      }
    }
  }

  static void WriteNow(String sOutputPath, String[] saToWrite) throws IOException {
    OutputPrimitives.writeOut(sOutputPath, saToWrite, true);
  }
}
