/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2012, J. M. Dieterich
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
import org.ogolem.core.CastException;
import org.ogolem.core.InitIOException;
import org.ogolem.generic.genericpool.GenericPool;

/**
 * Entry point for designing some new switches.
 *
 * @author Johannes Dieterich
 * @version 2012-06-12
 */
public class MainSwitches {

  public static void main(String[] argv) {

    int iNoOfThreads = 0;
    String sInputFile = "";

    if (argv.length < 1) {
      System.err.println("ERROR: We want both an input file and the number of threads. Aborting.");
      System.exit(22);
    } else if (argv.length != 2) {
      System.err.println("WARNING: No number of threads specified. Taking all cores now.");
      iNoOfThreads = Runtime.getRuntime().availableProcessors();
      sInputFile = argv[0];
    } else {
      // parse things
      try {
        iNoOfThreads = Integer.parseInt(argv[1]);
      } catch (Exception e) {
        System.err.println(
            "ERROR: Couldn't parse the number of threads "
                + e.toString()
                + ". "
                + "Taking all cores now.");
        iNoOfThreads = Runtime.getRuntime().availableProcessors();
      }
      sInputFile = argv[0];
    }

    // read the configuration file in
    SwitchesConfig swConfig = null;
    try {
      swConfig = SwitchesInput.readConfigIn(sInputFile);
    } catch (CastException e1) {
      System.err.println("ERROR: Can't configure myself. Aborting.");
      e1.printStackTrace(System.err);
      System.exit(111);
    } catch (InitIOException e2) {
      System.err.println("ERROR: Can't read/write somewhere. Aborting.");
      e2.printStackTrace(System.err);
      System.exit(112);
    } catch (Exception e3) {
      // and take care of all other cases
      System.err.println("ERROR: Something went wrong in configuring myself. Aborting.");
      e3.printStackTrace(System.err);
      System.exit(113);
    }

    // now do some "afterconfiguring"
    final ArrayList<Color> alColors = new ArrayList<>();

    final ColorPalette palette = ColorPalette.getReference();

    // yeah, this is pretty nasty. but this is no performance penalty
    for (int i = 0; i < SwitchesConfig.backbone.getConnectsCopy(true).length; i++) {
      alColors.add(palette.getRandomColor());
    }

    final Switch refSwitch = new Switch(new Backbone(SwitchesConfig.backbone), alColors);
    final GenericPool<Color, Switch> pool =
        GenericPool.getInstance(swConfig.getGenericConfig(), refSwitch);

    // do the initialization of the pool
    final ThreadingInits inits = new ThreadingInits(swConfig, iNoOfThreads, pool);
    inits.doTheInitialPoolFill(SwitchesConfig.iPoolSize, refSwitch);

    // let people know about intermediate results
    for (int i = 0; i < SwitchesConfig.iPoolSize; i++) {
      final Switch sw = pool.getIndividualAtPosition(i);

      final String[][] saCisTrans = sw.createPrintableCisTrans(swConfig.dBlowBondsFac);
      final String[] saColorFile = sw.createPrintableColors();

      final String sFileBase = "initrank" + i + "switch" + sw.getID();
      final String sColorFile = sFileBase + "-colors.col";
      final String sCisFile = sFileBase + "-cis.xyz";
      final String sTransFile = sFileBase + "-trans.xyz";

      // actually write it out
      try {
        Output.printMiscToFile(sColorFile, saColorFile);
        Output.printMiscToFile(sCisFile, saCisTrans[0]);
        Output.printMiscToFile(sTransFile, saCisTrans[1]);
      } catch (Exception e) {
        System.err.println(
            "ERROR: Failed to write results for initrank "
                + i
                + " continuing anyway."
                + e.toString());
      }
    }

    // do the global optimization
    final ThreadingGlobOpt globopt = new ThreadingGlobOpt(swConfig, iNoOfThreads, pool);
    globopt.doAllGlobOpts();

    // let people know about the results
    for (int i = 0; i < SwitchesConfig.iPoolSize; i++) {
      final Switch sw = pool.getIndividualAtPosition(i);

      final String[][] saCisTrans = sw.createPrintableCisTrans(swConfig.dBlowBondsFac);
      final String[] saColorFile = sw.createPrintableColors();

      final String sFileBase = "rank" + i + "switch" + sw.getID();
      final String sColorFile = sFileBase + "-colors.col";
      final String sCisFile = sFileBase + "-cis.xyz";
      final String sTransFile = sFileBase + "-trans.xyz";

      // actually write it out
      try {
        Output.printMiscToFile(sColorFile, saColorFile);
        Output.printMiscToFile(sCisFile, saCisTrans[0]);
        Output.printMiscToFile(sTransFile, saCisTrans[1]);
      } catch (Exception e) {
        System.err.println(
            "ERROR: Failed to write results for rank " + i + " continuing anyway. " + e.toString());
      }
    }

    // serialize the pool a last time
    try {
      Output.writeBinaryPool(pool, false);
    } catch (Exception e) {
      System.err.println("ERROR: Failes to serialize the final pool! " + e.toString());
    }
  }
}
