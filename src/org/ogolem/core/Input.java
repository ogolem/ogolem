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

import static org.ogolem.core.Constants.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericSanityCheck;
import org.ogolem.generic.IndividualWriter;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.helpers.Tuple;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.InquiryPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.md.MDConfig;
import org.ogolem.random.Lottery;
import org.ogolem.random.RNGenerator;
import org.ogolem.random.StandardRNG;

/**
 * This invokes the IOHandler and checks the resulting array of Strings for configuration options
 * and configures a GlobalConfig object.
 *
 * @author Johannes Dieterich
 * @version 2020-07-29
 */
public final class Input {

  public static GlobalConfig ConfigureMe(final String inputPath)
      throws InitIOException, CastException, Exception {
    return ConfigureMe(inputPath, true);
  }

  public static GlobalConfig ConfigureMe(final String inputPath, final boolean failOnMissingGeom)
      throws InitIOException, CastException, Exception {

    String[] data = null;
    try {
      data = InputPrimitives.readFileIn(inputPath);
    } catch (Exception e) {
      System.err.println("Failed to read input file!");
      throw new InitIOException(e);
    }

    if (!data[0].trim().equalsIgnoreCase("###OGOLEM###")) {
      throw new IOException(
          "Please specify a correct input file, denoted by ###OGOLEM### in the beginning!");
    }

    return configureFromBlock(data, inputPath, failOnMissingGeom);
  }

  public static GlobalConfig configureFromBlock(
      final String[] configInput, final String inputPath, final boolean failOnMissingGeom)
      throws InitIOException, CastException, Exception {

    final GlobalConfig globConf = new GlobalConfig();

    final Tuple3D<String, String, String> dirs =
        ManipulationPrimitives.outDirAndBaseName(inputPath);
    final String inputDir = dirs.getObject1();
    final String outputFolder = dirs.getObject2();
    final String baseName = dirs.getObject3();

    final String outputFile = outputFolder + File.separator + baseName + ".out";
    final String logFile = outputFolder + File.separator + baseName + ".log";

    globConf.outputFile = outputFile;
    globConf.outputFolder = outputFolder;

    globConf.stats = new GlobOptStatistics(logFile);

    String locoptString = "lbfgs:backend=xyz:universal";
    String initialLocOpt = null;
    String globOptString =
        "cluster{xover(sweden:cutstyle=1)mutation(germany:)}molecules{xover(germany:)mutation(germany:)}";
    String diversityString = "fitnessbased:1e-6";
    String rngString = "javarng:autoseed";
    int geomStart = -1;
    int geomEnd = -1;
    boolean hasChargeTag = false;
    int chargeStart = -1;
    int chargeEnd = -1;
    boolean hasSpinTag = false;
    int spinStart = -1;
    int spinEnd = -1;
    boolean hasDoFTag = false;
    int dofStart = -1;
    int dofEnd = -1;
    boolean hasConstrTag = false;
    int constrStart = -1;
    int constrEnd = -1;
    boolean hasEnvTag = false;
    int envStart = -1;
    int envEnd = -1;
    final List<Tuple<Integer, Integer>> backendDefTags = new LinkedList<>();

    MainParser:
    for (int i = 0; i < configInput.length; i++) {
      final String line = configInput[i];
      if (line.startsWith("//") || line.startsWith("#")) {
        // comment line: ignore it
        continue;
      } else if (line.equalsIgnoreCase("<GEOMETRY>")) {
        // found the begin of the geom environment, forward to the end of it
        geomStart = i;
        for (int j = i + 1; j < configInput.length; j++) {
          i++;
          if (configInput[j].trim().equalsIgnoreCase("</GEOMETRY>")) {
            geomEnd = j;
            continue MainParser;
          }
        }
        throw new Exception("No end </GEOMETRY> tag found.");
      } else if (line.equalsIgnoreCase("<CHARGES>")) {
        // found the begin of this environment, forward to the end of it
        hasChargeTag = true;
        chargeStart = i;
        for (int j = i + 1; j < configInput.length; j++) {
          i++;
          if (configInput[j].trim().equalsIgnoreCase("</CHARGES>")) {
            chargeEnd = j;
            continue MainParser;
          }
        }
        throw new Exception("No end </CHARGES> tag found.");
      } else if (line.equalsIgnoreCase("<SPINS>")) {
        // found the begin of this environment, forward to the end of it
        hasSpinTag = true;
        spinStart = i;
        for (int j = i + 1; j < configInput.length; j++) {
          i++;
          if (configInput[j].trim().equalsIgnoreCase("</SPINS>")) {
            spinEnd = j;
            continue MainParser;
          }
        }
        throw new Exception("No end </SPINS> tag found.");
      } else if (line.equalsIgnoreCase("<DOF>")) {
        // found the begin of this environment, forward to the end of it
        hasDoFTag = true;
        dofStart = i;
        for (int j = i + 1; j < configInput.length; j++) {
          i++;
          if (configInput[j].trim().equalsIgnoreCase("</DOF>")) {
            dofEnd = j;
            continue MainParser;
          }
        }
        throw new Exception("No end </DOF> tag found.");
      } else if (line.equalsIgnoreCase("<CONSTRAINTS>")) {
        // found the begin of this environment, forward to the end of it
        hasConstrTag = true;
        constrStart = i;
        for (int j = i + 1; j < configInput.length; j++) {
          i++;
          if (configInput[j].trim().equalsIgnoreCase("</CONSTRAINTS>")) {
            constrEnd = j;
            continue MainParser;
          }
        }
        throw new Exception("No end </CONSTRAINTS> tag found.");
      } else if (line.equalsIgnoreCase("<ENVIRONMENT>")) {
        // found the begin of this environment, forward to the end of it
        hasEnvTag = true;
        envStart = i;
        for (int j = i + 1; j < configInput.length; j++) {
          i++;
          if (configInput[j].trim().equalsIgnoreCase("</ENVIRONMENT>")) {
            envEnd = j;
            continue MainParser;
          }
        }
        throw new Exception("No end </ENVIRONMENT> tag found.");
      } else if (line.equalsIgnoreCase("<ADAPTIVE>")) {
        // found the begin of the adaptive environment, ONLY forward to the end of it to avoid
        // getting screamed at
        for (int j = i + 1; j < configInput.length; j++) {
          i++;
          if (configInput[j].trim().equalsIgnoreCase("</ADAPTIVE>")) {
            continue MainParser;
          }
        }
        throw new Exception("No end </ADAPTIVE> tag found.");
      } else if (line.equalsIgnoreCase("<MOLECULARDYNAMICS>")) {
        // found the begin of the MD environment, ONLY forward to the end of it to avoid getting
        // screamed at
        for (int j = i + 1; j < configInput.length; j++) {
          i++;
          if (configInput[j].trim().equalsIgnoreCase("</MOLECULARDYNAMICS>")) {
            continue MainParser;
          }
        }
        throw new Exception("No end </MOLECULARDYNAMICS> tag found.");
      } else if (line.startsWith("<CLUSTERBACKEND>")) {
        final int start = i;
        while (i < configInput.length) {
          if (configInput[i].trim().equalsIgnoreCase("</CLUSTERBACKEND>")) {
            final int end = i;
            final Tuple<Integer, Integer> tags = new Tuple<>(start, end);
            backendDefTags.add(tags);
            break;
          }
          i++;
        }
      } else if (line.startsWith("PoolSize=")) {
        final String sTemp2 = line.substring(9).trim();
        try {
          int iPoolSize = Integer.parseInt(sTemp2);
          if (iPoolSize < 0) {
            throw new Exception("No negative pool sizes allowed.");
          }
          globConf.poolSize = iPoolSize;
        } catch (Exception e) {
          System.err.println("Wrong input for PoolSize: " + e.toString() + "default of 100 used.");
        }
      } else if (line.startsWith("ClusterDetailedStats=")) {
        try {
          final boolean b = Boolean.parseBoolean(line.substring(21));
          globConf.enableDetailedStats = b;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for ClusterDetailedStats: " + e.toString() + ". NOT ENABLED!");
        }
      } else if (line.startsWith("DebugLevel=")) {
        final String s2 = line.substring(11).trim();
        try {
          final int deb = Integer.parseInt(s2);
          if (deb < 0) {
            throw new Exception("No negative debug levels allowed.");
          }
          GlobalConfig.DEBUGLEVEL = deb;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for DebugLevel: " + e.toString() + "default of 0 used. " + e.toString());
        }
      } else if (line.startsWith("GeneticRecordBufferSize=")) {
        final String s2 = line.substring(24).trim();
        try {
          final int recs = Integer.parseInt(s2);
          if (recs < 0) {
            throw new Exception("No negative buffer sizes for genetic records allowed.");
          }
          globConf.geneticRecordsToASCII = recs;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for GeneticRecordBufferSize: "
                  + e.toString()
                  + "default of 1000 used. "
                  + e.toString());
        }
      } else if (line.startsWith("restart=")) {
        final String sTemp2 = line.substring(8).trim();
        try {
          boolean bRestart = Boolean.parseBoolean(sTemp2);
          globConf.restart = bRestart;
        } catch (Exception e) {
          System.err.println("Wrong input for restart: " + e.toString() + " default false used.");
        }
      } else if (line.startsWith("GeneticRecordsToSerial=")) {
        final String sTemp2 = line.substring(23).trim();
        try {
          int iRecordsToSerial = Integer.parseInt(sTemp2);
          globConf.geneticRecordsToSerial = iRecordsToSerial;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for GeneticRecordsToSerial: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("GeometriesToSerial=")) {
        String sTemp2 = line.substring(19).trim();
        try {
          int iGeometriesToSerial = Integer.parseInt(sTemp2);
          globConf.geometriesToSerial = iGeometriesToSerial;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for GeometriesRecordsToSerial: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("SerializeAfterNewBest=")) {
        final String s2 = line.substring(22).trim();
        try {
          globConf.serializePoolAfterBest = Boolean.parseBoolean(s2);
        } catch (Exception e) {
          System.err.println(
              "Wrong input for SerializeAfterNewBest: " + e.toString() + ". Using default.");
        }
      } else if (line.startsWith("SilentMode=")) {
        final String s2 = line.substring(11).trim();
        try {
          globConf.silentMode = Boolean.parseBoolean(s2);
        } catch (Exception e) {
          System.err.println("Wrong input for SilentMode: " + e.toString() + ". Using default.");
        }
      } else if (line.startsWith("WriteEveryGeometry=")) {
        final String s2 = line.substring(19).trim();
        try {
          globConf.writeEveryGeom = Boolean.parseBoolean(s2);
        } catch (Exception e) {
          System.err.println(
              "Wrong input for WriteEveryGeometry: " + e.toString() + ". Using default.");
        }
      } else if (line.startsWith("AcceptableFitness=")) {
        final String sTemp2 = line.substring(18).trim();
        try {
          double dAcceptable = Double.parseDouble(sTemp2);
          globConf.acceptableFitness = dAcceptable;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for AcceptableFitness: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("DeleteOldGeoms=")) {
        final String sTemp2 = line.substring(15).trim();
        try {
          boolean bDeleteOldGeoms = Boolean.parseBoolean(sTemp2);
          globConf.deleteOldGeoms = bDeleteOldGeoms;
        } catch (Exception e) {
          System.err.println("Wrong input for DeleteOldGeoms: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("PrintGeomsBeforeLocOpt=")) {
        final String s2 = line.substring(23).trim();
        try {
          final boolean b = Boolean.parseBoolean(s2);
          globConf.printBeforeLocOpt = b;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for PrintGeomsBeforeLocOpt: " + e.toString() + ". Using default.");
        }
      } else if (line.startsWith("RandomNumberGenerator=")) {
        rngString = line.substring("RandomNumberGenerator=".length()).trim();
      } else if (line.startsWith("MolecularCoordinateMutation=")) {
        final String sTemp2 = line.substring(28).trim();
        switch (sTemp2) {
          case "some":
            globConf.whichMolecularMutation = GlobalConfig.MOLMUTCHOICE.SOME;
            break;
          case "one":
            globConf.whichMolecularMutation = GlobalConfig.MOLMUTCHOICE.ONE;
            break;
          default:
            throw new Exception(
                "Illegal choice " + sTemp2 + ". Please specify either some or one.");
        }
      } else if (line.startsWith("MutationPossibility=")) {
        final String s2 = line.substring(20).trim();
        try {
          final double d = Double.parseDouble(s2);
          globConf.mutatePossibility = d;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for MutationPossibility: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("MolMutationPossibility=")) {
        final String s2 = line.substring("MolMutationPossibility=".length()).trim();
        try {
          final double d = Double.parseDouble(s2);
          globConf.molMutatePoss = d;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for MolMutationPossibility: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("CrossoverPossibility=")) {
        final String s2 = line.substring(21).trim();
        try {
          final double d = Double.parseDouble(s2);
          globConf.crossPossibility = d;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for CrossoverPossibility: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("MolCrossoverPossibility=")) {
        final String s2 = line.substring("MolCrossoverPossibility=".length()).trim();
        try {
          final double d = Double.parseDouble(s2);
          globConf.molCrossPoss = d;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for MolCrossoverPossibility: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("GrowCell=")) {
        final String sTemp2 = line.substring(9).trim();
        try {
          boolean bGrowCell = Boolean.parseBoolean(sTemp2);
          globConf.growCell = bGrowCell;
        } catch (Exception e) {
          System.err.println("Wrong input for GrowCell: " + e.toString() + "default true used.");
        }
      } else if (line.startsWith("GeometryChoice=")) {
        globConf.whichGeomChoice = line.substring(15).trim();
      } else if (line.startsWith("NumberOfGlobIterations=")) {
        final String sTemp2 = line.substring(23).trim();
        try {
          int iNoOfIterations = Integer.parseInt(sTemp2);
          globConf.noOfGlobalSteps = iNoOfIterations;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for NumberOfGlobIterations: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("NamingOffset=")) {
        final String sTemp2 = line.substring(13).trim();
        try {
          int iNamingOffset = Integer.parseInt(sTemp2);
          globConf.namingOffsetGeoms = iNamingOffset;
        } catch (Exception e) {
          System.err.println("Wrong input for NamingOffset: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("DoGeometryNiching=")) {
        final String sTemp2 = line.substring(18).trim();
        globConf.doNiching = Boolean.parseBoolean(sTemp2);
      } else if (line.startsWith("IndividualsPerNicheAtMax=")) {
        final String sTemp2 = line.substring(25).trim();
        globConf.noOfIndividualsPerNicheAtMax = Integer.parseInt(sTemp2);
      } else if (line.startsWith("WhichClusterNicher=")) {
        final String sTemp2 = line.substring(19).trim();
        globConf.nicherString = sTemp2;
      } else if (line.startsWith("NichingAddsToStats=")) {
        final String sTemp2 = line.substring(19).trim();
        try {
          final int n = Integer.parseInt(sTemp2);
          globConf.addsToNicheStats = n;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Wrong input for NichingAddsToStats. " + e.toString() + " Using default.");
        }
      } else if (line.startsWith("DiversityCheck=")) {
        final String token = line.substring(15).trim();
        diversityString = token;
      } else if (line.startsWith("GlobOptAlgo=")) {
        String sTemp2 = line.substring(12).trim();
        globOptString = sTemp2;
      } else if (line.startsWith("GlobOptTries=")) {
        final String sTemp2 = line.substring(13).trim();
        try {
          int iGlobOptTries = Integer.parseInt(sTemp2);
          globConf.howManyTries = iGlobOptTries;
        } catch (Exception e) {
          System.err.println("Wrong input for GlobOptTries: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("LocOptAlgo=")) {

        // will be done later
        locoptString = line.substring(11).trim();

      } else if (line.startsWith("InitialLocOptAlgo=")) {

        // will be done later
        initialLocOpt = line.substring(18).trim();

      } else if (line.startsWith("CollisionDetection=")) {
        String sTemp2 = line.substring(19).trim();
        try {
          globConf.whichCollisionEngine = CollisionDetection.parseType(sTemp2);
        } catch (Exception e) {
          System.err.println(
              "Wrong input for CollisionDetection " + e.toString() + " default used.");
        }
      } else if (line.startsWith("DissociationDetection=")) {
        String sTemp2 = line.substring(22).trim();
        try {
          globConf.whichDissociationEngine = DissociationDetection.parseType(sTemp2);
        } catch (Exception e) {
          System.err.println(
              "Wrong input for DissociationDetection " + e.toString() + " default used.");
        }
      } else if (line.startsWith("InitialFillAlgo=")) {
        final String sTemp2 = line.substring(16).trim();
        try {
          globConf.whichInitialFill = GeometryInit.parseType(sTemp2);
        } catch (Exception e) {
          System.err.println(
              "Wrong input for InitialFillAlgo= " + e.toString() + " using default.");
        }
      } else if (line.startsWith("RatioExplDoFInit=")) {
        final String s2 = line.substring(17).trim();
        try {
          final float f = Float.parseFloat(s2);
          globConf.ratioExplDoFInit = f;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for RatioExplDoFInit= " + e.toString() + " using default.");
        }
      } else if (line.startsWith("MolecularCDInInit=")) {
        final String s2 = line.substring(18).trim();
        try {
          final boolean b = Boolean.parseBoolean(s2);
          globConf.doMolecularCD = b;
        } catch (Exception e) {
          System.err.println(
              "Wrong input for MolecularCDInInit= " + e.toString() + " using default.");
        }
      } else if (line.startsWith("CellSize=")) {
        final String sTemp2 = line.substring(9).trim();
        String[] saTemp = sTemp2.split("\\;", 3);
        for (int j = 0; j < 3; j++) {
          try {
            final double cell1D = Double.parseDouble(saTemp[j]);
            globConf.maxCellSize[j] = cell1D;
          } catch (Exception e) {
            System.err.println("Wrong input in CellSize: " + e.toString() + " default used.");
          }
        }
      } else if (line.startsWith("BlowBondDetect=")) {
        String sTemp2 = line.substring(15).trim();
        try {
          double dBlowBond = Double.parseDouble(sTemp2);
          globConf.blowFacBondDetect = dBlowBond;
        } catch (Exception e) {
          System.err.println("Wrong input in BlowBondDetect: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("MaxBondStretch=")) {
        String sTemp2 = line.substring(15).trim();
        try {
          double dBondStretch = Double.parseDouble(sTemp2);
          globConf.maxBondStretch = dBondStretch;
        } catch (Exception e) {
          System.err.println("Wrong input in MaxBondStretch: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("BlowInitialBonds=")) {
        String sTemp2 = line.substring(17).trim();
        try {
          double dBlowInitialBonds = Double.parseDouble(sTemp2);
          globConf.blowFacInitialBondDetect = dBlowInitialBonds;
        } catch (Exception e) {
          System.err.println("Wrong input in BlowInitialBonds: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("BlowEnvClusterClash=")) {
        String sTemp2 = line.substring("BlowEnvClusterClash=".length()).trim();
        try {
          final double d = Double.parseDouble(sTemp2);
          globConf.blowFacEnvClusterClashes = d;
        } catch (Exception e) {
          System.err.println(
              "Wrong input in BlowEnvClusterClash: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("OptimizationTarget=")) {
        globConf.fitnessFunctionConfig = line.substring(19).trim();
      } else if (line.startsWith("MaxIterLocOpt=")) {
        String sTemp2 = line.substring(14).trim();
        try {
          int iMaxIter = Integer.parseInt(sTemp2);
          globConf.maxIterLocOpt = iMaxIter;
        } catch (Exception e) {
          System.err.println("Wrong input in MaxIterLocOpt: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("SeedingPath=")) {
        final String sTemp2 = line.substring(12).trim();
        globConf.seedingFolder = sTemp2;
      } else if (line.startsWith("ThreshLocOptGradient=")) {
        final String sTemp2 = line.substring(21).trim();
        try {
          double dThresh = Double.parseDouble(sTemp2);
          globConf.threshLocOptGradient = dThresh;
        } catch (Exception e) {
          System.err.println(
              "Wrong input in ThreshLocOptGradient: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("ThreshLocOptCoord=")) {
        String sTemp2 = line.substring(18).trim();
        try {
          double dThresh = Double.parseDouble(sTemp2);
          globConf.threshLocOptCoord = dThresh;
        } catch (Exception e) {
          System.err.println(
              "Wrong input in ThreshLocOptCoord: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("TasksToSubmit=")) {
        final String sTemp2 = line.substring(14).trim();
        try {
          final int tasks = Integer.parseInt(sTemp2);
          globConf.tasksToSubmit = tasks;
        } catch (Exception e) {
          System.err.println("Wrong input in TasksToSubmit: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("IntLocOptMaxStep=")) {
        final String sTemp2 = line.substring(17).trim();
        try {
          double dMaxStep = Double.parseDouble(sTemp2);
          globConf.maxStepSizeLocOpt = dMaxStep;
        } catch (Exception e) {
          System.err.println("Wrong input in IntLocOptMaxStep: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("BlowFacDissoc=")) {
        String sTemp2 = line.substring(14);
        sTemp2 = sTemp2.trim();
        try {
          double dBlowDiss = Double.parseDouble(sTemp2);
          globConf.blowFacDissocDetect = dBlowDiss;
        } catch (Exception e) {
          System.err.println("Wrong input in BlowFacDissoc: " + e.toString() + " default used.");
        }
      } else if (line.startsWith("AdaptiveParameters=")) {
        final String sFile = line.substring(19).trim();
        try {
          final String[] saParams = Input.ReadFile(inputDir + File.separator + sFile);
          final org.ogolem.adaptive.AdaptiveParameters params =
              new org.ogolem.adaptive.AdaptiveParameters(saParams, -1);
          globConf.parameters = params;
        } catch (Exception e) {
          throw new CastException("Couldn't do anything with the parameters, exiting.", e);
        }
      } else if (line.startsWith("Penalty=")) {
        final String s = line.substring(8).trim();
        if (s.startsWith("firstincenter:")) {
          final String s2 = s.substring(14).trim();
          final String[] sa = s2.split("\\;");
          final double[] da = new double[3];
          try {
            for (int ii = 0; ii < 3; ii++) {
              da[ii] = Double.parseDouble(sa[ii].trim());
            }
          } catch (Exception e) {
            System.err.println(
                "ERROR: Problem in setting up firstincenter penalty function. "
                    + e.toString()
                    + ". Using no penalty now.");
          }
          globConf.penalty = new FirstInCenterPenalty(da[0], da[1], da[2]);
        } else {
          System.err.println("ERROR: Couldn't parse Penalty= " + s + ". Using no penalty.");
        }
      } else if (line.startsWith("PostSanityCheck=")) {
        final String sTemp2 = line.substring(16).trim();
        try {
          final boolean b = Boolean.parseBoolean(sTemp2);
          globConf.doPostSanityCheck = b;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Wrong input for PostSanityCheck. " + e.toString() + " Using default.");
        }
      } else if (line.startsWith("PostSanityCD=")) {
        final String sTemp2 = line.substring(13).trim();
        try {
          final boolean b = Boolean.parseBoolean(sTemp2);
          globConf.doPostCD = b;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Wrong input for PostSanityCD. " + e.toString() + " Using default.");
        }
      } else if (line.startsWith("PostSanityDD=")) {
        final String sTemp2 = line.substring(13).trim();
        try {
          final boolean b = Boolean.parseBoolean(sTemp2);
          globConf.doPostDD = b;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Wrong input for PostSanityDD. " + e.toString() + " Using default.");
        }
      } else if (line.startsWith("PreSanityCD=")) {
        final String sTemp2 = line.substring(12).trim();
        try {
          final boolean b = Boolean.parseBoolean(sTemp2);
          globConf.doPreCD = b;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Wrong input for PreSanityCD. " + e.toString() + " Using default.");
        }
      } else if (line.startsWith("PreSanityDD=")) {
        final String sTemp2 = line.substring(12).trim();
        try {
          final boolean b = Boolean.parseBoolean(sTemp2);
          globConf.doPreDD = b;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Wrong input for PreSanityDD. " + e.toString() + " Using default.");
        }
      } else {
        // unknown line: scream!
        throw new Exception("Unkown line " + line + " in global configuration parser.");
      }
    }

    // for adaptive methods, we search for <ADAPTIVE> tags
    boolean hasAdaptiveBlock = false;
    int startAdaptBlock = 0;
    int endAdaptBlock = 0;
    for (int i = 0; i < configInput.length; i++) {
      final String tmp = configInput[i].trim();
      if (tmp.equalsIgnoreCase("<ADAPTIVE>")) {
        startAdaptBlock = i;
        hasAdaptiveBlock = true;

      } else if (tmp.equalsIgnoreCase("</ADAPTIVE>")) {
        endAdaptBlock = i;
      }
    }
    if (hasAdaptiveBlock) {

      // save everything in between tags in another string array
      final int linesInBlock = endAdaptBlock - startAdaptBlock - 1;
      final String[] adaptiveBlock = new String[linesInBlock];
      int c = 0;
      for (int j = startAdaptBlock + 1; j < endAdaptBlock; j++) {
        adaptiveBlock[c] = configInput[j];
        c++;
      }

      // create a configuration object
      org.ogolem.adaptive.AdaptiveConf adaptiveConf =
          new org.ogolem.adaptive.AdaptiveConf(adaptiveBlock, globConf);
      // put the adaptive configuration into the global configuration
      globConf.adaptiveConf = adaptiveConf;
    } // end of adaptive input parsing

    // check for the presence of geometry tags and either exit (if we should not fail on missing
    // geometry definitions) or throw an exception
    if ((geomStart < 0 || geomEnd < 0) && failOnMissingGeom) {
      throw new Exception("Geometry tags somewhat corrupt or (partly) not existing.");
    } else if ((geomStart < 0 || geomEnd < 0) && !failOnMissingGeom) {
      return globConf;
    }

    // save everything in between tags in another string array
    final int geomBlockLength = geomEnd - geomStart;
    if (geomBlockLength <= 0) {
      throw new RuntimeException(
          "GEOMETRY block is of zero or negative length. Check tags for existence, proper order and spelling.");
    }
    final String[] saGeomContent = new String[geomBlockLength];
    int geomC = 0;
    for (int j = geomStart; j < geomEnd; j++) {
      saGeomContent[geomC] = configInput[j];
      geomC++;
    }

    // create the GeometryConfig
    final GeometryConfig gc = new GeometryConfig();
    for (final String saGeomContent1 : saGeomContent) {
      final String s = saGeomContent1.trim();
      if (s.startsWith("NumberOfParticles=")) {
        final String s2 = s.substring(18).trim();
        try {
          final int noParticles = Integer.parseInt(s2);
          gc.noOfParticles = noParticles;
        } catch (Exception e) {
          throw new Exception("Wrong input for NumberOfParticles: " + e.toString() + ".", e);
        }
      }
    }

    // search for various <MOLECULE> tags and put them in an array
    int molCount = 0;
    final int[] molStartLines = new int[gc.noOfParticles];
    final int[] molEndLines = new int[gc.noOfParticles];
    for (int i = 0; i < saGeomContent.length; i++) {
      final String s = saGeomContent[i].trim();
      if (s.equalsIgnoreCase("<MOLECULE>")) {
        molStartLines[molCount] = i;
      } else if (s.equalsIgnoreCase("</MOLECULE>")) {
        molEndLines[molCount] = i;
        molCount++;
      }
    }

    // Create the MoleculeConfig ArrayList
    final List<MoleculeConfig> molConfigs = new ArrayList<>(gc.noOfParticles);
    final List<BondInfo> bondMatrices = new ArrayList<>(gc.noOfParticles);
    for (int i = 0; i < gc.noOfParticles; i++) {
      bondMatrices.add(null);
    }

    // this is important, the total number of atoms
    final int[] molIDMapping = new int[gc.noOfParticles];
    int noAtomsTotal = 0;
    int tagPos = 0;
    for (int molId = 0; molId < gc.noOfParticles; molId++) {

      // create String array containing the information in between the <MOLECULE> tags
      final int noLines = molEndLines[tagPos] - molStartLines[tagPos] - 1;
      if (noLines <= 0) {
        throw new CastException(
            "Molecule tags w/o content or dangling tag. Make sure the NumberOfParticles matches the number of specfied molecules!");
      }
      final List<String> molBlock = new ArrayList<>(noLines);
      int noReps = 1;
      for (int j = molStartLines[tagPos] + 1; j < molEndLines[tagPos]; j++) {
        final String s = saGeomContent[j].trim();
        if (s.startsWith("MoleculeRepetitions=")) {
          final String s2 = s.substring(20).trim();
          noReps = Integer.parseInt(s2);
          // ignore this line
        } else {
          molBlock.add(s);
        }
      }

      for (int molRep = 0; molRep < noReps; molRep++) {

        // Configure the Molecules
        final MoleculeConfig mc = parseMoleculeBlock(molId, molBlock, bondMatrices, inputDir);
        noAtomsTotal += mc.noOfAtoms;

        // fix assertations caused by missing these in toCartesians()
        mc.spins = new short[mc.noOfAtoms];
        mc.charges = new float[mc.noOfAtoms];
        // assign molecular type if not done yet
        if (mc.sID.equalsIgnoreCase("N/A")) {
          mc.sID = mc.toCartesians().createDefaultMolecularType();
        }

        molConfigs.add(molId, mc);
        molIDMapping[molId] = tagPos;
        molId++;
      }

      molId--; // cure off by one
      tagPos++;
    }

    if (molConfigs.size() != gc.noOfParticles) {
      throw new Exception(
          "Mismatch of NumberOfParticles ("
              + gc.noOfParticles
              + ") and read number of particles ("
              + molConfigs.size()
              + ").");
    }

    /*
     * next thing: deal with other tags
     */

    /*
     * work on the charges
     */
    if (hasChargeTag) {

      if (chargeStart == -1 || chargeEnd == -1) {
        throw new CastException("Charge tags somewhat corrupt.");
      }

      // save everything in between tags in another string array
      final int diff = chargeEnd - chargeStart;
      if (diff < 0) {
        throw new RuntimeException("Dangling charge tag. Line difference is negative.");
      }
      final String[] saChargeContent = new String[diff - 1];
      int chargeC = 0;
      for (int j = chargeStart + 1; j < chargeEnd; j++) {
        saChargeContent[chargeC] = configInput[j];
        chargeC++;
      }

      final int[] atsPerMol = new int[gc.noOfParticles];
      for (int i = 0; i < gc.noOfParticles; i++) {
        atsPerMol[i] = molConfigs.get(i).noOfAtoms;
      }

      final float[][] charges =
          readChargesIn(saChargeContent, gc.noOfParticles, atsPerMol, molIDMapping);

      // put the charges in the right spot
      for (int i = 0; i < gc.noOfParticles; i++) {
        // there is a certain beauty in aggressive non-copying... ;-)
        molConfigs.get(i).charges = charges[i];
      }

    } else {
      // no charges, so why bothering acting on them? ;-)

      // initialize the charge arrays
      for (int i = 0; i < gc.noOfParticles; i++) {
        molConfigs.get(i).charges = new float[molConfigs.get(i).noOfAtoms];
      }
    }

    /*
     * next thing: spins
     */
    if (hasSpinTag) {

      if (spinStart == -1 || spinEnd == -1) {
        throw new CastException("Spin tags somewhat corrupt.");
      }

      // save everything in between tags in another string array
      final int iLineDiff = spinEnd - spinStart;
      if (iLineDiff < 0) {
        throw new RuntimeException("Dangling spin tag. Line difference is negative.");
      }
      final String[] saChargeContent = new String[iLineDiff - 1];
      int iSpinCounter = 0;
      for (int j = spinStart + 1; j < spinEnd; j++) {
        saChargeContent[iSpinCounter] = configInput[j];
        iSpinCounter++;
      }

      final int[] iaAtomsPerMol = new int[gc.noOfParticles];
      for (int i = 0; i < gc.noOfParticles; i++) {
        iaAtomsPerMol[i] = molConfigs.get(i).noOfAtoms;
      }

      final short[][] iaSpins =
          readSpinsIn(saChargeContent, gc.noOfParticles, iaAtomsPerMol, molIDMapping);

      // put the charges in the right spot
      for (int i = 0; i < gc.noOfParticles; i++) {
        // there is a certain beauty in aggressive non-copying... ;-)
        molConfigs.get(i).spins = iaSpins[i];
      }

    } else {
      // no spins, so why bothering acting on them? ;-)

      // initialize the spin arrays
      for (int i = 0; i < gc.noOfParticles; i++) {
        molConfigs.get(i).spins = new short[molConfigs.get(i).noOfAtoms];
      }
    }

    /*
     * next thing: degrees of freedom
     */
    if (hasDoFTag) {

      if (dofStart == -1 || dofEnd == -1) {
        throw new CastException("DoFs tags somewhat corrupt.");
      }

      // save everything in between tags in another string array
      int iLineDiff = dofEnd - dofStart;
      if (iLineDiff < 0) {
        throw new RuntimeException("Dangling DoF tag. Line difference is negative.");
      }
      String[] saDOFContent = new String[iLineDiff - 1];
      int iDOFCounter = 0;
      for (int j = dofStart + 1; j < dofEnd; j++) {
        saDOFContent[iDOFCounter] = configInput[j];
        iDOFCounter++;
      }

      for (final String dofLine : saDOFContent) {
        final String[] tmp = dofLine.split("\\;");
        // remove the whitespace in front and in the back
        for (int j = 0; j < tmp.length; j++) {
          tmp[j] = tmp[j].trim();
        }
        if (tmp.length == 2) {

          // get info which molecule is meant
          int iWhichMol = 0;
          try {
            iWhichMol = Integer.parseInt(tmp[0]);
          } catch (Exception e) {
            System.err.println("Failure in casting String to int in degrees of freedom.");
            throw new CastException(e);
          }

          if (tmp[1].equalsIgnoreCase("all")) {
            final boolean[][] dofs = molConfigs.get(iWhichMol).degreesOfFreedom;
            for (final boolean[] dof : dofs) {
              dof[0] = true;
              dof[1] = true;
              dof[2] = true;
            }
          } else {
            throw new InitIOException("Don't know what to do with " + tmp[1] + ".");
          }

        } else if (tmp.length == 3) {

          // get info which molecule, atom, coordinate is meant
          int iWhichMol = 0;
          int iWhichAtom = 0;
          int iCoord = 0;
          try {
            iWhichMol = Integer.parseInt(tmp[0]);
            iWhichAtom = Integer.parseInt(tmp[1]);
            switch (tmp[2]) {
              case "distance":
                iCoord = 0;
                break;
              case "angle":
                iCoord = 1;
                break;
              case "dihedral":
                iCoord = 2;
                break;
              default:
                throw new CastException(
                    "Illegal coordinate specification "
                        + tmp[2]
                        + " please specify one of: bond/angle/dihedral.");
            }
          } catch (Exception e) {
            System.err.println("Failure in casting String to int in degrees of freedom.");
            throw new CastException(e);
          }

          boolean found = false;
          for (int x = 0; x < molIDMapping.length; x++) {
            if (iWhichMol == molIDMapping[x]) {
              molConfigs.get(x).degreesOfFreedom[iWhichAtom][iCoord] = true;
              found = true;
            }
          }

          if (!found) {
            throw new CastException(
                "Charge parsing: Mol ID "
                    + iWhichMol
                    + " not existing. Range is 0 to "
                    + molIDMapping[molIDMapping.length - 1]);
          }
        } else {
          throw new InitIOException("Wrong number of identifiers for DOF input.");
        }
      }
    }

    /*
     * next thing: constraints
     */
    if (hasConstrTag) {

      if (constrStart == -1 || constrEnd == -1) {
        throw new CastException("Constraints tags somehwat corrupt.");
      }

      // save everything in between tags in another string array
      int iLineDiff = constrEnd - constrStart;
      if (iLineDiff < 0) {
        throw new RuntimeException("Dangling constraints tag. Line difference is negative.");
      }
      String[] saConstrContent = new String[iLineDiff - 1];
      int iConstrCounter = 0;
      for (int j = constrStart + 1; j < constrEnd; j++) {
        saConstrContent[iConstrCounter] = configInput[j];
        iConstrCounter++;
      }

      ConstrBlockLoop:
      for (final String saConstrContent1 : saConstrContent) {
        final String[] saTemp = saConstrContent1.split("\\;");
        // remove the whitespace in front and in the back
        for (int j = 0; j < saTemp.length; j++) {
          saTemp[j] = saTemp[j].trim();
        }
        // get info which molecule is meant
        int iWhichMol = 0;
        try {
          iWhichMol = Integer.parseInt(saTemp[0]);
        } catch (Exception e) {
          System.err.println("Failure in casting String to int in constraints.");
          throw new CastException(e);
        }
        // get the molecule config
        boolean found = false;
        MolMapLoop:
        for (int x = 0; x < molIDMapping.length; x++) {

          if (iWhichMol != molIDMapping[x]) {
            continue MolMapLoop;
          }
          found = true;

          final MoleculeConfig mc = molConfigs.get(iWhichMol);
          mc.constricted = true;
          if (mc.constraints == null) {
            mc.constraints = new boolean[3][mc.noOfAtoms];
          }

          if (saTemp[1].equalsIgnoreCase("all")) {
            for (int j = 0; j < mc.noOfAtoms; j++) {
              mc.constraints[0][j] = true;
              mc.constraints[1][j] = true;
              mc.constraints[2][j] = true;
            }
          } else if (saTemp[1].equalsIgnoreCase("allx")) {
            for (int j = 0; j < mc.noOfAtoms; j++) {
              mc.constraints[0][j] = true;
            }
          } else if (saTemp[1].equalsIgnoreCase("ally")) {
            for (int j = 0; j < mc.noOfAtoms; j++) {
              mc.constraints[1][j] = true;
            }
            continue;
          } else if (saTemp[1].equalsIgnoreCase("allz")) {
            for (int j = 0; j < mc.noOfAtoms; j++) {
              mc.constraints[2][j] = true;
            }
          } else {

            int iWhichAtom = 0;
            try {
              iWhichAtom = Integer.parseInt(saTemp[1]);
            } catch (Exception e) {
              System.err.println("Failure in casting String to int in constraints.");
              throw new CastException(e);
            }

            if (iWhichAtom >= mc.noOfAtoms) {
              throw new CastException("Atom choice insane in constraints.");
            }

            if (saTemp[2].equalsIgnoreCase("all")) {
              mc.constraints[0][iWhichAtom] = true;
              mc.constraints[1][iWhichAtom] = true;
              mc.constraints[2][iWhichAtom] = true;
              continue;
            } else {

              int iWhichCoord = 0;
              try {
                iWhichCoord = Integer.parseInt(saTemp[2]);
              } catch (Exception e) {
                System.err.println("Failure in casting String to int in constraints.");
                throw new CastException(e);
              }

              if (iWhichCoord > 2) {
                throw new CastException("Coordinate choice insane in constraints.");
              }

              mc.constraints[iWhichCoord][iWhichAtom] = true;
            }
          }
        }
        if (!found) {
          throw new CastException(
              "Constraints parsing: Mol ID "
                  + iWhichMol
                  + " not existing. Range is 0 to "
                  + molIDMapping[molIDMapping.length - 1]);
        }
      }
      /*            if(false){
      final boolean[][] constr = molConfigs.get(0).constraints;
      for(int i = 0; i< constr[0].length; i++){
      System.out.println("DEBUG: " + i + "\t" + constr[0][i] + "\t" + constr[1][i] + "\t" + constr[2][i]);
      }
      }*/
    }

    /*
     * INITIALIZE THE RANDOM NUMBER GENERATOR here,
     * since environments refer to it
     */
    RNGenerator rng = null;
    if (rngString.equalsIgnoreCase("javarng:autoseed")) {
      // get seed first
      final Random r = new Random();
      final long seed = r.nextLong();
      rng = new StandardRNG(seed);
    } else if (rngString.startsWith("javarng:seed=")) {
      final String s = rngString.substring("javarng:seed=".length()).trim();
      final long seed = Long.parseLong(s);
      rng = new StandardRNG(seed);
    } else {
      throw new RuntimeException("Illegal RNG configured: " + rngString);
    }

    System.out.println("INFO: RNG initialized: " + rng.getInformation());

    Lottery.setGenerator(rng);

    /*
     * work on environment tags
     */
    if (hasEnvTag) {

      if (envEnd == -1 || envStart == -1) {
        throw new CastException("Environment tags somehwat corrupt.");
      }

      // read the stuff in between the tags in
      final int iLineDiff = envEnd - envStart;
      if (iLineDiff < 0) {
        throw new RuntimeException("Dangling environment tag. Line difference is negative.");
      }
      final String[] saEnvContent = new String[iLineDiff - 1];
      int iEnvCounter = 0;
      for (int j = envStart + 1; j < envEnd; j++) {
        saEnvContent[iEnvCounter] = configInput[j];
        iEnvCounter++;
      }

      CartesianCoordinates envCartes = null;
      AllowedSpace space = null;
      EnvironmentFactory.KIND envKind = EnvironmentFactory.KIND.SIMPLE;
      List<Integer> listSecAtoms = new ArrayList<>();
      boolean bFlexyEnv = false;
      double dBlowFacInitEnvBonds = globConf.blowFacInitialBondDetect;
      int[] iaRefAtoms = new int[3];

      for (final String s : saEnvContent) {
        if (s.startsWith("//") || s.startsWith("#")) {
          // just go on, it's a comment
          continue;
        } else if (s.startsWith("EnvironmentCartes=")) {
          final String sCartes = s.substring(18).trim();

          // try to read the cartes in
          try {
            envCartes = readCartesFromFile(inputDir + File.separator + sCartes);
          } catch (Exception e) {
            System.err.println(
                "ERROR: Unsuccessful parsing of environment cartesians. " + e.toString());
            envCartes = null;
          }
        } else if (s.startsWith("EnvironmentKind=")) {
          final String s2 = s.substring(16).trim();

          if (s2.equalsIgnoreCase("simpleenvironment")) {
            envKind = EnvironmentFactory.KIND.SIMPLE;
          } else if (s2.equalsIgnoreCase("layerenvironment")) {
            envKind = EnvironmentFactory.KIND.LAYERONLY;
          } else if (s2.equalsIgnoreCase("surface")) {
            envKind = EnvironmentFactory.KIND.SURFACE;
          } else if (s2.equalsIgnoreCase("simplesurface")) {
            envKind = EnvironmentFactory.KIND.SIMPLESURFACE;
          } else {
            System.err.println(
                "ERROR: Can't recognize the environmental "
                    + "choice "
                    + s2
                    + " using SimpleEnvironment.");
            envKind = EnvironmentFactory.KIND.SIMPLE;
          }
        } else if (s.equalsIgnoreCase("flexyEnv")) {
          bFlexyEnv = true;
        } else if (s.startsWith("BlowFacEnvInit=")) {
          final String s2 = s.substring(15).trim();
          try {
            dBlowFacInitEnvBonds = Double.parseDouble(s2);
          } catch (Exception e) {
            System.err.println(
                "ERROR: Can't cast double choice for "
                    + "the initial environment blow factor. "
                    + e.toString()
                    + " Using default.");
          }
        } else if (s.startsWith("SecondaryAtoms=")) {
          final String s2 = s.substring(15).trim();
          final String[] saSecs = s2.split("\\;");
          for (String sSec : saSecs) {
            try {
              Integer iSec = Integer.parseInt(sSec);
              listSecAtoms.add(iSec);
            } catch (Exception e) {
              System.err.println(
                  "ERROR: Can't parse integer input" + " of secondary atoms. " + e.toString());
            }
          }
        } else if (s.startsWith("ReferencePoints=")) {
          final String s2 = s.substring(16).trim();
          final String[] sa = s2.split("\\;");
          try {
            for (int i = 0; i < 3; i++) {
              iaRefAtoms[i] = Integer.parseInt(sa[i]);
            }
          } catch (Exception e) {
            System.err.println("ERROR: Can't parse reference atoms. Using 0/0/0. " + e.toString());
            iaRefAtoms[0] = iaRefAtoms[1] = iaRefAtoms[2] = 0;
            continue;
          }

        } else if (s.startsWith("AllowedSpace=")) {
          final String s2 = s.substring(13);
          if (s2.startsWith("parallelepiped:")) {
            final String s3 = s2.substring(15);
            try {
              // parse the 4 corners, p.d. separated by an ;
              final String[] sa = s3.split("\\;");
              final double[][] corners = new double[4][3];
              for (int i = 0; i < 4; i++) {
                // now comma separated
                final String[] sa2 = sa[i].split("\\,");
                for (int j = 0; j < 3; j++) {
                  corners[i][j] = Double.parseDouble(sa2[j].trim()) * ANGTOBOHR;
                }
              }
              space = new ParallelepipedSpace(corners);
            } catch (Exception e) {
              System.err.println(
                  "ERROR: Exception occured when trying to set up parallelepiped. " + e.toString());
              space = null;
            }
          } else if (s2.startsWith("spherespace:")) {
            final String s3 = s2.substring(12).trim();
            try {
              // parse the middle and the radius
              final String[] sa = s3.split("\\;");
              final double[] middle = new double[3];
              // now comma separated
              final String[] sa2 = sa[0].split("\\,");
              for (int i = 0; i < 3; i++) {
                // now comma separated
                middle[i] = Double.parseDouble(sa2[i].trim()) * ANGTOBOHR;
              }
              double radius = Double.parseDouble(sa[1].trim()) * ANGTOBOHR;
              space = new SphereSpace(middle, radius);
            } catch (Exception e) {
              System.err.println(
                  "ERROR: Exception occured when trying to set up sphere. " + e.toString());
              space = null;
            }
          } else if (s2.startsWith("halfspherespace:")) {
            final String s3 = s2.substring(16).trim();
            try {
              // parse the middle and the radius
              final String[] sa = s3.split("\\;");
              final double[] middle = new double[3];
              // now comma separated
              final String[] sa2 = sa[0].split("\\,");
              for (int i = 0; i < 3; i++) {
                // now comma separated
                middle[i] = Double.parseDouble(sa2[i].trim()) * ANGTOBOHR;
              }
              double radius = Double.parseDouble(sa[1].trim()) * ANGTOBOHR;
              space = new HalfSphereSpace(middle, radius);
            } catch (Exception e) {
              System.err.println(
                  "ERROR: Exception occured when trying to set up half sphere. " + e.toString());
              space = null;
            }
          } else if (s2.startsWith("orbitspace:")) {
            final String s3 = s2.substring(11).trim();
            try {
              // parse the middle and the two orbs
              final String[] sa = s3.split("\\;");
              final double[] middle = new double[3];
              // now comma separated
              final String[] sa2 = sa[0].split("\\,");
              for (int i = 0; i < 3; i++) {
                middle[i] = Double.parseDouble(sa2[i].trim()) * ANGTOBOHR;
              }

              double lowOrb = Double.parseDouble(sa[1].trim()) * ANGTOBOHR;
              double highOrb = Double.parseDouble(sa[2].trim()) * ANGTOBOHR;
              space = new OrbitSpace(middle, lowOrb, highOrb);
            } catch (Exception e) {
              System.err.println(
                  "ERROR: Exception occured when trying to set up orbit. " + e.toString());
              space = null;
            }
          } else {
            System.err.println("ERROR: No such space: " + s2);
          }
        } else {
          // we do NOT recognize this input
          System.err.println("ERROR: Can't understand input " + s + ". Ignoring it.");
        }
      }

      // now just use the above collected information
      if (envCartes == null) {
        // not OK
        System.err.println(
            "ERROR: No cartesian coordinates for the "
                + "evironment construction exist. No environment can be "
                + "used therefore.");
      } else {
        if ((iaRefAtoms[0] == iaRefAtoms[1])
            || (iaRefAtoms[0] == iaRefAtoms[2])
            || (iaRefAtoms[1] == iaRefAtoms[2])) {
          System.err.println(
              "ERROR: Reference atoms for environment either " + "unset or bullshit");
          throw new CastException("Wrong references.");
        }

        if (space == null) {
          System.err.println("ERROR: Without allowed space, this is nonsense.");
          throw new CastException("No allowed space defined.");
        }

        // create the array
        final Atom[] references = new Atom[3];
        final AtomConfig ac1 = new AtomConfig();
        ac1.sID = envCartes.getAtomType(iaRefAtoms[0]);
        ac1.atomNo = AtomicProperties.giveAtomicNumber(ac1.sID);
        ac1.iID = iaRefAtoms[0];
        ac1.position = envCartes.getXYZCoordinatesOfAtom(iaRefAtoms[0]);

        final AtomConfig ac2 = new AtomConfig();
        ac2.sID = envCartes.getAtomType(iaRefAtoms[1]);
        ac2.atomNo = AtomicProperties.giveAtomicNumber(ac2.sID);
        ac2.iID = iaRefAtoms[1];
        ac2.position = envCartes.getXYZCoordinatesOfAtom(iaRefAtoms[1]);

        final AtomConfig ac3 = new AtomConfig();
        ac3.sID = envCartes.getAtomType(iaRefAtoms[2]);
        ac3.atomNo = AtomicProperties.giveAtomicNumber(ac3.sID);
        ac3.iID = iaRefAtoms[2];
        ac3.position = envCartes.getXYZCoordinatesOfAtom(iaRefAtoms[2]);

        references[0] = new Atom(ac1);
        references[1] = new Atom(ac2);
        references[2] = new Atom(ac3);

        // everything fine, let's construct the environment
        final Environment env =
            EnvironmentFactory.createEnvironment(
                envKind,
                envCartes,
                bFlexyEnv,
                listSecAtoms,
                dBlowFacInitEnvBonds,
                globConf.whichCollisionEngine,
                globConf.whichDissociationEngine,
                globConf.blowFacEnvClusterClashes,
                globConf.blowFacDissocDetect,
                references,
                space);

        gc.env = env;
      }
    }

    gc.geomMCs = molConfigs;
    if (molConfigs == null || molConfigs.isEmpty()) {
      System.err.println(
          "WARNING: No or an empty <GEOMETRY> definition is ONLY allowed if you are running adaptive.");
    }

    globConf.geoConf = gc;

    // get the bond information
    final BondInfo bonds = new SimpleBondInfo(noAtomsTotal);

    double dBlowBond = globConf.blowFacInitialBondDetect;

    int iEndSet = 0;

    for (int i = 0; i < gc.noOfParticles; i++) {

      final MoleculeConfig mc = gc.geomMCs.get(i);
      BondInfo molBonds = bondMatrices.get(i);
      if (molBonds == null) {
        molBonds = CoordTranslation.checkForBonds(mc.toCartesians(), dBlowBond);
      }

      final int iOffSet = iEndSet;
      iEndSet += mc.noOfAtoms;

      for (int j = iOffSet; j < iEndSet; j++) {
        for (int k = iOffSet; k < iEndSet; k++) {
          bonds.setBond(j, k, molBonds.bondType(j - iOffSet, k - iOffSet));
        }
      }
    }

    gc.bonds = bonds;

    // for molecular dynamics, we search for <MOLECULARDYNAMICS> tags
    boolean mdFlag = false;
    int startMDLine = -1;
    int endMDLine = -1;
    for (int i = 0; i < configInput.length; i++) {
      final String sTemp = configInput[i].trim();
      if (sTemp.equalsIgnoreCase("<MOLECULARDYNAMICS>")) {
        startMDLine = i;
        mdFlag = true;

      } else if (sTemp.equalsIgnoreCase("</MOLECULARDYNAMICS>")) {
        endMDLine = i;
      }
    }

    if (mdFlag) {
      // get the content in between the tags
      final String[] mdCont = new String[endMDLine - startMDLine - 1];
      System.arraycopy(configInput, startMDLine + 1, mdCont, 0, endMDLine - startMDLine - 1);

      // parse the stuff
      globConf.mdConf = MDConfig.parseConfig(mdCont, globConf.parameters, globConf.outputFolder);
    }

    HashMap<String, GenericBackend<Molecule, Geometry>> backendDefs = null; // this is totally fine!
    if (!backendDefTags.isEmpty()) {

      backendDefs = new HashMap<>(); // initialize

      // construct all the backend functions, since these are embedded into tags, it makes adding
      // more involved stuff easier, I hope
      for (final Tuple<Integer, Integer> t : backendDefTags) {

        String tag = null;
        String backend = null;
        String paramFile = null;

        for (int x = t.getObject1() + 1; x < t.getObject2(); x++) {
          final String line = configInput[x].trim();
          if (line.startsWith("#") || line.startsWith("//")) {
            continue;
          } else if (line.startsWith("BackendTag=")) {
            tag = line.substring(11).trim();
          } else if (line.startsWith("Backend=")) {
            backend = line.substring(8).trim();
          } else if (line.startsWith("Parameters=")) {
            paramFile = line.substring(11).trim();
          } else {
            throw new RuntimeException(
                "Unknown line: " + line + " in fitness function definition block.");
          }
        }

        if (tag == null) {
          throw new RuntimeException("No tag specified in backend definition block.");
        }
        if (backend == null) {
          throw new RuntimeException("No backend specified in backend definition block.");
        }

        AdaptiveParameters params;
        if (paramFile != null) {
          final String[] paramData = Input.ReadFile(paramFile);
          params = new org.ogolem.adaptive.AdaptiveParameters(paramData, -1);
        } else {
          params = globConf.parameters;
        }

        final BackendFactory backFactory = new BackendFactory(globConf, params);
        final GenericBackend<Molecule, Geometry> func = backFactory.parseBackend(backend);

        backendDefs.put(tag, func);
      }

      globConf.backendDefs = backendDefs;
    }

    /*
     * LAST THINGS: BUILD LOCOPT AND GLOBOPT OBJECTS
     */

    // now actually map the string to the local optimization method
    final LocOptFactory factory = new LocOptFactory(globConf, globConf.parameters, backendDefs);
    try {
      globConf.refNewton = factory.buildLocalOpt(locoptString);
    } catch (Exception e) {
      System.err.println("ERROR: Couldn't setup local optimization! " + e.toString());
      throw new InitIOException(e);
    }

    if (initialLocOpt != null) {
      try {
        for (final String s : globConf.getExample().makePrintableAbsoluteCoord(true)) {
          System.out.println(s);
        }
        globConf.initNewton = factory.buildLocalOpt(initialLocOpt);
      } catch (Exception e) {
        System.err.println("ERROR: Couldn't setup initial local optimization! " + e.toString());
        throw new InitIOException(e);
      }
    }

    final GenericSanityCheck<Molecule, Geometry> sanity =
        new GeometrySanityCheck(
            globConf.blowFacBondDetect,
            globConf.blowFacEnvClusterClashes,
            globConf.blowFacDissocDetect,
            globConf.doPreCD,
            globConf.doPreDD);
    final IndividualWriter<Geometry> writer = new GeometryWriter(globConf.geomDumpFolder);
    final double crossPoss = globConf.crossPossibility;
    final double mutPoss = globConf.mutatePossibility;
    final boolean printBeforeFitness = globConf.printBeforeLocOpt;
    final int noOfTries = globConf.howManyTries;
    final double molXOverProb = globConf.molCrossPoss;
    final double molMutProb = globConf.molMutatePoss;

    // map the string to the global optimization algorithm
    final GenericFitnessFunction<Molecule, Geometry> fitness =
        org.ogolem.core.FitnessFunctionFactory.build(
            globConf, globConf.refNewton, globConf.fitnessFunctionConfig);

    final GlobOptAlgoFactory globFactory =
        new GlobOptAlgoFactory(
            sanity,
            fitness,
            writer,
            crossPoss,
            mutPoss,
            printBeforeFitness,
            noOfTries,
            molXOverProb,
            molMutProb,
            globConf);
    globConf.opter = globFactory.translateToGlobOpt(globOptString);

    // map the string to the diversity checker
    globConf.diversityChecker = GlobalConfig.mapStringToDiversityCheck(diversityString);

    return globConf;
  }

  private static MoleculeConfig parseMoleculeBlock(
      final int id,
      final List<String> molBlock,
      final List<BondInfo> allBondInfos,
      final String inputDir)
      throws Exception {

    MoleculeConfig mc = new MoleculeConfig(true);
    mc.iID = id;
    int noAtoms = molBlock.size();
    if (molBlock.get(0).equalsIgnoreCase("flexy")) {

      mc.flexy = true;
      noAtoms--;

      // we expect a z-matrix
      ZMatrix zmat;

      if (molBlock.get(1).startsWith("MoleculePath=")) {
        // read the zmat from an auxiliary file in
        final String sPath = inputDir + File.separator + molBlock.get(1).substring(13).trim();

        if (sPath.endsWith(".zmat")) {
          zmat = readZmatFromFile(sPath);
        } else if (sPath.endsWith(".mopzmat")) {
          zmat = readMopZmatFromFile(sPath);
        } else {
          zmat = readZmatFromFile(sPath);
        }

        mc.noOfAtoms = zmat.getNoOfAtoms();

        if (molBlock.size() == 3 && molBlock.get(2).startsWith("BondMatrix=")) {
          // try to parse the bond matrix
          final String pathToBonds = molBlock.get(2).substring(11).trim();
          try {
            if (pathToBonds.endsWith(".list")) {
              final BondInfo bonds = readBondsListIn(pathToBonds, mc.noOfAtoms);
              allBondInfos.set(id, bonds);
            } else {
              final BondInfo bonds = readBondsIn(pathToBonds);
              allBondInfos.set(id, bonds);
            }
          } catch (Exception e) {
            System.err.println("ERROR: Couldn't read bonds in. " + e.toString());
            throw new InitIOException(e);
          }
        }

      } else {
        // first line is flexy
        mc.noOfAtoms = noAtoms;
        final String[] flexyCont = new String[noAtoms];
        for (int il = 0; il < mc.noOfAtoms; il++) {
          flexyCont[il] = molBlock.get(il + 1);
        }

        // if flexy is set, then we expect a z matrix as input.

        /*
         * The first number represents the atom. The second number is the bond length to the atom in the
         * third column. Then comes the bond angle and the atom to which the bond angle
         * is defined. Last the dihedral and the according atom.
         * For the first three atoms there is an exception. Some numbers are not populated.
         */
        zmat = new ZMatrix(mc.noOfAtoms);

        // THE FIRST LINE: Atom ID and nothing else
        zmat.setAZMatrixLine(0, flexyCont[0].trim(), 0.0, 0, 0.0, 0, 0.0, 0);

        // THE SECOND LINE: Atom ID, bond length, bond connect
        final String[] sa = flexyCont[1].split("\\;", 3);

        // trim the stuff
        for (int k = 0; k < 3; k++) {
          sa[k] = sa[k].trim();
        }
        try {
          final double bond = Double.parseDouble(sa[1]) * ANGTOBOHR;
          final int bondConn = Integer.parseInt(sa[2]) - 1;
          zmat.setAZMatrixLine(1, sa[0], bond, bondConn, 0.0, 0, 0.0, 0);
        } catch (Exception e) {
          throw new CastException("Failure to read in second z-Matrix line.", e);
        }

        // THE THIRD LINE: Atom ID, bond length, bond connect, bond angle, angle connect
        final String[] sa2 = flexyCont[2].split("\\;", 5);

        // trim the stuff
        for (int k = 0; k < 5; k++) {
          sa2[k] = sa2[k].trim();
        }
        try {
          final double bond = Double.parseDouble(sa2[1]) * ANGTOBOHR;
          final int bondConn = Integer.parseInt(sa2[2]) - 1;
          final double angle = Math.toRadians(Double.parseDouble(sa2[3]));
          final int angleConn = Integer.parseInt(sa2[4]) - 1;
          zmat.setAZMatrixLine(2, sa2[0], bond, bondConn, angle, angleConn, 0.0, 0);
        } catch (Exception e) {
          throw new CastException("Failure to read in third z-Matrix line.", e);
        }

        // FROM THE FOURTH LINE ON: EVERYTHING
        for (int k = 3; k < mc.noOfAtoms; k++) {
          final String[] sa3 = flexyCont[k].split("\\;", 7);

          // trim the stuff
          for (int kk = 0; kk < 7; kk++) {
            sa3[kk] = sa3[kk].trim();
          }
          try {
            final double bond = Double.parseDouble(sa3[1]) * ANGTOBOHR;
            final int bondConn = Integer.parseInt(sa3[2]) - 1;
            final double angle = Math.toRadians(Double.parseDouble(sa3[3]));
            final int angleConn = Integer.parseInt(sa3[4]) - 1;
            final double dihedral = Math.toRadians(Double.parseDouble(sa3[5]));
            final int dihedralConn = Integer.parseInt(sa3[6]) - 1;
            zmat.setAZMatrixLine(
                k, sa3[0], bond, bondConn, angle, angleConn, dihedral, dihedralConn);
          } catch (Exception e) {
            throw new CastException(e);
          }
        }
      }

      // creating the zmat is done, set it to the MoleculeConfig
      mc.zmat = zmat;

      // we create initial degrees of freedom matrixes here that we populate later
      final boolean[][] dofs = new boolean[mc.noOfAtoms][3];
      for (int t = 0; t < mc.noOfAtoms; t++) {
        for (int s = 0; s < 3; s++) {
          dofs[t][s] = false;
        }
      }
      mc.degreesOfFreedom = dofs;

      final CartesianCoordinates refCartes = CoordTranslation.zMatToCartesians(zmat);
      mc.atomNumbers = refCartes.getAllAtomNumbers();
      mc.atomTypes = refCartes.getAllAtomTypes();
      mc.refXYZ = refCartes.getAllXYZCoord();

    } else {

      // default is set to non-flexible
      if (molBlock.get(0).startsWith("MoleculePath=")) {
        // read the molecule from an auxiliary file in
        final String molPath = inputDir + File.separator + molBlock.get(0).substring(13).trim();

        CartesianCoordinates cartes;

        if (molPath.endsWith(".zmat")) {
          final ZMatrix zmat = readZmatFromFile(molPath);
          cartes = CoordTranslation.zMatToCartesians(zmat);
        } else if (molPath.endsWith(".regzmat")) {
          final ZMatrix zmat = readRegZmatFromFile(molPath);
          cartes = CoordTranslation.zMatToCartesians(zmat);
        } else if (molPath.endsWith(".xyz")) {
          cartes = readCartesFromFile(molPath);
        } else if (molPath.endsWith(".mopzmat")) {
          final ZMatrix zmat = readMopZmatFromFile(molPath);
          cartes = CoordTranslation.zMatToCartesians(zmat);
        } else {
          System.err.println(
              "WARNING: Unknown extension to file " + molPath + " assuming xyz format.");
          cartes = readCartesFromFile(molPath);
        }

        cartes.moveCoordsToCOM();
        // call the coordinate translation and get the MoleculeConfig
        mc =
            CoordTranslation.cartesianToBareMolConfig(
                cartes, id, mc.sID, false, null, mc.constricted, mc.constraints);

        if (molBlock.size() == 2 && molBlock.get(1).startsWith("BondMatrix=")) {
          // try to parse the bond matrix
          final String pathToBonds =
              inputDir + File.separator + molBlock.get(1).substring(11).trim();
          try {
            if (pathToBonds.endsWith(".list")) {
              final BondInfo bonds = readBondsListIn(pathToBonds, mc.noOfAtoms);
              allBondInfos.set(id, bonds);
            } else {
              final BondInfo bonds = readBondsIn(pathToBonds);
              allBondInfos.set(id, bonds);
            }
          } catch (Exception e) {
            System.err.println("ERROR: Couldn't read bonds in. " + e.toString());
            throw new InitIOException(e);
          }
        }
      } else {

        mc.noOfAtoms = noAtoms;

        final CartesianCoordinates cartes =
            new CartesianCoordinates(noAtoms, 1, new int[] {noAtoms});

        final String[] atomIDs = new String[mc.noOfAtoms];
        final double[][] xyz = cartes.getAllXYZCoord();
        for (int k = 0; k < noAtoms; k++) {
          // Fill the arrays needed for the cartesian
          final String[] saTemp = molBlock.get(k).split("\\;", 4);
          atomIDs[k] = saTemp[0].trim();
          try {
            xyz[0][k] = Double.parseDouble(saTemp[1].trim()) * ANGTOBOHR;
            xyz[1][k] = Double.parseDouble(saTemp[2].trim()) * ANGTOBOHR;
            xyz[2][k] = Double.parseDouble(saTemp[3].trim()) * ANGTOBOHR;
          } catch (Exception e) {
            System.err.println("Casting failed in reading xyz directly from input!");
            throw new CastException("Cast from String " + molBlock.get(k) + " unsuccessful.", e);
          }
        }

        // set all xyz coordinates
        cartes.setAllXYZ(xyz);
        // set the String ID
        cartes.setAllAtomTypes(atomIDs);
        cartes.moveCoordsToCOM();

        // call the coordinate translation and get the MoleculeConfig
        mc =
            CoordTranslation.cartesianToBareMolConfig(
                cartes, id, mc.sID, false, null, mc.constricted, mc.constraints);
      }
    }

    return mc;
  }

  /**
   * Reads the bonds in a matrix format in
   *
   * @param file the path to the bond matrix
   * @return the BondInfo, implementing SimpleBondInfo
   * @throws Exception in case something exceptional happens... ;-)
   */
  public static BondInfo readBondsIn(final String file) throws Exception {

    // read file
    final String[] data = ReadFile(file);
    final BondInfo bonds = new SimpleBondInfo(data.length);
    for (int i = 0; i < data.length; i++) {
      final String[] sa = data[i].trim().split("\\s+");
      for (int j = 0; j < data.length; j++) {
        bonds.setBond(i, j, Short.parseShort(sa[j]));
      }
    }

    return bonds;
  }

  /**
   * Reads the bonds in a list format in
   *
   * @param file the path to the list
   * @param noOfAtoms the number of atoms in the system
   * @return the BondInfo, implementing SimpleBondInfo
   * @throws Exception in case something exceptional happens... ;-)
   */
  public static BondInfo readBondsListIn(final String file, final int noOfAtoms) throws Exception {

    // read file
    final String[] data = ReadFile(file);

    // initialize with defaults
    final BondInfo bonds = new SimpleBondInfo(noOfAtoms);
    for (int i = 0; i < noOfAtoms; i++) {
      bonds.setBond(i, i, BondInfo.UNCERTAIN);
    }

    // parse data and manipulate bonds
    for (final String data1 : data) {
      if (data1.trim().isEmpty()) {
        continue;
      }
      final String[] sa = data1.trim().split("\\s+");
      final int i = Integer.parseInt(sa[0]);
      final int j = Integer.parseInt(sa[1]);
      final short b = Short.parseShort(sa[2]);
      // change the bonding matrix for those two
      bonds.setBond(i, j, b);
    }

    return bonds;
  }

  public static float[][] readChargesIn(
      final String[] chargeBlock,
      final int noOfMols,
      final int[] atsPerMol,
      final int[] molIDMapping)
      throws CastException {

    assert (chargeBlock != null); // may be 0 length
    assert (atsPerMol != null);
    assert (molIDMapping != null);
    assert (atsPerMol.length == molIDMapping.length);
    assert (noOfMols == molIDMapping.length);

    // initialize the charge arrays
    final float[][] faCharges = new float[noOfMols][];

    for (int i = 0; i < noOfMols; i++) {
      faCharges[i] = new float[atsPerMol[i]];
    }

    for (final String chargeBlock1 : chargeBlock) {
      final String[] saTemp = chargeBlock1.trim().split("\\;", 3);
      // get info which molecule, atom, charge is wanted
      int whichMol = 0;
      int whichAtom = 0;
      float charge = 0;
      try {
        whichMol = Integer.parseInt(saTemp[0].trim());
        whichAtom = Integer.parseInt(saTemp[1].trim());
        charge = Float.parseFloat(saTemp[2].trim());
      } catch (Exception e) {
        System.err.println("Failure in casting String to int in charges.");
        throw new CastException(e);
      }
      // set the charge in the appropriate place(s)
      boolean found = false;
      for (int x = 0; x < molIDMapping.length; x++) {
        if (whichMol == molIDMapping[x]) {
          faCharges[x][whichAtom] = charge;
          found = true;
        }
      }
      if (!found) {
        throw new CastException(
            "Charge parsing: Mol ID "
                + whichMol
                + " not existing. Range is 0 to "
                + molIDMapping[molIDMapping.length - 1]);
      }
    }

    return faCharges;
  }

  public static short[][] readSpinsIn(
      final String[] spinBlock, final int noOfMols, final int[] atsPerMol, final int[] molIDMapping)
      throws CastException {

    assert (spinBlock != null); // may be 0 length
    assert (atsPerMol != null);
    assert (molIDMapping != null);
    assert (atsPerMol.length == molIDMapping.length);
    assert (noOfMols == molIDMapping.length);

    // initialize the spin arrays
    final short[][] spins = new short[noOfMols][];

    for (int i = 0; i < noOfMols; i++) {
      spins[i] = new short[atsPerMol[i]];
    }

    for (final String spinBlock1 : spinBlock) {
      final String[] saTemp = spinBlock1.trim().split("\\;", 3);
      // get info which molecule, atom, spin is wanted
      int whichMol = 0;
      int whichAtom = 0;
      short spin = 0;
      try {
        whichMol = Integer.parseInt(saTemp[0].trim());
        whichAtom = Integer.parseInt(saTemp[1].trim());
        spin = Short.parseShort(saTemp[2].trim());
      } catch (Exception e) {
        System.err.println("Failure in casting String to int in spins.");
        throw new CastException(e);
      }
      // set the spin in the appropriate place(s)
      boolean found = false;
      for (int x = 0; x < molIDMapping.length; x++) {
        if (whichMol == molIDMapping[x]) {
          spins[x][whichAtom] = spin;
          found = true;
        }
      }
      if (!found) {
        throw new CastException(
            "Spin parsing: Mol ID "
                + whichMol
                + " not existing. Range is 0 to "
                + molIDMapping[molIDMapping.length - 1]);
      }
    }

    return spins;
  }

  static Geometry ReadSerializedGeometry(String sFilePath, String sFolderPath, boolean bDelete)
      throws SerialException, CastException, InitIOException {

    Object oObj = null;
    try {
      oObj =
          InputPrimitives.readBinInput(
              sFolderPath + System.getProperty("file.separator") + sFilePath);
    } catch (Exception e1) {
      System.err.println("Failed once to read geometry in! " + e1.toString());
      try {
        Thread.sleep(1000);
      } catch (InterruptedException enon) {
        // yeah, it is not necessarily too intelligent to have empty catch statements
        // but since we are already in a thrown exception, it should be fine IMHO
      }
      try {
        oObj =
            InputPrimitives.readBinInput(
                sFolderPath + System.getProperty("file.separator") + sFilePath);
      } catch (Exception e2) {
        System.err.println("Failed twice to read geometry in! " + e2.toString());
        try {
          Thread.sleep(1000);
        } catch (InterruptedException enon) {
          // yeah, it is not necessarily too intelligent to have empty catch statements
          // but since we are already in a thrown exception, it should be fine IMHO
        }
        try {
          oObj =
              InputPrimitives.readBinInput(
                  sFolderPath + System.getProperty("file.separator") + sFilePath);
        } catch (Exception e3) {
          System.err.println("Failed thrice to read geometry in! " + e3.toString());
          throw new SerialException(e3);
        }
      }
    }
    Geometry gGeo = null;
    try {
      gGeo = (Geometry) oObj;
    } catch (Exception e) {
      throw new CastException("Failure in cast from Object to Geometry!", e);
    }
    if (bDelete) {
      try {
        ManipulationPrimitives.remove(
            sFolderPath + System.getProperty("file.separator") + sFilePath);
      } catch (Exception e) {
        throw new InitIOException(e);
      }
    } else {
      // keep the old data
    }
    return gGeo;
  }

  @SuppressWarnings("unchecked")
  public static GenericPool<Molecule, Geometry> ReadSerializedPool(String sFolderPath, String sFile)
      throws SerialException, CastException {
    Object oObj = null;
    try {
      oObj =
          InputPrimitives.readBinInput(sFolderPath + System.getProperty("file.separator") + sFile);
    } catch (IOException | ClassNotFoundException e) {
      throw new SerialException(e);
    }
    GenericPool<Molecule, Geometry> pool;
    try {
      pool = (GenericPool<Molecule, Geometry>) oObj;
    } catch (Exception e) {
      throw new CastException("Failure in cast from object to pool!", e);
    }
    return pool;
  }

  public static GenericPool<Molecule, Geometry> ReadSerializedPool(String sFolderPath)
      throws SerialException, CastException {
    return ReadSerializedPool(sFolderPath, "pool.bin");
  }

  public static GlobalConfig ReadSerializedConfig(String sConfigPath)
      throws SerialException, CastException {
    Object oObj = null;
    try {
      oObj = InputPrimitives.readBinInput(sConfigPath);
    } catch (IOException | ClassNotFoundException e) {
      throw new SerialException(e);
    }
    GlobalConfig globconf = null;
    try {
      globconf = (GlobalConfig) oObj;
    } catch (Exception e) {
      throw new CastException("Failure in cast from Object to GlobalConfig!", e);
    }
    return globconf;
  }

  static double ReadEnergyMolproOutput(String sMolproOut) throws InitIOException, CastException {
    String[] saData;
    try {
      saData = InputPrimitives.readFileIn(sMolproOut);
    } catch (Exception e) {
      throw new InitIOException("Error in reading molpros output file.", e);
    }
    double dEnergy = 0.0;
    try {
      dEnergy = SearchMolproEnergy(saData);
    } catch (CastException e) {
      throw e;
    }
    return dEnergy;
  }

  static Gradient ReadGradientMolproOutput(String sMolproOut, int iNoOfAtoms)
      throws InitIOException, CastException {
    String[] saData;
    try {
      saData = InputPrimitives.readFileIn(sMolproOut);
    } catch (Exception e) {
      throw new InitIOException("Error in reading molpros output file.", e);
    }

    // find the beginning of the gradient
    int iGradientStart = 0;
    for (int i = 0; i < saData.length; i++) {
      if (saData[i].contains("GRADIENT FOR STATE 1.1")) {
        iGradientStart = i + 4;
        break;
      } else {
        // not found the start yet... next line
      }
    }

    // put the gradient data into a String array
    String[] saGradData = new String[iNoOfAtoms];
    for (int i = 0; i < iNoOfAtoms; i++) {
      saGradData[i] = saData[i + iGradientStart];
    }

    double[][] daGradient = new double[3][iNoOfAtoms];
    for (int i = 0; i < saGradData.length; i++) {
      String sTemp = saGradData[i].trim();
      // System.OUT.println("DEBUG: Dealing with line: " + tmp);
      // get rid of the atomic number
      String sTemp2 = sTemp.substring(sTemp.indexOf(" "));
      sTemp2 = sTemp2.trim();
      // first value: dE/dx
      int iEndOfNumber = sTemp2.indexOf(" ");
      sTemp = sTemp2.substring(0, iEndOfNumber);
      sTemp2 = sTemp2.substring(iEndOfNumber);
      // System.OUT.println("DEBUG: dE/dx " + tmp);
      try {
        final double dNumber = Double.parseDouble(sTemp);
        daGradient[0][i] = dNumber;
      } catch (Exception e) {
        throw new CastException(e);
      }
      // second value: dE/dy
      sTemp2 = sTemp2.trim();
      iEndOfNumber = sTemp2.indexOf(" ");
      sTemp = sTemp2.substring(0, iEndOfNumber);
      sTemp2 = sTemp2.substring(iEndOfNumber);
      // System.OUT.println("DEBUG: dE/dy " + tmp);
      try {
        final double dNumber = Double.parseDouble(sTemp);
        daGradient[1][i] = dNumber;
      } catch (Exception e) {
        throw new CastException(e);
      }
      // third value: dE/dz
      sTemp2 = sTemp2.trim();
      // System.OUT.println("DEBUG: dE/dz " + sTemp2);
      try {
        final double dNumber = Double.parseDouble(sTemp2);
        daGradient[2][i] = dNumber;
      } catch (Exception e) {
        throw new CastException(e);
      }
    }
    Gradient gradient = new Gradient();
    gradient.setGradientTotal(daGradient);

    double dEnergy = 0.0;
    try {
      dEnergy = SearchMolproEnergy(saData);
    } catch (CastException e) {
      throw e;
    }
    gradient.setTotalEnergy(dEnergy);
    return gradient;
  }

  private static double SearchMolproEnergy(String[] saData) throws CastException {
    int iEnergyLine = 0;
    for (int i = 0; i < saData.length; i++) {
      if (saData[i].contains("Variable memory released")) {
        iEnergyLine = i - 2;
        // the break is rather important since it can sometimes print the same line twice
        // (whyever...)
        break;
      }
    }
    String sTemp = saData[iEnergyLine].trim();

    double dEnergy = FixedValues.NONCONVERGEDENERGY;
    if (!sTemp.contains(" ")) {
      try {
        dEnergy = Double.parseDouble(sTemp);
      } catch (Exception e) {
        throw new CastException("Failure in casting molpros energy.", e);
      }
    } else {
      sTemp = sTemp.substring(0, sTemp.indexOf(" "));
      try {
        dEnergy = Double.parseDouble(sTemp);
      } catch (Exception e) {
        throw new CastException("Failure in casting molpros energy.", e);
      }
    }
    return dEnergy;
  }

  static CartesianCoordinates ReadXYZMolproOutput(
      String sMolproLog, int iNoOfAtoms, int iNoOfMolecules, int[] iaAtPerMol)
      throws InitIOException, CastException {

    CartesianCoordinates cartes = new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaAtPerMol);
    String[] saData;
    try {
      saData = InputPrimitives.readFileIn(sMolproLog);
    } catch (Exception e) {
      throw new InitIOException("Error in reading molpros log file.", e);
    }
    /*
     * get the "offset" for the optimized geometry
     */
    int iGeomOffset = saData.length - iNoOfAtoms - 3;

    /*
     * get the energy in a.u.
     */
    // System.OUT.println("DEBUG: The Energy line: " + saData[iGeomOffset]);
    String sEnergy = saData[iGeomOffset].substring(saData[iGeomOffset].indexOf("=") + 1);
    sEnergy = sEnergy.trim();
    // System.OUT.println("DEBUG: The Energy as String " + sEnergy);
    double dEnergy;
    try {
      dEnergy = Double.parseDouble(sEnergy);
      // System.OUT.println("DEBUG: Energy " + dEnergy);
    } catch (Exception e) {
      throw new CastException(e);
    }
    cartes.setEnergy(dEnergy);

    /*
     * get the coordinates  (molpro's output: angstrom)
     */
    String sTemp;
    String sTemp2;
    int iEndOfNumber;
    for (int i = 1 + iGeomOffset; i < iNoOfAtoms + 1 + iGeomOffset; i++) {
      sTemp = saData[i];
      // System.OUT.println("DEBUG: Dealing with: " + tmp);
      sTemp = sTemp.trim();
      // the atomic "name"
      sTemp2 = sTemp.substring(0, sTemp.indexOf(" "));
      cartes.setAtom(sTemp2, i - 1 - iGeomOffset);
      sTemp = sTemp.substring(sTemp.indexOf(" "));
      // System.OUT.println("DEBUG: Atom " + sTemp2);
      // the coordinates
      double[] daCoords = new double[3];

      // x
      sTemp = sTemp.trim();
      iEndOfNumber = sTemp.indexOf(" ");
      sTemp2 = sTemp.substring(0, iEndOfNumber);
      // System.OUT.println("DEBUG: Coordinate: " + sTemp2);
      try {
        daCoords[0] = Double.parseDouble(sTemp2) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException(e);
      }

      // y
      sTemp = sTemp.substring(sTemp.indexOf(" "));
      sTemp = sTemp.trim();
      iEndOfNumber = sTemp.indexOf(" ");
      sTemp2 = sTemp.substring(0, iEndOfNumber);
      // System.OUT.println("DEBUG: Coordinate: " + sTemp2);
      try {
        daCoords[1] = Double.parseDouble(sTemp2) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException(e);
      }

      // z
      sTemp = sTemp.substring(sTemp.indexOf(" "));
      sTemp = sTemp.trim();
      // System.OUT.println("DEBUG: Coordinate: " + tmp);
      try {
        daCoords[1] = Double.parseDouble(sTemp) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException(e);
      }

      cartes.setXYZCoordinatesOfAtom(daCoords, i - 1 - iGeomOffset);
    }

    return cartes;
  }

  static CartesianCoordinates ReadXYZDFTBPlusOutput(
      final String folder,
      final int noOfAtoms,
      final int noOfMolecules,
      final int[] noAtsPerMol,
      final short[] spins,
      final float[] charges)
      throws InitIOException, CastException {

    final String geomXYZFile = folder + System.getProperty("file.separator") + "geo_end.xyz";
    final CartesianCoordinates cartes =
        readCartesFromFile(geomXYZFile, noOfMolecules, noAtsPerMol, spins, charges);

    final String output = folder + System.getProperty("file.separator") + "detailed.out";
    try {
      final String[] data = InputPrimitives.readFileIn(output);
      for (int i = data.length - 1; i >= 0; i--) {
        final String line = data[i];
        if (line.contains("Total Electronic energy:")) {
          final String[] lineData = line.trim().split("\\s+");
          final double energy = Double.parseDouble(lineData[3]);
          cartes.setEnergy(energy);
        }
      }
    } catch (IOException | NumberFormatException e) {
      throw new CastException(e);
    }

    return cartes;
  }

  static CartesianCoordinates ReadXYZTinker(
      String sCartesianFile, int iNoOfAtoms, int iNoOfMolecules, int[] iaAtomsPerMol)
      throws IOException, CastException {

    final CartesianCoordinates cartes =
        new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaAtomsPerMol);
    final String[] saTinkerXYZOut = InputPrimitives.readFileIn(sCartesianFile);

    assert (saTinkerXYZOut != null);
    assert (saTinkerXYZOut.length >= iNoOfAtoms + 1);

    final double[][] xyz = cartes.getAllXYZCoord();

    // the chopping and casting games
    for (int i = 0; i < iNoOfAtoms; i++) {
      final String[] tmp = saTinkerXYZOut[i + 1].trim().split("\\s+");
      // the coordinates
      for (int j = 0; j < 3; j++) {
        try {
          xyz[j][i] = Double.parseDouble(tmp[j + 2]) * ANGTOBOHR;
        } catch (Exception e) {
          throw new CastException(e);
        }
      }
    }

    return cartes;
  }

  static void RemoveMolproFiles(String sMolproInput, boolean bLocOpt) throws InitIOException {
    String sBasis = sMolproInput.substring(0, sMolproInput.indexOf("."));
    try {
      ManipulationPrimitives.remove(sMolproInput);
    } catch (Exception e) {
      throw new InitIOException("Failure in deleting molpro input.", e);
    }
    try {
      ManipulationPrimitives.remove(sBasis + ".xml");
    } catch (Exception e) {
      throw new InitIOException("Failure in deleting molpro xml.", e);
    }
    try {
      ManipulationPrimitives.remove(sBasis + ".out");
    } catch (Exception e) {
      throw new InitIOException("Failure in deleting molpro output.", e);
    }
    if (bLocOpt) {
      try {
        ManipulationPrimitives.remove(sBasis + ".log");
      } catch (Exception e) {
        throw new InitIOException("Failure in deleting molpro log.", e);
      }
    }
  }

  public static void RemoveMopacFiles(String sMopacInput) throws InitIOException {
    String sBasis = sMopacInput.substring(0, sMopacInput.indexOf(".") + 1);
    String[] saListFiles;
    try {
      saListFiles = InquiryPrimitives.fileListWithPrefix(sBasis);
    } catch (Exception e) {
      throw new InitIOException(e);
    }

    for (final String saListFile : saListFiles) {
      try {
        ManipulationPrimitives.remove(saListFile);
      } catch (Exception e) {
        throw new InitIOException(e);
      }
    }
  }

  public static void RemoveOrcaFiles(String sOrcaBasis) throws InitIOException {
    String[] saListFiles;
    sOrcaBasis += ".";
    try {
      saListFiles = InquiryPrimitives.fileListWithPrefix(sOrcaBasis);
    } catch (Exception e) {
      throw new InitIOException(e);
    }

    for (final String saListFile : saListFiles) {
      try {
        ManipulationPrimitives.remove(saListFile);
      } catch (Exception e) {
        throw new InitIOException(e);
      }
    }
  }

  static void RemoveTinkerFiles(String sTinkerBasis) throws InitIOException {
    String[] saListFiles;
    sTinkerBasis += ".";
    try {
      saListFiles = InquiryPrimitives.fileListWithPrefix(sTinkerBasis);
    } catch (Exception e) {
      throw new InitIOException(e);
    }

    for (final String saListFile : saListFiles) {
      try {
        ManipulationPrimitives.remove(saListFile);
      } catch (Exception e) {
        throw new InitIOException(e);
      }
    }
  }

  public static String[] ReadFile(final String sFile) throws InitIOException {
    try {
      final String[] saFileContent = InputPrimitives.readFileIn(sFile);
      return saFileContent;
    } catch (Exception e) {
      throw new InitIOException("Couldn't read file " + sFile, e);
    }
  }

  /**
   * Returns a CartesianCoordinates object filed with everything but the energy.
   *
   * @param sFolder
   * @param saAtomTypes
   * @param iNoOfAtoms
   * @param iNoOfMolecules
   * @param iaAtomsPerMol
   * @return the cartesian coordinates
   * @throws ConvergenceException
   */
  static CartesianCoordinates ReadXYZNAMDOutput(
      String sFolder, String[] saAtomTypes, int iNoOfAtoms, int iNoOfMolecules, int[] iaAtomsPerMol)
      throws ConvergenceException {

    CartesianCoordinates cartes =
        new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaAtomsPerMol);
    cartes.setAllAtomTypes(saAtomTypes);

    /*
     * read the new coordinates in: opt.pdb
     */
    String sPDBFile = sFolder + System.getProperty("file.separator") + "opt.coor";
    double[][] daXYZ;
    try {
      daXYZ = ReadPDBToCartesians(sPDBFile, iNoOfAtoms);
    } catch (Exception e) {
      throw new ConvergenceException("Error in reading pdb file in.", e);
    }
    cartes.setAllXYZ(daXYZ);

    return cartes;
  }

  static double[][] ReadPDBToCartesians(String sWhichFile, int iNoOfAtoms) throws IOException {

    final double[][] daXYZ = new double[3][iNoOfAtoms];

    /*
     * read the complete pdb file in
     */
    final String[] saPDB = InputPrimitives.readFileIn(sWhichFile);

    for (int i = 0; i < iNoOfAtoms; i++) {
      final String[] tmp = saPDB[i + 1].substring(29).trim().split("\\s+");
      for (int j = 0; j < 3; j++) {
        // the actual cartesian coordinates
        daXYZ[j][i] = Double.parseDouble(tmp[j]) * ANGTOBOHR;
      }
    }
    return daXYZ;
  }

  public static CartesianCoordinates readCartesFromFile(
      final String sPath,
      final int noOfMolecules,
      final int[] atsPerMol,
      final short[] spins,
      final float[] charges)
      throws InitIOException, CastException {

    assert (spins.length == charges.length);
    assert (atsPerMol.length == noOfMolecules);

    // first try to read in the file
    final String[] fileCont = ReadFile(sPath);

    return parseCartesFromFileData(fileCont, noOfMolecules, atsPerMol, spins, charges);
  }

  public static CartesianCoordinates parseCartesFromFileData(
      final String[] fileCont,
      final int noOfMolecules,
      final int[] atsPerMol,
      final short[] spins,
      final float[] charges)
      throws InitIOException, CastException {

    assert (spins.length == charges.length);
    assert (atsPerMol.length == noOfMolecules);

    // the first line is our amounts of atoms
    int noOfAtoms = 0;
    try {
      noOfAtoms = Integer.parseInt(fileCont[0].trim());
    } catch (Exception e) {
      throw new CastException("Failure to parse the number of atoms.", e);
    }
    assert (spins.length == noOfAtoms);

    int calcNoAts = 0;
    for (int i = 0; i < atsPerMol.length; i++) {
      calcNoAts += atsPerMol[i];
    }
    if (calcNoAts != noOfAtoms) {
      throw new InitIOException(
          "ERROR: XYZ file does not contain the correct number of atoms! "
              + calcNoAts
              + " vs "
              + noOfAtoms);
    }

    /*
     * we consider (due to a lack of more information) the xyz to be a single
     * molecule
     */
    final CartesianCoordinates cartes =
        new CartesianCoordinates(noOfAtoms, noOfMolecules, atsPerMol);
    final double[][] xyz = cartes.getAllXYZCoord();

    // the second line is always a comment, so we start from the third
    for (int i = 0; i < noOfAtoms; i++) {

      final String[] tmp = fileCont[i + 2].trim().split("\\s+");
      cartes.setAtom(tmp[0], i);
      for (int coord = 0; coord < 3; coord++) {
        try {
          xyz[coord][i] = Double.parseDouble(tmp[coord + 1]) * ANGTOBOHR;
        } catch (Exception e) {
          throw new CastException("Failure to cast double to coordinate", e);
        }
      }
    }

    cartes.setAllCharges(charges);
    cartes.setAllSpins(spins);
    cartes.recalcAtomNumbers();

    return cartes;
  }

  /**
   * Returns a set of cartesian coordinates without spins or charges set and the information that it
   * is just one molecule spanning the whole cartesian!
   *
   * @param sPath
   * @return A cartesian coordinate object with the above mentioned restrictions.
   * @throws InitIOException
   * @throws CastException
   */
  public static CartesianCoordinates readCartesFromFile(String sPath)
      throws InitIOException, CastException {

    // first try to read in the file
    final String[] fileCont = ReadFile(sPath);

    // the first line is our amounts of atoms
    int noOfAtoms = 0;
    try {
      noOfAtoms = Integer.parseInt(fileCont[0].trim());
    } catch (Exception e) {
      throw new CastException("Failure to parse the number of atoms.", e);
    }

    /*
     * we consider (due to a lack of more information) the xyz to be a single
     * molecule
     */
    final int[] atsPerMol = {noOfAtoms};
    final CartesianCoordinates cartes = new CartesianCoordinates(noOfAtoms, 1, atsPerMol);

    final double[][] xyz = cartes.getAllXYZCoord();

    // the second line is always a comment, so we start from the third
    for (int i = 0; i < noOfAtoms; i++) {
      final String[] tmp = fileCont[i + 2].trim().split("\\s+");
      cartes.setAtom(tmp[0], i);
      for (int coord = 0; coord < 3; coord++) {
        try {
          xyz[coord][i] = Double.parseDouble(tmp[coord + 1]) * ANGTOBOHR;
        } catch (Exception e) {
          throw new CastException("Failure to cast double to coordinate", e);
        }
      }
    }
    cartes.recalcAtomNumbers();

    return cartes;
  }

  static CartesianCoordinates readCPMDCartesFromFile(
      final String sPath,
      final int iNoOfMolecules,
      final int[] iaAtsPerMol,
      final short[] iaSpins,
      final float[] faCharges)
      throws InitIOException, CastException {

    // first try to read in the file
    final String[] saFileCont = ReadFile(sPath);

    // the first line is our amounts of atoms
    int iNoOfAtoms = 0;
    try {
      iNoOfAtoms = Integer.parseInt(saFileCont[0].trim());
    } catch (Exception e) {
      throw new CastException("Failure to parse the number of atoms.", e);
    }

    /*
     * we consider (due to a lack of more information) the xyz to be a single
     * molecule
     */
    final CartesianCoordinates cartes =
        new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaAtsPerMol);

    // the second line is always a comment, so we start from the third
    final String[] atoms = cartes.getAllAtomTypes();
    final double[][] xyz = cartes.getAllXYZCoord();
    for (int i = 0; i < iNoOfAtoms; i++) {

      final String[] sa = saFileCont[i + 2].trim().split("\\s+");
      atoms[i] = sa[0];

      try {
        xyz[0][i] = Double.parseDouble(sa[1]) * ANGTOBOHR;
        xyz[1][i] = Double.parseDouble(sa[2]) * ANGTOBOHR;
        xyz[2][i] = Double.parseDouble(sa[3]) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException("Failure to cast double to coordinate", e);
      }
    }
    cartes.recalcAtomNumbers();

    return cartes;
  }

  public static ZMatrix readMopZmatFromFile(String path) throws Exception {
    // first try to read in the file
    final String[] fileCont = InputPrimitives.readFileIn(path);

    // the first line is our number of atoms
    int no = 0;
    try {
      no = Integer.parseInt(fileCont[0].trim());
    } catch (Exception e) {
      throw new CastException("ERROR: Failure to parse the number of atoms.", e);
    }

    final ZMatrix zmat = new ZMatrix(no);
    for (int i = 1; i < no + 1; i++) {
      final String[] sa = fileCont[i].trim().split("\\s+");
      final double bl = Double.parseDouble(sa[1]) * ANGTOBOHR;
      final double an = Math.toRadians(Double.parseDouble(sa[3]));
      final double di = Math.toRadians(Double.parseDouble(sa[5]));
      int bc = Integer.parseInt(sa[7]) - 1;
      if (bc < 0) {
        bc = 0;
      }
      int ac = Integer.parseInt(sa[8]) - 1;
      if (ac < 0) {
        ac = 0;
      }
      int dc = Integer.parseInt(sa[9]) - 1;
      if (dc < 0) {
        dc = 0;
      }
      zmat.setAZMatrixLine(i - 1, sa[0], bl, bc, an, ac, di, dc);
    }

    return zmat;
  }

  public static ZMatrix readZmatFromFile(final String path) throws Exception {

    // first try to read in the file
    final String[] fileCont = InputPrimitives.readFileIn(path);

    // the first line is our amounts of atoms
    final int noAtoms = Integer.parseInt(fileCont[0]);
    if (noAtoms < 2) {
      throw new Exception(
          "ERROR: OGOLEM assumes that you want to specify at least two atoms in the z-Matrix!");
    }

    final ZMatrix zmat = new ZMatrix(noAtoms);

    // THE SEC LINE: Atom ID and nothing else (in principle)
    final String[] line2 = fileCont[2].trim().split("\\s");
    zmat.setAZMatrixLine(0, line2[0], 0.0, 0, 0.0, 0, 0.0, 0);

    String sTempZMat;

    // THE THIRD LINE: Atom ID, bond length, bond connect
    final String[] line3 = fileCont[3].trim().split("\\s+");

    double dBondLength;
    try {
      dBondLength = Double.parseDouble(line3[1]) * ANGTOBOHR;
    } catch (Exception e) {
      throw new CastException(e);
    }
    int iBondConnect;
    try {
      iBondConnect = Integer.parseInt(line3[2]) - 1;
    } catch (Exception e) {
      throw new CastException(e);
    }
    zmat.setAZMatrixLine(1, line3[0], dBondLength, iBondConnect, 0.0, 0, 0.0, 0);

    if (noAtoms == 2) {
      return zmat;
    }

    // THE FOURTH LINE: Atom ID, bond length, bond connect, bond angle, angle connect
    sTempZMat = fileCont[4].trim();
    String[] saTemp = sTempZMat.split("\\s+");

    try {
      dBondLength = Double.parseDouble(saTemp[1]) * ANGTOBOHR;
    } catch (Exception e) {
      throw new CastException(e);
    }
    try {
      iBondConnect = Integer.parseInt(saTemp[2]) - 1;
    } catch (Exception e) {
      throw new CastException(e);
    }
    double dBondAngle;
    try {
      dBondAngle = Math.toRadians(Double.parseDouble(saTemp[3]));
    } catch (Exception e) {
      throw new CastException(e);
    }
    int iAngleConnect;
    try {
      iAngleConnect = Integer.parseInt(saTemp[4]) - 1;
    } catch (Exception e) {
      throw new CastException(e);
    }
    zmat.setAZMatrixLine(
        2, saTemp[0], dBondLength, iBondConnect, dBondAngle, iAngleConnect, 0.0, 0);

    // FROM THE FIFTH LINE ON: EVERYTHING
    double dDihedral;
    int iDihedralConnect;
    // three already done
    for (int k = 0; k < noAtoms - 3; k++) {
      sTempZMat = fileCont[k + 5].trim();

      saTemp = sTempZMat.split("\\s+");

      try {
        dBondLength = Double.parseDouble(saTemp[1]) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException(e);
      }
      try {
        iBondConnect = Integer.parseInt(saTemp[2]) - 1;
      } catch (Exception e) {
        throw new CastException(e);
      }
      try {
        dBondAngle = Math.toRadians(Double.parseDouble(saTemp[3]));
      } catch (Exception e) {
        throw new CastException(e);
      }
      try {
        iAngleConnect = Integer.parseInt(saTemp[4]) - 1;
      } catch (Exception e) {
        throw new CastException(e);
      }
      try {
        dDihedral = Math.toRadians(Double.parseDouble(saTemp[5]));
      } catch (Exception e) {
        throw new CastException(e);
      }
      try {
        iDihedralConnect = Integer.parseInt(saTemp[6]) - 1;
      } catch (Exception e) {
        throw new CastException(e);
      }
      zmat.setAZMatrixLine(
          k + 3,
          saTemp[0],
          dBondLength,
          iBondConnect,
          dBondAngle,
          iAngleConnect,
          dDihedral,
          iDihedralConnect);
    }

    return zmat;
  }

  public static ZMatrix readRegZmatFromFile(final String path) throws Exception {

    // first try to read in the file
    final String[] fileCont = InputPrimitives.readFileIn(path);

    // the first line is our amounts of atoms
    final int noAtoms = Integer.parseInt(fileCont[0]);
    if (noAtoms < 2) {
      throw new Exception(
          "ERROR: OGOLEM assumes that you want to specify at least two atoms in the z-Matrix!");
    }

    final ZMatrix zmat = new ZMatrix(noAtoms);

    // THE SEC LINE: Atom ID and nothing else (in principle)
    final String[] line2 = fileCont[2].trim().split("\\s");
    zmat.setAZMatrixLine(0, line2[0], 0.0, 0, 0.0, 0, 0.0, 0);

    String sTempZMat;

    // THE THIRD LINE: Atom ID, bond length, bond connect
    final String[] line3 = fileCont[3].trim().split("\\s+");

    try {
      final int iBondConnect = Integer.parseInt(line3[1]) - 1;
      final double dBondLength = Double.parseDouble(line3[2]) * ANGTOBOHR;
      zmat.setAZMatrixLine(1, line3[0], dBondLength, iBondConnect, 0.0, 0, 0.0, 0);
    } catch (Exception e) {
      throw new CastException(e);
    }

    if (noAtoms == 2) {
      return zmat;
    }

    // THE FOURTH LINE: Atom ID, bond length, bond connect, bond angle, angle connect
    sTempZMat = fileCont[4].trim();
    String[] saTemp = sTempZMat.split("\\s+");

    try {
      final int iBondConnect = Integer.parseInt(saTemp[1]) - 1;
      final double dBondLength = Double.parseDouble(saTemp[2]) * ANGTOBOHR;
      final int iAngleConnect = Integer.parseInt(saTemp[3]) - 1;
      final double dBondAngle = Math.toRadians(Double.parseDouble(saTemp[4]));
      zmat.setAZMatrixLine(
          2, saTemp[0], dBondLength, iBondConnect, dBondAngle, iAngleConnect, 0.0, 0);
    } catch (Exception e) {
      throw new CastException(e);
    }

    // FROM THE FIFTH LINE ON: EVERYTHING
    // three already done
    for (int k = 0; k < noAtoms - 3; k++) {
      sTempZMat = fileCont[k + 5].trim();

      saTemp = sTempZMat.split("\\s+");

      try {
        final int iBondConnect = Integer.parseInt(saTemp[1]) - 1;
        final double dBondLength = Double.parseDouble(saTemp[2]) * ANGTOBOHR;
        final int iAngleConnect = Integer.parseInt(saTemp[3]) - 1;
        final double dBondAngle = Math.toRadians(Double.parseDouble(saTemp[4]));
        final int iDihedralConnect = Integer.parseInt(saTemp[5]) - 1;
        final double dDihedral = Math.toRadians(Double.parseDouble(saTemp[5]));
        zmat.setAZMatrixLine(
            k + 3,
            saTemp[0],
            dBondLength,
            iBondConnect,
            dBondAngle,
            iAngleConnect,
            dDihedral,
            iDihedralConnect);
      } catch (Exception e) {
        throw new CastException(e);
      }
    }

    return zmat;
  }

  public static String[] listFilesInDirectory(final String sDir) {
    return InquiryPrimitives.fileList(sDir);
  }
}
