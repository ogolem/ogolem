/*
Copyright (c) 2010-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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
package org.ogolem.rmi;

import static org.ogolem.rmi.RMICodes.JOBSTATE.*;

import java.util.ArrayList;
import java.util.List;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.GenericPoolEntry;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.switches.Color;
import org.ogolem.switches.ColorPalette;
import org.ogolem.switches.Switch;
import org.ogolem.switches.SwitchesConfig;
import org.ogolem.switches.Taboos;

/**
 * A job for the global optimization of a switchable molecule.
 *
 * @author Johannes Dieterich
 * @version 2020-12-26
 */
final class SwitchGlobOptJob implements Job<Switch> {

  private static final long serialVersionUID = (long) 20101114;

  private final Switch refSwitch;
  private final SwitchesConfig config;
  private final ColorPalette palette;
  private final GenericPool<Color, Switch> pool;
  private final Taboos taboos;

  private final int poolSize;
  private final int globOptIter;

  private int countInitialSubs = 0;
  private int countInitialReturns = 0;
  private int countGlobOptSubs = 0;
  private int countGlobOptReturns = 0;

  private boolean initialTasksSubmitted = false;
  private boolean initialTasksComplete = false;
  private boolean allTasksSubmitted = false;
  private boolean allTasksComplete = false;

  SwitchGlobOptJob(final SwitchesConfig configuration, final Switch example) {
    this.refSwitch = example;
    this.config = configuration;
    this.palette = ColorPalette.getReference();
    this.pool = GenericPool.getInstance(configuration.getGenericConfig(), example);
    this.taboos = Taboos.getReference();
    this.poolSize = SwitchesConfig.iPoolSize;
    this.globOptIter = SwitchesConfig.iNoOfGlobIters;
  }

  @Override
  public Task<Switch> nextTask() {

    // any initial fill tasks left?
    if (!initialTasksSubmitted && !initialTasksComplete) {

      final Task<Switch> task = new SwitchInitTask(refSwitch, countInitialSubs, config, palette);

      countInitialSubs++;
      if (countInitialSubs == poolSize) {
        initialTasksSubmitted = true;
      }
      return task;
    } else if (initialTasksSubmitted && !initialTasksComplete) {
      // we need to wait ATM
      return null;
    }

    // any global optimization tasks left?
    if (!allTasksSubmitted && !allTasksComplete) {

      final int id = countGlobOptSubs + poolSize;
      final List<Switch> parents = pool.getParents();
      final Task<Switch> task = new SwitchGlobOptTask(id, parents, config, palette, taboos);

      countGlobOptSubs++;
      if (countGlobOptSubs == globOptIter) {
        allTasksSubmitted = true;
      }
      return task;
    } else if (allTasksSubmitted && !allTasksComplete) {
      // we need to wait ATM
      return null;
    } else {
      System.err.println(
          "ERROR: You shouldn't end up here in the RMI task getting. Contact author(s).");
      System.err.println(
          "       Please provide them with this: initTasksSub "
              + initialTasksSubmitted
              + " initTasksCompl "
              + initialTasksComplete
              + " allTasksSub "
              + allTasksSubmitted
              + " allTasksDone "
              + allTasksComplete);
      return null;
    }
  }

  @Override
  public RMICodes.JOBSTATE submitResult(Result<Switch> result) {
    if (!initialTasksComplete) {
      if (result.wasOK()) {
        pool.addIndividualForced(result.getResult(), result.getResult().getFitness());
        taboos.addTaboo(result.getResult());
      }

      countInitialReturns++;
      if (countInitialReturns == poolSize) {
        initialTasksComplete = true;
        writeInterData();
      }
    } else {

      if (result.wasOK()) {
        pool.addIndividual(result.getResult(), result.getResult().getFitness());
        taboos.addTaboo(result.getResult());
      }

      countGlobOptReturns++;
      if (countGlobOptReturns == globOptIter) {
        allTasksComplete = true;
        writeFinalData();
      }
    }

    if (jobWaiting()) {
      return WAITING;
    } else if (jobFinished()) {
      return FINISH;
    } else if (!jobWaiting() && !jobFinished()) {
      return CONTINUE;
    } else {
      // emergency
      System.err.println("ERROR: Unknown state in the OGOLEM RMI server.");
      return SNAFU;
    }
  }

  @Override
  public boolean jobWaiting() {
    return (initialTasksSubmitted && !initialTasksComplete)
        || (!allTasksSubmitted && !allTasksComplete);
  }

  @Override
  public boolean jobFinished() {
    return (initialTasksComplete && allTasksComplete);
  }

  private void writeInterData() {

    for (int i = 0; i < SwitchesConfig.iPoolSize; i++) {
      final Switch sw = pool.getIndividualAtPosition(i);

      final String[][] saCisTrans = sw.createPrintableCisTrans(config.dBlowBondsFac);
      final String[] saColorFile = sw.createPrintableColors();

      final String sFileBase = "initrank" + i + "switch" + sw.getID();
      final String sColorFile = sFileBase + "-colors.col";
      final String sCisFile = sFileBase + "-cis.xyz";
      final String sTransFile = sFileBase + "-trans.xyz";

      // actually write it out
      try {
        OutputPrimitives.writeOut(sColorFile, saColorFile, false);
        OutputPrimitives.writeOut(sCisFile, saCisTrans[0], false);
        OutputPrimitives.writeOut(sTransFile, saCisTrans[1], false);
      } catch (Exception e) {
        System.err.println(
            "ERROR: Failed to write results for initrank "
                + i
                + " continuing anyway."
                + e.toString());
      }
    }
  }

  private void writeFinalData() {

    for (int i = 0; i < SwitchesConfig.iPoolSize; i++) {
      final Switch sw = pool.getIndividualAtPosition(i);

      final String[][] saCisTrans = sw.createPrintableCisTrans(config.dBlowBondsFac);
      final String[] saColorFile = sw.createPrintableColors();

      final String sFileBase = "rank" + i + "switch" + sw.getID();
      final String sColorFile = sFileBase + "-colors.col";
      final String sCisFile = sFileBase + "-cis.xyz";
      final String sTransFile = sFileBase + "-trans.xyz";

      // actually write it out
      try {
        OutputPrimitives.writeOut(sColorFile, saColorFile, false);
        OutputPrimitives.writeOut(sCisFile, saCisTrans[0], false);
        OutputPrimitives.writeOut(sTransFile, saCisTrans[1], false);
      } catch (Exception e) {
        System.err.println(
            "ERROR: Failed to write results for initrank "
                + i
                + " continuing anyway."
                + e.toString());
      }
    }
  }

  @Override
  public synchronized List<Task<Switch>> nextInitTasks(final int maxTasks, final int noProxies) {

    // this function shall only be called for the initial tasks!
    if (initialTasksComplete) {
      return new ArrayList<>();
    }

    final int noTasks =
        Math.min(maxTasks, (int) Math.ceil((double) this.poolSize / (double) noProxies));

    final List<Task<Switch>> jobs = new ArrayList<>();
    for (int i = 0; i < noTasks; i++) {
      final Task<Switch> t = nextTask();
      if (t == null) {
        break;
      } // tasks exhausted for the time being
      jobs.add(t);
    }

    return jobs;
  }

  @Override
  public synchronized Tuple3D<Boolean, Long, Integer> nextGlobOptChunk(final int maxTasks) {

    if (!initialTasksComplete) {
      return new Tuple3D<>(false, 0l, 0);
    }

    final int tasks = Math.min(maxTasks, globOptIter - countGlobOptSubs);
    final long futureStart = countGlobOptSubs + poolSize;
    if (tasks <= 0) {
      return new Tuple3D<>(false, futureStart, 0);
    } else {
      countGlobOptSubs += tasks;
      return new Tuple3D<>(true, futureStart, tasks);
    }
  }

  @Override
  public synchronized List<Switch> mergePools(
      int noOfAssocResults,
      List<Switch> clientPool,
      final int maxStructsBack,
      final long lastStart) {

    if (!initialTasksComplete) {

      clientPool.forEach(
          (t) -> {
            pool.addIndividualForced(t, t.getFitness());
          });
      countInitialReturns += noOfAssocResults;
      if (countInitialReturns == poolSize) {
        initialTasksComplete = true;
        writeInterData();
      }
    } else {

      clientPool.forEach(
          (t) -> {
            final long id = t.getID();
            if (!(id < lastStart)) {
              pool.addIndividual(t, t.getFitness());
            }
          });
      countGlobOptReturns += noOfAssocResults;
      if (countGlobOptReturns == globOptIter) {
        allTasksComplete = true;
        writeFinalData();
      }
    }

    final List<Switch> mergedPool = new ArrayList<>(maxStructsBack);
    int c = 0;
    for (final GenericPoolEntry<Color, Switch> entry : pool.getAllIndividuals()) {
      c++;
      if (entry != null && entry.getIndividual() != null) {
        mergedPool.add(entry.getIndividual());
      } else {
        break; // end of pool
      }
      if (c >= maxStructsBack) {
        break;
      }
    }

    return mergedPool;
  }
}
