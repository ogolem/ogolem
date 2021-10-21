/*
Copyright (c) 2014, J. M. Dieterich
              2015-2021, J. M. Dieterich and B. Hartke
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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import org.ogolem.core.Output;
import org.ogolem.generic.Configuration;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericGlobalOptimization;
import org.ogolem.generic.GenericInitializer;
import org.ogolem.generic.IndividualReader;
import org.ogolem.generic.IndividualWriter;
import org.ogolem.generic.Optimizable;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.GenericPoolEntry;
import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.rmi.RMICodes.JOBSTATE;

/**
 * A generic global optimization job.
 *
 * @author Johannes Dieterich
 * @version 2021-10-15
 */
public class GenericGlobOptJob<E, T extends Optimizable<E>> implements Job<T> {

  private static final long serialVersionUID = (long) 20211015;
  private static final boolean DEBUG = false;

  private final String outFolder;
  private final String outFile;
  private final String poolTmpFile;

  private final Configuration<E, T> config;
  private final GenericPool<E, T> pool;
  private final GenericHistory<E, T> history;

  private final int poolSize;
  private final int globOptIter;

  private final List<String> seeds;

  private final IndividualWriter<T> writer;
  private final IndividualReader<T> reader;
  private final GenericFitnessFunction<E, T> fitFunc;
  private final GenericInitializer<E, T> initer;
  private final GenericGlobalOptimization<E, T> globOpter;

  private final boolean doNiching;
  private final NicheComputer<E, T> nicheComp;

  private final ReentrantReadWriteLock stateLock;
  private final ReentrantReadWriteLock poolLock;

  private int countInitialSubs = 0;
  private int countInitialReturns = 0;
  private int countGlobOptSubs = 0;
  private int countGlobOptReturns = 0;

  private boolean initialTasksSubmitted = false;
  private boolean initialTasksComplete = false;
  private boolean allTasksSubmitted = false;
  private boolean allTasksComplete = false;
  private boolean resultsWritten = false;

  public GenericGlobOptJob(
      final Configuration<E, T> config,
      final GenericPool<E, T> pool,
      final GenericHistory<E, T> history,
      final IndividualReader<T> reader,
      final IndividualWriter<T> writer,
      final String outputFolder,
      final String outputFile)
      throws Exception {

    assert (pool != null);
    assert (history != null);
    assert (config != null);
    this.config = config;
    this.pool = pool;
    this.history = history;
    this.poolSize = pool.getPoolSize();
    this.reader = reader;
    this.writer = writer;
    this.outFile = outputFile;
    this.outFolder = outputFolder;

    this.globOptIter = config.getNumberOfGlobalSteps();
    this.initer = config.getInitializer();
    this.fitFunc = config.getFitnessFunction();
    this.globOpter = config.getGlobalOptimization();

    this.nicheComp = config.getNicheComputer();
    this.doNiching = (nicheComp != null);

    // take care of seeding
    final String seedFolder = config.seedFolder();
    this.seeds = new ArrayList<>();
    if (seedFolder != null && !seedFolder.isEmpty()) {
      final File f = new File(seedFolder);
      if (!f.exists() && !f.isDirectory()) {
        throw new RuntimeException(
            "Seed folder " + seedFolder + " either not existing or not a directory.");
      }

      for (final String s : f.list()) {
        seeds.add(s);
      }
    }

    this.poolTmpFile = config.intermediatePoolFile();

    writeStartData();

    this.stateLock = new ReentrantReadWriteLock();
    this.poolLock = new ReentrantReadWriteLock();
  }

  @Override
  public Task<T> nextTask() {

    final ReentrantReadWriteLock.ReadLock rl = stateLock.readLock();
    rl.lock();

    // any initial fill tasks left?
    if (!initialTasksSubmitted && !initialTasksComplete) {

      rl.unlock();

      Task<T> task;
      if (!seeds.isEmpty()) {
        // local optimization of a seed
        task =
            new GenericSeedTask<>(
                reader,
                countInitialSubs,
                pool.getExample(),
                fitFunc,
                seeds.remove(seeds.size() - 1));
      } else {
        // real init
        final T work = pool.getExample();
        task = new GenericInitTask<>(initer, countInitialSubs, work, fitFunc);
      }

      stateLock.writeLock().lock();

      countInitialSubs++;
      if (countInitialSubs == poolSize) {
        initialTasksSubmitted = true;
      }
      stateLock.writeLock().unlock();
      return task;
    } else if (initialTasksSubmitted && !initialTasksComplete) {
      // we need to wait ATM
      rl.unlock();
      return null;
    }

    // any global optimization tasks left?
    if (!allTasksSubmitted && !allTasksComplete) {
      rl.unlock();

      final long id = countGlobOptSubs + poolSize;
      final List<T> parents = pool.getParents();

      final Task<T> task = new GenericGlobOptTask<>(globOpter, parents.get(0), parents.get(1), id);

      // add to the history something reasonable to catch the case when the globopt is unsuccessful
      history.addFamily(parents.get(0).getID(), parents.get(1).getID(), (int) id, false, false);

      stateLock.writeLock().lock();

      countGlobOptSubs++;
      if (countGlobOptSubs == globOptIter) {
        allTasksSubmitted = true;
      }
      stateLock.writeLock().unlock();
      return task;
    } else if (allTasksSubmitted && !allTasksComplete) {
      // we need to wait ATM
      rl.unlock();
      return null;
    } else {
      rl.unlock();
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
  public synchronized List<Task<T>> nextInitTasks(final int maxTasks, final int noProxies) {
    // the synchronized here makes sense, as that nextTask() is locking everything

    assert (maxTasks > 0);
    assert (noProxies > 0);

    if (DEBUG) {
      System.out.println("DEBUG: Asking for " + maxTasks + " tasks with " + noProxies);
    }

    // this function shall only be called for the initial tasks!
    if (initialTasksComplete) {
      if (DEBUG) {
        System.out.println("DEBUG: Done with initial tasks.");
      }
      return new ArrayList<>();
    }

    final int noTasks =
        Math.min(maxTasks, (int) Math.ceil((double) this.poolSize / (double) noProxies));

    if (DEBUG) {
      System.out.println(
          "DEBUG: Handing out "
              + noTasks
              + " to client who asked for "
              + maxTasks
              + " we have "
              + noProxies
              + " proxies and a pool size of "
              + poolSize);
    }

    final List<Task<T>> jobs = new ArrayList<>();
    for (int i = 0; i < noTasks; i++) {
      final Task<T> t = nextTask();
      if (t == null) {
        break;
      } // tasks exhausted for the time being
      jobs.add(t);
    }

    if (DEBUG) {
      System.out.println("DEBUG: Returning " + jobs.size() + " initial tasks to do.");
    }

    return jobs;
  }

  @Override
  public Tuple3D<Boolean, Long, Integer> nextGlobOptChunk(final int maxTasks) {

    if (DEBUG) {
      System.out.println(
          "DEBUG: Being asked for new chunk. Status: "
              + initialTasksComplete
              + " maxTasks "
              + maxTasks);
      System.out.println(
          "DEBUG: globOptIter "
              + globOptIter
              + " countGlobOptSubs "
              + countGlobOptSubs
              + " poolSize "
              + poolSize);
    }

    final ReentrantReadWriteLock.ReadLock rl = stateLock.readLock();
    int tasks = -42;
    long futureStart = -42;
    rl.lock();
    try {
      if (!initialTasksComplete) {
        return new Tuple3D<>(false, 0l, 0);
      }
      tasks = Math.min(maxTasks, globOptIter - countGlobOptSubs);
      futureStart = countGlobOptSubs + poolSize;

      if (tasks <= 0) {
        if (DEBUG) {
          System.out.println("DEBUG: Globpt tasks less than 0, returning no tasks.");
        }

        return new Tuple3D<>(false, futureStart, 0);
      }
    } finally {
      rl.unlock();
    }

    stateLock.writeLock().lock();
    try {
      countGlobOptSubs += tasks;
      if (DEBUG) {
        System.out.println("DEBUG: Returning " + tasks + " globopt tasks to do.");
      }

      return new Tuple3D<>(true, futureStart, tasks);
    } finally {
      stateLock.writeLock().unlock();
    }
  }

  @Override
  public RMICodes.JOBSTATE submitResult(final Result<T> result) {

    assert (result != null);

    final ReentrantReadWriteLock.ReadLock rl = stateLock.readLock();
    rl.lock();

    if (!initialTasksComplete) {
      rl.unlock();

      if (result.wasOK()) {
        if (doNiching) {
          final Niche niche = nicheComp.computeNiche(result.getResult());
          poolLock.writeLock().lock();
          pool.addIndividualForcedUnsync(
              result.getResult(), niche, result.getResult().getFitness());
          poolLock.writeLock().unlock();
        } else {
          poolLock.writeLock().lock();
          pool.addIndividualForcedUnsync(result.getResult(), null, result.getResult().getFitness());
          poolLock.writeLock().unlock();
        }
      }

      stateLock.writeLock().lock();

      countInitialReturns++;
      if (countInitialReturns == poolSize) {
        initialTasksComplete = true;
        writeInterData();
      }
      stateLock.writeLock().unlock();
    } else {

      rl.unlock();

      if (result.wasOK()) {

        if (DEBUG)
          System.out.println(
              "DEBUG: Result about to be put in from client " + result.getClientID());

        final T res = result.getResult();

        boolean accepted;
        boolean wasNull;
        if (res == null) {
          System.err.println("ERROR: Returned result was null. Contact the author(s) please.");
        } else {
          if (doNiching) {
            final Niche niche = nicheComp.computeNiche(res);
            poolLock.writeLock().lock();
            accepted = pool.addIndividualUnsync(res, niche, res.getFitness());
            poolLock.writeLock().unlock();
          } else {
            poolLock.writeLock().lock();
            accepted = pool.addIndividualUnsync(res, null, res.getFitness());
            poolLock.writeLock().unlock();
          }
          wasNull = false;
          history.addFamily(
              res.getFatherID(), res.getMotherID(), (int) res.getID(), accepted, wasNull, res);
          if (DEBUG)
            System.out.println(
                "DEBUG: Result "
                    + res.getID()
                    + " from "
                    + result.getClientID()
                    + " was added to pool and history. Mother ID: "
                    + res.getMotherID()
                    + ", father ID: "
                    + res.getFatherID());
        }
      } else {
        final long[] la = result.getFamily();
        assert (la != null);
        assert (history != null);
        history.addFamily(la[0], la[1], (int) la[2], false, true);
        if (DEBUG) System.out.println("DEBUG: Result " + la[2] + " not OK but added to history.");
      }

      stateLock.writeLock().lock();

      countGlobOptReturns++;
      if (countGlobOptReturns == globOptIter) {
        allTasksComplete = true;
        writeFinalData();
        this.resultsWritten = true;
      }

      stateLock.writeLock().unlock();
    }

    final boolean waiting = jobWaiting();
    final boolean finished = jobFinished();

    if (waiting) {
      if (DEBUG) {
        System.out.println("DEBUG: Telling to wait.");
      }
      return JOBSTATE.WAITING;
    } else if (finished) {
      if (DEBUG) {
        System.out.println("DEBUG: Telling to finish.");
      }
      return JOBSTATE.FINISH;
    } else if (!waiting && !finished) {
      if (DEBUG) {
        System.out.println("DEBUG: Telling to continue.");
      }
      return JOBSTATE.CONTINUE;
    } else {
      // emergency
      System.err.println("ERROR: Unknown state in the OGOLEM RMI server.");
      return JOBSTATE.SNAFU;
    }
  }

  @Override
  public List<T> mergePools(
      final int noOfAssocResults,
      final List<T> clientPool,
      final int maxStructsBack,
      final long lastStart) {

    if (DEBUG) {
      System.out.println(
          "DEBUG: Attempting to merge pools " + noOfAssocResults + "  " + maxStructsBack);
    }

    final List<T> mergedPool = new ArrayList<>(maxStructsBack);

    final ReentrantReadWriteLock.ReadLock rl = stateLock.readLock();
    rl.lock();
    if (!initialTasksComplete) {
      rl.unlock();

      if (clientPool.size() != noOfAssocResults) {
        // should raise a flag
        System.err.println(
            "WARNING: Client got "
                + noOfAssocResults
                + " inits but is only returning "
                + clientPool.size()
                + " results.");
      }

      poolLock.writeLock().lock();

      clientPool.forEach(
          (t) -> {
            final long id = t.getID();
            // should be even easier: if ID is below our previous offset -> discard
            if (!(id < lastStart)) {
              if (doNiching) {
                final Niche niche = nicheComp.computeNiche(t);
                pool.addIndividualForcedUnsync(t, niche, t.getFitness());
              } else {
                pool.addIndividualForcedUnsync(t, null, t.getFitness());
              }
            }
          });

      int c = 0;
      for (final GenericPoolEntry<E, T> entry : pool) {
        c++;
        if (entry != null && entry.individual() != null) {
          mergedPool.add(entry.individual());
        } else {
          break; // end of pool
        }
        if (c >= maxStructsBack) {
          break;
        }
      }

      poolLock.writeLock().unlock();
      stateLock.writeLock().lock();

      countInitialReturns += noOfAssocResults;
      if (countInitialReturns == poolSize) {
        initialTasksComplete = true;
        writeInterData();
        if (DEBUG) {
          System.out.println("DEBUG: Initial tasks all complete.");
        }
      }

      stateLock.writeLock().unlock();
    } else {
      rl.unlock();
      poolLock.writeLock().lock();

      clientPool.forEach(
          (t) -> {
            final long id = t.getID();
            // should be even easier: if ID is below our previous offset -> discard
            if (!(id < lastStart)) {
              if (doNiching) {
                final Niche niche = nicheComp.computeNiche(t);
                pool.addIndividualUnsync(t, niche, t.getFitness());
              } else {
                pool.addIndividualUnsync(t, null, t.getFitness());
              }
            }
          });

      int c = 0;
      for (final GenericPoolEntry<E, T> entry : pool) {
        c++;
        if (entry != null && entry.individual() != null) {
          mergedPool.add(entry.individual());
        } else {
          break; // end of pool
        }
        if (c >= maxStructsBack) {
          break;
        }
      }

      poolLock.writeLock().unlock();
      stateLock.writeLock().lock();

      countGlobOptReturns += noOfAssocResults;

      if (DEBUG) {
        System.out.println(
            "DEBUG: Associated results: "
                + noOfAssocResults
                + " returned globopts: "
                + countGlobOptReturns
                + " vs globoptiters: "
                + globOptIter);
      }

      if (countGlobOptReturns >= globOptIter) {
        allTasksComplete = true;
        writeFinalData();
        this.resultsWritten = true;
        if (DEBUG) {
          System.out.println("DEBUG: Globopt tasks all complete.");
        }
      }

      stateLock.writeLock().unlock();
    }

    if (DEBUG) {
      System.out.println(
          "DEBUG: Done merging pools "
              + noOfAssocResults
              + "  "
              + maxStructsBack
              + "  "
              + mergedPool.size());
    }

    return mergedPool;
  }

  @Override
  public boolean jobWaiting() {

    stateLock.readLock().lock();
    final boolean waiting =
        (initialTasksSubmitted && !initialTasksComplete)
            || (allTasksSubmitted && !allTasksComplete);
    if (DEBUG) {
      System.out.println(
          "DEBUG: Job waiting? "
              + initialTasksSubmitted
              + " "
              + initialTasksComplete
              + " "
              + allTasksSubmitted
              + " "
              + allTasksComplete);
      System.out.println(
          "DEBUG: raw watiing: "
              + countInitialSubs
              + " "
              + countInitialReturns
              + " "
              + countGlobOptSubs
              + " "
              + countGlobOptReturns);
    }
    stateLock.readLock().unlock();

    return waiting;
  }

  @Override
  public boolean jobFinished() {

    boolean accFitReached = false;
    poolLock.readLock().lock();
    try {
      accFitReached =
          pool.acceptableFitnessReached(); // this can be outside the synchronized block of states!
    } finally {
      poolLock.readLock().unlock();
    }

    boolean finished = false;
    boolean resWritten = true;
    stateLock.readLock().lock();
    try {
      finished = (initialTasksComplete && allTasksComplete) || accFitReached;
      if (DEBUG) {
        System.out.println(
            "DEBUG: Job finished? "
                + initialTasksComplete
                + " "
                + allTasksComplete
                + " "
                + accFitReached);
        System.out.println(
            "DEBUG: raw watiing: "
                + countInitialSubs
                + " "
                + countInitialReturns
                + " "
                + countGlobOptSubs
                + " "
                + countGlobOptReturns);
      }
      resWritten = this.resultsWritten;
    } finally {
      stateLock.readLock().unlock();
    }

    if (finished && !resWritten) {
      stateLock.writeLock().lock();
      try {
        writeFinalData();
        this.resultsWritten = true;
      } finally {
        stateLock.writeLock().unlock();
      }
    }

    return finished;
  }

  private void writeStartData() {

    final String[] header = Output.getHeader();
    final String[] generific =
        new String[] {
          "####################################################################",
          "####################################################################",
          "                  AND OUR PROBLEM OF THE DAY IS:",
          "                  " + pool.getExample().getClass().getCanonicalName(),
          "####################################################################",
          "####################################################################",
        };
    final List<String> confDat = config.getFormattedConfig();
    try {
      OutputPrimitives.writeOut(outFile, header, true);
      OutputPrimitives.writeOut(outFile, generific, true);
      OutputPrimitives.writeOut(outFile, confDat, true);
    } catch (IOException e) {
      System.err.println("Failed to write header or config to output file. Ignoring...");
      e.printStackTrace(System.err);
    }
  }

  private void writeInterData() {

    final List<String> initFitnesses = new ArrayList<>();
    initFitnesses.add("#-----------------------------------------------------------");
    initFitnesses.add("#");
    initFitnesses.add("# Individuals in initial pool coming, fitness units unknown.");
    initFitnesses.add("#");
    initFitnesses.add("# pool position        individual id                 fitness");

    final List<String> poolCont = pool.getFormattedPool();
    poolCont.forEach(
        (s) -> {
          initFitnesses.add(s);
        });

    initFitnesses.add("#");
    initFitnesses.add("#-----------------------------------------------------------");

    poolLock.readLock().lock();
    try {
      OutputPrimitives.writeObjToBinFile(poolTmpFile, pool);
      OutputPrimitives.writeOut(outFile, initFitnesses, true);

      // print the individuals out as well
      int c = 0;
      for (final GenericPoolEntry<E, T> entry : pool) {
        final T ind = entry.individual();
        final String file =
            outFolder + File.separator + "initrank" + c + "individual" + ind.getID();
        writer.writeIndividual(ind, file);
        c++;
      }
    } catch (Exception e) {
      System.err.println("Failure in writing initial pool fitnesses out. Ignoring.");
      e.printStackTrace(System.err);
    }

    poolLock.readLock().unlock();
  }

  private void writeFinalData() {

    history.flushRecords();
    history.writeTotalStats();

    final List<String> finalFitnesses = new ArrayList<>();
    finalFitnesses.add("#-----------------------------------------------------------");
    finalFitnesses.add("#");
    finalFitnesses.add("# Individuals in final pool coming, fitness units unknown.");
    finalFitnesses.add("#");
    finalFitnesses.add("# pool position        individual id                 fitness");

    final List<String> poolCont2 = pool.getFormattedPool();
    poolCont2.forEach(
        (s) -> {
          finalFitnesses.add(s);
        });

    finalFitnesses.add("#");
    finalFitnesses.add("#-----------------------------------------------------------");
    poolLock.readLock().lock();
    try {
      OutputPrimitives.writeObjToBinFile(outFolder + File.separator + "finalPool.bin", pool);
      OutputPrimitives.writeOut(outFile, finalFitnesses, true);

      // print the individuals out as well
      int c = 0;
      for (final GenericPoolEntry<E, T> entry : pool) {
        final T ind = entry.individual();
        final String file = outFolder + File.separator + "rank" + c + "individual" + ind.getID();
        writer.writeIndividual(ind, file);
        c++;
      }
    } catch (Exception e) {
      System.err.println("Failure in writing initial pool fitnesses out. Ignoring.");
      e.printStackTrace(System.err);
    }
    poolLock.readLock().unlock();
  }
}
