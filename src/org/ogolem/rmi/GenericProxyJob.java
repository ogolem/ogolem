/*
Copyright (c) 2014, J. M. Dieterich
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
package org.ogolem.rmi;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import org.ogolem.generic.Configuration;
import org.ogolem.generic.GenericGlobalOptimization;
import org.ogolem.generic.Optimizable;
import org.ogolem.generic.generichistory.GenericHistory;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.GenericPoolEntry;
import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.rmi.RMICodes.JOBSTATE;

/**
 * A generic global optimization job.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class GenericProxyJob<E, T extends Optimizable<E>> implements Job<T> {

  private static final long serialVersionUID = (long) 20150804;
  private static final boolean DEBUG = false;

  private final RMICommunication<T> serverComm;
  private final int myID;
  private List<Task<T>> initTasks = null;
  private final int maxTasks;
  private final boolean doMaxStructsExchange;
  private final int maxStructsExchange;

  private final GenericPool<E, T> pool;
  private final GenericHistory<E, T> history;
  private final GenericGlobalOptimization<E, T> globOpter;

  private final AliveCheck aliveCheck;

  private final boolean doNiching;
  private final NicheComputer<E, T> nicheComp;

  private final ReentrantReadWriteLock stateLock;

  private boolean globOptComplete = false;

  private int initTasksGotten = 0;
  private int initTasksHandedOut = 0;
  private int initTasksReturned = 0;
  private boolean initialTasksSubmitted = false;
  private boolean initialTasksComplete = false;

  private long currIDStart = 0;
  private int tasksGotten = 0;
  private int tasksHandedOut = 0;
  private int tasksReturned = 0;
  private boolean allTasksSubmitted = false;

  private RMICodes.JOBSTATE lastServerComm = JOBSTATE.CONTINUE;

  public GenericProxyJob(
      final Configuration<E, T> config,
      final GenericPool<E, T> pool,
      final GenericHistory<E, T> history,
      final RMICommunication<T> comm,
      final int maxTasks,
      final int maxStructsExchange,
      final AliveCheck alive,
      final int myID)
      throws Exception {

    assert (pool != null);
    assert (history != null);
    assert (comm != null);
    assert (maxTasks > 0);
    this.pool = pool;
    this.history = history;
    this.serverComm = comm;
    this.globOpter = config.getGlobalOptimization();

    this.nicheComp = config.getNicheComputer();
    this.doNiching = (nicheComp != null);

    this.maxTasks = maxTasks;
    this.maxStructsExchange = maxStructsExchange;
    this.doMaxStructsExchange =
        !(maxStructsExchange < 0 || maxStructsExchange >= pool.getPoolSize());
    this.aliveCheck = alive;
    this.myID = myID;
    this.stateLock = new ReentrantReadWriteLock();
  }

  @Override
  public Task<T> nextTask() {

    if (DEBUG) {
      System.out.println(
          "DEBUG: State of proxy: init tasks submitted? "
              + initialTasksSubmitted
              + " init tasks complete? "
              + initialTasksComplete
              + " curr ID offset "
              + currIDStart
              + " #tasks gotten "
              + tasksGotten
              + " #tasks handed out "
              + tasksHandedOut
              + " #tasks returned "
              + tasksReturned
              + " #tasks submitted "
              + allTasksSubmitted);
    }

    // initialize initTasks, if necessary
    final ReentrantReadWriteLock.ReadLock rl = stateLock.readLock();
    rl.lock();

    if (initTasks == null) {

      rl.unlock();

      stateLock.writeLock().lock();
      try {
        initTasks = serverComm.getInitialTasks(maxTasks, myID);
        if (DEBUG) {
          System.out.println("DEBUG: Initial chunk request: " + initTasks.size());
        }
        initTasksGotten = initTasks.size();
      } catch (Exception e) {
        System.err.println("ERROR: Exception in initial proxy to server communication.");
        e.printStackTrace(System.err);
        System.exit(28);
      } finally {
        stateLock.writeLock().unlock();
      }

      rl.lock(); // to be locked before the next check
    }

    // any initial fill tasks left?
    if (!initialTasksSubmitted && !initialTasksComplete) {

      rl.unlock();

      stateLock.writeLock().lock();
      initTasksHandedOut++;
      if (initTasksHandedOut == initTasksGotten) {
        initialTasksSubmitted = true;
      }
      final Task<T> task = initTasks.remove(0); // remove topmost, slightly slower but cleaner
      stateLock.writeLock().unlock();

      return task;
    } else if (initialTasksSubmitted && !initialTasksComplete) {
      // we need to wait ATM

      // yes, this is slow. however, to avoid doing something to lastServerComm in the write region
      // which is not picked up, we need to do this here. bummer.
      rl.unlock();
      stateLock.writeLock().lock();

      // here, we need to account for the situation that we were told to wait, ergo DID NOT get a
      // chunk of init tasks
      // and now have nothing to hand over, which means the client has nothing to return, so no new
      // chunk request gets triggered
      if (lastServerComm == JOBSTATE.WAITING) {

        // ask for a new chunk
        try {
          final List<Task<T>> moreTasks = serverComm.getInitialTasks(maxTasks, myID);
          if (DEBUG) {
            System.out.println("DEBUG: Gotten " + moreTasks.size() + " more initial tasks.");
          }
          if (!moreTasks.isEmpty()) {
            // more tasks
            this.initTasksGotten += moreTasks.size();
            for (int i = 0; i < moreTasks.size(); i++) {
              this.initTasks.add(moreTasks.remove(i));
            }
          }
          lastServerComm = JOBSTATE.CONTINUE;
        } catch (Exception e) {
          System.err.println("Exception when trying to get a ");
          e.printStackTrace(System.err);
        }

        stateLock.writeLock().unlock();

        return null; // we tell it to wait either way. the client will ask us again and can then get
        // a proper task.
      }

      stateLock.writeLock().unlock();
      if (DEBUG) {
        System.out.println("DEBUG: initial tasks submitted, but not complete. Waiting.");
      }

      return null;
    }

    // any global optimization tasks left?
    if (!allTasksSubmitted && !globOptComplete) {

      rl.unlock();
      stateLock.writeLock().lock();

      final long id = currIDStart + tasksHandedOut;

      final List<T> parents = pool.getParents();

      if (DEBUG) {
        System.out.println(
            "DEBUG: ID "
                + id
                + " mother "
                + parents.get(0).getID()
                + " father "
                + parents.get(1).getID());
      }
      final Task<T> task = new GenericGlobOptTask<>(globOpter, parents.get(0), parents.get(1), id);

      // add to the history something reasonable to catch the case when the globopt is unsuccessful
      history.addFamily(parents.get(0).getID(), parents.get(1).getID(), (int) id, false, false);

      tasksHandedOut++;
      if (tasksHandedOut >= tasksGotten) {
        allTasksSubmitted = true;
      }
      stateLock.writeLock().unlock();

      return task;
    } else if (allTasksSubmitted && !globOptComplete) {

      // yes, this is slow. however, to avoid doing something to lastServerComm in the write region
      // which is not picked up, we need to do this here. bummer.
      rl.unlock();
      stateLock.writeLock().lock();

      // here, we need to account for the situation that we were told to wait, ergo DID NOT get a
      // chunk of globopt tasks
      // and now have nothing to hand over, which means the client has nothing to return, so no new
      // chunk request gets triggered
      if (lastServerComm == JOBSTATE.WAITING) {

        // ask for a new chunk
        try {
          final Tuple3D<RMICodes.JOBSTATE, Long, Integer> tup =
              serverComm.getGlobOptChunk(maxTasks, myID);
          if (DEBUG) {
            System.out.println(
                "DEBUG: Got a chunk of global optimization tasks "
                    + tup.getObject1()
                    + " "
                    + tup.getObject2()
                    + " "
                    + tup.getObject3());
          }
          if (null != tup.getObject1())
            switch (tup.getObject1()) {
              case CONTINUE:
                // all fine :-)
                this.currIDStart = tup.getObject2();
                this.tasksGotten = tup.getObject3();
                this.tasksHandedOut = 0;
                this.tasksReturned = 0;
                this.allTasksSubmitted = false;
                if (DEBUG) {
                  System.out.println(
                      "DEBUG: Chunk request: all fine. " + currIDStart + " " + tasksGotten);
                }
                lastServerComm = JOBSTATE.CONTINUE;
                break;
              case FINISH:
                // all exhausted
                if (DEBUG) {
                  System.out.println("DEBUG: Chunk request: empty.");
                }
                globOptComplete = true;
                lastServerComm = JOBSTATE.FINISH;
                break;
              case WAITING:
                if (DEBUG) {
                  System.out.println("DEBUG: Chunk request: wait.");
                }
                // do wait
                lastServerComm = JOBSTATE.WAITING;
                break;
              default:
                // unknown status for server-proxy communication
                lastServerComm = JOBSTATE.SNAFU;
                throw new Exception(
                    "Unknown status "
                        + tup.getObject1()
                        + " "
                        + tup.getObject2()
                        + "  "
                        + tup.getObject3()
                        + " for server-proxy communication.");
            }
        } catch (Exception e) {
          System.err.println("Exception when trying to get a ");
          e.printStackTrace(System.err);
        }

        stateLock.writeLock().unlock();

        return null; // we tell it to wait either way. the client will ask us again and can then get
        // a proper task.
      }

      // we need to wait ATM
      stateLock.writeLock().unlock();
      if (DEBUG) {
        System.out.println("DEBUG: globopt tasks submitted, but not complete. Waiting.");
      }

      return null;
    } else {
      rl.unlock();
      System.err.println(
          "ERROR: You shouldn't end up here in the RMI task getting. Contact author(s).");
      System.err.println(
          "       Please provide them with this state of proxy: init tasks submitted? "
              + initialTasksSubmitted
              + " init tasks complete? "
              + initialTasksComplete
              + " curr ID offset "
              + currIDStart
              + " #tasks gotten "
              + tasksGotten
              + " #tasks handed out "
              + tasksHandedOut
              + " #tasks returned "
              + tasksReturned
              + " #tasks submitted "
              + allTasksSubmitted);
      return null;
    }
  }

  @Override
  public List<Task<T>> nextInitTasks(final int maxTasks, final int noProxies) {
    throw new Error("Proxy should not be called from Proxy! (I)");
  }

  @Override
  public Tuple3D<Boolean, Long, Integer> nextGlobOptChunk(final int noProxies) {
    throw new Error("Proxy should not be called from Proxy! (III)");
  }

  @Override
  public RMICodes.JOBSTATE submitResult(final Result<T> result) {

    final ReentrantReadWriteLock.ReadLock rl = stateLock.readLock();
    rl.lock();

    if (!initialTasksComplete) {
      rl.unlock();

      if (result.wasOK()) {
        if (doNiching) {
          final Niche niche = nicheComp.computeNiche(result.getResult());
          pool.addIndividualForced(result.getResult(), niche, result.getResult().getFitness());
        } else {
          pool.addIndividualForced(result.getResult(), result.getResult().getFitness());
        }
      }

      stateLock.writeLock().lock();
      initTasksReturned++;
      if (initTasksReturned == initTasksGotten) {

        try {
          // see if there are more, only merge once there are not more to do
          final List<Task<T>> moreTasks = serverComm.getInitialTasks(maxTasks, myID);
          if (DEBUG) {
            System.out.println("DEBUG: Gotten " + moreTasks.size() + " more initial tasks.");
          }
          if (!moreTasks.isEmpty()) {
            // more tasks
            this.initTasksGotten += moreTasks.size();
            for (int i = 0; i < moreTasks.size(); i++) {
              this.initTasks.add(moreTasks.remove(i));
            }

            stateLock.writeLock().unlock(); // SUPERIMPORTANT
            lastServerComm = JOBSTATE.CONTINUE;
            return JOBSTATE.CONTINUE;
          }
        } catch (Exception e) {
          System.err.println("ERROR: Couldn't ask for more initial tasks.");
          e.printStackTrace(System.err);
          System.exit(997);
        }

        initialTasksComplete = true;
        try {
          if (doMaxStructsExchange) {
            final List<T> myPool = new ArrayList<>();
            int c = 0;
            int offset = 0;
            // always synchronize our top, NEW individuals with the server
            while (c < maxStructsExchange
                && pool.getCurrentPoolSize() > 0
                && offset < pool.getCurrentPoolSize()) {
              final T ind =
                  pool.getIndividualAtPosition(
                      offset); // we remove from the top (ignoring known IDs)
              final long id = ind.getID();
              if (id < currIDStart) {
                offset++; // increment the offset
                continue;
              }
              myPool.add(ind);
              pool.removeIndividualAtPos(offset);
              c++;
            }
            final List<T> merged =
                serverComm.synchronizePool(myPool, myID, maxStructsExchange, currIDStart);
            if (merged.size() > maxStructsExchange) {
              throw new Exception("ERROR: Getting more individuals back than I thought!");
            }
            merged.forEach(
                (t) -> {
                  // augment pool... we however will need to check if the IDs are the same!
                  // (i.e., the main pool may give something back that actually originated from us
                  // but then, we also want to get "older" IDs that originated from other proxies
                  if (t.getID() < currIDStart) {
                    // as this is "expensive", only do this for the IDs that are older
                    if (doNiching) {
                      final Niche niche = nicheComp.computeNiche(t);
                      pool.addIndividualForcedUnsyncCheckID(t, niche, t.getFitness(), 5);
                    } else {
                      pool.addIndividualForcedUnsyncCheckID(t, null, t.getFitness(), 5);
                    }
                  } else {
                    if (doNiching) {
                      final Niche niche = nicheComp.computeNiche(t);
                      pool.addIndividualForced(t, niche, t.getFitness());
                    } else {
                      pool.addIndividualForced(t, t.getFitness());
                    }
                  }
                });
            if (DEBUG) {
              System.out.println("DEBUG: Initial pool merging with max struct exchange done.");
            }
          } else {
            final List<T> myPool = new ArrayList<>();
            pool.getAllIndividuals()
                .forEach(
                    (entry) -> {
                      myPool.add(entry.individual());
                    });
            final List<T> merged =
                serverComm.synchronizePool(myPool, myID, pool.getPoolSize(), currIDStart);
            if (merged.size() > pool.getPoolSize()) {
              throw new Exception("WARNING: Getting more individuals back than I thought!");
            }
            pool.emptyPool();
            merged.forEach(
                (t) -> {
                  // replace pool...
                  if (doNiching) {
                    final Niche niche = nicheComp.computeNiche(t);
                    pool.addIndividualForced(t, niche, t.getFitness());
                  } else {
                    pool.addIndividualForced(t, t.getFitness());
                  }
                });
            if (DEBUG) {
              System.out.println("DEBUG: Initial full pool merging done.");
            }
          }
          // and get a chunk of global optimization!
          final Tuple3D<RMICodes.JOBSTATE, Long, Integer> tup =
              serverComm.getGlobOptChunk(maxTasks, myID);
          if (DEBUG) {
            System.out.println(
                "DEBUG: Got a chunk of global optimization tasks "
                    + tup.getObject1()
                    + " "
                    + tup.getObject2()
                    + " "
                    + tup.getObject3());
          }
          if (null != tup.getObject1())
            switch (tup.getObject1()) {
              case CONTINUE:
                // all fine :-)
                this.currIDStart = tup.getObject2();
                this.tasksGotten = tup.getObject3();
                this.tasksHandedOut = 0;
                this.tasksReturned = 0;
                this.allTasksSubmitted = false;
                if (DEBUG) {
                  System.out.println(
                      "DEBUG: Chunk request: all fine. " + currIDStart + " " + tasksGotten);
                }
                lastServerComm = JOBSTATE.CONTINUE;
                break;
              case FINISH:
                // all exhausted
                if (DEBUG) {
                  System.out.println("DEBUG: Chunk request: empty.");
                }
                globOptComplete = true;
                lastServerComm = JOBSTATE.FINISH;
                break;
              case WAITING:
                if (DEBUG) {
                  System.out.println("DEBUG: Chunk request: wait.");
                }
                // do wait
                lastServerComm = JOBSTATE.WAITING;
                return JOBSTATE.WAITING;
              default:
                // unknown status for server-proxy communication
                lastServerComm = JOBSTATE.SNAFU;
                throw new Exception(
                    "Unknown status "
                        + tup.getObject1()
                        + " "
                        + tup.getObject2()
                        + "  "
                        + tup.getObject3()
                        + " for server-proxy communication.");
            }
        } catch (Exception e) {
          System.err.println("ERROR: Couldn't merge initial pool.");
          e.printStackTrace(System.err);
          System.exit(993);
        }
      }
      stateLock.writeLock().unlock();
    } else {
      rl.unlock();

      if (result.wasOK()) {

        if (DEBUG)
          System.out.println(
              "DEBUG: Result about to be put in from client "
                  + result.getClientID()
                  + " with ID "
                  + result.getResult().getID());

        final T res = result.getResult();

        boolean accepted;
        boolean wasNull;
        if (res == null) {
          System.err.println("ERROR: Returned result was null. Contact the author(s) please.");
        } else {
          if (doNiching) {
            final Niche niche = nicheComp.computeNiche(res);
            accepted = pool.addIndividual(res, niche, res.getFitness());
          } else {
            accepted = pool.addIndividual(res, res.getFitness());
          }
          wasNull = false;
          history.addFamily(
              res.getMotherID(), res.getFatherID(), (int) res.getID(), accepted, wasNull, res);
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
        history.addFamily(la[0], la[1], (int) la[2], false, true);
        if (DEBUG) System.out.println("DEBUG: Result " + la[2] + " not OK but added to history.");
      }

      stateLock.writeLock().lock();
      tasksReturned++;
      if (tasksReturned == tasksGotten) {
        // that was that chunk. merge and try to get a new chunk
        try {
          if (DEBUG) {
            System.out.println("DEBUG: Trying to merge and get a new global optimization chunk.");
          }
          if (doMaxStructsExchange) {
            final List<T> myPool = new ArrayList<>();
            int c = 0;
            // always synchronize our top individuals with the server
            for (final GenericPoolEntry<E, T> entry : pool.getAllIndividuals()) {
              myPool.add(entry.individual());
              pool.removeIndividualAtPos(c);
              c++;
              if (c >= maxStructsExchange) {
                break;
              }
            }
            final List<T> merged =
                serverComm.synchronizePool(myPool, myID, maxStructsExchange, currIDStart);
            if (merged.size() > maxStructsExchange) {
              throw new Exception("Getting more individuals back than I thought!");
            }
            merged.forEach(
                (t) -> {
                  // augment pool...
                  if (doNiching) {
                    final Niche niche = nicheComp.computeNiche(t);
                    pool.addIndividualForced(t, niche, t.getFitness());
                  } else {
                    pool.addIndividualForced(t, t.getFitness());
                  }
                });
            if (DEBUG) {
              System.out.println("DEBUG: Pool merging with max struct exchange done.");
            }
          } else {
            final List<T> myPool = new ArrayList<>();
            pool.getAllIndividuals()
                .forEach(
                    (entry) -> {
                      myPool.add(entry.individual());
                    });
            final List<T> merged =
                serverComm.synchronizePool(myPool, myID, pool.getPoolSize(), currIDStart);
            pool.emptyPool();
            merged.forEach(
                (t) -> {
                  // replace pool...
                  if (doNiching) {
                    final Niche niche = nicheComp.computeNiche(t);
                    pool.addIndividualForced(t, niche, t.getFitness());
                  } else {
                    pool.addIndividualForced(t, t.getFitness());
                  }
                });
            if (DEBUG) {
              System.out.println("DEBUG: Full pool merging done.");
            }
          }

          if (DEBUG) {
            System.out.println("DEBUG: Merging done. Getting new chunk of glob opt tasks.");
          }

          final Tuple3D<RMICodes.JOBSTATE, Long, Integer> tup =
              serverComm.getGlobOptChunk(maxTasks, myID);
          if (DEBUG) {
            System.out.println(
                "DEBUG: Got another chunk of global optimization tasks "
                    + tup.getObject1()
                    + " "
                    + tup.getObject2()
                    + " "
                    + tup.getObject3());
          }
          if (null != tup.getObject1())
            switch (tup.getObject1()) {
              case CONTINUE:
                // all fine :-)
                this.currIDStart = tup.getObject2();
                this.tasksGotten = tup.getObject3();
                this.tasksHandedOut = 0;
                this.tasksReturned = 0;
                this.allTasksSubmitted = false;
                if (DEBUG) {
                  System.out.println(
                      "DEBUG: Chunk request (II): all fine. " + currIDStart + " " + tasksGotten);
                }
                lastServerComm = JOBSTATE.CONTINUE;
                break;
              case FINISH:
                // all exhausted
                if (DEBUG) {
                  System.out.println("DEBUG: Chunk request (II): empty.");
                }
                globOptComplete = true;
                lastServerComm = JOBSTATE.FINISH;
                break;
              case WAITING:
                // do wait
                if (DEBUG) {
                  System.out.println("DEBUG: Chunk request (II): wait.");
                }
                lastServerComm = JOBSTATE.WAITING;
                return JOBSTATE.WAITING;
              default:
                // unknown status for server-proxy communication
                lastServerComm = JOBSTATE.SNAFU;
                throw new Exception(
                    "Unknown status "
                        + tup.getObject1()
                        + " "
                        + tup.getObject2()
                        + "  "
                        + tup.getObject3()
                        + " for server-proxy communication.");
            }

        } catch (Exception e) {
          System.err.println("ERROR: Couldn't merge intermediate pool or get new chunk.");
          e.printStackTrace(System.err);
          System.exit(994);
        }
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
      System.err.println("ERROR: Unknown state in the OGOLEM RMI proxy.");
      return JOBSTATE.SNAFU;
    }
  }

  @Override
  public List<T> mergePools(
      final int noOfAssocResults,
      final List<T> clientPool,
      final int maxStructsBack,
      final long lastStart) {

    throw new Error("Proxy should not be called from Proxy! (II)");
  }

  @Override
  public boolean jobWaiting() {

    stateLock.readLock().lock();
    try {
      final boolean waiting =
          (initialTasksSubmitted && !initialTasksComplete)
              || (tasksHandedOut == tasksGotten) && (tasksReturned < tasksGotten);
      if (DEBUG) {
        System.out.println(
            "DEBUG: Job waiting? "
                + initialTasksSubmitted
                + " "
                + initialTasksComplete
                + " "
                + tasksHandedOut
                + " "
                + tasksGotten
                + " "
                + tasksReturned);
      }

      return waiting;
    } finally {
      stateLock.readLock().unlock();
    }
  }

  @Override
  public boolean jobFinished() {

    stateLock.readLock().lock();
    try {
      final boolean isFinished = (initialTasksComplete && globOptComplete);
      if (isFinished && !aliveCheck.wasMarkedDone()) {
        try {
          if (DEBUG) {
            System.out.println("DEBUG: Telling live check to shut down.");
          }
          aliveCheck.done();
        } catch (Exception e) {
          System.err.println(
              "WARNING: Couldn't properly shut the "
                  + "AliveCheck down. Exiting now. "
                  + e.toString()
                  + " Proxy ID is "
                  + myID
                  + ".");
          e.printStackTrace(System.err);
        }
      }
      if (DEBUG) {
        System.out.println("DEBUG: Job finished? " + initialTasksComplete + " " + globOptComplete);
      }

      return isFinished;
    } finally {
      stateLock.readLock().unlock();
    }
  }
}
