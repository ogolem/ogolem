/*
Copyright (c) 2014, J. M. Dieterich
              2020-2021, J. M. Dieterich and B. Hartke
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
package org.ogolem.generic.mpi;

import java.util.Random;
import org.ogolem.generic.Optimizable;
import org.ogolem.helpers.Tuple;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.random.Lottery;
import org.ogolem.random.RNGenerator;
import org.ogolem.random.StandardRNG;
import org.ogolem.rmi.Job;
import org.ogolem.rmi.Result;
import org.ogolem.rmi.Task;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A generic class for MPI based parallelization of our optimization problems. Please note: you must
 * call MPI_INIT from OUTSIDE this and decide already, if the process is queen (rank 0) or drone
 * (rank 1 - N)! Note that the queen does NOT do any quantum of work herself here.
 *
 * @author Johannes Dieterich
 * @version 2021-07-28
 */
public class GenericMPIOptimization<E, T extends Optimizable<E>> {

  private static final Logger LOG = LoggerFactory.getLogger(GenericMPIOptimization.class);
  private static final int NEXTTASK = 0;
  private static final int EXITDONE = 1;
  private static final int WAITFOR = 2;
  private static final int DUMMY = 3;
  private static final int RESULT = 4;

  private GenericMPIOptimization() {}

  @SuppressWarnings("unchecked")
  public static <T> void runAsQueen(final MPIInterface mpi, final Job<T> job) throws Exception {

    final boolean debug = LOG.isDebugEnabled();
    if (debug) LOG.debug("Entering generic MPI globopt as queen. BRACE YOURSELF BIG TIMES!");

    final int noProcs = mpi.mpiCommSize();
    if (noProcs <= 1) {
      throw new RuntimeException("Trying to actually work all by myself is not in my nature. Bye.");
    }

    /*
     * initial broadcast
     */
    final char[] initMessage = "Hello from MPI queen, all is well!".toCharArray();
    final Tuple<Integer, char[]> bcastRet = mpi.mpiBcast(initMessage, 0, 0);
    if (bcastRet.getObject1() != 0)
      throw new RuntimeException("Error in sending initial bcast. " + bcastRet);

    /*
     * fill up all drones once
     */
    int taskCounter = 0;
    for (int proc = 1; proc < noProcs; proc++) {
      final Task<T> task = job.nextTask();
      if (task == null) {
        // this case may happen if pool size < #drones
        if (debug) LOG.debug("Queen: Sending WAITFOR to " + proc);
        final int sendRes = mpi.mpiSend(new byte[1], 0, proc, WAITFOR);
        if (sendRes != 0)
          throw new RuntimeException("Failure to send wait to rank " + proc + "\t" + sendRes);
      } else {
        final byte[] data = OutputPrimitives.writeObjToByteArray(task);
        final int retSend = mpi.mpiSend(data, 0, proc, 0);
        if (retSend != 0)
          throw new RuntimeException(
              "Error in sending initial task "
                  + taskCounter
                  + " to "
                  + proc
                  + " with error "
                  + retSend);

        taskCounter++;
        if (debug) LOG.debug("Queen: Send task out to " + proc + " counter is " + taskCounter);
      }
    }

    while (!job.jobFinished()) {

      // receive something
      if (debug) LOG.debug("Queen: Waiting to receive, counter is " + taskCounter);
      final MPIInterface.MPIStatus rcvStat =
          mpi.mpiRecvBytes(0, MPIInterface.ANY_SOURCE, MPIInterface.ANY_TAG);
      if (rcvStat.err != 0) {
        throw new RuntimeException(
            "Error in receiving w/ counter "
                + taskCounter
                + " tag "
                + rcvStat.tag
                + " from "
                + rcvStat.source
                + " error "
                + rcvStat.err);
      }

      if (debug) LOG.debug("Queen: Received, tag " + rcvStat.tag + " from " + rcvStat.source);
      if (rcvStat.tag != DUMMY) {
        final Result<T> result = (Result<T>) InputPrimitives.readByteInput(rcvStat.msg);
        taskCounter--;

        // submit it
        job.submitResult(result);
      }

      if (job.jobFinished()) { // no need to get a new task out
        if (debug) LOG.debug("Queen: Sending WAITFOR (DONE) to " + rcvStat.source);
        final int sendRes = mpi.mpiSend(new byte[1], 0, rcvStat.source, WAITFOR);
        if (sendRes != 0)
          throw new RuntimeException(
              "Failure to send wait to rank " + rcvStat.source + "\t" + sendRes);
        break;
      }

      // and give something out
      final Task<T> task = job.nextTask();
      if (task == null) {
        if (debug) LOG.debug("Queen: Sending WAITFOR to " + rcvStat.source);
        final int sendRes = mpi.mpiSend(new byte[1], 0, rcvStat.source, WAITFOR);
        if (sendRes != 0)
          throw new RuntimeException(
              "Failure to send wait to rank " + rcvStat.source + "\t" + sendRes);
      } else {
        if (debug) LOG.debug("Queen: Sending task to " + rcvStat.source);
        final byte[] data = OutputPrimitives.writeObjToByteArray(task);
        final int sendRes = mpi.mpiSend(data, 0, rcvStat.source, NEXTTASK);
        if (sendRes != 0)
          throw new RuntimeException(
              "Failure to send next task to rank " + rcvStat.source + "\t" + sendRes);

        taskCounter++;
      }
    }

    /*
     * receive any stragglers
     */
    if (taskCounter > 0) {
      LOG.info("Queen: Job finished. Waiting for " + taskCounter + " drones to report.");
    }

    while (taskCounter > 0) {
      if (debug) LOG.debug("Queen: Waiting for stragglers " + taskCounter);
      // receive outstanding
      final MPIInterface.MPIStatus rcvStat =
          mpi.mpiRecvBytes(0, MPIInterface.ANY_SOURCE, MPIInterface.ANY_TAG);
      taskCounter--;
    }

    LOG.info("Queen: Finalizing across all drones.");

    /*
     * tell everybody to go kill themselves
     */
    for (int proc = 1; proc < noProcs; proc++) {
      if (debug) LOG.debug("Queen: Sending EXITDONE to " + proc);
      final int sendRes = mpi.mpiSend(new byte[1], 0, proc, EXITDONE);
      if (sendRes != 0)
        throw new RuntimeException("Failure to send exit done to rank " + proc + "\t" + sendRes);
    }
  }

  @SuppressWarnings("unchecked")
  public static <X, Y extends Optimizable<X>> void runAsDrone(
      final MPIInterface mpi, final long waittime) throws Exception {

    // init Lottery - it is unreasonable to assume a distributed job could synchronize Lottery,
    // hence default init
    final Random r = new Random();
    final long seed = r.nextLong();
    final RNGenerator rng = new StandardRNG(seed);
    Lottery.setGenerator(rng);

    final int myRank = mpi.mpiCommRank();

    final boolean debug = LOG.isDebugEnabled();
    if (debug)
      LOG.debug(
          "Entering generic MPI globopt as drone ( " + myRank + " ). BRACE YOURSELF BIG TIMES!");

    // try to receive the initial broadcast
    final char[] initMessage = new char[34];
    final Tuple<Integer, char[]> bcastRet = mpi.mpiBcast(initMessage, myRank, 0);
    if (bcastRet.getObject1() != 0)
      throw new RuntimeException("Error in receiving initial bcast. " + bcastRet);
    final String sInit = new String(bcastRet.getObject2()).trim();

    if (!sInit.equalsIgnoreCase("Hello from MPI queen, all is well!")) {
      throw new RuntimeException(
          "Initial message from queen was not what "
              + myRank
              + " expected. Exiting. Received message: "
              + sInit);
    }

    // do stuff as long as there is something to do
    int taskCounter = 0;
    for (; ; ) {

      // always the same idea: ask the queen for a new task
      if (debug) LOG.debug("Drone " + myRank + " asking for new task...");
      final MPIInterface.MPIStatus status = mpi.mpiRecvBytes(myRank, 0, MPIInterface.ANY_TAG);
      if (status.err != 0) {
        throw new RuntimeException("Failure to receive on " + myRank + " with error " + status.err);
      }

      if (debug) LOG.debug("Drone " + myRank + " received " + status.tag);

      final int messTag = status.tag;
      if (messTag == EXITDONE || (messTag != WAITFOR && messTag != NEXTTASK)) {
        if (debug) LOG.debug("Queen tells me to quit: doing so! Tag was " + messTag);
        break;
      }

      if (messTag == WAITFOR) {
        // sleep some to not congest the queen
        Thread.sleep(waittime);

        // dummy send
        if (debug) LOG.debug("Drone " + myRank + " dummy sending... ");
        final int retSend = mpi.mpiSend(new byte[1], myRank, 0, DUMMY);
        if (retSend != 0)
          throw new RuntimeException(
              "Error in sending result back from " + myRank + " with error " + retSend);
        if (debug) LOG.debug("Drone " + myRank + " sent dummy.");
      } else {

        // there is a task: execute it
        taskCounter++;
        if (debug) LOG.debug("There is a new task! No: " + taskCounter);

        final Task<Y> task = (Task<Y>) InputPrimitives.readByteInput(status.msg);
        final Result<Y> result = task.executeTask(myRank);

        // report back
        final byte[] resData = OutputPrimitives.writeObjToByteArray(result);
        if (debug) LOG.debug("Drone " + myRank + " sending result...");
        final int retSend = mpi.mpiSend(resData, myRank, 0, RESULT);
        if (retSend != 0)
          throw new RuntimeException(
              "Error in sending result back from " + myRank + " with error " + retSend);
        if (debug) LOG.debug("Drone " + myRank + " sent result.");
      }
    }
  }
}
