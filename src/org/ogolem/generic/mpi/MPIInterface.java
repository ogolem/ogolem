/*
Copyright (c) 2020-2023, J. M. Dieterich and B. Hartke
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

import java.io.Serializable;
import java.lang.foreign.*;
import java.lang.invoke.MethodHandle;
import org.ogolem.helpers.Tuple;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Interface through Java's foreign interface to an MPI library. Assumes a COMM_WORLD communicator
 * everywhere. Tested with OpenMPI.
 *
 * @author Johannes M Dieterich
 * @version 2023-09-03
 */
public class MPIInterface {

  private static final Logger LOG = LoggerFactory.getLogger(MPIInterface.class);

  private static boolean initialized = false;
  private static boolean finalized = false;

  /* the following constants are our internal choices for ANY_SOURCE, ANY_TAG which the wrapper translates into the platform-dependent ones */
  public static final int ANY_SOURCE = -42;
  public static final int ANY_TAG = -42;

  private final boolean debug;
  private final String soName;
  private final Arena arena;
  private final MethodHandle mpiInitH;
  private final MethodHandle mpiAbortH;
  private final MethodHandle mpiFinalizeH;
  private final MethodHandle mpiWtimeH;
  private final MethodHandle mpiRankH;
  private final MethodHandle mpiSizeH;
  private final MethodHandle mpiProbeH;
  private final MethodHandle mpiSendH;
  private final MethodHandle mpiRecvH;
  private final MethodHandle mpiBcastH;

  public MPIInterface(final String soName, final boolean debug) throws Exception {
    assert (soName != null);
    assert (!soName.isEmpty());

    this.debug = debug;
    this.soName = soName;
    this.arena = Arena.ofShared();

    Linker linker = Linker.nativeLinker();

    try (Arena arena = Arena.ofConfined()) {

      System.loadLibrary(soName);
      final SymbolLookup libmpi = SymbolLookup.loaderLookup();
      final MemorySegment mpiInitSeg = libmpi.find("ogo_MPI_Init").orElseThrow();
      final MemorySegment mpiAbortSeg = libmpi.find("ogo_MPI_Abort").orElseThrow();
      final MemorySegment mpiFinalizeSeg = libmpi.find("ogo_MPI_Finalize").orElseThrow();
      final MemorySegment mpiWtimeSeg = libmpi.find("ogo_MPI_Wtime").orElseThrow();
      final MemorySegment mpiRankSeg = libmpi.find("ogo_MPI_Comm_rank").orElseThrow();
      final MemorySegment mpiSizeSeg = libmpi.find("ogo_MPI_Comm_size").orElseThrow();
      final MemorySegment mpiProbeSeg = libmpi.find("ogo_MPI_Probe_bytes").orElseThrow();
      final MemorySegment mpiSendSeg = libmpi.find("ogo_MPI_Send_bytes").orElseThrow();
      final MemorySegment mpiRecvSeg = libmpi.find("ogo_MPI_Recv_bytes").orElseThrow();
      final MemorySegment mpiBcastSeg = libmpi.find("ogo_MPI_Bcast_char").orElseThrow();

      final var mpiInitF = FunctionDescriptor.of(ValueLayout.JAVA_INT);
      this.mpiInitH = linker.downcallHandle(mpiInitSeg, mpiInitF);

      final var mpiAbortF = FunctionDescriptor.of(ValueLayout.JAVA_INT, ValueLayout.JAVA_INT);
      this.mpiAbortH = linker.downcallHandle(mpiAbortSeg, mpiAbortF);

      final var mpiFinalizeF = FunctionDescriptor.of(ValueLayout.JAVA_INT);
      this.mpiFinalizeH = linker.downcallHandle(mpiFinalizeSeg, mpiFinalizeF);

      final var mpiWtimeF = FunctionDescriptor.of(ValueLayout.JAVA_DOUBLE);
      this.mpiWtimeH = linker.downcallHandle(mpiWtimeSeg, mpiWtimeF, Linker.Option.critical(true));

      final var mpiRankF = FunctionDescriptor.of(ValueLayout.JAVA_INT, ValueLayout.ADDRESS);
      this.mpiRankH = linker.downcallHandle(mpiRankSeg, mpiRankF);

      final var mpiSizeF = FunctionDescriptor.of(ValueLayout.JAVA_INT, ValueLayout.ADDRESS);
      this.mpiSizeH = linker.downcallHandle(mpiSizeSeg, mpiSizeF);

      final var mpiProbeF = FunctionDescriptor.of(ValueLayout.JAVA_INT, ValueLayout.ADDRESS);
      this.mpiProbeH = linker.downcallHandle(mpiProbeSeg, mpiProbeF);

      final var mpiSendF =
          FunctionDescriptor.of(
              ValueLayout.JAVA_INT,
              ValueLayout.ADDRESS,
              ValueLayout.JAVA_INT,
              ValueLayout.JAVA_INT,
              ValueLayout.JAVA_INT);
      this.mpiSendH = linker.downcallHandle(mpiSendSeg, mpiSendF);

      final var mpiRecvF =
          FunctionDescriptor.of(
              ValueLayout.JAVA_INT,
              ValueLayout.ADDRESS,
              ValueLayout.JAVA_INT,
              ValueLayout.JAVA_INT,
              ValueLayout.JAVA_INT);
      this.mpiRecvH = linker.downcallHandle(mpiRecvSeg, mpiRecvF);

      final var mpiBcastF =
          FunctionDescriptor.of(
              ValueLayout.JAVA_INT,
              ValueLayout.ADDRESS,
              ValueLayout.JAVA_INT,
              ValueLayout.JAVA_INT);
      this.mpiBcastH = linker.downcallHandle(mpiBcastSeg, mpiBcastF);
    }
  }

  public int mpiInit(final String[] args) {

    if (initialized) {
      LOG.warn("MPI_Init() already called. Ignoring.");
      return 0;
    }

    int initRet = 0;
    try {
      initRet = (int) mpiInitH.invoke();
      initialized = true;
      Runtime.getRuntime()
          .addShutdownHook(
              new Thread(
                  () -> {
                    if (initialized && !finalized) {
                      LOG.info("Shutdown hook: finalizing MPI.");
                      mpiFinalize();
                    }
                  }));
    } catch (Exception e) {
      e.printStackTrace(System.err);
      return 42;
    } catch (Throwable t) {
      t.printStackTrace(System.err);
      return 43;
    }

    if (debug) LOG.info("Successful MPI_Init().");

    return initRet;
  }

  public int mpiAbort(final int errCode) {
    try {
      final int ret = (int) mpiAbortH.invoke(errCode);
      if (debug) LOG.info("Successful MPI_Abort() w/ " + errCode);
      return ret;
    } catch (Exception e) {
      e.printStackTrace();
      return 42;
    } catch (Throwable t) {
      t.printStackTrace(System.err);
      return 43;
    }
  }

  public int mpiFinalize() {

    if (finalized) {
      LOG.warn("MPI_Finalize() already called. Ignoring.");
      return 0;
    }

    try {
      final int ret = (int) mpiFinalizeH.invoke();
      finalized = true;
      if (debug) LOG.info("Successful MPI_Finalize().");
      return ret;
    } catch (Exception e) {
      e.printStackTrace();
      return 42;
    } catch (Throwable t) {
      t.printStackTrace(System.err);
      return 43;
    }
  }

  public double mpiWtime() {
    try {
      final double wtime = (double) mpiWtimeH.invoke();
      if (debug) LOG.info("Successful MPI_Wtime() w/ " + wtime);
      return wtime;
    } catch (Exception e) {
      e.printStackTrace();
      return 42.0;
    } catch (Throwable t) {
      t.printStackTrace(System.err);
      return 43.0;
    }
  }

  public int mpiCommRank() throws Exception {

    int rank = -999;
    int ret = 0;
    MemorySegment rankSeg = arena.allocate(ValueLayout.JAVA_INT);
    try {
      ret = (int) mpiRankH.invoke(rankSeg);
      rank = rankSeg.getAtIndex(ValueLayout.JAVA_INT, 0);
    } catch (Exception e) {
      throw e;
    } catch (Throwable t) {
      throw new Exception(t);
    }

    if (ret != 0 || rank < 0) {
      throw new RuntimeException("Failure in MPI_Comm_rank() " + ret + "\t, rank: " + rank);
    }

    if (debug) LOG.info("Successful MPI_Comm_rank() w/ " + rank);

    return rank;
  }

  public int mpiCommSize() throws Exception {

    int size = -999;
    int ret = 0;
    MemorySegment sizeSeg = arena.allocate(ValueLayout.JAVA_INT);
    try {
      ret = (int) mpiSizeH.invoke(sizeSeg);
      size = sizeSeg.getAtIndex(ValueLayout.JAVA_INT, 0);
    } catch (Exception e) {
      throw e;
    } catch (Throwable t) {
      throw new Exception(t);
    }

    if (ret != 0 || size <= 0) {
      throw new RuntimeException("Failure in MPI_Comm_size() " + ret + "\t, size: " + size);
    }

    if (debug) LOG.info("Successful MPI_Comm_size() w/ " + size);

    return size;
  }

  public Tuple<Integer, char[]> mpiBcast(char[] message, final int myRank, final int root) {

    final int count = message.length;

    int ret = 0;
    try {
      final byte[] bytes = new String(message).getBytes(java.nio.charset.StandardCharsets.UTF_8);
      MemorySegment msgSeg = arena.allocateFrom(ValueLayout.JAVA_BYTE, bytes);
      ret = (int) mpiBcastH.invoke(msgSeg, bytes.length, root);

      if (myRank != root) {
        final byte[] rcvBytes = msgSeg.toArray(ValueLayout.JAVA_BYTE);
        final String s = new String(rcvBytes, java.nio.charset.StandardCharsets.UTF_8);
        final char[] c = s.toCharArray();
        System.arraycopy(c, 0, message, 0, Math.min(c.length, message.length));
      }

    } catch (Exception e) {
      e.printStackTrace();
      return new Tuple<>(42, message);
    } catch (Throwable t) {
      t.printStackTrace();
      return new Tuple<>(43, message);
    }

    if (debug)
      LOG.info(
          "Successful MPI_Bcast() w/ "
              + new String(message)
              + "(count:"
              + count
              + ") on rank "
              + myRank);

    return new Tuple<>(ret, message);
  }

  public int mpiSend(final byte[] message, final int myRank, final int toRank, final int tag) {

    final int count = message.length;

    int ret = 0;
    try {
      MemorySegment msgSeg = arena.allocateFrom(ValueLayout.JAVA_BYTE, message);
      ret = (int) mpiSendH.invoke(msgSeg, count, toRank, tag);

    } catch (Exception e) {
      e.printStackTrace();
      return 42;
    } catch (Throwable t) {
      t.printStackTrace();
      return 43;
    }

    if (debug)
      LOG.info(
          "Successful MPI_Send() w/ (count:"
              + count
              + ") from rank "
              + myRank
              + " to "
              + toRank
              + " with tag "
              + tag);

    return ret;
  }

  public MPIStatus mpiRecvBytes(final int myRank, final int sourceRank, final int tag) {

    // first MPI_Probe to figure out (if applicable) source, tag, buffer size
    final int[] inputData = new int[] {sourceRank, tag, 0};
    int[] outputData;

    int ret = 0;
    try {
      MemorySegment dataSeg = arena.allocateFrom(ValueLayout.JAVA_INT, inputData);
      ret = (int) mpiProbeH.invoke(dataSeg);
      outputData = dataSeg.toArray(ValueLayout.JAVA_INT);

      if (ret != 0) return new MPIStatus(null, sourceRank, tag, ret);

      if (debug)
        LOG.info(
            "Successful MPI_Probe() on "
                + myRank
                + " from "
                + sourceRank
                + " with tag "
                + tag
                + ", probe results: "
                + outputData[0]
                + " "
                + outputData[1]
                + " "
                + outputData[2]);

      // then actual Recv after allocating a fittingly sized buffer
      byte[] data;
      MemorySegment recvSeg = arena.allocate(outputData[2]);
      ret = (int) mpiRecvH.invoke(recvSeg, outputData[2], outputData[0], outputData[1]);
      data = recvSeg.toArray(ValueLayout.JAVA_BYTE);

      if (debug)
        LOG.info("Successful MPI_Recv() from " + sourceRank + " on " + myRank + " with tag " + tag);

      // assemble the status object to return
      return new MPIStatus(data, outputData[0], outputData[1], ret);

    } catch (Exception e) {
      e.printStackTrace();
      return new MPIStatus(null, sourceRank, tag, 42);
    } catch (Throwable t) {
      t.printStackTrace();
      return new MPIStatus(null, sourceRank, tag, 43);
    }
  }

  public static class MPIStatus implements Serializable {

    private static final long serialVersionUID = 20201223;

    public final int source;
    public final int tag;
    public final int err;
    public final byte[] msg;

    MPIStatus(final byte[] msg, final int source, final int tag, final int err) {
      this.source = source;
      this.tag = tag;
      this.err = err;
      this.msg = msg;
    }
  }
}
