/*
Copyright (c) 2020-2021, J. M. Dieterich and B. Hartke
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

import static jdk.incubator.foreign.CLinker.*;

import java.io.Serializable;
import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodType;
import jdk.incubator.foreign.*;
import org.ogolem.helpers.Tuple;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Interface through Java's foreign interface to an MPI library. Assumes a COMM_WORLD communicator
 * everywhere. Tested with OpenMPI.
 *
 * @author Johannes M Dieterich
 * @version 2021-07-27
 */
public class MPIInterface implements Serializable {

  private static final long serialVersionUID = 20210727;
  private static final Logger LOG = LoggerFactory.getLogger(MPIInterface.class);

  /* the following constants are our internal choices for ANY_SOURCE, ANY_TAG which the wrapper translates into the platform-dependent ones */
  public static final int ANY_SOURCE = -42;
  public static final int ANY_TAG = -42;

  private final boolean debug;
  private final String soName;
  private final LibraryLookup libmpi;
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
    this.libmpi = LibraryLookup.ofLibrary(soName);

    final LibraryLookup.Symbol mpiInitSym = libmpi.lookup("ogo_MPI_Init").get();
    final LibraryLookup.Symbol mpiAbortSym = libmpi.lookup("ogo_MPI_Abort").get();
    final LibraryLookup.Symbol mpiFinalizeSym = libmpi.lookup("ogo_MPI_Finalize").get();
    final LibraryLookup.Symbol mpiWtimeSym = libmpi.lookup("ogo_MPI_Wtime").get();
    final LibraryLookup.Symbol mpiRankSym = libmpi.lookup("ogo_MPI_Comm_rank").get();
    final LibraryLookup.Symbol mpiSizeSym = libmpi.lookup("ogo_MPI_Comm_size").get();
    final LibraryLookup.Symbol mpiProbeSym = libmpi.lookup("ogo_MPI_Probe_bytes").get();
    final LibraryLookup.Symbol mpiSendSym = libmpi.lookup("ogo_MPI_Send_bytes").get();
    final LibraryLookup.Symbol mpiRecvSym = libmpi.lookup("ogo_MPI_Recv_bytes").get();
    final LibraryLookup.Symbol mpiBcastSym = libmpi.lookup("ogo_MPI_Bcast_char").get();

    this.mpiInitH =
        CLinker.getInstance()
            .downcallHandle(
                mpiInitSym, MethodType.methodType(int.class), FunctionDescriptor.of(C_INT));

    this.mpiAbortH =
        CLinker.getInstance()
            .downcallHandle(
                mpiAbortSym,
                MethodType.methodType(int.class, int.class),
                FunctionDescriptor.of(C_INT, C_INT));

    this.mpiFinalizeH =
        CLinker.getInstance()
            .downcallHandle(
                mpiFinalizeSym, MethodType.methodType(int.class), FunctionDescriptor.of(C_INT));

    this.mpiWtimeH =
        CLinker.getInstance()
            .downcallHandle(
                mpiWtimeSym, MethodType.methodType(double.class), FunctionDescriptor.of(C_DOUBLE));

    this.mpiRankH =
        CLinker.getInstance()
            .downcallHandle(
                mpiRankSym,
                MethodType.methodType(int.class, MemoryAddress.class),
                FunctionDescriptor.of(C_INT, C_POINTER));

    this.mpiSizeH =
        CLinker.getInstance()
            .downcallHandle(
                mpiSizeSym,
                MethodType.methodType(int.class, MemoryAddress.class),
                FunctionDescriptor.of(C_INT, C_POINTER));

    this.mpiProbeH =
        CLinker.getInstance()
            .downcallHandle(
                mpiProbeSym,
                MethodType.methodType(int.class, MemoryAddress.class),
                FunctionDescriptor.of(C_INT, C_POINTER));

    this.mpiSendH =
        CLinker.getInstance()
            .downcallHandle(
                mpiSendSym,
                MethodType.methodType(
                    int.class, MemoryAddress.class, int.class, int.class, int.class),
                FunctionDescriptor.of(C_INT, C_POINTER, C_INT, C_INT, C_INT));

    this.mpiRecvH =
        CLinker.getInstance()
            .downcallHandle(
                mpiRecvSym,
                MethodType.methodType(
                    int.class, MemoryAddress.class, int.class, int.class, int.class),
                FunctionDescriptor.of(C_INT, C_POINTER, C_INT, C_INT, C_INT));

    this.mpiBcastH =
        CLinker.getInstance()
            .downcallHandle(
                mpiBcastSym,
                MethodType.methodType(int.class, MemoryAddress.class, int.class, int.class),
                FunctionDescriptor.of(C_INT, C_POINTER, C_INT, C_INT));
  }

  public int mpiInit(final String[] args) {

    int initRet = 0;
    try (MemorySegment argc = MemorySegment.allocateNative(4);
        MemorySegment argv = MemorySegment.ofArray(new char[] {'o'})) {
      argc.copyFrom(MemorySegment.ofArray(new int[] {0}));
      initRet = (int) mpiInitH.invokeExact();
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
      final int ret = (int) mpiAbortH.invokeExact(errCode);
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
    try {
      final int ret = (int) mpiFinalizeH.invokeExact();
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
      final double wtime = (double) mpiWtimeH.invokeExact();
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
    try (MemorySegment rankSeg = MemorySegment.allocateNative(4)) {
      ret = (int) mpiRankH.invokeExact(rankSeg.address());
      rank = rankSeg.toIntArray()[0];
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
    try (MemorySegment sizeSeg = MemorySegment.allocateNative(4)) {
      ret = (int) mpiSizeH.invokeExact(sizeSeg.address());
      size = sizeSeg.toIntArray()[0];
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

    final String s = new String(message);
    final int count = s.length();

    int ret = 0;
    try (MemorySegment msgSeg =
        MemorySegment.allocateNative(count + 1)) { // +1 for NULL termination
      msgSeg.copyFrom(CLinker.toCString(s));
      ret = (int) mpiBcastH.invokeExact(msgSeg.address(), count, root);

      if (myRank != root) {
        message = CLinker.toJavaString(msgSeg).toCharArray();
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
    try (MemorySegment msgSeg = MemorySegment.allocateNative(count)) {
      msgSeg.copyFrom(MemorySegment.ofArray(message));
      ret = (int) mpiSendH.invokeExact(msgSeg.address(), count, toRank, tag);

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
    try (MemorySegment dataSeg = MemorySegment.allocateNative(3 * 4)) {
      dataSeg.copyFrom(MemorySegment.ofArray(inputData));
      ret = (int) mpiProbeH.invokeExact(dataSeg.address());
      outputData = dataSeg.toIntArray();

    } catch (Exception e) {
      e.printStackTrace();
      return new MPIStatus(null, sourceRank, tag, 42);
    } catch (Throwable t) {
      t.printStackTrace();
      return new MPIStatus(null, sourceRank, tag, 43);
    }

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
    try (MemorySegment dataSeg = MemorySegment.allocateNative(outputData[2])) {
      ret =
          (int)
              mpiRecvH.invokeExact(dataSeg.address(), outputData[2], outputData[0], outputData[1]);
      data = dataSeg.toByteArray();

    } catch (Exception e) {
      e.printStackTrace();
      return new MPIStatus(null, sourceRank, tag, 42);
    } catch (Throwable t) {
      t.printStackTrace();
      return new MPIStatus(null, sourceRank, tag, 43);
    }

    if (debug)
      LOG.info("Successful MPI_Recv() from " + sourceRank + " on " + myRank + " with tag " + tag);

    // assemble the status object to return
    return new MPIStatus(data, outputData[0], outputData[1], ret);
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
