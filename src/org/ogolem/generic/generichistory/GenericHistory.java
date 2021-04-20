/*
Copyright (c) 2010, J. M. Dieterich and B. Hartke
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
package org.ogolem.generic.generichistory;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.io.Serializable;
import java.nio.charset.Charset;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import org.ogolem.generic.Optimizable;
import org.ogolem.generic.genericpool.GenericPoolEntry;
import org.ogolem.io.OutputPrimitives;

/**
 * This holds kind of a genetic record on all families that were created during the program's run.
 * It also tries to gain information from previous global optimization steps.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public final class GenericHistory<E, T extends Optimizable<E>> implements Serializable {

  private static final long serialVersionUID = (long) 20150710;

  private long totalNullCounter = 0;
  private long totalAccCounter = 0;
  private long totalCounter = 0;
  private final String whereToWriteSerial;
  private final boolean beSilent;
  private String whereToAppendASCII;
  private int counter = 0;
  private int serialCounter = 0;
  private final int addsToSerial;
  private final int addsToASCII;
  private boolean firstASCII = true;
  private List<GeneticRecord<E, T>> records;

  // the lock for the history
  private final ReentrantReadWriteLock lock;
  private final Lock roLock;
  private final Lock rwLock;

  /*
   * Constructor
   */
  /*
   * This is specific for the singleton design pattern, since the
   * constructor is private, we need to construct the object ourselves.
   */
  @SuppressWarnings("rawtypes")
  private static GenericHistory history;

  private GenericHistory(
      final int recordsToASCII,
      final int recordsToSerial,
      final String asciiAppend,
      final String binOut,
      final boolean beSilent) {

    this.records = new LinkedList<>();
    // the out can be set to this
    this.whereToWriteSerial = binOut;
    this.addsToASCII = recordsToASCII;
    this.addsToSerial = recordsToSerial;
    this.whereToAppendASCII = asciiAppend;

    // lock
    lock = new ReentrantReadWriteLock();
    roLock = lock.readLock();
    rwLock = lock.writeLock();
    this.beSilent = beSilent;
  }

  @SuppressWarnings({"rawtypes", "unchecked"})
  public static synchronized <E, T extends Optimizable<E>> GenericHistory<E, T> getReference(
      final GenericHistoryConfig config) {
    if (history == null)
      history =
          new GenericHistory(
              config.recordsToASCII,
              config.recordsToSerial,
              config.asciiAppend,
              config.binOut,
              config.silentMode);
    return history;
  }

  /*
   * Methods
   */

  /**
   * Sets the ASCII out to a specified file. Should be done before using the history.
   *
   * @param path
   */
  public void setASCIIOut(final String path) {
    rwLock.lock();
    try {
      this.whereToAppendASCII = path;
    } finally {
      rwLock.unlock();
    }
  }

  /**
   * The easier version of the overloaded addFamily method. This is should be used when no valuable
   * information was gained in the global optimization step (meaning: the geometry was returned as
   * being null).
   *
   * @param mother the id of the mother of this individual
   * @param father the id of the father of this individual
   * @param child the id of the child
   * @param accepted whether it got accepted to the pool
   * @param wasNullIndividual whether it was null
   */
  public void addFamily(
      final long mother,
      final long father,
      final int child,
      final boolean accepted,
      final boolean wasNullIndividual) {
    // this routine does not need to be synchronized as the other function is
    addFamily(mother, father, child, accepted, wasNullIndividual, null);
  }

  /**
   * This is the other overloaded addFamily method. In contrary to the previous one, one also hands
   * over the resulting geometry. Using this methods also triggers the serialization of the
   * GeneticHistory after a certain number of adds. This is a stub for cultural algorithms.
   *
   * @param mother the id of the mother of this individual
   * @param father the id of the father of this individual
   * @param child the id of the child
   * @param accepted whether it got accepted to the pool
   * @param wasNullIndividual whether it was null
   * @param individual the actual individual
   */
  public void addFamily(
      final long mother,
      final long father,
      final long child,
      final boolean accepted,
      final boolean wasNullIndividual,
      final T individual) {

    rwLock.lock();
    try {
      if (beSilent) {
        return;
      }

      // this is the overloaded version for the cultural algorithms... ;-)
      counter++;
      totalCounter++;
      if (wasNullIndividual) {
        totalNullCounter++;
      }
      if (accepted) {
        totalAccCounter++;
      }

      if (individual != null) {
        final GeneticRecord<E, T> record =
            new GeneticRecord<>(child, mother, father, accepted, wasNullIndividual);
        records.add(record);
        serialCounter++;

        if (serialCounter >= addsToSerial) {
          // serialize it and set it to zero again
          serializeMe();
          serialCounter = 0;
        }
      }

      if (counter >= addsToASCII) {
        writeASCII();
        counter = 0;
      }
    } finally {
      rwLock.unlock();
    }
  }

  private void writeASCII() {
    writeASCII2();
    // XXX this needs to be undone for cultural algorithms!
    records.clear();
  }

  public void flushRecords() {
    rwLock.lock();
    try {
      if (counter != 0) writeASCII();
    } finally {
      rwLock.unlock();
    }
  }

  private void writeASCII2() {

    final String sep = System.getProperty("line.separator");
    BufferedWriter buffwriter = null;
    try {
      buffwriter =
          new BufferedWriter(
              new OutputStreamWriter(
                  new FileOutputStream(whereToAppendASCII, true), Charset.forName("UTF-8")));

      if (firstASCII) {
        buffwriter.write("The following genetic history was created during the run.");
        buffwriter.write(sep);
        buffwriter.write(
            "\t"
                + "Child"
                + "\t"
                + "Mother"
                + "\t"
                + "Father"
                + "\t"
                + "got accepted?"
                + "\t"
                + "unsucessful genetic operations?");
        buffwriter.write(sep);
        firstASCII = false;
      }

      for (final GeneticRecord<E, T> record : records) {
        buffwriter.write(record.formattedHistoryLine() + sep);
      }

      buffwriter.close();
    } catch (IOException e) {
      System.err.println("WARNING: Couldn't write ASCII history. Continuing...");
      e.printStackTrace(System.err);
    } finally {
      if (buffwriter != null) {
        try {
          buffwriter.close();
        } catch (IOException e) {
          System.err.println("WARNING: Couldn't write ASCII history. Continuing...");
          e.printStackTrace(System.err);
        }
      }
    }
  }

  public void writeTotalStats() {

    roLock.lock();
    try {
      final String sep = System.getProperty("line.separator");
      BufferedWriter buffwriter = null;
      try {
        buffwriter =
            new BufferedWriter(
                new OutputStreamWriter(
                    new FileOutputStream(whereToAppendASCII, true), Charset.forName("UTF-8")));
        buffwriter.write("-----------------------------------------------------------" + sep);
        buffwriter.write("Overall genetic history" + sep);
        buffwriter.write("\ttotal number of trial genetic pool entries:    " + totalCounter + sep);
        buffwriter.write(
            "\tpercentage of unsuccessful genetic operations: "
                + String.format(
                    Locale.US, "%3.1f", (double) totalNullCounter / (double) totalCounter * 100)
                + "%"
                + sep);
        buffwriter.write(
            "\tpercentage of accepted pool entries:           "
                + String.format(
                    Locale.US, "%3.1f", (double) totalAccCounter / (double) totalCounter * 100)
                + "%"
                + sep);
        buffwriter.write("-----------------------------------------------------------" + sep);

        buffwriter.close();
      } catch (IOException e) {
        System.err.println("WARNING: Couldn't write total stats in history. Continuing...");
        e.printStackTrace(System.err);
      } finally {
        if (buffwriter != null) {
          try {
            buffwriter.close();
          } catch (IOException e) {
            System.err.println("WARNING: Couldn't write total stats in history. Continuing...");
            e.printStackTrace(System.err);
          }
        }
      }
    } finally {
      roLock.unlock();
    }
  }

  public void addInitialGeometries(List<T> initialGeoms) {
    // XXX this is a stub for the cultural algorithms
  }

  public void addInitialGeometries(Iterator<GenericPoolEntry<E, T>> initialGeoms) {
    // XXX this is a stub for the cultural algorithms
  }

  private void serializeMe() {

    try {
      OutputPrimitives.writeObjToBinFile(whereToWriteSerial, this);
    } catch (Exception e) {
      // some error occured during the serialization attempt, try again
      e.printStackTrace(System.err);
      try {
        OutputPrimitives.writeObjToBinFile(whereToWriteSerial, this);
      } catch (Exception e2) {
        // ok, it failed a second time. too bad...
        e2.printStackTrace(System.err);
        System.err.println(
            "Double failure in trying to write the genetic history out. Continuing without guarantee now!");
      }
    }
  }

  // TODO cultural algorithms.
  /**
   * The specific method for writing the GeneticHistory out (serialized). Needed due to some static
   * fields that should be conserved.
   */
  private void writeObject(ObjectOutputStream oos) throws IOException {
    oos.writeObject(records);
  }

  /**
   * The specific method for reading the GeneticHistory in, when serialized. Needed due to some
   * static fields that should be conserved. Suppresses the unchecked warning for an conversion to
   * the LinkedList of GeneticRecord objects.
   */
  @SuppressWarnings(value = "unchecked")
  private void readObject(ObjectInputStream ois) throws IOException, ClassNotFoundException {
    records = (LinkedList<GeneticRecord<E, T>>) ois.readObject();
  }
}
