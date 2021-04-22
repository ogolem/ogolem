/*
Copyright (c) 2010-2013, J. M. Dieterich
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
package org.ogolem.io;

import java.io.*;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

/**
 * All output wrapped up.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public final class OutputPrimitives {

  /**
   * Write an object as a binary (assuming that it supports serialization).
   *
   * @param fileName
   * @param obj
   * @throws IOException
   */
  public static void writeObjToBinFile(final String fileName, final Serializable obj)
      throws IOException {

    try (final ObjectOutputStream outStream =
        new ObjectOutputStream(new FileOutputStream(fileName))) {
      outStream.writeObject(obj);
    } catch (IOException e) {
      throw e;
    }
  }

  /**
   * Writes an array of data to a file.
   *
   * @param path
   * @param data
   * @param noOverride
   * @throws IOException
   */
  public static void writeOut(final String path, final String[] data, final boolean noOverride)
      throws IOException {

    final String sep = System.lineSeparator();
    try (final BufferedWriter buffwriter =
        new BufferedWriter(
            new OutputStreamWriter(
                new FileOutputStream(path, noOverride), Charset.forName("UTF-8")))) {
      for (final String line : data) {
        buffwriter.write(line);
        buffwriter.write(sep);
      }
    } catch (IOException e) {
      throw e;
    }
  }

  /**
   * Writes an arraylist of data to a file.
   *
   * @param path
   * @param data
   * @param noOverride
   * @throws IOException
   */
  public static void writeOut(final String path, final List<String> data, final boolean noOverride)
      throws IOException {

    final String sep = System.lineSeparator();
    try (final BufferedWriter buffwriter =
        new BufferedWriter(
            new OutputStreamWriter(
                new FileOutputStream(path, noOverride), Charset.forName("UTF-8")))) {
      for (final String line : data) {
        buffwriter.write(line);
        buffwriter.write(sep);
      }
    } catch (IOException e) {
      throw e;
    }
  }

  /**
   * Writes a String to a file.
   *
   * @param path
   * @param data
   * @param noOverride
   * @throws IOException
   */
  public static void writeOut(final String path, final String data, final boolean noOverride)
      throws IOException {

    try (final BufferedWriter buffwriter =
        new BufferedWriter(
            new OutputStreamWriter(
                new FileOutputStream(path, noOverride), Charset.forName("UTF-8")))) {
      buffwriter.write(data);
    } catch (IOException e) {
      throw e;
    }
  }

  /**
   * Creates a symbolic link.
   *
   * @param orig the original file/folder
   * @param link the link name
   * @throws IOException
   */
  public static void createLink(final String orig, final String link) throws IOException {

    final Path target = Paths.get(orig).toAbsolutePath();
    final Path newLink = Paths.get(link).toAbsolutePath();

    final File linkF = newLink.toFile();
    if (linkF.exists()) {
      System.err.println(
          "WARNING: Link " + link + " exists. Assuming previous linking attempt and returning.");
      return;
    }

    try {
      Files.createSymbolicLink(newLink, target);
    } catch (UnsupportedOperationException x) {
      throw new IOException(x);
    }
  }

  public static void copyFile(final String orig, final String copy) throws IOException {

    try (final FileChannel inChannel = new FileInputStream(new File(orig)).getChannel();
        final FileChannel outChannel = new FileOutputStream(new File(copy)).getChannel()) {
      inChannel.transferTo(0, inChannel.size(), outChannel);
    } catch (IOException e) {
      throw new IOException("Failure to copy file.", e);
    }
  }

  public static void createAFolder(final String folderName) throws IOException {

    final File f = new File(folderName);

    // try to figure whether the file exists
    if (f.exists() && f.isDirectory()) {
      return;
    } else if (f.exists() && !f.isDirectory()) {
      System.err.println(
          "ERROR: Failure to create folder " + folderName + ". A file of this name exists.");
      return;
    }

    if (!f.mkdirs()) {
      System.err.println("ERROR: Failure to create folder " + folderName + ". Reason unknown.");
    }
  }

  public static void createAFolder(final String folderName, final boolean tmpDir)
      throws IOException {

    String parentFolder =
        (tmpDir) ? System.getProperty("java.io.tmpdir") : System.getProperty("user.dir");
    if (!parentFolder.endsWith(System.getProperty("file.separator"))) {
      parentFolder += System.getProperty("file.separator");
    }

    final String totalPath = parentFolder + folderName;
    createAFolder(totalPath);
  }

  // TODO add more features and replace the various i/o stuff in the rest of the framework, also
  // remove all custom IOhandlers, Input and Output classes as far as possible
}
