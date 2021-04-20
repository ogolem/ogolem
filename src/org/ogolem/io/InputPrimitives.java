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
import java.nio.charset.Charset;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * All input wrapped up.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public final class InputPrimitives {

  /**
   * Reads a serialized object.
   *
   * @param fileName
   * @return an object for the serialized data
   * @throws IOException
   * @throws ClassNotFoundException
   */
  public static Object readBinInput(final String fileName)
      throws IOException, ClassNotFoundException {

    final File file = new File(fileName);
    if (!file.exists()) {
      throw new IOException("File does not exist.");
    } else if (file.length() > FixedValues.MAXFILESIZEREADING) {
      throw new IOException("File way bigger than expected. Must be garbage.");
    }

    Object obj = null;
    try (final ObjectInputStream objectStream =
        new ObjectInputStream(new FileInputStream(fileName))) {
      obj = objectStream.readObject();
    } catch (IOException | ClassNotFoundException e) {
      throw e;
    }

    return obj;
  }

  /**
   * Reads a file into memory.
   *
   * @param fileName
   * @return the complete file content as a String array
   * @throws IOException
   */
  public static String[] readFileIn(final String fileName) throws IOException {

    final File file = new File(fileName);
    if (!file.exists()) {
      throw new IOException("File " + fileName + " does not exist.");
    } else if (file.length() > FixedValues.MAXFILESIZEREADING) {
      throw new IOException("File " + fileName + " way bigger than expected. Must be garbage.");
    }

    String line;
    final LinkedList<String> ll = new LinkedList<>();
    try (final BufferedReader buffreader =
        new BufferedReader(
            new InputStreamReader(new FileInputStream(fileName), Charset.forName("UTF-8")))) {
      while ((line = buffreader.readLine()) != null) {
        ll.add(line);
      }
    } catch (IOException e) {
      throw e;
    }

    final Iterator<String> it = ll.iterator();
    final String[] data = new String[ll.size()];

    int counter = 0;
    while (it.hasNext()) {
      data[counter] = it.next();
      counter++;
    }

    return data;
  }
}
