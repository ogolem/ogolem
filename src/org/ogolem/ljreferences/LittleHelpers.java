/**
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
package org.ogolem.ljreferences;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author Johannes Dieterich
 * @version 2010-03-08
 */
final class LittleHelpers {

    static String[] readFileIn(final String sFileName) throws IOException{
        
        BufferedReader buffreader = null;
        String line;
        int iNumberOfRows = 0;
        try {
            buffreader = new BufferedReader(new FileReader(sFileName));
            while ((line = buffreader.readLine()) != null) {
                iNumberOfRows++;
            }
        } catch (IOException e) {
            throw new IOException("Error occured during reading in file!", e);
        } finally {
            if (buffreader != null) {
                try {
                    buffreader.close();
                } catch (IOException e) {
                    throw new IOException("Error occured during closing file!", e);
                }
            }
        }
        String[] saInputData = new String[iNumberOfRows];
        int i = 0;
        try {
            buffreader = new BufferedReader(new FileReader(sFileName));
            while ((line = buffreader.readLine()) != null) {
                saInputData[i] = line;
                i++;
            }
        } catch (IOException e) {
            throw new IOException("Error occured during reading in file!", e);
        } finally {
            if (buffreader != null) {
                try {
                    buffreader.close();
                } catch (IOException e) {
                    throw new IOException("Error occured during closing file!", e);
                }
            }
        }
        return saInputData;
    }


    static void writeNow(String sOutputPath, String[] saToWrite) throws IOException {

        BufferedWriter buffwriter = null;
        int iLength = saToWrite.length;

        try {
            buffwriter = new BufferedWriter(new FileWriter(sOutputPath, true));
            for (int i = 0; i < iLength; i++) {
                buffwriter.write(saToWrite[i]);
                buffwriter.write(System.getProperty("line.separator"));
            }
        } catch (IOException e) {
            throw new IOException("Error occured during writing!", e);
        } finally {
            if (buffwriter != null) {
                try {
                    buffwriter.close();
                } catch (IOException e) {
                    throw new IOException("Error occured after writing, when trying to close files!", e);
                }
            }
        }
    }
}
