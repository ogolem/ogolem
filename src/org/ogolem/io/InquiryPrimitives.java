/**
Copyright (c) 2013, J. M. Dieterich
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

import java.io.File;
import java.io.IOException;

/**
 * Primitive functions for inquiring the sate of files/folders/etc.
 * @author Johannes Dieterich
 * @version 2013-01-05
 */
public class InquiryPrimitives {
    
    public static boolean doesFileExist(final String path){
        final File f = new File(path);
        
        return f.exists();
    }
    
   public static String[] fileListWithPrefix(final String prefix) throws IOException {
        final FileFilter filter = new FileFilter(prefix, true);
        final File f = new File(System.getProperty("user.dir"));
        final String[] saList = f.list(filter);
        
        return saList;
    }

    public static String[] fileListWithSuffix(final String suffix) throws IOException {
        final FileFilter filter = new FileFilter(suffix, false);
        final File f = new File(System.getProperty("user.dir"));
        final String[] saList = f.list(filter);

        return saList;
    }

    public static String[] fileListWithPrefix(final String prefix, final String dir) throws IOException {
        final FileFilter filter = new FileFilter(prefix, true);
        final File f = new File(dir);
        final String[] saList = f.list(filter);

        return saList;
    }

    public static String[] fileListWithSuffix(final String suffix, final String dir) throws IOException {
        final FileFilter filter = new FileFilter(suffix, false);
        final File f = new File(dir);
        final String[] saList = f.list(filter);

        return saList;
    }

    public static String[] fileList(final String path){
        final File f = new File(path);
        final String[] saList = f.list();

        return saList;
    }
}
