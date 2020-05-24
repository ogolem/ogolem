/**
Copyright (c) 2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.macrobenchmarks;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Some helpers solely for the macrobenchmarks
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
class Helpers {
    
    private static final Logger LOG = LoggerFactory.getLogger(Helpers.class);    
    private static final boolean DEBUG = false;

    static int executeJavaProcess(final String workDir, final String[] args) throws Exception {
        
        final String javaHome = System.getProperty("java.home");
        final String javaBin = javaHome + File.separator + "bin" + File.separator + "java";
        final String pathToOurJar = new File(System.getProperty("java.class.path")).getAbsolutePath();
 
        if(DEBUG){
            final Properties p = System.getProperties();
            p.list(System.out);
        }

        final List<String> command = new ArrayList<>();
        command.add(javaBin);
        command.add("-ea");
        command.add("-jar");
        command.add(pathToOurJar);
        
        for(final String arg : args){
            command.add(arg);
        }
        
        LOG.debug("Executing " + javaBin + " with jar " + pathToOurJar);
        
        final ProcessBuilder builder = new ProcessBuilder(command);
        final Process proc = builder.directory(new File(workDir)).redirectOutput(new File(workDir + File.separator + "bench.out")).redirectError(new File(workDir + File.separator + "bench.err")).start();
        proc.waitFor();
                
        return proc.exitValue();
    }
    
    static class Rank0Filter implements FilenameFilter {

        @Override
        public boolean accept(File arg0, String arg1) {
            return (arg1.startsWith("rank0individual"));
        }
    }
}
