/**
Copyright (c) 2010     , J. M. Dieterich
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
package org.ogolem.core;

import java.io.Serializable;
import java.util.Calendar;
import java.text.SimpleDateFormat;
import org.ogolem.io.OutputPrimitives;

/**
 * This is a really simple stub for a statistics object.
 * @author Johannes Dieterich
 * @version 2010-11-24
 */
public final class GlobOptStatistics implements Serializable{

    private static final long serialVersionUID = (long) 20101124;

    private static final String DATE_FORMAT_NOW = "yyyy-MM-dd HH:mm:ss";

    private final String logFile;

    /**
     * Constructs the statistics object.
     * @param log The path to the log file.
     */
    public GlobOptStatistics(final String log){
        this.logFile = log;
    }

    /**
     * If the individual is the new best individual, a line is printed to the
     * log file with the current date and time.
     * @param individual the ID of the individual
     * @param position the position in the pool
     */
    public void individualAddedToPool(long individual, int position, double fitness){
        if(position == 0){
            // print a a line to file
            final Calendar cal = Calendar.getInstance();
            final SimpleDateFormat sdf = new SimpleDateFormat(DATE_FORMAT_NOW);
            final String date = sdf.format(cal.getTime());
            final String line = "[" + date + "]" + "   new best individual " + individual + "\t fitness " + fitness;
            final String[] data = {line};
            
            try{
                OutputPrimitives.writeOut(logFile, data, true);
            } catch(Exception e){
                System.err.println("WARNING: Couldn't write to log file. " + e.toString());
            }
        }
    }
}
