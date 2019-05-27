/**
Copyright (c) 2010     , J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2016, J. M. Dieterich and B. Hartke
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

import java.io.Serializable;

/**
 * A configuration object for the generic history.
 * @author Johannes Dieterich
 * @version 2016-04-19
 */
public class GenericHistoryConfig implements Serializable {
    
    private static final long serialVersionUID = (long) 20131231;
    public int offset = 0;
    public int recordsToSerial = 1000;
    public int recordsToASCII = 1000;
    public String binOut = "genetic-history.bin";
    public String asciiAppend = "default-ascii-genhistory";
    public boolean silentMode;
    
    public String getMyConfig(){
        
        if(silentMode){
            return "silent mode is ON";
        }
        
        String s = "";
        s += "genetic history log file is:   " + asciiAppend + "\n"
           + "genetic history binary is:     " + binOut + "\n"
           + "records to accumulate for log: " + recordsToASCII + "\n"
           + "records to accumulate for bin: " + recordsToSerial + "\n"
           + "internal records offset:       " + offset + "\n";
        
        return s;
    }
}
