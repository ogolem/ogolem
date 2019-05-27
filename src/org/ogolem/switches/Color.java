/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010     , J. M. Dieterich
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
package org.ogolem.switches;

import java.io.Serializable;

/**
 * Value object for the colors.
 * @author Johannes Dieterich
 * @version 2010-11-14
 */
public final class Color implements Serializable{

    private static final long serialVersionUID = (long) 20100121;

    private final int iThisColor;

    private final String sThisColor;

    Color(final String sColor){
        if(sColor.equalsIgnoreCase("black")){
            iThisColor = 0;
        } else if(sColor.equalsIgnoreCase("white")){
            iThisColor = 1;
        } else if(sColor.equalsIgnoreCase("red")){
            iThisColor = 2;
        } else if(sColor.equalsIgnoreCase("blue")){
            iThisColor = 3;
        } else if(sColor.equalsIgnoreCase("gray")){
            iThisColor = 4;
        } else if(sColor.equalsIgnoreCase("magenta")){
            iThisColor = 5;
        } else if(sColor.equalsIgnoreCase("green")){
            iThisColor = 6;
        } else if(sColor.equalsIgnoreCase("yellow")){
            iThisColor = 7;
        } else if(sColor.equalsIgnoreCase("brown")){
            iThisColor = 8;
        } else {
            System.err.println("WARNING: Didn't know the color " + sColor + " assuming black now.");
            iThisColor = 0;
        }
        this.sThisColor = sColor;
    }

    Color(final Color refColor){
        iThisColor = refColor.getThisColor();
        sThisColor = refColor.getThisColorString();
    }


    int getThisColor(){
        return iThisColor;
    }

    String getThisColorString(){
        return sThisColor;
    }
}
