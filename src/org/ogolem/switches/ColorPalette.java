/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2012, J. M. Dieterich
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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;

/**
 * Describes what colors are even available for usage.
 * @author Johannes Dieterich
 * @version 2012-06-24
 */
public final class ColorPalette implements Serializable{

    private static final long serialVersionUID = (long) 2010114;

    private final ArrayList<Color> alAvailColors;

    private static final ColorPalette palette = new ColorPalette();

    private ColorPalette(){
        this.alAvailColors = new ArrayList<>();
    }

    /**
     * Hand over a reference to the singleton.
     * @return a reference to the palette
     */
    public static ColorPalette getReference(){
        return palette;
    }

    /**
     * Set all available colors into the ColorPalette.
     * @param availColors
     */
    synchronized void initializeColors(final ArrayList<Color> availColors){
        final Iterator<Color> itAvails = availColors.iterator();

        while(itAvails.hasNext()){
            alAvailColors.add(itAvails.next());
        }
    }

    /**
     * Returns a random color of the set of allowed ones.
     * @return a random color
     */
    public Color getRandomColor(){

        if(alAvailColors.isEmpty() || alAvailColors.size() == 1){
            System.err.println("WARNING: Doesn't make much sense to have 0/1 available colors. Returning null.");
            return null;
        }

        final Random random = new Random();
        final int iWhich = random.nextInt(alAvailColors.size());
        final Color randCol = new Color(alAvailColors.get(iWhich));

        return randCol;
    }

    /**
     * Returns which colors are possible.
     * @return an array of possible colors
     */
    String[] getPossibleColors(){

        final String[] saPossCols = new String[9];
        saPossCols[0] = "white";
        saPossCols[1] = "black";
        saPossCols[2] = "red";
        saPossCols[3] = "blue";
        saPossCols[4] = "gray";
        saPossCols[5] = "magenta";
        saPossCols[6] = "green";
        saPossCols[7] = "yellow";
        saPossCols[8] = "brown";

        return saPossCols;
    }
}
