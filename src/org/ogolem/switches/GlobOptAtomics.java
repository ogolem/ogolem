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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;

/**
 * Provides atomic functions for the global optimization.
 * @author Johannes Dieterich
 * @version 2012-06-12
 */
final class GlobOptAtomics {

    static ArrayList<ArrayList<Color>> genotypeCross(final ArrayList<Color> alMother, final ArrayList<Color> alFather){

        final Random random = new Random();
        final int iLength = alFather.size();

        final ArrayList<Color> alChildOne = new ArrayList<>(iLength);
        final ArrayList<Color> alChildTwo = new ArrayList<>(iLength);

        // first we want the crossover point
        final int iCrossPos = random.nextInt(iLength);

        // now we cross the colors
        for(int i = 0; i < iCrossPos; i++){
            alChildOne.add(alFather.get(i));
            alChildTwo.add(alMother.get(i));
        }

        for(int i = iCrossPos; i < iLength; i++){
            alChildOne.add(alMother.get(i));
            alChildTwo.add(alFather.get(i));
        }

        final ArrayList<ArrayList<Color>> alResult = new ArrayList<>(2);
        alResult.add(alChildOne);
        alResult.add(alChildTwo);

        return alResult;
    }

    static ArrayList<Color> genotypeMutation(final ArrayList<Color> alStart, final boolean bMoreMutation){

        final Random random = new Random();

        // first figure out whether we even want to mutate
        int iRandom = random.nextInt(20);

        // probability is 5.0 percent, so any number is as good as another one
        if (iRandom == 1) {
            // do mutation
            final ArrayList<Color> alEnd = new ArrayList<>(alStart.size());

            final ColorPalette palette = ColorPalette.getReference();

            // copy all colors
            final Iterator<Color> itColors = alStart.iterator();
            while(itColors.hasNext()){
                alEnd.add(new Color(itColors.next()));
            }

            // now we figure out where we want to mutate
            int iPosition = random.nextInt(alEnd.size());

            // and then we just mutate there (and after if wanted)
            if(!bMoreMutation){
                // really just at this spot
                final Color randCol = palette.getRandomColor();
                alEnd.set(iPosition, randCol);
            } else {
                // we mutate from this spot on
                for(int i = iPosition; i < alEnd.size(); i++){
                    final Color randCol = palette.getRandomColor();
                    alEnd.set(i, randCol);
                }
            }

            return alEnd;
        } else {
            // no mutation needed
            return alStart;
        }
    }
}
