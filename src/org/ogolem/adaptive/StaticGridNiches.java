/**
Copyright (c) 2010, J. M. Dieterich and B. Hartke
              2015-2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive;

import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;

/**
 * Follows a niching approach to obtain a higher diversity in the pool of parameters using a static
 * N-D grid over the search space.
 * @author Johannes Dieterich
 * @author Mark Dittner
 * @version 2020-04-29
 */
final class StaticGridNiches implements NicheComputer<Double,AdaptiveParameters>{

    private static final long serialVersionUID = (long) 20150526;
    
    private final int iNoOfNichesPerDim;
    private final double[] daBoxSize;
    private final double[][] daBoxes;
    private final boolean bUseRatios;

    StaticGridNiches(final int iNichesPerDim, final double[] minBorders, final double[] maxBorders, final boolean useRatios){
        this.iNoOfNichesPerDim = iNichesPerDim;
        this.daBoxSize = new double[minBorders.length];
        this.daBoxes = new double[minBorders.length][iNichesPerDim];
        for(int i = 0; i < daBoxSize.length; i++){
            daBoxSize[i] = Math.abs(minBorders[i] - maxBorders[i])/iNoOfNichesPerDim;

            for(int j = 0 ; j < iNichesPerDim; j++){
                if(j == 0){
                    daBoxes[i][j] = minBorders[i] + daBoxSize[i];
                    continue;
                }
                daBoxes[i][j] = daBoxSize[i] + daBoxes[i][j-1];
            }
        }
        this.bUseRatios = useRatios;
    }
    
    private StaticGridNiches(final StaticGridNiches orig){
        this.iNoOfNichesPerDim = orig.iNoOfNichesPerDim;
        this.bUseRatios = orig.bUseRatios;
        this.daBoxSize = orig.daBoxSize.clone();
        this.daBoxes = new double[orig.daBoxes.length][];
        for(int i = 0; i < daBoxes.length; i++){
            this.daBoxes[i] = orig.daBoxes[i].clone();
        }
    }
    
    @Override
    public StaticGridNiches copy(){
        return new StaticGridNiches(this);
    }

    @Override
    public Niche computeNiche(AdaptiveParameters params){

        final double[] daParams = params.getAllParamters();
        final int[] iaBins = new int[daParams.length];

        DimsLoop:
        for(int i = 0; i < daParams.length; i++){

            BoxesLoop:
            for(int j = 0; j < daBoxes[i].length; j++){
                if(daParams[i] < daBoxes[i][j]){
                    // found the box
                    iaBins[i] = j;

                    break BoxesLoop;
                }

                // might be (due to locopt) bigger than the biggest bin
                if(j == daBoxes[i].length -1){
                    iaBins[i] = j;
                    break BoxesLoop;
                }
            }
        }

        String sID = "";
        if(bUseRatios){
            final int[] iaRatios = new int[iNoOfNichesPerDim];
            for(int i = 0; i < iaBins.length; i++){
                iaRatios[iaBins[i]]++;
            }

            // assemble the niche ID
            for(int i = 0; i < iaRatios.length; i++){
                sID += "r" + i + "_" + iaRatios[i];
            }

        } else{
            // assemble the niche ID
            for(int i = 0; i < iaBins.length; i++){
                sID += "b" + iaBins[i];
            }
        }

        // create the niche object
        final Niche niche = new Niche(sID);

        return niche;
    }
}
