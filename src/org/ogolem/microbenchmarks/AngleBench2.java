/**
Copyright (c) 2019, J. M. Dieterich and B. Hartke
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
package org.ogolem.microbenchmarks;

import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CoordTranslation;

/**
 *Benchmark the speed of computing an angle from three vectors of Cartesian coordinates
 * @author Johannes Dieterich
 * @version 2019-12-29
 */
class AngleBench2 implements SingleMicroBenchmark {
    
    private final double[] coord1 = new double[3];
    private final double[] coord2 = new double[3];
    private final double[] coord3 = new double[3];
    
    AngleBench2(){
        final CartesianCoordinates kana = CartesianCoordinatesLibrary.getKanamycinAMP2Opt();
        final double[][] xyz = kana.getAllXYZCoordsCopy();
        
        coord1[0] = xyz[0][10];
        coord1[1] = xyz[1][10];
        coord1[2] = xyz[2][10];
        
        coord2[0] = xyz[0][30];
        coord2[1] = xyz[1][30];
        coord2[2] = xyz[2][30];
        
        coord3[0] = xyz[0][68];
        coord3[1] = xyz[1][68];
        coord3[2] = xyz[2][68];
    }
    
    @Override
    public double runSingle() throws Exception {
        return CoordTranslation.calcAngle(coord1, coord2, coord3);
    }

    @Override
    public String name() {
        return "angle benchmark II";
    }
}
