/**
Copyright (c) 2010, J. M. Dieterich
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
package org.ogolem.core;

import contrib.jama.*;
import org.ogolem.random.Lottery;

/**
 * A cuboidal space.
 * @author Johannes Dieterich
 * @version 2016-12-18
 */
final class ParallelepipedSpace implements AllowedSpace{

    private final static long serialVersionUID = (long) 20100921;

    private final Lottery random;

    /**
     * 4 corners required to define the parallelpiped.
     */
    private final double[][] corners;

    /**
     * 3 cell vectors.
     */
    private final double[][] cellVectors;

    private final Matrix invTrans;
    
    ParallelepipedSpace(final double[][] corners){
        this.random = Lottery.getInstance();
        this.corners = corners;
        this.cellVectors = new double[3][3];

        for(int j = 0; j < 3; j++){
            for(int i = 0; i < 3; i++){
                cellVectors[j][i] = this.corners[j+1][i] - this.corners[0][i];
            }
        }
        
        final Matrix matTrans = new Matrix(3,3);
        final double[][] trans = matTrans.getArray();
        trans[0][0] = cellVectors[0][0];
        trans[1][0] = cellVectors[0][1];
        trans[2][0] = cellVectors[0][2];
        trans[0][1] = cellVectors[1][0];
        trans[1][1] = cellVectors[1][1];
        trans[2][1] = cellVectors[1][2];
        trans[0][2] = cellVectors[2][0];
        trans[1][2] = cellVectors[2][1];
        trans[2][2] = cellVectors[2][2];
        
        invTrans = matTrans.inverse();
    }
    
    private ParallelepipedSpace(final ParallelepipedSpace orig){
        
        // easy as that, no deep copies needed since these are read-only fields anyway
        this.random = Lottery.getInstance();
        this.corners = orig.corners.clone();
        this.cellVectors = orig.cellVectors.clone();
        this.invTrans = orig.invTrans.copy();        
    }

    @Override
    public ParallelepipedSpace clone(){
        return new ParallelepipedSpace(this);
    }

    @Override
    public double[] getPointInSpace(){

        final double[] randVec = new double[3];
        for(int i= 0; i < 3; i++){
             randVec[i] = random.nextDouble();
        }

        final double[] point = new double[3];
        for(int i = 0; i < 3; i++){
            point[i] = randVec[0]*cellVectors[0][i]
                    + randVec[1]*cellVectors[1][i]
                    + randVec[2]*cellVectors[2][i];
        }

        // now we have a point in the "unit cell", we just need to "move it"
        for(int i = 0; i < 3; i++){
            point[i] += corners[0][i];
        }
        
        assert(isPointInSpace(point));
        
        return point;
    }

    @Override
    public boolean isPointInSpace(final double[] point){

        // move inside of the cell (perhaps)
        final double[] tmpPoint = new double[3];
        for(int i = 0; i < 3; i++){
            tmpPoint[i] = point[i] - corners[0][i];
        }

        // tranform into the coordinate system of the parallelepiped
        final Matrix matPoint = new Matrix(3,1);
        final double[][] pointArray = matPoint.getArray();
        pointArray[0][0] = tmpPoint[0];
        pointArray[1][0] = tmpPoint[1];
        pointArray[2][0] = tmpPoint[2];

        // invert the matrix to achieve the transformation matrix
        final Matrix matTransformed = invTrans.times(matPoint);
        final double[][] transformed = matTransformed.getArray();

        for(int i = 0; i < 3; i++){
            if(transformed[i][0] > 1.0){
                // outside of our new coordinate system
                return false;
            }
        }

        // if we arrive here, everything is fine
        return true;
    }
}
