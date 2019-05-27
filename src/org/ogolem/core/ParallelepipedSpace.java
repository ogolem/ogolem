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

import contrib.jama.*;
import java.util.Random;

/**
 * A cuboidal space.
 * @author Johannes Dieterich
 * @version 2010-09-22
 */
final class ParallelepipedSpace implements AllowedSpace{

    private final static long serialVersionUID = (long) 20100921;

    private final Random random;

    /**
     * 4 corners required to define the parallelpiped.
     */
    private final double[][] daCorners;

    /**
     * 3 cell vectors.
     */
    private final double[][] daCellVectors;

    ParallelepipedSpace(final double[][] corners){
        this.random = new Random();
        this.daCorners = corners;
        this.daCellVectors = new double[3][3];

        for(int j = 0; j < 3; j++){
            for(int i = 0; i < 3; i++){
                daCellVectors[j][i] = daCorners[j+1][i] - daCorners[0][i];
            }
        }

    }

    @Override
    public ParallelepipedSpace clone(){

        double[][] daCornersClone = new double[4][3];
        for(int i = 0; i < 4; i++){
            daCornersClone[i] = this.daCorners[i].clone();
        }

        return new ParallelepipedSpace(daCornersClone);
    }

    @Override
    public double[] getPointInSpace(){

        double[] daRand = new double[3];
        for(int i= 0; i < 3; i++){
            daRand[i] = random.nextDouble();
        }

        double[] daPoint = new double[3];
        for(int i = 0; i < 3; i++){
            daPoint[i] = daRand[0]*daCellVectors[0][i]
                    + daRand[1]*daCellVectors[1][i]
                    + daRand[2]*daCellVectors[2][i];
        }

        // now we have a point in the "unit cell", we just need to "move it"
        for(int i = 0; i < 3; i++){
            daPoint[i] += daCorners[0][i];
        }

        return daPoint;
    }

    @Override
    public boolean isPointInSpace(double[] point){

        // move inside of the cell (perhaps)
        double[] daTempPoint = new double[3];
        for(int i = 0; i < 3; i++){
            daTempPoint[i] = point[i] - daCorners[0][i];
        }

        // tranform into the coordinate system of the parallelepiped
        Matrix matTrans = new Matrix(3,3);
        double[][] daTrans = matTrans.getArray();
        daTrans[0][0] = daCellVectors[0][0];
        daTrans[1][0] = daCellVectors[0][1];
        daTrans[2][0] = daCellVectors[0][2];
        daTrans[0][1] = daCellVectors[1][0];
        daTrans[1][1] = daCellVectors[1][1];
        daTrans[2][1] = daCellVectors[1][2];
        daTrans[0][2] = daCellVectors[2][0];
        daTrans[1][2] = daCellVectors[2][1];
        daTrans[2][2] = daCellVectors[2][2];

        Matrix matPoint = new Matrix(3,1);
        double[][] daPoint = matPoint.getArray();
        daPoint[0][0] = daTempPoint[0];
        daPoint[1][0] = daTempPoint[1];
        daPoint[2][0] = daTempPoint[2];

        // invert the matrix to achieve the transformation matrix
        Matrix matTransformed = matTrans.inverse().times(matPoint);
        double[][] daTransformed = matTransformed.getArray();

        for(int i = 0; i < 3; i++){
            if(daTransformed[i][0] > 1.0){
                // outside of our new coordinate system
                return false;
            }
        }

        // if we arrive here, everything is fine
        return true;
    }
}
