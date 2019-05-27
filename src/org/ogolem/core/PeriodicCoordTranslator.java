/**
Copyright (c) 2013, J. M. Dieterich
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

import contrib.jama.Matrix;

/**
 * Ensures coordinates to be in the same periodic cell w/ the same offset.
 * Also (in the future!) might take care of carrying out symmetry operations.
 * @author Johannes Dieterich
 * @version 2013-01-07
 */
public class PeriodicCoordTranslator {
    
    /*
     * Best for symmetry: interface Symmetrizer and dependency injection here?
     */
    
    private static final boolean DEBUG = false;
    private static final double MOVEINCR = 10000.0;
    
    /**
     * Translates cartesian coordinates into fractional periodic ones given a cell.
     * @param cartes A cartesian coordinates object.
     * @param cell the periodic cell.
     * @param anker the position of our anker (which is the atom with the "smallest" set of coordinates).
     * @param correctMultiples if all periodic coordinates should be scaled to be in the interval [0,1].
     * @return Fractional periodic coordinates (may be a multiple of 1.0).
     */
    public static double[][] coordsToPeriodic(final CartesianCoordinates cartes, final double[][] cell, final double[] anker,
            final boolean correctMultiples){
        
        // for the time being no symmetry so the number of atoms is correct
        final int noAtomsPer = cartes.getNoOfAtoms();
        final double[][] xyz = cartes.getAllXYZCoordsCopy();
        
        // first find the atom with the smallest coordinates
        // as asmall hack: move all atoms somewhere really positive (+MOVEINCR)
        // and check which atom then has the smallest distance to 0/0/0
        int minAtomID = -1;
        double minDistSq = Double.MAX_VALUE;
        for(int i = 0; i < noAtomsPer; i++){
            final double x = xyz[0][i]+MOVEINCR;
            final double y = xyz[1][i]+MOVEINCR;
            final double z = xyz[2][i]+MOVEINCR;
            final double distSq = x*x+y*y+z*z;
            if(distSq < minDistSq){
                minAtomID = i;
                minDistSq = distSq;
            }
        }
        
        // now translate all atoms so the found atom is on the anker position
        final double[] trans2Ank = new double[3];
        trans2Ank[0] = xyz[0][minAtomID]-anker[0];
        trans2Ank[1] = xyz[1][minAtomID]-anker[1];
        trans2Ank[2] = xyz[2][minAtomID]-anker[2];
        for(int c = 0; c < 3; c++){
            for(int i = 0; i < noAtomsPer; i++){
                xyz[c][i] -= trans2Ank[c];
            }
        }
        
        if(DEBUG){
            for(int i = 0; i < 3; i++){
                final double dX = xyz[0][minAtomID]-anker[0];
                final double dY = xyz[1][minAtomID]-anker[1];
                final double dZ = xyz[2][minAtomID]-anker[2];
                final double distSq = dX*dX + dY*dY + dZ*dZ;
                System.out.println("DEBUG: After translation anker atom " + minAtomID + " is " + Math.sqrt(distSq) + " away from anker position.");
            }
        }
        
        // tranform into the coordinate system of the parallelepiped
        final Matrix matTrans = new Matrix(3,3);
        final double[][] transform = matTrans.getArray();
        transform[0][0] = cell[0][0];
        transform[1][0] = cell[0][1];
        transform[2][0] = cell[0][2];
        transform[0][1] = cell[1][0];
        transform[1][1] = cell[1][1];
        transform[2][1] = cell[1][2];
        transform[0][2] = cell[2][0];
        transform[1][2] = cell[2][1];
        transform[2][2] = cell[2][2];

        final Matrix matPoint = new Matrix(xyz);

        // invert the matrix to get the transformation matrix
        final Matrix matTransformed = matTrans.inverse().times(matPoint);
        // get the fractional coordinates
        final double[][] periodic = matTransformed.getArray();
        
        if(correctMultiples){
            for(int c = 0; c < 3; c++){
                for(int at = 0; at < noAtomsPer; at++){
                    periodic[c][at] %= 1.0;
                }
            }
        }
        
        return periodic;
    }
    
    public static double[][] coordsFromPeriodic(final double[][] periodic, final double[][] cell){
        
        // again we assume that the number of atoms is constant (aka no symmetry)
        final int noAtoms = periodic[0].length;
        final double[][] xyz = new double[3][noAtoms];
        
        for(int c = 0; c < 3; c++){
            for(int at = 0; at < noAtoms; at++){
                xyz[c][at] = periodic[0][at]*cell[0][c]
                           + periodic[1][at]*cell[1][c]
                           + periodic[2][at]*cell[2][c];
            }
        }
        
        return xyz;
    }
}
