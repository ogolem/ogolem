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

import java.util.Random;

/**
 * Defines a sperical space.
 * @author Johannes Dieterich
 * @version 2010-09-22
 */
final class SphereSpace implements AllowedSpace{

    private final static long serialVersionUID = (long) 20100922;

    private final Random random;

    private final double dRadius;

    private final double[] daMiddle;

    SphereSpace(double[] middle, double radius){
        this.daMiddle = middle;
        this.dRadius = radius;
        this.random = new Random();
    }

    @Override
    public SphereSpace clone(){
        double[] daTemp = daMiddle.clone();

        return new SphereSpace(daTemp, this.dRadius);
    }

    @Override
    public double[] getPointInSpace(){

        // we need some randoms
        double[][] daSpherical = new double[3][1];
        // radius
        daSpherical[0][0] = dRadius * random.nextDouble();

        // phi and omega
        daSpherical[1][0] = random.nextDouble() * 2.0 * Math.PI;
        daSpherical[2][0] = random.nextDouble() * Math.PI;

        // translate to cartesian
        double[][] daCartes = CoordTranslation.sphericalToCartesianCoord(daSpherical);

        // move with respect to middle
        double[] daPoint = new double[3];
        for(int i = 0; i < 3; i++){
            daPoint[i] = daCartes[i][0] + daMiddle[i];
        }

        return daPoint;

    }

    @Override
    public boolean isPointInSpace(double[] point){


        // move with respect to middle of sphere
        final double[][] daTemp = new double[3][1];
        for(int i = 0; i < 3; i++){
            daTemp[i][0] = point[i] - daMiddle[i];
        }

        // translate to spherical coords
        final double[][] daSpherical = CoordTranslation.cartesianToSphericalCoord(daTemp);

        // check
        return (daSpherical[0][0] <= dRadius);
    }

}
