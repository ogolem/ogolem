/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010, J. M. Dieterich
              2016-2020, J. M. Dieterich and B. Hartke
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

import java.io.Serializable;
import java.util.List;
import org.ogolem.helpers.Tuple;

/**
 * Defines what all potential types of environments need to be capable of.
 * @author Johannes Dieterich
 * @version 2017-01-01
 */
public interface Environment extends Cloneable, Serializable{

    public static enum ENVIRONMENTTYPE {SURFACE, CAVITY};
    
    Environment clone();

    boolean doesItFit(final CartesianCoordinates clusterCartes, final BondInfo clusterBonds);

    Tuple<CartesianCoordinates, BondInfo> marryThem(final CartesianCoordinates clusterCartes, final BondInfo clusterBonds);
    
    CartesianCoordinates marryThem(final CartesianCoordinates clusterCartes);

    CartesianCoordinates divorceThem(final CartesianCoordinates completeCartes);

    void initializeConnections(final CartesianCoordinates clusterCartes, final BondInfo clusterBonds);

    void mutateConnections(final CartesianCoordinates clusterCartes);

    List<Environment> createOffspring(final Environment father);

    List<Integer> whichSecondaryAtoms();

    boolean isEnvironmentRigid();

    double[] returnMyGenom();

    void putGenomIn(double[] genom);

    int atomsInEnv();
    
    void moveCluster(final int direction, final double move);
    
    CartesianCoordinates getEnvironmentCartes();
    
    ENVIRONMENTTYPE getEnvironmentType();
}
