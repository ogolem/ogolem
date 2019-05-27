/**
Copyright (c) 2010-2012, J. M. Dieterich
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

import java.io.Serializable;
import java.util.ArrayList;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.Topology;
import org.ogolem.helpers.Tuple3D;

/**
 * Defines an adaptive interaction term.
 * @author Johannes Dieterich
 * @version 2012-03-08
 */
public interface AdaptiveInteractionTerm extends Serializable, Cloneable {

    AdaptiveInteractionTerm clone();

    double partialInteraction(Topology topology, AdaptiveParameters params);

    double partialParamGradient(Topology topology, AdaptiveParameters params, double[] grad);

    Gradient partialCartesianGradient(Topology topo, AdaptiveParameters params);

    /**
     * Reports which parameters are required for this set of geometries.
     * @param cartesians
     * @param topologies List of topologies in the reference data. Size may be (depending on the implementation) different then the one of the list of cartesians.
     * @param sMethod
     * @return A 3D-Tupel, which parameters are required (names, number per key, total number).
     */
    Tuple3D<String[],int[],Integer> requiredParams(ArrayList<CartesianCoordinates> cartesians, ArrayList<Topology> topologies, String sMethod);

    double[][] bordersForMyParams(AdaptiveParameters params);
}
