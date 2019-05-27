/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
              2015-2016, J. M. Dieterich and B. Hartke
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
import java.util.List;
import org.ogolem.adaptive.genericfitness.GenericReferencePoint;
import org.ogolem.adaptive.genericfitness.PropertyCalculator;
import org.ogolem.adaptive.genericfitness.ReferenceInputData;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.properties.Property;

/**
 * Kind of atomics for ALL kinds of reparametrizable methods.
 * @author Johannes Dieterich
 * @version 2016-07-16
 */
public interface Adaptivable extends Cloneable, Serializable {
    
    Adaptivable clone();

    @Deprecated
    double energyOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, int geomID, final BondInfo bonds);

    @Deprecated
    double gradientOfStructWithParams(final CartesianCoordinates cartes,
            final AdaptiveParameters params, int geomID, final BondInfo bonds,
            final double[] grad);

    double[][] minMaxBordersForParams(final AdaptiveParameters params);

    AdaptiveParameters createInitialParameterStub(final ArrayList<CartesianCoordinates> refCartes,
            final String sMethod);

    <T extends Property, V extends ReferenceInputData<T>> PropertyCalculator<T,V> getCalculatorForProperty(final T property, final V data);
    
    List<? extends Property> runAllPropertyCalcs(final AdaptiveParameters params, final List<GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>>> referencePoints);
}
