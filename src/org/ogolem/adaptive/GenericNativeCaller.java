/*
Copyright (c) 2018-2020, J. M. Dieterich and B. Hartke
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

import java.util.ArrayList;
import java.util.List;
import org.ogolem.adaptive.genericfitness.GenericReferencePoint;
import org.ogolem.adaptive.genericfitness.PropertyCalculator;
import org.ogolem.adaptive.genericfitness.ReferenceGenericScalarData;
import org.ogolem.adaptive.genericfitness.ReferenceInputData;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.io.InputPrimitives;
import org.ogolem.properties.GenericScalarProperty;
import org.ogolem.properties.Property;

/**
 * An interface to call a native library for evaluations.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class GenericNativeCaller implements Adaptivable {

  private static final long serialVersionUID = (long) 20180617;

  private final int noParams;
  private final double[][] borders;
  private long handle;

  private static native long initialize(final int noParams);

  private static native double fitnessScalar(
      final long handle,
      final int noParams,
      final double[] parameters,
      final int typeID,
      final long pointID);

  /*
   * Loading of external library libNOPFL.so (for UNIX) with a JNI interface
   */
  static {
    System.loadLibrary("NOPFL"); // native OGOLEM parameter fitting library
  }

  GenericNativeCaller(final int noParams) throws Exception {
    assert (noParams > 0);
    this.noParams = noParams;
    this.borders = new double[2][noParams];

    final String[] dat = InputPrimitives.readFileIn("genericnative.borders");
    for (int i = 0; i < noParams; i++) {
      final String[] line = dat[i].trim().split("\\s+");
      borders[0][i] = Double.parseDouble(line[0]);
      borders[1][i] = Double.parseDouble(line[1]);
    }

    handle = initialize(noParams);
  }

  private GenericNativeCaller(final GenericNativeCaller orig) {
    this.noParams = orig.noParams;
    this.borders = orig.borders;
  }

  @Override
  public GenericNativeCaller copy() {
    return new GenericNativeCaller(this);
  }

  @Override
  public double energyOfStructWithParams(
      CartesianCoordinates cartes, AdaptiveParameters params, int geomID, BondInfo bonds) {
    throw new UnsupportedOperationException(
        "Not supported. Don't use generic native caller for this!");
  }

  @Override
  public double gradientOfStructWithParams(
      CartesianCoordinates cartes,
      AdaptiveParameters params,
      int geomID,
      BondInfo bonds,
      double[] grad) {
    throw new UnsupportedOperationException(
        "Not supported. Don't use generic native caller for this!");
  }

  @Override
  public double[][] minMaxBordersForParams(AdaptiveParameters params) {
    return borders;
  }

  @Override
  public AdaptiveParameters createInitialParameterStub(
      ArrayList<CartesianCoordinates> refCartes, String sMethod) {

    final String[] atoms = new String[] {"hartkenium"};
    final int[] paramsPerAt = new int[] {noParams};
    final AdaptiveParameters paramStub =
        new AdaptiveParameters(noParams, -1, atoms, paramsPerAt, sMethod);

    // add parameter descriptions
    final String[] descr = paramStub.getAllDescriptions();
    for (int i = 0; i < descr.length; i++) {
      descr[i] = "generic native parameter " + i;
    }

    return paramStub;
  }

  @Override
  public <T extends Property, V extends ReferenceInputData<T>>
      PropertyCalculator<T, V> getCalculatorForProperty(T property, V data) {

    // instantiate nested classes
    if (property instanceof GenericScalarProperty) {
      if (data instanceof ReferenceGenericScalarData) {
        return (PropertyCalculator<T, V>) new ScalarCalculator();
      }
    }

    System.err.println(
        "No calculator for property "
            + property.name()
            + " with data "
            + data.getClass().getName()
            + " in GenericNativeCaller.");

    return null;
  }

  @Override
  public List<? extends Property> runAllPropertyCalcs(
      AdaptiveParameters params,
      List<GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>>>
          referencePoints) {

    final List<Property> allProps = new ArrayList<>();
    for (final GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>> refPoint :
        referencePoints) {
      final Property p = refPoint.getReferenceProperty();
      if (p instanceof GenericScalarProperty) {
        final GenericScalarProperty ref = (GenericScalarProperty) p;
        final ReferenceGenericScalarData dat =
            (ReferenceGenericScalarData) refPoint.getReferenceInputData();
        final PropertyCalculator<GenericScalarProperty, ReferenceGenericScalarData> calc =
            getCalculatorForProperty(ref, dat);

        final GenericScalarProperty prop = calc.calculateProperty(params, dat);
        allProps.add(prop);
      } else {
        // error
        throw new RuntimeException("Unknown property type " + p.name() + ".");
      }
    }

    return allProps;
  }

  class ScalarCalculator
      implements PropertyCalculator<GenericScalarProperty, ReferenceGenericScalarData> {

    private static final long serialVersionUID = (long) 20180123;

    ScalarCalculator() {}

    private ScalarCalculator(final ScalarCalculator orig) {}

    @Override
    public ScalarCalculator copy() {
      return new ScalarCalculator(this);
    }

    @Override
    public GenericScalarProperty calculateProperty(
        AdaptiveParameters p, ReferenceGenericScalarData data) {

      final long pointID = data.getPointID();
      final int typeID = data.getTypeID();
      final double[] params = p.getAllParamters();

      final double scalar = fitnessScalar(handle, noParams, params, typeID, pointID);

      return new GenericScalarProperty(scalar, typeID);
    }

    @Override
    public GenericScalarProperty calculatePropertyGradient(
        AdaptiveParameters p, ReferenceGenericScalarData data, double[] grad) {
      return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
    }
  }
}
