/*
Copyright (c) 2016-2021, J. M. Dieterich and B. Hartke
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

import static org.junit.jupiter.api.Assertions.*;
import static org.ogolem.core.Constants.ANGTOBOHR;
import static org.ogolem.core.Constants.EVTOHARTREE;

import org.junit.jupiter.api.Test;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Gradient;
import org.ogolem.core.SimpleBondInfo;
import org.ogolem.core.Topology;

/**
 * @author Johannes Dieterich
 * @version 2021-07-17
 */
public class AdaptiveSWG2BodyTermTest {

  /** Test of partialInteraction method, of class AdaptiveSWG3BodyTerm. */
  @Test
  public void testPartialInteraction1NC() {

    System.out.println("partialInteraction topo 1 no caching");
    final Topology topology = setupTopo1();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(false, 0.2);
    final double expResult = 0.579173;
    final double result = instance.partialInteraction(topology, params);
    assertEquals(expResult, result, 1e-6);
  }

  @Test
  public void testPartialInteraction1WC() {

    System.out.println("partialInteraction topo 1 with caching");
    final Topology topology = setupTopo1();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(true, 0.2);
    final double expResult = 0.579173;
    final double result = instance.partialInteraction(topology, params);
    assertEquals(expResult, result, 1e-6);
  }

  @Test
  public void testPartialInteraction2NC() {

    System.out.println("partialInteraction topo 2 no caching");
    final Topology topology = setupTopo2();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(false, 0.2);
    final double expResult = 0.266691;
    final double result = instance.partialInteraction(topology, params);
    assertEquals(expResult, result, 1e-6);
  }

  @Test
  public void testPartialInteraction2WC() {

    System.out.println("partialInteraction topo 2 with caching");
    final Topology topology = setupTopo2();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(true, 0.2);
    final double expResult = 0.266691;
    final double result = instance.partialInteraction(topology, params);
    assertEquals(expResult, result, 1e-6);
  }

  @Test
  public void testPartialInteraction3NC() {

    System.out.println("partialInteraction topo 3 no caching");
    final Topology topology = setupTopo3();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(false, 0.2);
    final double expResult = 1.967588;
    final double result = instance.partialInteraction(topology, params);
    assertEquals(expResult, result, 1e-6);
  }

  @Test
  public void testPartialInteraction3WC() {

    System.out.println("partialInteraction topo 3 with caching");
    final Topology topology = setupTopo3();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(true, 0.2);
    final double expResult = 1.967588;
    final double result = instance.partialInteraction(topology, params);
    assertEquals(expResult, result, 1e-6);
  }

  @Test
  public void testPartialCartesianGradient1NC() {

    System.out.println("partialCartesianGradient topo 1 no caching");
    final Topology topology = setupTopo1();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(false, 0.2);
    final double expResult = 0.579173;
    final Gradient g = instance.partialCartesianGradient(topology, params);
    final double result = g.getTotalEnergy();
    assertEquals(expResult, result, 1e-6);

    final Gradient numGrad = NumericalGradients.calculateGrad(params, instance, topology);

    final boolean same = g.compare(numGrad, 1e-6, true);
    assertTrue(same);
  }

  @Test
  public void testPartialCartesianGradient1WC() {

    System.out.println("partialCartesianGradient topo 1 with caching");
    final Topology topology = setupTopo1();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(true, 0.2);
    final double expResult = 0.579173;
    final Gradient g = instance.partialCartesianGradient(topology, params);
    final double result = g.getTotalEnergy();
    assertEquals(expResult, result, 1e-6);

    final Gradient numGrad = NumericalGradients.calculateGrad(params, instance, topology);

    final boolean same = g.compare(numGrad, 1e-6, true);
    assertTrue(same);
  }

  @Test
  public void testPartialCartesianGradient2NC() {

    System.out.println("partialCartesianGradient topo 2 no caching");
    final Topology topology = setupTopo2();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(false, 0.2);
    final double expResult = 0.266691;
    final Gradient g = instance.partialCartesianGradient(topology, params);
    final double result = g.getTotalEnergy();
    assertEquals(expResult, result, 1e-6);

    final Gradient numGrad = NumericalGradients.calculateGrad(params, instance, topology);

    final boolean same = g.compare(numGrad, 1e-6, true);
    assertTrue(same);
  }

  @Test
  public void testPartialCartesianGradient2WC() {

    System.out.println("partialCartesianGradient topo 2 with caching");
    final Topology topology = setupTopo2();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(true, 0.2);
    final double expResult = 0.266691;
    final Gradient g = instance.partialCartesianGradient(topology, params);
    final double result = g.getTotalEnergy();
    assertEquals(expResult, result, 1e-6);

    final Gradient numGrad = NumericalGradients.calculateGrad(params, instance, topology);

    final boolean same = g.compare(numGrad, 1e-6, true);
    assertTrue(same);
  }

  @Test
  public void testPartialCartesianGradient3NC() {

    System.out.println("partialCartesianGradient topo 3 no caching");
    final Topology topology = setupTopo3();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(false, 0.2);
    final double expResult = 1.967588;
    final Gradient g = instance.partialCartesianGradient(topology, params);
    final double result = g.getTotalEnergy();
    assertEquals(expResult, result, 1e-6);

    final Gradient numGrad = NumericalGradients.calculateGrad(params, instance, topology);

    final boolean same = g.compare(numGrad, 1e-6, true);
    assertTrue(same);
  }

  @Test
  public void testPartialCartesianGradient3WC() {

    System.out.println("partialCartesianGradient topo 3 with caching");
    final Topology topology = setupTopo3();
    final AdaptiveParameters params = setupParam1();
    final AdaptiveSWG2BodyTerm instance = new AdaptiveSWG2BodyTerm(true, 0.2);
    final double expResult = 1.967588;
    final Gradient g = instance.partialCartesianGradient(topology, params);
    final double result = g.getTotalEnergy();
    assertEquals(expResult, result, 1e-6);

    final Gradient numGrad = NumericalGradients.calculateGrad(params, instance, topology);

    final boolean same = g.compare(numGrad, 1e-6, true);
    assertTrue(same);
  }

  private Topology setupTopo1() {

    final CartesianCoordinates cartes = new CartesianCoordinates(3, 3, new int[] {1, 1, 1});
    final String[] atoms = cartes.getAllAtomTypes();
    atoms[0] = "Si";
    atoms[1] = "Si";
    atoms[2] = "Si";
    cartes.recalcAtomNumbers();

    final double[][] xyz = cartes.getAllXYZCoord();
    xyz[0][0] = 0.0;
    xyz[1][0] = 0.0;
    xyz[2][0] = 1.5 * ANGTOBOHR;
    xyz[0][1] = 0.0;
    xyz[1][1] = 0.0;
    xyz[2][1] = 0.0;
    xyz[0][2] = 0.0;
    xyz[1][2] = 0.0;
    xyz[2][2] = -1.5 * ANGTOBOHR;

    final BondInfo bonds = new SimpleBondInfo(3);

    final Topology topo = new Topology(cartes, bonds);

    return topo;
  }

  private Topology setupTopo2() {

    final CartesianCoordinates cartes = new CartesianCoordinates(3, 3, new int[] {1, 1, 1});
    final String[] atoms = cartes.getAllAtomTypes();
    atoms[0] = "Si";
    atoms[1] = "Si";
    atoms[2] = "Si";
    cartes.recalcAtomNumbers();

    final double[][] xyz = cartes.getAllXYZCoord();
    xyz[0][0] = 0.0;
    xyz[1][0] = 0.0;
    xyz[2][0] = 1.6 * ANGTOBOHR;
    xyz[0][1] = 0.0;
    xyz[1][1] = 0.0;
    xyz[2][1] = 0.0;
    xyz[0][2] = 0.0;
    xyz[1][2] = 1.6 * ANGTOBOHR;
    xyz[2][2] = 0.0;

    final BondInfo bonds = new SimpleBondInfo(3);

    final Topology topo = new Topology(cartes, bonds);

    return topo;
  }

  private Topology setupTopo3() {

    final CartesianCoordinates cartes = new CartesianCoordinates(3, 3, new int[] {1, 1, 1});
    final String[] atoms = cartes.getAllAtomTypes();
    atoms[0] = "Si";
    atoms[1] = "Si";
    atoms[2] = "Si";
    cartes.recalcAtomNumbers();

    final double[][] xyz = cartes.getAllXYZCoord();
    xyz[0][0] = 0.0;
    xyz[1][0] = 0.0;
    xyz[2][0] = 1.5 * ANGTOBOHR;
    xyz[0][1] = 0.0;
    xyz[1][1] = 0.0;
    xyz[2][1] = 0.0;
    xyz[0][2] = 1.0 * ANGTOBOHR;
    xyz[1][2] = 0.0;
    xyz[2][2] = 2.0 * ANGTOBOHR;

    final BondInfo bonds = new SimpleBondInfo(3);

    final Topology topo = new Topology(cartes, bonds);

    return topo;
  }

  private AdaptiveParameters setupParam1() {

    final int numberOfParameters = 4;
    final int id = -1;

    final String[] keys = new String[1];
    keys[0] = "adaptiveswg2b:SiSi";

    final int[] noOfParamsPerKey = new int[] {4};

    final String paramsForMethod = "ADAPTIVESWGTESTONLY";

    final AdaptiveParameters params =
        new AdaptiveParameters(numberOfParameters, id, keys, noOfParamsPerKey, paramsForMethod);

    final double[] p = params.getAllParamters();
    p[0] = 3.77 * ANGTOBOHR;
    p[1] = 16.3 * EVTOHARTREE;
    p[2] = 11.581 * ANGTOBOHR * ANGTOBOHR * ANGTOBOHR * ANGTOBOHR;
    p[3] = 2.095 * ANGTOBOHR;

    return params;
  }
}
