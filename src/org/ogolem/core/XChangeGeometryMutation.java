/*
Copyright (c) 2014, J. M. Dieterich
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

import java.util.List;
import org.ogolem.generic.GenericMutation;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * A mutation operator exchanging particles.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class XChangeGeometryMutation implements GenericMutation<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20140327;
  public static final double DEFAULTGAUSSWIDTH = 0.4;
  public static final int SINGLEEXCHANGEMODE = 0;
  public static final int MULTIPLEEXCHANGEMODE = 1;

  private final Lottery r;
  private final double gaussWidth;
  private final int mode;

  XChangeGeometryMutation(final int mode, final double gaussWidth) {
    assert (gaussWidth > 0.0);
    assert (mode == 1 || mode == 0);
    this.gaussWidth = gaussWidth;
    this.mode = mode;
    this.r = Lottery.getInstance();
  }

  XChangeGeometryMutation(final XChangeGeometryMutation orig) {
    this.gaussWidth = orig.gaussWidth;
    this.mode = orig.mode;
    this.r = Lottery.getInstance();
  }

  @Override
  public XChangeGeometryMutation copy() {
    return new XChangeGeometryMutation(this);
  }

  @Override
  public String getMyID() {
    return "XChange mutation\n\tmode: " + mode + "\n\t Gaussian width " + gaussWidth;
  }

  @Override
  public Geometry mutate(final Geometry orig) {

    assert (orig.getNumberOfIndieParticles() >= 2);

    final Geometry mutatedGeom = new Geometry(orig);

    // find out how many molecules we will exchange
    int noToXChange;
    if (mode == SINGLEEXCHANGEMODE) {
      noToXChange = 2;
    } else if (mode == MULTIPLEEXCHANGEMODE) {
      final double rand = RandomUtils.halfgaussDouble(0.0, 1.0, gaussWidth);
      // since we want at least one pair,  choices are reduced by two
      final int noOfChoices = mutatedGeom.getNumberOfIndieParticles() - 1;
      final double oneSlot = 1. / (double) noOfChoices;
      final int choice = (int) Math.floor(rand / oneSlot);

      noToXChange = choice + 1;
    } else {
      throw new RuntimeException("Illegal mode " + mode + " in XChange operator");
    }

    /*
     * we want to know now, WHICH molecules we put into our exchange pool
     */
    final int noIndies = mutatedGeom.getNumberOfIndieParticles();

    final List<Integer> toBeXChanged = RandomUtils.rndListOfPoints(noToXChange, noIndies);

    /*
     * the actual xchanging: always pairwise, of course... ;-)
     */
    for (int i = 0; i < noToXChange; i++) {

      // get the pair
      final int start = toBeXChanged.get(i);
      final int partner = r.nextInt(noIndies);

      final double[] com1 =
          mutatedGeom.getMoleculeAtPosition(start).getExternalCenterOfMass().clone();
      final double[] euler1 = mutatedGeom.getMoleculeAtPosition(start).getOrientation().clone();

      final double[] com2 =
          mutatedGeom.getMoleculeAtPosition(partner).getExternalCenterOfMass().clone();
      final double[] euler2 = mutatedGeom.getMoleculeAtPosition(partner).getOrientation().clone();

      // xchange
      mutatedGeom.getMoleculeAtPosition(start).setExternalCenterOfMass(com2);
      mutatedGeom.getMoleculeAtPosition(start).setOrientation(euler2);

      mutatedGeom.getMoleculeAtPosition(partner).setExternalCenterOfMass(com1);
      mutatedGeom.getMoleculeAtPosition(partner).setOrientation(euler1);
    }

    return mutatedGeom;
  }
}
