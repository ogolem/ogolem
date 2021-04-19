/*
Copyright (c) 2013-2014, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Hartke
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

import static org.ogolem.core.GlobOptAtomics.randomRotation;

import java.util.List;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;

/**
 * A wrapper around the merging xover.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class IcelandGeometryXOver implements GenericCrossover<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20130327;
  private static final boolean DEBUG = false;

  private final MergingPhenoXOver merger;

  IcelandGeometryXOver(
      final CartesianFullBackend back,
      final double blowColl,
      final double blowDiss,
      final CollisionDetectionEngine colldetect,
      final MergingPhenoXOver.MergingPhenoConfig config,
      final boolean doCDDD,
      final boolean doRandomTrans) {
    this.merger =
        new MergingPhenoXOver(back, blowColl, blowDiss, colldetect, config, doCDDD, doRandomTrans);
  }

  IcelandGeometryXOver(final IcelandGeometryXOver orig) {
    this.merger = orig.merger.copy();
  }

  @Override
  public IcelandGeometryXOver copy() {
    return new IcelandGeometryXOver(this);
  }

  @Override
  public String getMyID() {
    return "ICELAND GEOMETRY GLOBOPT\n\tmerger: " + merger.getMyID();
  }

  @Override
  public Tuple<Geometry, Geometry> crossover(
      final Geometry mother, final Geometry father, final long futureID) {

    if (mother.getNumberOfIndieParticles() == 1) {

      final Geometry gChild1 = new Geometry(mother);
      final Geometry gChild2 = new Geometry(father);

      return new Tuple<>(gChild1, gChild2);
    }

    // create the children
    Geometry gChild1 = new Geometry(father);
    gChild1.setID(futureID);
    Geometry gChild2 = new Geometry(mother);
    gChild2.setID(futureID);

    final CartesianCoordinates cartesMother = gChild2.getCartesians();
    final CartesianCoordinates cartesFather = gChild1.getCartesians();

    cartesMother.moveCoordsToCOM();
    cartesFather.moveCoordsToCOM();

    if (DEBUG) {
      final String[] printCart1 = cartesMother.createPrintableCartesians();
      final String[] printCart2 = cartesFather.createPrintableCartesians();
      System.out.println("DEBUG: Cartesians prior to rotation: mother");
      for (final String s : printCart1) {
        System.out.println(s);
      }
      System.out.println("DEBUG: Cartesians prior to rotation: father");
      for (final String s : printCart2) {
        System.out.println(s);
      }
    }

    // rotate the two
    randomRotation(cartesMother);
    randomRotation(cartesFather);

    if (DEBUG) {
      final String[] printCart1 = cartesMother.createPrintableCartesians();
      final String[] printCart2 = cartesFather.createPrintableCartesians();
      System.out.println("DEBUG: Cartesians after rotation: mother");
      for (final String s : printCart1) {
        System.out.println(s);
      }
      System.out.println("DEBUG: Cartesians after rotation: father");
      for (final String s : printCart2) {
        System.out.println(s);
      }
    }

    // translate it to geometries again
    CoordTranslation.updateGeometryFromCartesian(cartesFather, gChild1);
    CoordTranslation.updateGeometryFromCartesian(cartesMother, gChild2);

    // set the IDs of mother and father
    gChild1.setMother(mother.getID());
    gChild1.setFather(father.getID());
    gChild2.setMother(mother.getID());
    gChild2.setFather(father.getID());

    gChild1.setID(futureID);
    gChild2.setID(futureID);

    // in case of environments, add it in at this point
    if (mother.containsEnvironment()) {
      // the environment
      final List<Environment> envChilds =
          mother.getEnvironment().createOffspring(father.getEnvironment());
      gChild1.setEnvironment(envChilds.get(0));
      gChild2.setEnvironment(envChilds.get(1));
    }

    return merger.crossover(gChild1, gChild2, futureID);
  }

  @Override
  public short hasPriority() {
    return -1;
  }
}
