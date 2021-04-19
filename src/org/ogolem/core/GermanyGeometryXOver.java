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

import java.util.ArrayList;
import java.util.List;
import org.ogolem.generic.GenericCrossover;
import org.ogolem.helpers.Tuple;
import org.ogolem.random.RandomUtils;

/**
 * Implements the most fundamental crossover: genotype. Used in the "germany" algorithms. It sucks.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public class GermanyGeometryXOver implements GenericCrossover<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20140327;
  private final double gaussWidth;

  GermanyGeometryXOver(final double gaussWidth) {
    this.gaussWidth = gaussWidth;
  }

  GermanyGeometryXOver(final GermanyGeometryXOver orig) {
    this.gaussWidth = orig.gaussWidth;
  }

  @Override
  public GermanyGeometryXOver copy() {
    return new GermanyGeometryXOver(this);
  }

  @Override
  public String getMyID() {
    return "GERMANY\ngeometry genotype crossover\n\t Gaussian width: " + gaussWidth;
  }

  @Override
  public Tuple<Geometry, Geometry> crossover(
      Geometry mother, Geometry father, final long futureID) {

    // the children geometries
    Geometry gChildOne;
    Geometry gChildTwo;

    int iNoOfIndies = mother.getNumberOfIndieParticles();
    // System.out.println("The number of Indies is reported to be " + iNoOfIndies);

    final double dRanGauss = RandomUtils.gaussDouble(-1.0, 1.0, gaussWidth);

    /*
     * now we have a Gaussian distributes random number in between -1.0 and 1.0.
     * mapping of that to the actual cutting position takes place now
     */

    // the point where the cut through the vectors will take place
    int iCut = 0;

    // the "size" of each bin
    double dBinSize = 2.0 / (double) iNoOfIndies;
    // System.out.println("GermanyGlobOpt: the binsize:" + dBinSize);
    // the initial start is at -1.0
    double dOffset = -1.0 + dBinSize;

    // loop till we find the right spot to cut
    while (dOffset <= 1.0) {
      if (dRanGauss < dOffset) {
        // the random number is in the bin
        break;
      } else {
        dOffset += dBinSize;
        iCut++;
      }
    }

    // the geometry configurations
    GeometryConfig gcChildOne = new GeometryConfig();
    GeometryConfig gcChildTwo = new GeometryConfig();

    GeometryConfig gcMother = mother.returnMyConfig();
    GeometryConfig gcFather = father.returnMyConfig();

    // set some small values
    gcChildOne.motherID = gcMother.lID;
    gcChildTwo.motherID = gcMother.lID;
    gcChildOne.fatherID = gcFather.lID;
    gcChildTwo.fatherID = gcFather.lID;
    gcChildOne.noOfParticles = iNoOfIndies;
    gcChildTwo.noOfParticles = iNoOfIndies;

    // the actual chopping of MC's...
    ArrayList<MoleculeConfig> alMCOne = new ArrayList<>(iNoOfIndies);
    ArrayList<MoleculeConfig> alMCTwo = new ArrayList<>(iNoOfIndies);

    // System.out.println("GermanyGlobOpt: the cutting number is: " + iCut);

    for (int i = 0; i < iCut; i++) {
      alMCOne.add(i, gcMother.geomMCs.get(i));
      alMCTwo.add(i, gcFather.geomMCs.get(i));
    }

    for (int i = iCut; i < iNoOfIndies; i++) {
      alMCOne.add(i, gcFather.geomMCs.get(i));
      alMCTwo.add(i, gcMother.geomMCs.get(i));
    }

    // put the MCs into the GC
    gcChildOne.geomMCs = alMCOne;
    gcChildTwo.geomMCs = alMCTwo;

    if (mother.containsEnvironment()) {
      // the environment
      final List<Environment> envChilds = gcChildOne.env.createOffspring(gcChildTwo.env);
      gcChildOne.env = envChilds.get(0);
      gcChildTwo.env = envChilds.get(1);
    }

    // create the two children
    gChildOne = new Geometry(gcChildOne);
    gChildTwo = new Geometry(gcChildTwo);

    return new Tuple<>(gChildOne, gChildTwo);
  }

  @Override
  public short hasPriority() {
    return -1;
  }
}
