/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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

/**
 * A decorator for the geometry initialization routines.
 *
 * @author Johannes Dieterich
 * @version 2020-08-09
 */
public final class GeometryInit implements GeometryInitialization {

  private static final long serialVersionUID = (long) 20160402;

  public static enum INITSTYLE {
    PACKASCENDING,
    PACKBYSIZE,
    PACKRANDOMLY,
    LAYERPACKINGASCENDING,
    LAYERPACKINGBYSIZE,
    LAYERPACKINGRANDOMLY,
    RANDOMWITHDD,
    RANDOMWITHOUTDD
  };

  public static final INITSTYLE DEFAULTINIT = INITSTYLE.PACKASCENDING;

  private final GeometryInitialization init;

  public GeometryInit(final INITSTYLE whichInit) {

    switch (whichInit) {
      case PACKBYSIZE:
        init =
            new PackingInit(
                NorwayGeometryMutation.MUTMODE.BYSIZE, NorwayGeometryMutation.PACKDIM.THREED);
        break;
      case PACKRANDOMLY:
        init =
            new PackingInit(
                NorwayGeometryMutation.MUTMODE.RANDOM, NorwayGeometryMutation.PACKDIM.THREED);
        break;
      case PACKASCENDING:
        init =
            new PackingInit(
                NorwayGeometryMutation.MUTMODE.ASCENDING, NorwayGeometryMutation.PACKDIM.THREED);
        break;
      case LAYERPACKINGBYSIZE:
        init =
            new PackingInit(
                NorwayGeometryMutation.MUTMODE.BYSIZE, NorwayGeometryMutation.PACKDIM.TWOD);
        break;
      case LAYERPACKINGRANDOMLY:
        init =
            new PackingInit(
                NorwayGeometryMutation.MUTMODE.RANDOM, NorwayGeometryMutation.PACKDIM.TWOD);
        break;
      case LAYERPACKINGASCENDING:
        init =
            new PackingInit(
                NorwayGeometryMutation.MUTMODE.ASCENDING, NorwayGeometryMutation.PACKDIM.TWOD);
        break;
      case RANDOMWITHDD:
        init = new RandomizedGeomInit(true);
        break;
      case RANDOMWITHOUTDD:
        init = new RandomizedGeomInit(false);
        break;
      default:
        System.err.println(
            "ERROR: No such initialization routine representing "
                + whichInit.name()
                + " using packing without size consideration now.");
        init =
            new PackingInit(
                NorwayGeometryMutation.MUTMODE.ASCENDING, NorwayGeometryMutation.PACKDIM.THREED);
        break;
    }
  }

  @Override
  public Geometry initTheGeometry(
      Geometry geom,
      final CollisionDetection.CDTYPE whichCollisionDetection,
      final DissociationDetection.DDTYPE whichDissociationDetection,
      final double[] cellSize,
      final double blowDissDetect,
      final double blowCollDetect,
      final float explDoFRatio,
      final boolean molecularCD) {

    final Geometry geomNew = new Geometry(geom);

    return init.initTheGeometry(
        geomNew,
        whichCollisionDetection,
        whichDissociationDetection,
        cellSize,
        blowDissDetect,
        blowCollDetect,
        explDoFRatio,
        molecularCD);
  }

  public static INITSTYLE parseType(final String type) throws Exception {

    if (type.equalsIgnoreCase("packbysize")) {
      return INITSTYLE.PACKBYSIZE;
    } else if (type.equalsIgnoreCase("packrandomly")) {
      return INITSTYLE.PACKRANDOMLY;
    } else if (type.equalsIgnoreCase("packascending")) {
      return INITSTYLE.PACKASCENDING;
    } else if (type.equalsIgnoreCase("packlayerbysize")) {
      return INITSTYLE.LAYERPACKINGBYSIZE;
    } else if (type.equalsIgnoreCase("packlayerrandomly")) {
      return INITSTYLE.LAYERPACKINGRANDOMLY;
    } else if (type.equalsIgnoreCase("packlayerascending")) {
      return INITSTYLE.LAYERPACKINGASCENDING;
    } else if (type.equalsIgnoreCase("randomwithdd")) {
      return INITSTYLE.RANDOMWITHDD;
    } else if (type.equalsIgnoreCase("randomwithoutdd")) {
      return INITSTYLE.RANDOMWITHOUTDD;
    }

    throw new Exception("Illegal geometry initialization " + type + ".");
  }
}
