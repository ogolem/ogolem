/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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
 * Packs the geometry for initialization.
 *
 * @author Johannes Dieterich
 * @version 2020-08-09
 */
public final class PackingInit implements GeometryInitialization {

  private static final long serialVersionUID = (long) 20200809;

  private final NorwayGeometryMutation.MUTMODE packingMode;
  private final NorwayGeometryMutation.PACKDIM packDim;

  public PackingInit(
      final NorwayGeometryMutation.MUTMODE mode, final NorwayGeometryMutation.PACKDIM packDim) {
    this.packingMode = mode;
    this.packDim = packDim;
  }

  @Override
  public Geometry initTheGeometry(
      final Geometry geom,
      final CollisionDetection.CDTYPE whichCollDetect,
      final DissociationDetection.DDTYPE whichDissocDetect,
      final double[] cellSize,
      final double blowDiss,
      final double blowColl,
      final float explDoFRatio,
      final boolean molecularCD) {

    final NorwayGeometryMutation norway =
        new NorwayGeometryMutation(
            packDim, whichCollDetect, blowColl, blowDiss, whichDissocDetect, packingMode);

    final Geometry packed = norway.mutate(geom);

    return packed;
  }
}
