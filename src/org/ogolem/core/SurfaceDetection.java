/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2011, J. M. Dieterich
              2015-2020, J. M. Dieterich and B. Harkte
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
 * Works as a decorator to all surface dections implemented.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
class SurfaceDetection implements SurfaceDetectionEngine {

  private static final long serialVersionUID = (long) 20150428;

  public static enum SURFDETECTTYPE {
    HARTKESURFACEDETECTTESS,
    HARTKESURFACEDETECTLONGLAT
  };

  public static final SURFDETECTTYPE DEFAUTLSURFDETECT = SURFDETECTTYPE.HARTKESURFACEDETECTTESS;

  protected final SURFDETECTTYPE whichSurfDetect;
  protected final SurfaceDetectionEngine surfaceDetect;

  SurfaceDetection(final SURFDETECTTYPE whichSurfaceDetection) {
    throw new RuntimeException("Surface detection not a stable feature yet, sorry!");
  }

  @Override
  public Surface detectTheSurface(final CartesianCoordinates cartes) {
    return surfaceDetect.detectTheSurface(cartes);
  }

  @Override
  public SurfaceDetection copy() {
    return new SurfaceDetection(whichSurfDetect);
  }

  static SURFDETECTTYPE parseType(final String type) throws Exception {
    if (type.equalsIgnoreCase("hartketess")) {
      return SURFDETECTTYPE.HARTKESURFACEDETECTTESS;
    } else if (type.equalsIgnoreCase("hartkelonglat")) {
      return SURFDETECTTYPE.HARTKESURFACEDETECTLONGLAT;
    }

    throw new Exception("Illegal surface detection " + type + ".");
  }
}
