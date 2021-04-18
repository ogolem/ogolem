/*
Copyright (c) 2013, J. M. Dieterich
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

import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;

/**
 * Nicher working based on the number of interior molecules
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
class InteriorMoleculeNicheComp implements NicheComputer<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20200429;
  private static final boolean DEBUG = true;
  private final SurfaceDetectionEngine surfer;

  InteriorMoleculeNicheComp(final SurfaceDetectionEngine surfer) {
    this.surfer = surfer;
  }

  private InteriorMoleculeNicheComp(final InteriorMoleculeNicheComp orig) {
    this.surfer = orig.surfer.copy();
  }

  @Override
  public InteriorMoleculeNicheComp copy() {
    return new InteriorMoleculeNicheComp(this);
  }

  @Override
  public Niche computeNiche(final Geometry g) {

    final Surface surf = surfer.detectTheSurface(g.getCartesians());
    final int[] molSurf = surf.giveAllSurfaceMolecules();

    int numInternal = 0;
    for (int i = 0; i < g.getNumberOfIndieParticles(); i++) {
      boolean found = false;
      for (int j = 0; j < molSurf.length; j++) {
        if (molSurf[j] == i) {
          found = true;
          break;
        }
      }
      if (DEBUG) {
        System.out.println("DEBUG: Molecule " + i + " found in surface? " + found);
      }
      if (!found) {
        numInternal++;
      }
    }

    if (DEBUG) {
      System.out.println("DEBUG: Number of interior molecules " + numInternal);
    }

    return new Niche("interiormolecule" + numInternal);
  }
}
