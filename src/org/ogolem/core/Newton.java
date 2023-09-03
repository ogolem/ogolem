/*
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2012, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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

import java.io.Serializable;
import org.ogolem.generic.Copyable;

/**
 * Defining an interface for local optimization methods.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
public interface Newton extends Copyable, Serializable {

  Newton copy();

  /**
   * Locally optimize one molecule/building block of the cluster.
   *
   * @param mStartMolecule The unoptimized molecule
   * @return An optimized molecule.
   */
  Molecule localOptimization(Molecule mStartMolecule);

  Geometry localOptimization(Geometry gStartGeometry);

  /**
   * It depends on the implementation, whether the input cartesian stays unchanged or not. The usage
   * of this method might be problematic in case that e.g. the implementation relies on the length
   * of certain auxiliary information. USE WITH CARE AND IMPLEMENT PROPERLY!
   *
   * @param id An identifier.
   * @param cartes A set of cartesian coordinates.
   * @param constraints The constraints of these coordinates.
   * @param isConstricted if the structure is at all constricted.
   * @return An optimized set of cartesian coordinates.
   * @throws Exception if something during the optimization went wrong.
   */
  CartesianCoordinates cartesToCartes(
      long id,
      CartesianCoordinates cartes,
      boolean[][] constraints,
      boolean isConstricted,
      final BondInfo bonds)
      throws Exception;

  String myIDandMethod();

  /**
   * @return The Backend, if available. Else null.
   */
  CartesianFullBackend getBackend();

  long getNumberOfGeomLocalOpts();

  long getNumberOfMolLocalOpts();
}
