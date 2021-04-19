/*
Copyright (c) 2016, J. M. Dieterich
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

      This product includes software developed at the Universities of
      Kiel, Goettingen (Germany) and Princeton University (USA) by its
      contributors: J. M. Dieterich and B. Hartke.

    * Neither the name of the University of Kiel, the University of Goettingen,
      Princeton University nor the names of its contributors may be used to
      endorse or promote products derived from this software without specific
      prior written permission.

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
 * An interface describing structural data value types. Really just a very high-level description!
 *
 * @author Johannes M Dieterich
 * @version 2020-12-30
 */
public interface StructuralData extends Serializable, Copyable {

  public StructuralData copy();

  public CartesianCoordinates getCartesianCoordinates();

  public int getNoOfAtoms();

  public int getNoOfMolecules();

  public short[] getAllAtomNumbers();

  public short[] getAllSpins();

  public float[] getAllCharges();

  public void setAllCharges(final float[] newCharges);

  public void setAllSpins(final short[] newSpins);
}
