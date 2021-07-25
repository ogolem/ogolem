/*
Copyright (c) 2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.math;

/**
 * A matrix implementation for a symmetric matrix without diagonal.
 *
 * @author Johannes Dieterich
 * @version 2020-12-21
 */
public class SymmetricMatrixNoDiag implements Matrix {

  private static final long serialVersionUID = (long) 20201221;

  private final int noRowsCols;
  private final double[] buffer;

  public SymmetricMatrixNoDiag(final int noRowsCols) {
    assert (noRowsCols >= 0);
    this.noRowsCols = noRowsCols;
    this.buffer = new double[noRowsCols * (noRowsCols - 1) / 2];
  }

  @Override
  public int noRows() {
    return noRowsCols;
  }

  @Override
  public int noCols() {
    return noRowsCols;
  }

  @Override
  public void setElement(final int i, final int j, final double val) {

    // we always assume j > i
    final int row = Math.min(j, i);
    final int col = Math.max(j, i);

    final int idx =
        (noRowsCols * (noRowsCols - 1) / 2)
            - (noRowsCols - row) * ((noRowsCols - row) - 1) / 2
            + col
            - row
            - 1;

    buffer[idx] = val;
  }

  @Override
  public double getElement(final int i, final int j) {

    // we always assume j > i
    final int row = Math.min(j, i);
    final int col = Math.max(j, i);

    final int idx =
        (noRowsCols * (noRowsCols - 1) / 2)
            - (noRowsCols - row) * ((noRowsCols - row) - 1) / 2
            + col
            - row
            - 1;

    return buffer[idx];
  }

  /** Returns the n*(n-1)/2 storage buffer where j > i always. */
  @Override
  public double[] underlyingStorageBuffer() {
    return buffer;
  }
}
