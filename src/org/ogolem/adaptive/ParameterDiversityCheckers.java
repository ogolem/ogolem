/*
Copyright (c) 2012-2014, J. M. Dieterich
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
package org.ogolem.adaptive;

import org.ogolem.generic.genericpool.DiversityChecker;
import org.ogolem.generic.genericpool.GenericPoolEntry;

/**
 * @author Johannes Dieterich
 * @version 2014-01-20
 */
public class ParameterDiversityCheckers {

  /** Checks the diversity of two structures based on the actual parameters */
  public static class ParamsDiversityChecker
      implements DiversityChecker<Double, AdaptiveParameters> {

    private static final long serialVersionUID = (long) 20131231;
    private final double divThresh;

    public ParamsDiversityChecker(final double thresh) {
      this.divThresh = thresh;
    }

    @Override
    public boolean areDiverse(
        final GenericPoolEntry<Double, AdaptiveParameters> individuum1,
        final GenericPoolEntry<Double, AdaptiveParameters> individuum2) {

      final double[] daParams1 = individuum1.getIndividual().getAllParamters();
      final double[] daParams2 = individuum2.getIndividual().getAllParamters();

      for (int i = 0; i < daParams1.length; i++) {
        if (Math.abs(daParams1[i] - daParams2[i]) > divThresh) {
          // diverse enough in (at least) one coordinate
          return true;
        }
      }

      return false;
    }

    @Override
    public String getMyName() {
      return "parameter diversity checker\n\tthreshold " + divThresh;
    }
  }
}
