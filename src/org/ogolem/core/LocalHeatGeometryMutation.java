/*
Copyright (c) 2017-2020, J. M. Dieterich and B. Hartke
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

import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.GenericMutation;
import org.ogolem.heat.LocalHeatPulses;

/**
 * A mutation using local heat pulses.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
class LocalHeatGeometryMutation implements GenericMutation<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20170305;

  private final GenericFitnessFunction<Molecule, Geometry> locopt;
  private final CollisionDetectionEngine cd;
  private final LocalHeatPulses.Configuration heatConfig;
  private final double blowBonds;
  private final double blowDissoc;
  private final double acceptableFitness;

  LocalHeatGeometryMutation(
      final GenericFitnessFunction<Molecule, Geometry> refLoc,
      final LocalHeatPulses.Configuration heatConfig,
      final CollisionDetection.CDTYPE cdType,
      final double blowBonds,
      final double blowDissoc,
      final double acceptableFitness) {

    assert (refLoc != null);
    assert (heatConfig != null);
    this.locopt = refLoc;
    this.heatConfig = heatConfig;
    this.cd = new CollisionDetection(cdType);
    this.blowBonds = blowBonds;
    this.blowDissoc = blowDissoc;
    this.acceptableFitness = acceptableFitness;

    // make sure we are NOT trying to reset to random (otherwise below handing over null will give
    // NPE)
    this.heatConfig.resetToRandom = false;
  }

  private LocalHeatGeometryMutation(final LocalHeatGeometryMutation orig) {
    this.cd = orig.cd.copy();
    this.heatConfig = orig.heatConfig;
    this.locopt = orig.locopt.copy();
    this.blowBonds = orig.blowBonds;
    this.blowDissoc = orig.blowDissoc;
    this.acceptableFitness = orig.acceptableFitness;
  }

  @Override
  public LocalHeatGeometryMutation copy() {
    return new LocalHeatGeometryMutation(this);
  }

  @Override
  public String getMyID() {

    String id = "Local heat pulses: " + locopt.getMyID();
    id += "\n" + heatConfig.printConfig();

    return id;
  }

  @Override
  public Geometry mutate(Geometry orig) {
    return LocalHeatPulses.cycle(
        orig,
        locopt,
        cd,
        blowBonds,
        blowDissoc,
        orig.getBondInfo(),
        heatConfig,
        acceptableFitness,
        null);
  }
}
