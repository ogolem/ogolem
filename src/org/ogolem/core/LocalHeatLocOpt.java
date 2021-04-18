/*
Copyright (c) 2012-2013, J. M. Dieterich
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

import org.ogolem.generic.GenericAbstractLocOpt;
import org.ogolem.generic.GenericInitializer;
import org.ogolem.generic.GenericLocOpt;
import org.ogolem.heat.LocalHeatPulses;

/**
 * Uses Moebius et al local heat pulses (which are meant as a global optimization) as a local
 * optimization backend.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class LocalHeatLocOpt extends GenericAbstractLocOpt<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20170625;

  private final GenericLocOpt<Molecule, Geometry> locopt;
  private final CollisionDetectionEngine cd;
  private final LocalHeatPulses.Configuration heatConfig;
  private final double acceptableFitness;
  private final double blowBonds;
  private final double blowDissoc;
  private final GlobalConfig config;

  LocalHeatLocOpt(
      final GlobalConfig config,
      final GenericLocOpt<Molecule, Geometry> refLoc,
      final LocalHeatPulses.Configuration heatConfig)
      throws Exception {
    super(refLoc.getBackend());
    this.locopt = refLoc;
    this.cd = new CollisionDetection(config.whichCollisionEngine);
    this.heatConfig = heatConfig;
    this.acceptableFitness = config.acceptableFitness;
    this.blowBonds = config.blowFacBondDetect;
    this.blowDissoc = config.blowFacDissocDetect;
    this.config = config;
  }

  private LocalHeatLocOpt(final LocalHeatLocOpt orig) {
    super(orig);
    this.locopt = orig.locopt.copy();
    this.cd = orig.cd.copy();
    this.heatConfig = new LocalHeatPulses.Configuration(orig.heatConfig);
    this.acceptableFitness = orig.acceptableFitness;
    this.blowBonds = orig.blowBonds;
    this.blowDissoc = orig.blowDissoc;
    this.config = orig.config;
  }

  @Override
  public LocalHeatLocOpt copy() {
    return new LocalHeatLocOpt(this);
  }

  @Override
  public String getMyID() {
    String id = "Local heat pulses: " + locopt.getMyID();
    id += "\n" + heatConfig.printConfig();

    return id;
  }

  @Override
  protected Geometry optimize(Geometry gStartGeometry) {

    // we pull a new initer here, as we otherwise above run into an internal cyclic dependency (the
    // initer depends on the backend, which depends on the locopt, which depends on the... initer).
    GenericInitializer<Molecule, Geometry> initer = null;
    try {
      initer = config.getInitializer();
    } catch (Exception e) {
      throw new RuntimeException("Getting the initializer not working. This should not happen!", e);
    }

    assert (initer != null);

    return LocalHeatPulses.cycle(
        gStartGeometry,
        locopt,
        cd,
        blowBonds,
        blowDissoc,
        gStartGeometry.getBondInfo(),
        heatConfig,
        acceptableFitness,
        initer);
  }
}
