/*
Copyright (c) 2012-2014, J. M. Dieterich
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

import java.io.File;
import java.util.Iterator;
import org.ogolem.generic.IndividualWriter;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.properties.Property;

/**
 * Writes one geometry out.
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class GeometryWriter implements IndividualWriter<Geometry> {

  private static final long serialVersionUID = (long) 20140330;
  private final String outFolder;

  public GeometryWriter(final String outFolder) {
    this.outFolder = outFolder;
  }

  @Override
  public GeometryWriter copy() {
    return new GeometryWriter(this.outFolder);
  }

  @Override
  public void writeIndividual(final Geometry g) {
    // write xyz out
    try {
      Output.printXYZSingleGeometry(g, true, outFolder);
    } catch (Exception e) {
      System.err.println("WARNING: Couldn't write geometry " + g.getID() + ". " + e.toString());
    }

    // see about properties
    final Iterator<Property> it = g.getPropertyIterator();
    if (!it.hasNext()) {
      // empty iterator, probably the more usual case
      return;
    }

    final String propOut = outFolder + File.separator + "geometry" + g.getID() + ".properties";
    String s = "";
    while (it.hasNext()) {
      final Property p = it.next();
      s += p.name() + "\n\t" + p.printableProperty() + "\n\n";
    }

    try {
      OutputPrimitives.writeOut(propOut, s, true);
    } catch (Exception e) {
      System.err.println(
          "WARNING: Couldn't write properties of geometry " + g.getID() + ". " + e.toString());
    }
  }

  @Override
  public void writeIndividual(final Geometry g, final String toFile) {
    // write xyz out
    try {
      final String[] dat = g.makePrintableAbsoluteCoord(true);
      OutputPrimitives.writeOut(toFile + ".xyz", dat, true);
    } catch (Exception e) {
      System.err.println("WARNING: Couldn't write geometry " + g.getID() + ". " + e.toString());
    }

    // see about properties
    final Iterator<Property> it = g.getPropertyIterator();
    if (!it.hasNext()) {
      // empty iterator, probably the more usual case
      return;
    }

    final String propOut = toFile + ".properties";
    String s = "";
    while (it.hasNext()) {
      final Property p = it.next();
      s += p.name() + "\n\t" + p.printableProperty() + "\n\n";
    }

    try {
      OutputPrimitives.writeOut(propOut, s, true);
    } catch (Exception e) {
      System.err.println(
          "WARNING: Couldn't write properties of geometry " + g.getID() + ". " + e.toString());
    }
  }
}
