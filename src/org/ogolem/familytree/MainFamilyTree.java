/*
Copyright (c) 2015i-2021, J. M. Dieterich and B. Hartke
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
package org.ogolem.familytree;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.GenericPoolEntry;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.OutputPrimitives;

/**
 * Builds a gv/dot format file containing the directed acyclic graph of the family tree.
 *
 * @author Johannes Dieterich
 * @version 2021-10-15
 */
public class MainFamilyTree {

  public static void run(final String args[]) {

    if (args[0].equalsIgnoreCase("help")) {
      System.out.println("This is the genetic history visualizer: it builds a family tree.");
      System.out.println("Mandatory input:");
      System.out.println(" * the path to an ogolem output file");
      System.out.println("Optional input:");
      System.out.println(" * the path to a configuration file");
      return;
    }

    final String outFile = args[0];
    VisualizationConfig config = null;
    if (args.length > 1) {
      try {
        final String[] data = InputPrimitives.readFileIn(args[1]);
        config = new VisualizationConfig(data);
      } catch (Exception e) {
        System.err.println("Failure to configure visualization printer.");
        e.printStackTrace(System.err);
        System.exit(1);
      }
    } else {
      config = new VisualizationConfig();
    }

    assert (config != null);

    // parse the output file
    String[] outData = null;
    try {
      outData = InputPrimitives.readFileIn(outFile);
    } catch (Exception e) {
      System.err.println("Failure to read the output file in.");
      e.printStackTrace(System.err);
      System.exit(2);
    }

    assert (outData != null);

    final List<GeneticRecord> allRecords = new ArrayList<>();
    boolean startFound = false;
    boolean firstAfterFound = true;
    long maxNoChildren = 0;
    long minIDInHistory = Long.MAX_VALUE;
    for (final String outLine : outData) {
      if (outLine.trim().startsWith("The following genetic history was created during the run.")) {
        startFound = true;
      } else if (startFound) {
        if (firstAfterFound) {
          firstAfterFound = false;
          continue;
        }
        if (outLine.trim().startsWith("---------------------")) {
          // end reached
          break;
        }
        // parse line to find minimum id
        final String[] info = outLine.trim().split("\\s+");
        try {
          final int id = Integer.parseInt(info[0]);
          minIDInHistory = Math.min(minIDInHistory, id);

          // add this genetic record
          final int mother = Integer.parseInt(info[1]);
          final int father = Integer.parseInt(info[2]);
          final boolean wasAccepted = Boolean.parseBoolean(info[3]);
          final boolean wasNull = Boolean.parseBoolean(info[4]);

          final GeneticRecord record = new GeneticRecord(id, mother, father, wasAccepted, wasNull);
          allRecords.add(record);

        } catch (Exception e) {
          System.err.println("Failure to parse line information (I).");
          e.printStackTrace(System.err);
          System.exit(3);
        }
      }
    }

    // add also the 0 - (minIDInHistory-1) records!
    for (int oid = 0; oid < minIDInHistory; oid++) {
      final GeneticRecord oRec = new GeneticRecord(oid, -1, -1, true, false);
      allRecords.add(oRec);
    }

    // sort by ID
    Collections.sort(allRecords, new GeneticRecordComparator());

    // now mark the parenthood
    for (final GeneticRecord rec : allRecords) {

      final int father = rec.father();
      final int mother = rec.mother();

      if (father >= 0) {
        final GeneticRecord faRec = allRecords.get(father);
        faRec.incrementChildren();
        maxNoChildren = Math.max(faRec.noChildren(), maxNoChildren);
      }

      if (mother >= 0) {
        final GeneticRecord moRec = allRecords.get(mother);
        moRec.incrementChildren();
        maxNoChildren = Math.max(moRec.noChildren(), maxNoChildren);
      }
    }

    if (config.addAllFinalPoolIndividuals) {
      // read the pool in and mark the individuals
      try {
        final GenericPool<?, ?> pool =
            (GenericPool<?, ?>) InputPrimitives.readBinInput(config.poolFile);
        for (int pos = 0; pos < pool.getCurrentPoolSize(); pos++) {
          final GenericPoolEntry<?, ?> entry = pool.getEntryAtPosition(pos);
          final int id = (int) entry.individual().getID();
          allRecords.get(id).setInFinalPool(true);
        }
      } catch (Exception e) {
        System.err.println("Failed to read the pool file " + config.poolFile);
        e.printStackTrace(System.err);
        System.exit(31);
      }
    }

    final Map<Integer, Integer> map = new HashMap<>();
    if (config.removeAllNonChildren) {
      // weed out all the entries which have not reproduced and save the position of those that have
      int countMultiChildren = 0;
      final List<GeneticRecord> toDelete = new ArrayList<>();
      for (int i = 0; i < allRecords.size(); i++) {
        final GeneticRecord rec = allRecords.get(i);
        boolean secondaryCond = true;
        if (config.addAllFinalPoolIndividuals) {
          secondaryCond = !rec.isInFinalPool();
        }
        if (rec.noChildren() == 0 && secondaryCond) {
          // remove
          toDelete.add(rec);
        } else {
          // keep
          map.put(rec.id(), countMultiChildren);
          countMultiChildren++;
        }
      }

      // remove all
      toDelete.forEach(
          (rec) -> {
            allRecords.remove(rec);
          });
    } else {
      // weed out all the entries which were null and/or were not accepted
      int countAccepted = 0;
      final List<GeneticRecord> toDelete = new ArrayList<>();
      for (int i = 0; i < allRecords.size(); i++) {
        final GeneticRecord rec = allRecords.get(i);
        boolean secondaryCond = true;
        if (config.addAllFinalPoolIndividuals) {
          secondaryCond = !rec.isInFinalPool();
        }
        if ((rec.wasNull() && secondaryCond) || (!rec.wasAccepted() && secondaryCond)) {
          // remove
          toDelete.add(rec);
        } else {
          // keep
          map.put(rec.id(), countAccepted);
          countAccepted++;
        }
      }

      // remove all
      toDelete.forEach(
          (rec) -> {
            allRecords.remove(rec);
          });
    }

    // setup the dot file
    final List<String> datOut = new LinkedList<>();
    datOut.add("digraph FamilyTree {");
    datOut.add("   bgcolor = \"" + config.backGround + "\";");

    // first all the vertices
    for (int vertID = 0; vertID < allRecords.size(); vertID++) {

      final GeneticRecord rec = allRecords.get(vertID);

      String hex = null;
      if (rec.isInFinalPool()) {
        hex = colorToHex(config.vertexColorFinalPool, config.transparency);
      } else {
        hex = colorToHex(config.vertexColor, config.transparency);
      }

      if (config.vertexSizeChildrenBased) {
        // TODO: change size depending on #children
        System.err.println(
            "No support yet to adjust the vertex size based on the number of children, sorry...");
        System.exit(42);
      }
      String label = "id=" + rec.id();
      if (config.verboseLabel) {
        // add the number of children to the label
        label += " / children=" + rec.noChildren();
      }
      datOut.add(
          "   "
              + (vertID + 1)
              + " [label=\""
              + label
              + "\", shape = \""
              + config.vertexShape
              + "\", style = \""
              + config.vertexStyle
              + "\", color = \""
              + config.circleRing
              + "\", fillcolor = \""
              + hex
              + "\"];");
    }

    // now all the edges
    for (int vertID = 0; vertID < allRecords.size(); vertID++) {

      final GeneticRecord rec = allRecords.get(vertID);

      // the edge to the mother
      final String hexMum = colorToHex(config.edgeColorAsMum, config.transparency);
      final int mumID = rec.mother();
      if (mumID >= 0) {
        final Integer mumVert = map.get(mumID);
        if (mumVert == null) {
          System.err.println(
              "Vertex (mum) " + mumID + " not left in map. Child id " + rec.id() + ". Exiting.");
          System.exit(10);
        }
        final int mumVertex = mumVert;
        datOut.add(
            "   "
                + (mumVertex + 1)
                + " -> "
                + (vertID + 1)
                + " [style = \""
                + config.edgeStyle
                + "\", color = \""
                + hexMum
                + "\", penwidth = \""
                + config.penWidth
                + "\"];");
      }

      // the edge to the father
      final String hexDad = colorToHex(config.edgeColorAsDad, config.transparency);
      final int dadID = rec.father();
      if (dadID >= 0) {
        final Integer dadVert = map.get(dadID);
        if (dadVert == null) {
          System.err.println(
              "Vertex (dad) " + mumID + " not left in map. Child id " + rec.id() + ". Exiting.");
          System.exit(10);
        }
        final int dadVertex = dadVert;
        datOut.add(
            "   "
                + (dadVertex + 1)
                + " -> "
                + (vertID + 1)
                + " [style = \""
                + config.edgeStyle
                + "\", color = \""
                + hexDad
                + "\", penwidth = \""
                + config.penWidth
                + "\"];");
      }
    }

    // print to file
    datOut.add("}");
    try {
      OutputPrimitives.writeOut(config.dotOut, datOut, true);
    } catch (Exception e) {
      System.err.println("ERROR: Couldn't write final dot file out.");
      e.printStackTrace(System.err);
    }
  }

  private static String colorToHex(final float[] colDef, final float trans) {
    assert (colDef.length == 3);
    final Color col = new Color(colDef[0], colDef[1], colDef[2], trans);
    return "#"
        + toBrowserHexValue(col.getRed())
        + toBrowserHexValue(col.getGreen())
        + toBrowserHexValue(col.getBlue())
        + toBrowserHexValue(col.getAlpha());
  }

  private static String toBrowserHexValue(int number) {
    final StringBuilder builder = new StringBuilder(Integer.toHexString(number & 0xff));
    while (builder.length() < 2) {
      builder.append("0");
    }
    return builder.toString().toUpperCase();
  }
}
