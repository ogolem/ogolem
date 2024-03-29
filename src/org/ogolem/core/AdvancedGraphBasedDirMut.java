/*
Copyright (c) 2013-2015, J. M. Dieterich and B. Hartke
              2017-2020, J. M. Dieterich and B. Hartke
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

import static org.ogolem.core.DirMutOptStrategies.PointOptStrategy;
import static org.ogolem.core.DirMutPointProviders.PointProvider;

import java.awt.Color;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import org.ogolem.generic.GenericFitnessBackend;
import org.ogolem.generic.GenericLocOpt;
import org.ogolem.generic.GenericMutation;
import org.ogolem.helpers.Tuple;
import org.ogolem.math.BoolSymmetricMatrixNoDiag;

/**
 * A directed mutation using a graph-based analysis of the cluster in question.
 *
 * @author Johannes Dieterich
 * @author Bernd Hartke
 * @version 2020-12-30
 */
public class AdvancedGraphBasedDirMut implements GenericMutation<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20200627;

  private final boolean DEBUG;
  private final double blowBonds;
  private final CollisionDetectionEngine collDetect;
  private final double blowColl;
  private final GenericLocOpt<Molecule, Geometry> locopt;
  private final GenericFitnessBackend<Molecule, Geometry> back;
  private final boolean doLocOpt;
  private final boolean doCDFirst;
  private final double realCD;
  private final double realDD;
  private final PointProvider pointProv;
  private final PointOptStrategy pointOpt;
  private final boolean fullyRelaxed;
  private final boolean markMovedMolsUnmovable;
  private final boolean doCDForEveryTrial;
  private final int noMoved;
  private final boolean envAware;

  /**
   * A GDM constructor.
   *
   * @param config the GDM configuration
   * @param collDetect the collision detection engine to be used
   * @param locopt the local optimization engine (including a backend!) to be used
   * @param realCDBlow the CD blow factor used in the rest of the program, for debugging purposes
   * @param realDDBlow the DD blow factor used in the rest of the program, for debugging purposes
   * @param doLocOptFirst do a local optimization first, prior to the GDM
   * @param doCDFirst do a CD first, AFTER the local optimization (if defined) and the GDM
   * @param provider the point provider (e.g. a grid)
   * @param strategy the strategy what to do on each grid point (e.g. an optimization or only a
   *     rigid body opt or...)
   * @param envAware if an environment is present in the system, should the algorithm be aware of
   *     this to compute collisions and bonds?
   */
  public AdvancedGraphBasedDirMut(
      final AdvancedGraphBasedDirMut.GDMConfiguration config,
      final CollisionDetectionEngine collDetect,
      final GenericLocOpt<Molecule, Geometry> locopt,
      final double realCDBlow,
      final double realDDBlow,
      final boolean doLocOptFirst,
      final boolean doCDFirst,
      final PointProvider provider,
      final PointOptStrategy strategy,
      final boolean envAware) {
    this.pointProv = provider;
    this.pointOpt = strategy;
    this.DEBUG = (GlobalConfig.DEBUGLEVEL > 1);
    this.blowBonds = config.blowBonds;
    this.collDetect = collDetect;
    this.blowColl = config.blowCD;
    this.doLocOpt = doLocOptFirst;
    this.doCDFirst = doCDFirst;
    if (locopt.getBackend() == null) {
      throw new RuntimeException("GDM must have a non-null Backend from Newton!");
    }
    this.locopt = locopt;
    this.back = locopt.getBackend();
    this.realCD = realCDBlow;
    this.realDD = realDDBlow;
    this.fullyRelaxed = config.fullOptEachMove;
    this.markMovedMolsUnmovable = config.markMovedMolsUnmoveable;
    this.noMoved = config.noOfMolsToMove;
    this.doCDForEveryTrial = config.doCDEveryTrial;
    this.envAware = envAware;

    assert (blowBonds > 0.0);
    assert (blowColl > 0.0);
  }

  private AdvancedGraphBasedDirMut(final AdvancedGraphBasedDirMut orig) {
    this.DEBUG = (GlobalConfig.DEBUGLEVEL > 1);
    this.blowBonds = orig.blowBonds;
    this.collDetect = orig.collDetect.copy();
    this.blowColl = orig.blowColl;
    this.back = orig.back.copy();
    this.locopt = orig.locopt.copy();
    this.realCD = orig.realCD;
    this.realDD = orig.realDD;
    this.doLocOpt = orig.doLocOpt;
    this.fullyRelaxed = orig.fullyRelaxed;
    this.doCDFirst = orig.doCDFirst;
    this.pointProv = orig.pointProv.copy();
    this.pointOpt = orig.pointOpt.copy();
    this.markMovedMolsUnmovable = orig.markMovedMolsUnmovable;
    this.noMoved = orig.noMoved;
    this.doCDForEveryTrial = orig.doCDForEveryTrial;
    this.envAware = orig.envAware;
  }

  @Override
  public AdvancedGraphBasedDirMut copy() {
    return new AdvancedGraphBasedDirMut(this);
  }

  @Override
  public String getMyID() {
    return "advanced graph based directed mutation: \n\t point provider: "
        + pointProv.getMyID()
        + "\n\t opt strategy? "
        + pointOpt.getMyID()
        + "\n\t blow factor bonds "
        + blowBonds
        + "\n\t blow factor collision "
        + blowColl
        + "\n\t CD first? "
        + doCDFirst
        + "\n\t opt first? "
        + doLocOpt
        + "\n\t number of moved mols: "
        + noMoved
        + "\n\t mark moved mols unmoveable? "
        + markMovedMolsUnmovable
        + "\n\t fully relaxed? "
        + fullyRelaxed
        + "\n\t do CD for every trial? "
        + doCDForEveryTrial;
  }

  @Override
  public Geometry mutate(final Geometry geom) {

    if (DEBUG) {
      System.out.println("DEBUG: Running our fancy ADVANCED graph based directed mutation.");
    }

    Geometry work = new Geometry(geom);

    if (doLocOpt) {
      if (DEBUG) {
        System.out.println("DEBUG: doing a local optimization first.");
      }
      work = locopt.fitness(work, false);
    }

    if (doCDFirst) {
      if (DEBUG) {
        System.out.println("DEBUG: doing a CD now.");
      }
      final boolean hasColl =
          collDetect.checkOnlyForCollision(work.getCartesians(), blowColl, work.getBondInfo());
      if (hasColl) {
        if (DEBUG) {
          System.out.println("DEBUG: found a collision, returning null.");
        }
        return null;
      }
    }

    final int noMols = work.getNumberOfIndieParticles();
    final int noAtoms = work.getNumberOfAtoms();
    final int noEnvAtoms = (work.containsEnvironment()) ? work.getEnvironment().atomsInEnv() : 0;

    // get the energy before we manipulate
    final double[] currCoords = back.getActiveCoordinates(work.copy());
    final double eBefore = back.fitness(currCoords, 42);

    final List<Integer> movedMols = new ArrayList<>(noMoved);
    for (int moveCounter = 0; moveCounter < noMoved; moveCounter++) {

      if (fullyRelaxed) {
        if (DEBUG) {
          System.out.println("DEBUG: Optimizing at step " + moveCounter);
        }
        work = locopt.fitness(work, false);
        if (DEBUG) {
          System.out.println("DEBUG: Post optimization energy: " + work.getFitness());
        }
      }

      // create the connectivity matrix for the geometry (so in between molecules)
      // 1) get the full one
      final SimpleBondInfo fullBonds =
          (work.containsEnvironment() && envAware)
              ? CoordTranslation.checkForBondsIncludingEnvironment(work, blowBonds)
              : CoordTranslation.checkForBonds(work, blowBonds);

      final BoolSymmetricMatrixNoDiag bondMat = fullBonds.getFullBondMatrix();

      if (DEBUG) {
        System.out.println("Bonds with blow factor " + blowBonds);
        final boolean[] bondsBuffer = bondMat.underlyingStorageBuffer();
        System.out.println(Arrays.toString(bondsBuffer));
      }

      // figure out the number of connections from each molecule to another one (do not count
      // internal bonds)
      final int[] connections = new int[noMols];
      int[][] connCounter = null;
      if (DEBUG) {
        connCounter = new int[noMols][noMols];
      }
      int atCounter = 0;
      for (int i = 0; i < noMols; i++) {
        final int ats = work.getMoleculeAtPosition(i).getNumberOfAtoms();
        for (int at = atCounter; at < atCounter + ats; at++) {
          // go over the bonds (remove self interactions!)
          // we DO count if multiple atoms of the same molecule connect to another molecule
          // not only is this easier, but it makes more sense for the application at hand
          // which is to find the least connected molecule
          for (int x = 0; x < noAtoms; x++) {
            if (at == x) continue;
            if (bondMat.getElement(at, x)) {
              if (x < atCounter || x >= atCounter + ats) {
                if (DEBUG) {

                  // figure out the molecule this atom belongs to
                  int otherMol = -1;
                  int count = 0;
                  for (int y = 0; y < noAtoms; y++) {
                    count += work.getMoleculeAtPosition(y).getNumberOfAtoms();
                    if (count > x) {
                      otherMol = y;
                      break;
                    }
                  }

                  System.out.println(
                      "Connection found for molecule " + i + " " + "to molecule " + otherMol);
                  connCounter[i][otherMol]++;
                  connCounter[otherMol][i]++;
                }
                connections[i]++;
              }
            }
          }

          if (work.containsEnvironment() && envAware) {
            // count connections with the environment
            for (int x = 0; x < noEnvAtoms; x++) {
              if (at == x + noAtoms) continue;
              if (bondMat.getElement(at, x + noAtoms)) {
                connections[i]++;
              }
            }
          }
        }

        atCounter += ats;
      }

      if (DEBUG) {
        System.out.println("DEBUG: Following molecular connection information obtained:");
        for (int i = 0; i < noMols; i++) {
          System.out.println(" " + i + "\t" + connections[i]);
        }
      }

      // sort the molecules by connectivity in ascending order
      // XXX not too efficient, could be fused with above loop but whatever
      final List<Tuple<Integer, Integer>> sortedMolsByConn = new ArrayList<>(noMols);
      Outer:
      for (int i = 0; i < noMols; i++) {
        final Tuple<Integer, Integer> tup = new Tuple<>(i, connections[i]);
        sortedMolsByConn.add(tup);
      }

      Collections.sort(sortedMolsByConn, new MolConnComparator());

      if (DEBUG) {
        int count = 0;
        for (final Tuple<Integer, Integer> tup : sortedMolsByConn) {
          System.out.println(
              "DEBUG: At " + count + " is " + tup.getObject1() + " with " + tup.getObject2());
          count++;
        }
      }

      int moveMol = -1;
      if (markMovedMolsUnmovable) {
        // System.out.println("Trying to find one that was unmarked...");
        // find the least connected one that has not been moved yet (so to avoid moving the same
        // molecule around and around)
        for (final Tuple<Integer, Integer> tup : sortedMolsByConn) {
          // System.out.println("Having " + tup.getObject1() + " with " + tup.getObject2());
          if (!movedMols.contains(
              (int) tup.getObject1())) { // downcast so that the object reference changes!
            // System.out.println("Found " + tup.getObject1());
            moveMol = tup.getObject1();
            movedMols.add(moveMol);
            break;
          }
        }
      } else {
        // simple, just the least connected
        moveMol = sortedMolsByConn.get(0).getObject1();
        movedMols.add(moveMol);
      }

      // find the move partner (something not so connected, no self-move)
      int movePartner = -1;
      for (final Tuple<Integer, Integer> tup : sortedMolsByConn) {
        if (tup.getObject1() != moveMol) {
          movePartner = tup.getObject1();
          break;
        }
      }

      assert (moveMol >= 0);
      assert (movePartner >= 0);

      if (DEBUG) {
        System.out.println(
            "DEBUG: Moving molecule " + moveMol + " with moving partner " + movePartner);
      }

      if (DEBUG) {

        final float[] colEnOne = new float[] {0.0f, 0.0f, 1.0f};
        final float[] colEnTwo = new float[] {1.0f, 0.0f, 0.0f};
        final float[] colTrOne = new float[] {0.0f, 0.0f, 1.0f};
        final float[] colTrTwo = new float[] {1.0f, 0.0f, 0.0f};
        final float transparency = 1.0f;
        final String backGround = "white";
        final String circleRing = "black";
        final String vertexShape = "circle";
        final String vertexStyle = "filled";
        final String edgeStyle = "filled";
        final int penWidth = 3;

        // dot file output :-)
        System.out.println("graph G {");
        System.out.println("   bgcolor = \"" + backGround + "\";");

        // first the vertices
        int maxConns = 0;
        for (int i = 0; i < connections.length; i++) {
          maxConns = Math.max(maxConns, connections[i] / 2);
        }
        final float[] colDef = new float[3];
        for (int i = 0; i < connections.length; i++) {
          // calculate our color in RGB space
          /*final double p = connections[i]/(2*maxConns);
          colDef[0] = (float) (colEnTwo[0] * p + colEnOne[0] * (1 - p));
          colDef[1] = (float) (colEnTwo[1] * p + colEnOne[1] * (1 - p));
          colDef[2] = (float) (colEnTwo[2] * p + colEnOne[2] * (1 - p));
          if(colDef[0] > 1.0f){colDef[0] = 1.0f;}
          if(colDef[1] > 1.0f){colDef[1] = 1.0f;}
          if(colDef[2] > 1.0f){colDef[2] = 1.0f;}

          final String hex = colorToHex(colDef, transparency);*/
          String hex;
          if (i == moveMol) {
            hex = "#FF0000";
          } else if (i == movePartner) {
            hex = "#003366";
          } else {
            hex = "#FFFFFF";
          }

          final String labelStr = "" + connections[i];
          System.out.println(
              "   "
                  + (i + 1)
                  + " [label=\""
                  + labelStr
                  + "\", shape = \""
                  + vertexShape
                  + "\", style = \""
                  + vertexStyle
                  + "\", color = \""
                  + circleRing
                  + "\", fillcolor = \""
                  + hex
                  + "\"];");
        }

        // now the edges
        for (int i = 0; i < connections.length - 1; i++) {
          for (int j = i + 1; j < connections.length; j++) {
            if (connCounter[i][j] <= 0) {
              continue;
            }

            /*final double p = connCounter[i][j]/(2*maxConns);
            colDef[0] = (float) (colTrTwo[0] * p + colTrOne[0] * (1 - p));
            colDef[1] = (float) (colTrTwo[1] * p + colTrOne[1] * (1 - p));
            colDef[2] = (float) (colTrTwo[2] * p + colTrOne[2] * (1 - p));
            if(colDef[0] > 1.0f){colDef[0] = 1.0f;}
            if(colDef[1] > 1.0f){colDef[1] = 1.0f;}
            if(colDef[2] > 1.0f){colDef[2] = 1.0f;}

            final String hex = colorToHex(colDef, transparency);*/
            final String hex = "#000000";

            System.out.println(
                "   "
                    + (i + 1)
                    + " -- "
                    + (j + 1)
                    + " [style = \""
                    + edgeStyle
                    + "\", color = \""
                    + hex
                    + "\", penwidth = \""
                    + penWidth
                    + "\"];");
          }
        }

        System.out.println("}");
      }

      // hand this info over to the point provider in exchange for a list of trial points
      Geometry bestFoundGeom = null;
      double bestFoundFitness = Double.MAX_VALUE;
      assert (moveMol >= 0);
      assert (movePartner >= 0);

      double zUpperEnv = -Double.MAX_VALUE;
      Environment env = (work.containsEnvironment()) ? work.getEnvironment() : null;
      if (work.containsEnvironment()
          && envAware
          && env.getEnvironmentType() == Environment.ENVIRONMENTTYPE.SURFACE) {
        // this is surface only code: figure out where the upper z-component of the surface is
        // assumes the surface to be in x-y plane and the cluster on top of it
        final CartesianCoordinates cW = work.getCartesiansWithEnvironment();
        final double[][] cartes = cW.getAllXYZCoord();
        final int atomsInEnv = env.atomsInEnv();
        final int clusterAtoms = work.getNumberOfAtoms();
        final int totalAtoms = cartes[0].length;
        assert (totalAtoms - atomsInEnv == clusterAtoms);
        int highestEnvAtom = -1;
        for (int atom = clusterAtoms; atom < totalAtoms; atom++) {
          if (cartes[2][atom] > zUpperEnv) {
            zUpperEnv = cartes[2][atom];
            highestEnvAtom = atom;
          }
        }

        // add the atomic radius of that environment atom and some for the moved mols atoms
        final short atomNoHighestEnv = cW.getAtomNumberOfAtom(highestEnvAtom);
        final double radiusEnvAtom = AtomicProperties.giveRadius(atomNoHighestEnv);
        zUpperEnv += radiusEnvAtom;

        final short[] atomNosMoved = work.getMoleculeAtPosition(moveMol).getAtomNumbers();
        double avgMolRadius = 0.0;
        for (int i = 0; i < atomNosMoved.length; i++) {
          avgMolRadius += AtomicProperties.giveRadius(atomNosMoved[i]);
        }
        avgMolRadius /= atomNosMoved.length;

        // fudge a bit
        avgMolRadius *= 1.1;

        zUpperEnv += avgMolRadius;
      }

      while (pointProv.hasNextPoint()) {
        if (DEBUG) {
          System.out.println("DEBUG: Working on new trial point.");
        }

        final int moved = pointProv.putNextPointIn(work, sortedMolsByConn, moveMol, movePartner);
        if (work.containsEnvironment()
            && envAware
            && env.getEnvironmentType() == Environment.ENVIRONMENTTYPE.SURFACE) {
          // check if COM of new point is in the surface, if so: do not bother
          final double zNew = work.getCOM(moved)[2];
          if (zNew <= zUpperEnv) {
            continue;
          }
        }

        if (doCDForEveryTrial) {
          final boolean hasColl = CollisionDetection.checkOnlyForCollision(work, blowColl);
          if (hasColl) {
            if (DEBUG) {
              System.out.println("DEBUG: Found a collision for this trial point.");
            }
            continue;
          }
        }

        // now so "something" whatever the optimization strategy is, really
        final Tuple<Double, Geometry> res = pointOpt.pointOpt(work, locopt, moved, movePartner);
        if (bestFoundGeom == null) {
          if (DEBUG) {
            System.out.println("DEBUG: Unconditionally making this best found.");
          }
          bestFoundGeom = res.getObject2().copy();
          bestFoundFitness = res.getObject1();
          continue;
        }
        if (res.getObject1() < bestFoundFitness) {
          if (DEBUG) {
            System.out.println("DEBUG: Conditionally making this best found.");
          }
          bestFoundGeom = res.getObject2().copy();
          bestFoundFitness = res.getObject1();
        }
      }

      pointProv.reset(); // reset so that we can get a new set of trial points :-)

      if (bestFoundGeom != null) {

        // well, just continue with the best we found :-)
        final Geometry mutationResult = bestFoundGeom.copy();
        mutationResult.setFitness(bestFoundFitness);

        work = mutationResult;
      }
      // continue the loop over molecules :-)
    }

    if (work.getFitness() >= (eBefore - 1E-8)) {
      System.out.println(
          "INFO: Advanced GDM tried hard but could not find a substantially better configuration than "
              + eBefore
              + " vs best found "
              + work.getFitness()
              + " continuing but complaining.");
      /*
       * the usual BXH disclaimer applies that this makes sense anyways because it adds new genome etc pp.
       */
    } else {
      final double percImprove = (eBefore - work.getFitness()) * 100 / eBefore;
      System.out.println(
          "INFO: Advanced GDM did find a better configuration: "
              + String.format("%6.2f", percImprove)
              + " % better. Initial: "
              + eBefore
              + " now: "
              + work.getFitness());
    }

    return work;
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

  static class GDMConfiguration implements Serializable {

    private static final long serialVersionUID = (long) 20140526;

    int noOfMolsToMove = 1;
    boolean markMovedMolsUnmoveable = true;
    boolean fullOptEachMove = false;
    double blowBonds = -1;
    double blowCD = -1;
    boolean doCDEveryTrial = true;
  }

  private static class MolConnComparator
      implements Comparator<Tuple<Integer, Integer>>, Serializable {

    private static final long serialVersionUID = (long) 20201230;

    @Override
    public int compare(final Tuple<Integer, Integer> t, final Tuple<Integer, Integer> t1) {
      final int conn1 = t.getObject2();
      final int conn2 = t1.getObject2();
      if (conn1 == conn2) {
        return 0;
      } else if (conn1 < conn2) {
        return -1;
      } else if (conn1 > conn2) {
        return 1;
      } else {
        throw new RuntimeException("GDM comperator: WTF?!");
      }
    }
  }
}
