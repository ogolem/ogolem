/*
Copyright (c) 2012-2013, J. M. Dieterich and B. Hartke
              2015, J. M. Dieterich and B. Hartke
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

import static org.ogolem.core.Constants.ANGTOBOHR;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.ogolem.generic.GenericMutation;

/**
 * Directed mutation using knowledge of which atoms are on the surface of the cluster.
 *
 * @author Bernd Hartke
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
public class SurfaceDirectedMutation implements GenericMutation<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 20130923;
  private static final boolean DEBUG = false;

  private static final double MINDIFFOLDNEWSQ = 1.0 * ANGTOBOHR * ANGTOBOHR;
  private static final double SCALEFACMOLSIZE = 1.1;
  public static final int SIMPLESURFMODE = 0;
  public static final int TRIANGSURFMODE = 1;

  private final SurfaceDetectionEngine surfDetect;
  private final int surfMode;
  private final double blowColl;
  private final CartesianFullBackend back;
  private final CollisionDetectionEngine collDetect;

  SurfaceDirectedMutation(
      final SurfaceDetectionEngine surfDetect,
      final int surfMode,
      final double blowColl,
      final CartesianFullBackend back,
      final CollisionDetectionEngine collDetect) {
    this.surfDetect = surfDetect;
    this.back = back;
    this.blowColl = blowColl;
    this.surfMode = surfMode;
    this.collDetect = collDetect;
  }

  SurfaceDirectedMutation(final SurfaceDirectedMutation orig) {
    this.back = orig.back.copy();
    this.blowColl = orig.blowColl;
    this.collDetect = orig.collDetect.copy();
    this.surfDetect = orig.surfDetect.copy();
    this.surfMode = orig.surfMode;
  }

  @Override
  public SurfaceDirectedMutation copy() {
    return new SurfaceDirectedMutation(this);
  }

  @Override
  public String getMyID() {
    throw new UnsupportedOperationException(
        "Not supported yet."); // To change body of generated methods, choose Tools | Templates.
  }

  @Override
  public Geometry mutate(Geometry geom) {

    // please note that the backend probably will NOT support to evaluate anything but the *same*
    // system as always.
    // so no mixing, erasing of molecules, adding of molecules,...

    // nonsense message just to see if we get here at all...
    if (DEBUG) {
      System.out.println("DEBUG: hi there, this is SurfaceDirMut!");
    }

    final Geometry cluster = new Geometry(geom);
    final CartesianCoordinates c = cluster.getCartesians();
    final double[] energyParts = new double[c.getNoOfMolecules()];
    final double e =
        back.energyCalculation(
            cluster.getID(),
            -1,
            c.getAll1DCartes(),
            c.getAllAtomTypes(),
            c.getAllAtomNumbers(),
            c.getAllAtomsPerMol(),
            energyParts,
            c.getNoOfAtoms(),
            c.getAllCharges(),
            c.getAllSpins(),
            cluster.getBondInfo(),
            c.containedEnvType() == CartesianCoordinates.ENVTYPE.RIGID);
    if (DEBUG) {
      System.out.println("DEBUG: e is " + e);
    }
    // TODO expand DirectedMutation to do its stuff for n>1 molecules.
    // for the moment, just pretend that we only want to move one molecule...
    // (to move more, worstIndex could simply return an integer array of the N worst molecules,
    // and then the remainder could be looped over N -- with only some slight adjustments
    // in a few locations, hopefully ;-) ;-)
    // The more tricky considerations are: If the 1st molecule is moved, this should result
    // in a changed geometry, which would then effectively block the found vacancy for further
    // placement attempts (hopefully). But does it really make sense to proceed further w/o
    // re-calculating the energy contributions? And w/o doing a locopt? But if all that is done,
    // this definitely becomes too expensive for something that is done at random for rather
    // arbitrary pool members.)
    final int indexMolMove = worstIndex(geom, energyParts);
    if (DEBUG) {
      System.out.println("DEBUG: index of molecule with worst Econtrib: " + indexMolMove);
    }

    // remove this molecule from our work geometry
    final Molecule movMol = cluster.removeMolecule(indexMolMove);

    // here we save the "old" position of this molecule; this is used later to avoid
    // rediscovering it in the grid search for a new position -- which would render the
    // whole DirMut rather pointless... (on the other hand, one could think about NOT doing
    // this check, which would be equivalent to assuming that there _are_ better vacancies
    // than the current/old point, with the fallback option of rediscovering the old point
    // if there are _no_ better vacancies... Not sure if this really makes sense, though :-( )
    final double[] oldPoint = movMol.getExternalCenterOfMass().clone();

    // find the surface of the remaining cluster
    final Surface surf = surfDetect.detectTheSurface(cluster.getCartesians());

    // estimate the diameter of the molecule
    final double movMolDiam = estimateDiameter(movMol);

    // create a list of possible positions
    List<double[]> possNewPos;
    if (surfMode == SIMPLESURFMODE) {
      possNewPos =
          listOfPosOver(
              surf,
              cluster,
              movMol,
              indexMolMove,
              movMolDiam,
              collDetect,
              blowColl,
              cluster.getBondInfo(),
              oldPoint);
    } else if (surfMode == TRIANGSURFMODE) {
      possNewPos =
          listOfPosTriag(
              surf,
              cluster,
              movMol,
              indexMolMove,
              movMolDiam,
              collDetect,
              blowColl,
              cluster.getBondInfo(),
              oldPoint);
    } else {
      System.err.println(
          "ERROR: Surface mode " + surfMode + " unsupported. Returning non-mutated geometry.");
      return geom;
    }

    if (possNewPos.isEmpty()) {
      if (DEBUG) {
        System.out.println("DEBUG: No fitting new position found, returning previous geometry");
      }
      return new Geometry(geom);
    }

    // now we have the most promising vacancy => move the worst molecule there:
    // are the following two lines...
    final double[] bestPoint =
        findBestPoint(
            back, possNewPos, movMol, indexMolMove, cluster, energyParts, cluster.getBondInfo());

    movMol.setExternalCenterOfMass(bestPoint);
    cluster.addMolecule(movMol, indexMolMove);
    // equivalent to the next few lines (except for the issue of the Eulers, of course)
    //        final double[] best6Point = new double[6];
    //        best6Point[0] = bestPoint[0];
    //        best6Point[1] = bestPoint[1];
    //        best6Point[2] = bestPoint[2];
    //        best6Point[3] = 0.0; // here, some random Euler needed
    //        best6Point[4] = 0.0; // here, some random Euler needed
    //        best6Point[5] = 0.0; // here, some random Euler needed
    //        work.setExtCoordMolecule(best6Point,indexMolMove);
    // TODO: of course, this Euler stuff needs a _real_ solution: naively, one should try placements
    // with many different orientations for each vacancy test point, but obviously this would be
    // a big effort in computer time which would not really pay off.
    // The cheapest "solution" is to ignore the Eulers (what I am probably doing here),
    // the next less cheap "solution" would be to try a _very_ limited number of different Eulers
    // (which could be a standard set or a random set).

    return cluster;
    // REMARK please note that this *will* be afterwards optimized (if it passes CD/DD)
    // so no previous extra optimization if possible
  }

  private static int worstIndex(final Geometry geom, final double[] energyParts) {

    // returns index of molecule with worst (largest = "most positive") energy contribution to a
    // geometry
    // (actually, this could be a specific implementation of a far more general type which would
    // return
    // the element that is best/worst by any criterion, from any collection (or also the best/worst
    // N entries)
    int index = 0;
    for (int i = 0; i < energyParts.length; i++) {
      if (DEBUG) {
        System.out.println("DEBUG:" + i + ": " + energyParts[i]);
      }
      if (energyParts[i] > energyParts[index]) {
        index = i;
      }
    }

    if (DEBUG) {
      System.out.println("DEBUG: Worst index is " + index);
    }

    return index;
  }

  private static double[] findBestPoint(
      final CartesianFullBackend back,
      final List<double[]> alGoodPoints,
      final Molecule mol,
      final int molPos,
      final Geometry cluster,
      final double[] energyParts,
      final BondInfo bonds) {

    if (DEBUG) {
      String methodID = back.getMethodID();
      System.out.println("DEBUG: Method in DirectedMutation " + methodID);
    } // just for debugging, to verify backend ID

    final int numberOfAtoms = cluster.getNumberOfAtoms();
    double[] bestPoint = new double[3];
    double eBest = Double.POSITIVE_INFINITY;
    int iter = 0;
    for (final double[] aPoint : alGoodPoints) {
      if (DEBUG) {
        System.out.println("DEBUG: Testing point " + java.util.Arrays.toString(aPoint));
      }
      mol.setExternalCenterOfMass(aPoint);
      cluster.addMolecule(mol, molPos);
      final CartesianCoordinates c = cluster.getCartesians();
      final double[] xyz1D = c.getAll1DCartes();
      final double e =
          back.energyCalculation(
              cluster.getID(),
              iter,
              xyz1D,
              c.getAllAtomTypes(),
              c.getAllAtomNumbers(),
              c.getAllAtomsPerMol(),
              energyParts,
              numberOfAtoms,
              c.getAllCharges(),
              c.getAllSpins(),
              bonds,
              c.containedEnvType() == CartesianCoordinates.ENVTYPE.RIGID);
      if (DEBUG) {
        System.out.println("DEBUG:    has energy: " + e);
      }
      if (e < eBest) {
        eBest = e;
        bestPoint = aPoint;
      }
      cluster.removeMolecule(molPos);
      iter++;
    }
    if (DEBUG) {
      System.out.println("DEBUG: best point was: " + Arrays.toString(bestPoint));
    }
    if (DEBUG) {
      System.out.println("DEBUG:    with energy: " + eBest);
    }

    return bestPoint;
  }

  private static List<double[]> listOfPosTriag(
      final Surface surf,
      final Geometry cluster,
      final Molecule mol,
      final int molPos,
      final double molSize,
      final CollisionDetectionEngine collDetect,
      final double blowColl,
      final BondInfo bonds,
      final double[] oldPoint) {

    final int[] molIDs = surf.giveAllSurfaceMolecules();
    final double[][] extCOMs = cluster.getAllCOMs();
    final double cutoffSq = 2.5 * 2.5 * molSize * molSize;

    final List<double[]> possPoints = new ArrayList<>();
    for (final int i : molIDs) {

      for (final int j : molIDs) {

        if (i == j) {
          continue;
        }
        final double dXIJ = extCOMs[0][i] - extCOMs[0][j];
        final double dYIJ = extCOMs[1][i] - extCOMs[1][j];
        final double dZIJ = extCOMs[2][i] - extCOMs[2][j];
        final double distIJSq = dXIJ * dXIJ + dYIJ * dYIJ + dZIJ * dZIJ;
        if (distIJSq > cutoffSq) {
          continue;
        }

        for (final int k : molIDs) {
          if (i == k || j == k) {
            continue;
          }
          final double dXIK = extCOMs[0][i] - extCOMs[0][k];
          final double dYIK = extCOMs[1][i] - extCOMs[1][k];
          final double dZIK = extCOMs[2][i] - extCOMs[2][k];
          final double distIKSq = dXIK * dXIK + dYIK * dYIK + dZIK * dZIK;
          if (distIKSq > cutoffSq) {
            continue;
          }

          final double dXJK = extCOMs[0][j] - extCOMs[0][k];
          final double dYJK = extCOMs[1][j] - extCOMs[1][k];
          final double dZJK = extCOMs[2][j] - extCOMs[2][k];
          final double distJKSq = dXJK * dXJK + dYJK * dYJK + dZJK * dZJK;
          if (distJKSq > cutoffSq) {
            continue;
          }

          // ok, the triangle of the points i/j/k is not too big
          // find the center of this triangle
          final double[] center = new double[3];
          center[0] = (extCOMs[0][i] + extCOMs[0][j] + extCOMs[0][k]) / 3;
          center[1] = (extCOMs[1][i] + extCOMs[1][j] + extCOMs[1][k]) / 3;
          center[2] = (extCOMs[2][i] + extCOMs[2][j] + extCOMs[2][k]) / 3;

          // perhaps this already works?
          final double dX = center[0] - oldPoint[0];
          final double dY = center[1] - oldPoint[1];
          final double dZ = center[2] - oldPoint[2];
          final double distCOMCOMSq = dX * dX + dY * dY + dZ * dZ;
          if (distCOMCOMSq > MINDIFFOLDNEWSQ) {
            // far enough away from the old position
            final double[] comMol = mol.getExternalCenterOfMass();
            comMol[0] = center[0];
            comMol[1] = center[1];
            comMol[2] = center[2];
            cluster.addMolecule(mol, molPos);
            final CartesianCoordinates cart = cluster.getCartesians();

            final CollisionInfo collinfo = collDetect.checkForCollision(cart, blowColl, bonds);

            if (!collinfo.hasCollision()) {
              final double[] point = comMol.clone();
              if (DEBUG) {
                System.out.println("DEBUG: added nonColl point " + Arrays.toString(point));
              }
              possPoints.add(point);
              cluster.removeMolecule(molPos);
              continue; // continue so not to add another point further away from the surface!
            }
            cluster.removeMolecule(molPos);
          }

          // try a second time further away from the surface (at least, if the triangle points
          // outwards)
          final double dist =
              Math.sqrt(center[0] * center[0] + center[1] * center[1] + center[2] * center[2]);
          final double idealDist = dist + SCALEFACMOLSIZE * molSize;
          final double scale = idealDist / dist;

          final double[] comMol = mol.getExternalCenterOfMass();
          comMol[0] = center[0] * scale;
          comMol[1] = center[1] * scale;
          comMol[2] = center[2] * scale;

          final double dX2 = comMol[0] - oldPoint[0];
          final double dY2 = comMol[1] - oldPoint[1];
          final double dZ2 = comMol[2] - oldPoint[2];
          final double distCOMCOMSq2 = dX2 * dX2 + dY2 * dY2 + dZ2 * dZ2;
          if (distCOMCOMSq2 > MINDIFFOLDNEWSQ) {
            // far enough away from the old position
            cluster.addMolecule(mol, molPos);
            final CartesianCoordinates cart = cluster.getCartesians();

            final CollisionInfo collinfo = collDetect.checkForCollision(cart, blowColl, bonds);

            if (!collinfo.hasCollision()) {
              final double[] point = comMol.clone();
              if (DEBUG) {
                System.out.println("DEBUG: added nonColl point " + Arrays.toString(point));
              }
              possPoints.add(point);
            }
            cluster.removeMolecule(molPos);
          }
        }
      }
    }

    return possPoints;
  }

  private static List<double[]> listOfPosOver(
      final Surface surf,
      final Geometry cluster,
      final Molecule mol,
      final int molPos,
      final double molSize,
      final CollisionDetectionEngine collDetect,
      final double blowColl,
      final BondInfo bonds,
      final double[] oldPoint) {

    final int[] molIDs = surf.giveAllSurfaceMolecules();

    final List<double[]> possPoints = new ArrayList<>();
    for (final int molID : molIDs) {

      if (DEBUG) {
        System.out.println("DEBUG: Working on surface molecule " + molID);
      }

      // get the COM of this surface mol and scale it to be the new COM of the to be set molecule
      final double[] com = cluster.getCOM(molID);
      final double[] comMol = mol.getExternalCenterOfMass();

      final double dist = Math.sqrt(com[0] * com[0] + com[1] * com[1] + com[2] * com[2]);
      final double idealDist = dist + SCALEFACMOLSIZE * molSize;
      final double scale = idealDist / dist;

      comMol[0] = com[0] * scale;
      comMol[1] = com[1] * scale;
      comMol[2] = com[2] * scale;

      final double dX = comMol[0] - oldPoint[0];
      final double dY = comMol[1] - oldPoint[1];
      final double dZ = comMol[2] - oldPoint[2];
      final double distCOMCOMSq = dX * dX + dY * dY + dZ * dZ;
      if (distCOMCOMSq > MINDIFFOLDNEWSQ) {
        // far enough away from the old position
        cluster.addMolecule(mol, molPos);
        final CartesianCoordinates cart = cluster.getCartesians();

        final CollisionInfo collinfo = collDetect.checkForCollision(cart, blowColl, bonds);

        if (!collinfo.hasCollision()) {
          final double[] point = comMol.clone();
          if (DEBUG) {
            System.out.println("DEBUG: added nonColl point " + Arrays.toString(point));
          }
          possPoints.add(point);
        }
        cluster.removeMolecule(molPos);
      }
    }

    if (DEBUG) {
      System.out.println("DEBUG: number of collision-free positions: " + possPoints.size());
    }

    return possPoints;
  }

  private static double estimateDiameter(final Molecule mol) {

    final CartesianCoordinates molC = mol.getCartesians();
    molC.moveCoordsToCOM();
    final double[][] xyz = molC.getAllXYZCoord();

    // measure the maximal distance to COM
    double maxDistSq = 0.0;
    for (int i = 0; i < molC.getNoOfAtoms(); i++) {

      final double distSq = xyz[0][i] * xyz[0][i] + xyz[1][i] * xyz[1][i] + xyz[2][i] * xyz[2][i];
      if (distSq > maxDistSq) {
        maxDistSq = distSq;
      }
    }

    if (DEBUG) {
      System.out.println(
          "DEBUG: Estimating to be moved molecule to have a diameter of " + Math.sqrt(maxDistSq));
    }

    return Math.sqrt(maxDistSq) * 2;
  }
}
