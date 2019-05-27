/**
Copyright (c) 2009, J. M. Dieterich and B. Hartke
              2010-2013, J. M. Dieterich
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

import java.util.ArrayList;

/**
 * This calculates distances between points in space using different
 * (and also different scaling!) algorithms.
 * @author Johannes Dieterich
 * @version 2013-12-04
 */
public class DistanceCalc{
    
    private static final boolean DEBUG = false;
    
    /**
     * This follows Edgar Dijkstra's algorithm to calculate the length
     * of the shortest path to each reachable destination starting from
     * a given source. Scales with O**2.
     * For a reference implementation see, e.g.:
     * Timothy Budd, Classic Data Structures in Java, 2001
     * @param G The WeightedGraph.
     * @param s Starting atom.
     * @return The shortest path from the specified atom to all others.
     */
    static int[] dijkstraDistance(WeightedGraph G, int s) {
        // shortest known distance from "s"
        final double[] daDist = new double[G.size()];

        // preceeding node in path
        final int[] iaPredNodes = new int[G.size()];

        // all false in the beginning since nothing has been visited yet
        final boolean[] baVisited = new boolean[G.size()];

        // assign maximal values to all but the "self" distance
        for (int i = 0; i < daDist.length; i++) {
            daDist[i] = Double.MAX_VALUE;
        }
        daDist[s] = 0;

        for (int i = 0; i < daDist.length; i++) {
            final int iNext = minVertex(daDist, baVisited);
            baVisited[iNext] = true;

            // the shortest path to next is dist[next] and via pred[next].

            final int[] n = G.neighbors(iNext);
            for (int j = 0; j < n.length; j++) {
                final int v = n[j];
                final double d = daDist[iNext] + G.getWeight(iNext, v);
                if (daDist[v] > d) {
                    daDist[v] = d;
                    iaPredNodes[v] = iNext;
                }
            }
        }

        // pred[s]==0 does not contain anything useful
        return iaPredNodes;
    }

    private static int minVertex(double[] dist, boolean[] v) {
        double x = Double.MAX_VALUE;
        int y = -1;   // graph not connected, or no unvisited vertices
        for (int i = 0; i < dist.length; i++) {
            if (!v[i] && dist[i] < x) {
                y = i;
                x = dist[i];
            }
        }
        return y;
    }

    static void printPath(WeightedGraph G, int[] pred, int s, int e) {
        final ArrayList<String> path = new ArrayList<>();
        int x = e;
        while (x != s) {
            path.add(0, G.getLabel(x));
            x = pred[x];
        }
        path.add(0, G.getLabel(s));
        System.out.println(path);
    }

    /**
     * Compares two paths.
     * @param G1 The first weighted graph for comparision.
     * @param G2 The second weighted graph for comparision.
     * @param pred1 Shortest Dijkstra path for G1.
     * @param pred2 Shortest Dijkstra path for G2.
     * @return true if they are the same.
     */
    static boolean compareTwoPaths(WeightedGraph G1, WeightedGraph G2, int[] pred1, int[] pred2) {
        boolean bSame = true;
        for (int i = 0; i < pred1.length; i++) {
            if (pred1[i] != pred2[i]) {
                bSame = false;
            }
        }
        if (bSame) {
            // paths are completely identical
            return bSame;
        }

        // check further if they are perhaps just different numbered
        // TODO implement better comparision
        return bSame;
    }

    /**
     * Follows Warshalls Algorithm to solve the all-pair reachability problem
     * using the adjacency-matrix representation of a graph.
     * Converts the adjacency matrix into a reachability matrix.
     * Scales with O**3.
     *
     * For a reference implementation and explanation see, e.g.:
     * Timothy Budd, Classic Data Structures in Java, 2001
     * @param adjacency The adjacency matrix.
     * @return The matrix transformed to a reachability matrix.
     */
    static int[][] warshallDistances(final int[][] adjacency) {
        final int length = adjacency.length;
        for (int k = 0; k < length; k++) {
            for (int i = 0; i < length; i++) {
                for (int j = 0; j < length; j++) {
                    /*
                     * if there is a path from vertex i to vertex j that does
                     * not go through any vertex higher than k then value is 1.
                     * |= : bitwise OR and ASSIGN
                     * & : bitwise AND
                     */
                    adjacency[i][j] |= adjacency[i][k] & adjacency[k][j];
                }
            }
        }
        // adjacency is now actually the reachability matrix
        return adjacency;
    }

    /**
     * Based on DFS with a loop check, this calculates a reachability.
     * @param adjacencyMatrix an adjacency matrix of the graph
     * @return true if everything is reachable from the first vertex
     */
    public static boolean dfsReachability(final boolean[][] adjacencyMatrix) {
        return dfsReachability(adjacencyMatrix, adjacencyMatrix.length);
    }
        
    /**
     * Based on DFS with a loop check, this calculates a reachability.
     * @param adjacencyMatrix an adjacency matrix of the graph
     * @param noNodes number of nodes to consider in the adjacency matrix.
     * @return true if everything is reachable from the first vertex
     */
    public static boolean dfsReachability(final boolean[][] adjacencyMatrix, final int noNodes) {

        // first all marked 0 by allocation
        final int[] states = new int[noNodes];
        runRecursiveDFS(0, states, adjacencyMatrix,noNodes);

        /*
         * now our states should (in case that they are all reachable from the first one) be marked as
         * finished (2): check this
         */
        for (int i = 0; i < noNodes; i++) {
            final int state = states[i];
            if (state != 2){
                if(DEBUG) System.err.println("DEBUG: Unconnected for " + i + " is " + state);
                return false;
            }
        }

        // if we end up here, everything is OK and everything is reachable from the 0'th vertex
        return true;
    }

    private static void runRecursiveDFS(final int state, final int[] states,
            final boolean[][] adjacencyMatrix, final int noNodes) {

        // 1: visited, not yet finished
        states[state] = 1;

        for (int i = 0; i < noNodes; i++) {
            // call recursively if not yet visited and connected
            if (adjacencyMatrix[state][i] && states[i] == 0) {
                runRecursiveDFS(i, states, adjacencyMatrix,noNodes);
            }
        }

        // 2: now we are done with this vertex
        states[state] = 2;
    }
}
