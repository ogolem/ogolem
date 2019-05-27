/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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

/**
 * Dissociation detection within in a cluster of molecules.
 * @author Johannes Dieterich
 * @version 2015-04-28
 */
public class DissociationDetection {
    
    public static enum DDTYPE{DFS,WARSHALL};
    public static final DDTYPE DEFAULTDD = DDTYPE.DFS;
    
    /**
     * @param pairWiseDists The pairwise distances between two atoms.
     * @param saAtoms An array of Strings containing the atomic IDs.
     * @param nos the atomic numbers
     * @param dissBlow the blow factor for the atomic radii
     * @param whichDetection which detection algorithm
     * @return True if dissociation has been detected.
     */
    public static boolean checkForDissociation(final double[][] pairWiseDists,
            final String[] saAtoms, final short[] nos, final double dissBlow,
            final DDTYPE whichDetection) {
        
        if(whichDetection == DDTYPE.DFS){
            // use DFS
            final int noOfAtoms = nos.length;

            // prepare an adjacency matrix
            final boolean[][] adjacency = new boolean[noOfAtoms][noOfAtoms];

            for (int i = 0; i < noOfAtoms; i++) {
                final double rad1 = AtomicProperties.giveRadius(nos[i]);
                for (int j = i; j < noOfAtoms; j++) { // on purpose from i as this sets the diagonal to true
                    final double rad2 = AtomicProperties.giveRadius(nos[j]);
                    final double summedRadii = (rad1 + rad2) * dissBlow;
                    if (summedRadii >= pairWiseDists[i][j]) {
                        adjacency[i][j] = true;
                        adjacency[j][i] = true;
                    } else{
                        adjacency[i][j] = false;
                        adjacency[j][i] = false;
                    }
                }
            }

            // call DFS
            final boolean isConnected = DistanceCalc.dfsReachability(adjacency);

            return !isConnected;

        } else if (whichDetection == DDTYPE.WARSHALL) {
            // use Warshall

            final int noOfAtoms = nos.length;

            // the adjacency matrix
            final int[][] adjacency = new int[noOfAtoms][noOfAtoms];

            for (int i = 0; i < noOfAtoms; i++) {
                final double rad1 = AtomicProperties.giveRadius(nos[i]);
                for (int j = i; j < noOfAtoms; j++) {
                    final double rad2 = AtomicProperties.giveRadius(nos[j]);
                    final double summedRadii = (rad1 + rad2) * dissBlow;
                    if (summedRadii >= pairWiseDists[i][j]) {
                        // atoms within reach
                        adjacency[i][j] = 1;
                        adjacency[j][i] = 1;
                    }
                }
            }

            // get the reachability matrix
            final int[][] reachabilities = DistanceCalc.warshallDistances(adjacency);
            final int length = reachabilities.length;
            /*
             * loop over the elements to check for zeros
             * since our bonds are not "directed"and the matrix is symmetric,
             * we can assume that it is sufficient to check just the upper half
             * of it.
             */
            for (int i = 0; i < length - 1; i++) {
                for (int j = i + 1; j < length; j++) {
                    if (reachabilities[i][j] == 0) {
                        // dissociation found
                        return true;
                    }
                }
            }
            return false;
        } else {
            System.err.println("ERROR: Can't handle dissociation detection " + whichDetection + " assuming everything to be fine now.");
            return false;
        }
    }

    /* TODO we would like a non O(N^2) scaling DD. this should be possible
     * considering a grid-like approach putting molecules into cells and
     * trying to check whether there are cells w/o molecules. Problem is:
     * the cells need to be big enough to yield a performance bettering and
     * at the same time small enough to not render the DD useless. This
     * means that there needs to be user input and the Warshall as a fallback
     * alternative. Meaning this needs to be turned into a gateway again.
     * I love them. :-)
     * So first really get the gridlike CD to work properly and then use
     * the information gained there for this. :-)
     */
    
    public static DDTYPE parseType(final String type) throws Exception {
        
        if(type.equalsIgnoreCase("dfs")){
            return DDTYPE.DFS;
        } else if(type.equalsIgnoreCase("warshall")){
            return DDTYPE.WARSHALL;
        }
        
        throw new Exception("Illegal dissociation detection " + type + ".");
    }
}
