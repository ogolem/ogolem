/**
Copyright (c) 2009, J. M. Dieterich and B. Hartke
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
 * Heavily adapted for our particular problem.
 * @author Johannes Dieterich
 * @version 2009-05-20
 */
class WeightedGraph {

    // the adjacency matrix
    private final double[][] daEdges;

    private final String[] saLabels;

    private final int iNoOfAtoms;


    /**
     * This is the first constructor for the WeightedGraph. It requires pairwise
     * distances for the specific cartesian coordinates.
     * @param cartes The cartesian coordinates we want to create the WeightedGraph from.
     * @param dBlowDiss A blow factor for dissociation detection.
     * @param daPairwiseDists The pairwise distances between all atoms.
     */
    WeightedGraph (CartesianCoordinates cartes, double dBlowDiss, double[][] daPairwiseDists){
        this.iNoOfAtoms = cartes.getNoOfAtoms();
        this.saLabels = cartes.getAllAtomTypes();
        this.daEdges = daPairwiseDists;

        // now assign negative values to all non-connected atoms

        for(int i = 0; i < daEdges.length - 1; i++){
            for(int j = i+1; j < daEdges.length; j++){
                double dRadiiAdded = dBlowDiss* (AtomicProperties.giveRadius(saLabels[i])
                        + AtomicProperties.giveRadius(saLabels[j]));
                if(daEdges[i][j] > dRadiiAdded){
                    daEdges[i][j] = -1;
                }
                daEdges[j][i] = daEdges[i][j];
            }
        }
    }

    /**
     * This is the second constructor for the WeightedGraph. It does not require
     * pairwise distances but instead calculates them itself.
     * @param cartes The cartesian coordinates we want to create the WeightedGraph from.
     * @param dBlowDiss A blow factor for dissociation detection.
     */
    WeightedGraph (CartesianCoordinates cartes, double dBlowDiss){
        this.iNoOfAtoms = cartes.getNoOfAtoms();
        this.saLabels = cartes.getAllAtomTypes();
        this.daEdges = setEdgesUp(cartes.getAllXYZCoord(), dBlowDiss);
    }

    int size() {
        return saLabels.length;
    }

    void setLabel(int iVertex, String sLabel) {
        saLabels[iVertex] = sLabel;
    }

    String getLabel(int iVertex) {
        return saLabels[iVertex];
    }

    void addEdge(int iSource, int iTarget, int iWeight) {
        daEdges[iSource][iTarget] = iWeight;
    }

    boolean isEdge(int iSource, int iTarget) {
        return daEdges[iSource][iTarget] > 0;
    }

    void removeEdge(int iSource, int iTarget) {
        daEdges[iSource][iTarget] = 0;
    }

    double getWeight(int iSource, int iTarget) {
        return daEdges[iSource][iTarget];
    }

    int[] neighbors(int iVertex) {
        int iCounter = 0;
        for (int i = 0; i < daEdges[iVertex].length; i++) {
            if (daEdges[iVertex][i] > 0) {
                iCounter++;
            }
        }
        final int[] iaAnswer = new int[iCounter];
        iCounter = 0;
        for (int i = 0; i < daEdges[iVertex].length; i++) {
            if (daEdges[iVertex][i] > 0) {
                iaAnswer[iCounter++] = i;
            }
        }
        return iaAnswer;
    }

    void printToSystemOut() {
        for (int j = 0; j < daEdges.length; j++) {
            System.out.print(saLabels[j] + ": ");
            for (int i = 0; i < daEdges[j].length; i++) {
                if (daEdges[j][i] > 0) {
                    System.out.print(saLabels[i] + ":" + daEdges[j][i] + " ");
                }
            }
            System.out.println();
        }
    }

    /**
     * Setting up the edges as an adjacency matrix and putting negative values
     * into all non-connected pairs.
     * @param daCartesianCoordinates The cartesian coordinates.
     * @param dBlowDiss A blow factor for dissociation detection.
     * @return The edges as an adjacency matrix.
     */
    private double[][] setEdgesUp(double[][] daCartesianCoordinates, double dBlowDiss){
        double[][] daSetUpEdges = new double[daCartesianCoordinates[0].length]
                [daCartesianCoordinates[0].length];

        for(int i = 0; i < daCartesianCoordinates[0].length-1; i++){
            for(int j = i + 1; j < daCartesianCoordinates[0].length; j++){
                // calculate the value of the edge
                daSetUpEdges[i][j] = Math.sqrt(
                        Math.pow((daCartesianCoordinates[0][j] - daCartesianCoordinates[0][i]), 2.0)
                        + Math.pow((daCartesianCoordinates[1][j] - daCartesianCoordinates[1][i]), 2.0)
                        + Math.pow((daCartesianCoordinates[2][j] - daCartesianCoordinates[2][i]), 2.0));
                double dRadiiAdded = dBlowDiss* (AtomicProperties.giveRadius(saLabels[i])
                        + AtomicProperties.giveRadius(saLabels[j]));
                if(daSetUpEdges[i][j] > dRadiiAdded){
                    daSetUpEdges[i][j] = -1;
                }
                daSetUpEdges[j][i] = daSetUpEdges[i][j];
            }
        }
        
        return daSetUpEdges;
    }
}
