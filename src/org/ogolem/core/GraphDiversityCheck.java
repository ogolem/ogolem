/**
Copyright (c) 2010     , J. M. Dieterich and B. Hartke
              2010-2012, J. M. Dieterich
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

import org.jgrapht.EdgeFactory;
import org.jgrapht.WeightedGraph;
import org.jgrapht.graph.ClassBasedEdgeFactory;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.WeightedPseudograph;

/**
 * Checks the diversity of two geometries by means of graph based operations.
 * @author Johannes Dieterich
 * @version 2012-06-24
 */
final class GraphDiversityCheck {

    //XXX: feature enhancement required
    private final boolean bUseSubgraphs = false;

    private final double dBlowDissoc;

    GraphDiversityCheck(final double blowFacDissoc){
        this.dBlowDissoc = blowFacDissoc;
    }

    boolean checkDiversity(final Geometry geom1, final Geometry geom2){

        // first assign graphs to both geometries
        final WeightedGraph<String,MyWeightedEdge> graph1 = constructGraph(geom1);
        final WeightedGraph<String,MyWeightedEdge> graph2 = constructGraph(geom2);

        // compare them


        return true;
    }

    WeightedGraph<String,MyWeightedEdge> constructGraph(final Geometry geom){

        final EdgeFactory<String,MyWeightedEdge> edgeFac =
                new ClassBasedEdgeFactory<>(MyWeightedEdge.class);
        final WeightedGraph<String,MyWeightedEdge> graph =
                new WeightedPseudograph<>(edgeFac);

        final CartesianCoordinates cartes = geom.getCartesians();
        final int[] iaAtsPerMol = cartes.getAllAtomsPerMol();
        final double[][] daXYZ = cartes.getAllXYZCoord();
        final String[] saAtoms = cartes.getAllAtomTypes();

        /*
         * populate the vertices
         */
        final String[] saVertices = new String[iaAtsPerMol.length];
        int iOffset = 0;
        int iEndset = 0;
        int[] iaAtomOffsets = new int[iaAtsPerMol.length];
        for(int i = 0; i < iaAtsPerMol.length; i++){
            iEndset += iaAtsPerMol[i];
            if(i != 0){
                iaAtomOffsets[i] = iaAtomOffsets[i-1] + iaAtsPerMol[i];
            }

            for(int j = iOffset; j < iEndset; j++){
                saVertices[i] += saAtoms[j];
            }
            graph.addVertex(saVertices[i]);
            iOffset += iaAtsPerMol[i];
        }

        /*
         * populate the edges
         */
        int iAtomOffset = 0;
        int iCurrMolecule = -1;
        String sFirstVertex = "";
        String sSecondVertex = "";
        for(int i = 0; i < saAtoms.length-iaAtsPerMol[iaAtsPerMol.length-1]; i++){
            final double dRadius1 = AtomicProperties.giveRadius(saAtoms[i]);

            if(i == iAtomOffset){
                // we are at the next "molecule"
                iCurrMolecule++;
                iAtomOffset += iaAtsPerMol[iCurrMolecule];
                sFirstVertex = saVertices[iCurrMolecule];
            }

            for(int j = iAtomOffset; j < saAtoms.length; j++){

                //TODO second vertex: fix
                //iaAtomsOffsets;

                final double dRadius2 = AtomicProperties.giveRadius(saAtoms[j]);
                final double dAddedRadii = dRadius1 + dRadius2;
                final double dAtomicDist = Math.sqrt(Math.pow((daXYZ[0][i]-daXYZ[0][j]),2)
                        + Math.pow((daXYZ[1][i]-daXYZ[1][j]),2) + Math.pow((daXYZ[2][i]-daXYZ[2][j]),2));

                if(dAtomicDist < (dBlowDissoc*dAddedRadii)){
                    // we have an edge
                    final MyWeightedEdge edge = edgeFac.createEdge(sFirstVertex, sSecondVertex);
                    edge.setWeight(dAtomicDist);

                    final boolean bAdded = graph.addEdge(sFirstVertex, sSecondVertex, edge);

                    assert bAdded;
                    // TODO probably, we have a problem with the redundancy of the edges
                    
                }
            }
        }

        return graph;
    }

    /*
     * how about making a building block the vertex and measure how many connections
     * each has with each. if we assign "molecular types" we can interchange
     * additionally, we should put weight (the distance) on the graph to be able to
     * check that as well
     * type them with String CHOCArHeNe whatever ;-)
     * next version should then allow for subgraphs (inside the molecule)
     */

     private class MyWeightedEdge extends DefaultEdge{

         private static final long serialVersionUID = (long) 2010727;

         protected double weight = WeightedGraph.DEFAULT_EDGE_WEIGHT;

         void setWeight(double weight){
             this.weight = weight;
         }

         double getWeight(){
             return weight;
         }
     }
}
