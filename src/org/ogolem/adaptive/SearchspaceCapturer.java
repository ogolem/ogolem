/**
Copyright (c) 2011-2012, J. M. Dieterich
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
package org.ogolem.adaptive;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;

/**
 *
 * @author Johannes Dieterich
 * @version 2012-06-24
 */
final class SearchspaceCapturer {

    static final int GATYPE = 1;
    static final int MCMTYPE = 2;
    static final int HITTYPE = 3;

    private static final int[] grid = {150,150};
    private static final int dumpEach = 10;
    private static final LinkedList<SearchspacePoint> pointsList = new LinkedList<>();

    static synchronized void addPoint(final double[] point, final int type, final double val){
        final SearchspacePoint p = new SearchspacePoint(point, type, val);
        pointsList.add(p);
    }

    static void dumpData(final double[][] borders, final double accFit) throws Exception{
        
        if(borders[0].length > 2 || grid.length > 2) throw new UnsupportedOperationException("Sorry, no more then 2D ATM please.");

        // allocate the matrices
        final int[][] matrix = new int[grid[0]+2][grid[1]+2];
        final int[][] mini = new int[3][3];

        int which = 0;
        int count = 0;
        int count2 = 0;
        double bestFit = Double.MAX_VALUE;
        for(final SearchspacePoint p : pointsList){
            which++;
            final double[] po = p.getPoint();
            final int ty = p.getType();

            if(bestFit <= accFit) continue;

            bestFit = Math.min(p.getValue(), bestFit);
            boolean plotSurely;
            if(bestFit == p.getType()) plotSurely = true;
            else plotSurely = false;

            // bin
            int x;
            int y;
            if(po[0] < borders[0][0]){
                x = 1;
            } else if(po[0] > borders[1][0]){
                x = grid[0];
            } else{
                x = (int) Math.round((po[0]-borders[0][0])/(borders[1][0]-borders[0][0])*grid[0]) + 1;
            }
            if(po[1] < borders[0][1]){
                y = 1;
            } else if(po[1] > borders[1][1]){
                y = grid[1];
            } else{
                y = (int) Math.round((po[1]-borders[0][1])/(borders[1][1]-borders[0][1])*grid[1]) + 1;
            }

            if(x < 1) x = 1;
            else if(x > grid[0]) x = grid[0];

            if(y < 1) y = 1;
            else if(y > grid[1]) y = grid[1];

            // save & mark
            mini[0][0] = matrix[x-1][y-1];
            mini[1][0] = matrix[x  ][y-1];
            mini[2][0] = matrix[x+1][y-1];
            mini[0][1] = matrix[x-1][y  ];
            mini[1][1] = matrix[x  ][y  ];
            mini[2][1] = matrix[x+1][y  ];
            mini[0][2] = matrix[x-1][y+1];
            mini[1][2] = matrix[x  ][y+1];
            mini[2][2] = matrix[x+1][y+1];

            matrix[x-1][y-1] = HITTYPE;
            matrix[x  ][y-1] = HITTYPE;
            matrix[x+1][y-1] = HITTYPE;
            matrix[x-1][y  ] = HITTYPE;
            matrix[x  ][y  ] = HITTYPE;
            matrix[x+1][y  ] = HITTYPE;
            matrix[x-1][y+1] = HITTYPE;
            matrix[x  ][y+1] = HITTYPE;
            matrix[x+1][y+1] = HITTYPE;


            // dump matrix?
            if(count2 >= dumpEach || which == pointsList.size()-1 || plotSurely){
                dumpMatrix(count, matrix);
                count++;
                count2 = 0;
            } else{
                count2++;
            }

            // normalize
            matrix[x-1][y-1] = mini[0][0];
            matrix[x  ][y-1] = mini[1][0];
            matrix[x+1][y-1] = mini[2][0];
            matrix[x-1][y  ] = mini[0][1];
            if(mini[1][1] == 0) matrix[x  ][y  ] = ty;
            else matrix[x  ][y  ] = mini[1][1];
            matrix[x+1][y  ] = mini[2][1];
            matrix[x-1][y+1] = mini[0][2];
            matrix[x  ][y+1] = mini[1][2];
            matrix[x+1][y+1] = mini[2][2];
        }

    }

    //TODO hook up, test and plotting stuff

    private static class SearchspacePoint{

        private final double[] p;
        private final int t;
        private final double v;

        SearchspacePoint(final double[] point, final int type, final double val){
            this.p = point;
            this.t = type;
            this.v = val;
        }

        int getType(){
            return t;
        }

        double[] getPoint(){
            return p;
        }

        double getValue(){
            return v;
        }
    }

    private static void dumpMatrix(int count, int[][] matrix){

        final String path = "matrix" + count + ".mat";

        final String sep = System.getProperty("line.separator");
        BufferedWriter buffwriter = null;
        try {
            buffwriter = new BufferedWriter(new FileWriter(path, false));
            for (int i = 1; i < matrix.length-1; i++) {
                String line = "";
                for(int j = 1; j < matrix[i].length-1; j++){
                    line += "  " + matrix[i][j];
                }
                buffwriter.write(line);
                buffwriter.write(sep);
            }
            buffwriter.close();
        } catch (IOException e) {
            e.printStackTrace(System.err);
        } finally {
            if (buffwriter != null) {
                try {
                    buffwriter.close();
                } catch (IOException e) {
                    e.printStackTrace(System.err);
                }
            }
        }
    }
}
