/**
Copyright (c) 2013, J. M. Dieterich
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

import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;

/**
 * Based off Bernd's old LJ GA paper and a private communication. Project
 * the cluster into xy. Although not limited to LJ clusters (obviously...)
 * most likely not useful for much else (according to BXH).
 * @author Johannes Dieterich
 * @author Bernd Hartke
 * @version 2015-05-27
 */
class LJProjNicheComp implements NicheComputer<Molecule,Geometry> {
    
    private static final long serialVersionUID = (long) 20140107;
    public static final double defXIncr = Math.toRadians(10.0);
    public static final double defYIncr = Math.toRadians(10.0);
    public static final double defGridIncr = 2.0*Constants.ANGTOBOHR;
    public static final int defBinStart = 0;
    public static final int defBinWidth = 1;

    private final double[][] scr = new double[3][3];
    private final double xIncr;
    private final double yIncr;
    private final double gridIncr;
    private final int binStart;
    private final int binWidth;
    
    LJProjNicheComp(final double xIncr, final double yIncr, final double gridIncr, final int binStart, int binWidth){
        
        if(Math.abs(xIncr-yIncr)>1E-6){
            System.err.println("WARNING: LJProjNicheComp: different xIncr and yIncr: " + xIncr + " vs " + yIncr + ". Allowed but discouraged.");
        }
        
        if(binWidth<1){
            System.err.println("WARNING: LJProjNicheComp: binWidth = " + binWidth + " makes no sense. Assuming binWidth = 1");
            binWidth = 1;
        }
        
        this.xIncr = xIncr;
        this.yIncr = yIncr;
        this.gridIncr = gridIncr;
        this.binStart = binStart;
        this.binWidth = binWidth;
    }
    
    private LJProjNicheComp(final LJProjNicheComp orig){
        this.gridIncr = orig.gridIncr;
        this.xIncr = orig.xIncr;
        this.yIncr = orig.yIncr;
        this.binStart = orig.binStart;
        this.binWidth = orig.binWidth;
    }
    
    @Override
    public LJProjNicheComp clone(){
        return new LJProjNicheComp(this);
    }
    
    /**
     * Projects the cluster into xy and extracts a niche out of a couple of attempts.
     * @param g
     * @return the niche this geometry can be associated to 
     */
    @Override
    public Niche computeNiche(final Geometry g){
        
        final CartesianCoordinates c = g.getCartesians();
        c.moveCoordsToCOM();
        final int noAtoms = c.getNoOfAtoms();
        final double[][] xyz = c.getAllXYZCoordsCopy();
        final double[][] xyzSave = new double[3][noAtoms];
        final double[][] rotated = new double[3][noAtoms];
        
        int minNonEmpty = Integer.MAX_VALUE;
        int maxPerCellForMin = -1;
        
        int maxPerCell = 0;
        int nonEmptyForMax = -1;
        
        int[][] binsCache = new int[100][100]; // good starting allocation.
        final int[] currBinCacheSize = new int[]{100,100};
        final int[] thisBinSize = new int[]{0,0};
        double xRot = 0.0;
        while(xRot <= (Math.PI+0.00001)){
            
            // save this state
            System.arraycopy(xyz[0], 0, xyzSave[0], 0, noAtoms);
            System.arraycopy(xyz[1], 0, xyzSave[1], 0, noAtoms);
            System.arraycopy(xyz[2], 0, xyzSave[2], 0, noAtoms);
            
            double yRot = 0.0;
            while(yRot <= (Math.PI+0.00001)){
                
                // check grid for this rotation
                doNiching(xyz,gridIncr,binsCache,currBinCacheSize,thisBinSize);
                final int numberOfBins1 = thisBinSize[0];
                final int numberOfBins2 = thisBinSize[1];
                final int numberOfBins = numberOfBins1 * numberOfBins2;
                
                int nonEmpty = 0;
                int maxPerBin = 0;
                for (int i = 0; i < thisBinSize[0]; i++) {
                    for (int j = 0; j < thisBinSize[1]; j++) {
                        if (binsCache[i][j] > 0) {
                            nonEmpty++;
                            maxPerBin = Math.max(maxPerBin, binsCache[i][j]);
                        }
                    }
                }
                
                // according to BXH: 1) find the orientation that has the minimum number
                // of occupied bins and store the max number of atoms in one bin
                // 2) find the orientation that has the maximum number of atoms in
                // one bin and store the number of occupied bins for this
                
                nonEmpty = (int) Math.round(nonEmpty / (double) numberOfBins * 100.0);
                
                if(nonEmpty < minNonEmpty){
                    minNonEmpty = nonEmpty;
                    maxPerCellForMin = maxPerBin;
                }
                
                if(maxPerBin > maxPerCell){
                    maxPerCell = maxPerBin;
                    nonEmptyForMax = nonEmpty;
                }
                
                // rotate by the incr one further
                CoordTranslation.rotateXYZAroundX(xyz, yIncr, rotated, scr);
                
                // copy rotated over
                System.arraycopy(rotated[0], 0, xyz[0], 0, noAtoms);
                System.arraycopy(rotated[1], 0, xyz[1], 0, noAtoms);
                System.arraycopy(rotated[2], 0, xyz[2], 0, noAtoms);
                
                yRot += yIncr;
            }
            
            // rotate one incr further (using the stored xyz coordinates)
            CoordTranslation.rotateXYZAroundX(xyzSave, xIncr, rotated, scr);
            
            // copy rotated over
            System.arraycopy(rotated[0], 0, xyz[0], 0, noAtoms);
            System.arraycopy(rotated[1], 0, xyz[1], 0, noAtoms);
            System.arraycopy(rotated[2], 0, xyz[2], 0, noAtoms);
            
            xRot += xIncr;
        }
        
        
        //return new Niche("ljprojniche: " + minNonEmpty + "nonempty" + maxPerCellForMin + "/" + maxPerCell + "maxpercell" + nonEmptyForMax);
        nonEmptyForMax -= binStart; // shift to start of niche bins
        final int nicheNumber = nonEmptyForMax / binWidth; // intentional integer division, to get niche bin number
        final int nicheLow = binStart + nicheNumber * binWidth;
        final int nicheHigh = nicheLow + binWidth;
        
        return new Niche("ljprojniche: " + nicheLow + "-" + nicheHigh);
    }
    
    private static void doNiching(final double[][] xyz, final double gridIncr, int[][] binsCache,
            final int[] currCacheSize, final int[] noXY){
        
        // figure out how big our grid should be
        double maxX = Double.NEGATIVE_INFINITY;
        double minX = Double.POSITIVE_INFINITY;
        double maxY = Double.NEGATIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY;
        for(int i = 0; i < xyz[0].length; i++){
            maxX = Math.max(maxX, xyz[0][i]);
            minX = Math.min(minX, xyz[0][i]);
        }
        for(int i = 0; i < xyz[0].length; i++){
            maxY = Math.max(maxY, xyz[1][i]);
            minY = Math.min(minY, xyz[1][i]);
        }
        
        final int noX = (int) Math.ceil((maxX-minX)/gridIncr);
        final int noY = (int) Math.ceil((maxY-minY)/gridIncr);
        
        // start the binning ("project in xy plane" -> remove z component)
        noXY[0] = noX;
        noXY[1] = noY;
        if(noX <= currCacheSize[0] && noY <= currCacheSize[1]){
            // zero
            for(int x = 0; x < noX; x++){
                for(int y = 0; y < noY; y++){
                    binsCache[x][y] = 0;
                }
            }
        } else {
            final int allocX = Math.max(currCacheSize[0],noX);
            final int allocY = Math.max(currCacheSize[1],noY);
            currCacheSize[0] = allocX;
            currCacheSize[1] = allocY;
            // reallocate
            binsCache = new int[noX][noY];
        }
        
        for(int i = 0; i < xyz[0].length; i++){
            
            final int binX = (int) Math.floor((xyz[0][i]-minX)/gridIncr);
            final int binY = (int) Math.floor((xyz[1][i]-minY)/gridIncr);
            binsCache[binX][binY]++;
        }
    }
}
