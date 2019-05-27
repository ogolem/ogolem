/**
Copyright (c) 2014, J. M. Dieterich
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
package org.ogolem.dimerizer;

import contrib.edu.princeton.eac.Grid;
import org.ogolem.core.Geometry;
import org.ogolem.core.GeometryConfig;
import org.ogolem.core.GlobalConfig;
import org.ogolem.core.MoleculeConfig;
import org.ogolem.core.Newton;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * The configuration object with dimerizer specific stuff.
 * @author Johannes Dieterich
 * @version 2013-03-07
 */
public class DimerizerConfig {
    
    private static final Logger log = LoggerFactory.getLogger(DimerizerConfig.class);
    
    Geometry geom;
    Grid<Geometry> grid;
    Newton locopt;
    
    public DimerizerConfig(final String[] inputData, final GlobalConfig globConf) throws Exception {
        this(inputData,globConf,6);
    }
    
    public DimerizerConfig(final String[] inputData, final GlobalConfig globConf,
            final int firstXDim) throws Exception{
        
        assert(firstXDim > 0 && firstXDim < 7);
        
        if(!inputData[0].trim().equalsIgnoreCase("###OGOLEMDIMERS###")){
            throw new Exception("Configuration file must begin with ###OGOLEMDIMERS###!");
        }
        for(int i = 1; i < inputData.length; i++){
            final String line = inputData[i].trim();
            if(line.startsWith("//") || line.startsWith("#")) {
                continue;
            } else if(line.startsWith("Grid=")){
                
                final double[] starts = new double[firstXDim];
                final double[] ends = new double[firstXDim];
                final double[] space = new double[firstXDim];
                
                final String[] sa = line.substring(5).trim().split("\\/");
                for(int x = 0; x < firstXDim; x++){
                    final String s = sa[x];
                    final String[] sax = s.trim().split("\\:");
                    starts[x] = Double.parseDouble(sax[0].trim());
                    ends[x] = Double.parseDouble(sax[1].trim());
                    space[x] = Double.parseDouble(sax[2].trim());
                    if(space[x] < 0 || ends[x] < starts[x]){
                        throw new Exception("Something wrong with grid dim " + x);
                    }
                }
                
                if(firstXDim == 6){
                
                    // check the sanity of the Euler angles (yes, not all fuckups...)
                    if(starts[3] < -Math.PI){starts[3] = -Math.PI;}
                    if(ends[3] >  Math.PI){ends[3] =  Math.PI;}
        
                    if(starts[4] < -0.5*Math.PI){starts[4] = -0.5*Math.PI;}
                    if(ends[4] >  0.5*Math.PI){ends[4] =  0.5*Math.PI;}
        
                    if(starts[5] < -Math.PI){starts[5] = -Math.PI;}
                    if(ends[5] >  Math.PI){ends[5] =  Math.PI;}
                }
                
                log.info("Trying to setup grid.");
                
                this.grid = new Grid<>(starts,ends,space);
                
                log.info("Grid setup complete.");
            } else{
                throw new Exception("Unknown configuration option: " + line);
            }
        }
        
        this.locopt = globConf.getRefNewton().clone();
        
        // setup constraints (we will need them)
        final GeometryConfig gc = globConf.geoConfCopy();
        final MoleculeConfig mc1 = gc.geomMCs.get(0);
        mc1.constricted = true;
        mc1.constraints = new boolean[3][mc1.noOfAtoms];
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < mc1.noOfAtoms; j++){
                mc1.constraints[i][j] = false;
            }
        }
        
        final MoleculeConfig mc2 = gc.geomMCs.get(1);
        mc2.constricted = true;
        mc2.constraints = new boolean[3][mc2.noOfAtoms];
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < mc2.noOfAtoms; j++){
                mc2.constraints[i][j] = false;
            }
        }
        
        this.geom = new Geometry(gc);
        if(geom.getNumberOfIndieParticles() != 2){
            throw new Exception("Dimerizer wants 2 (!) molecules in geometry. Not more, not less.");
        }
        
        if(grid == null || locopt == null){
            throw new Exception("Something missing in config. Need a grid and a locopt.");
        }
    }
    
    public Grid<Geometry> getGrid(){
        return this.grid;
    }
}
