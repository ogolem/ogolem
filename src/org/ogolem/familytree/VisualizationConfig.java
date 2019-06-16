/**
Copyright (c) 2015, J. M. Dieterich and B. Hartke
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

import java.io.Serializable;

/**
 * The configuration. Mainly things such as colors etc.
 * @author Johannes Dieterich
 * @version 2015-07-11
 */
class VisualizationConfig implements Serializable {
    
    private static final long serialVersionUID = (long) 20150710;
    
    float[] vertexColor = new float[]{1.0f,1.0f,1.0f};
    float[] edgeColorAsMum = new float[]{1.0f,0.0f,0.0f};
    float[] edgeColorAsDad = new float[]{0.0f,0.0f,1.0f};
    float[] vertexColorFinalPool = new float[]{0.75f,0.0f,0.0f};
    float transparency = 1.0f;

    int penWidth = 3;

    String backGround = "black";
    String circleRing = "white";
    
    boolean verboseLabel = false;
    
    boolean vertexSizeChildrenBased = false;
    
    boolean removeAllNonChildren = false;
    
    boolean addAllFinalPoolIndividuals = false; // only if above true makes sense!

    String poolFile = "pool.bin";
    
    String vertexShape = "circle";
    String vertexStyle = "filled";
    String edgeStyle = "filled";
    
    String dotOut = "familytree.gv";
    
    /**
     * Construct a default configuration
     */
    VisualizationConfig(){};
    
    /**
     * Construct a custom configuration from a file.
     * @param config an array of strings containing the configuration file, one entry per line
     * @throws Exception if some parsing error occurs
     */
    VisualizationConfig(final String[] inp) throws Exception{
        
        for (final String inp1 : inp) {
            final String s = inp1.trim();
            if(s.startsWith("#") || s.startsWith("//")){
                continue;
            } else if(s.startsWith("VertexShape=")){
                vertexShape = s.substring(12).trim();
            } else if(s.startsWith("VertexStyle=")){
                vertexStyle = s.substring(12).trim();
            } else if(s.startsWith("EdgeStyle=")){
                edgeStyle = s.substring(10).trim();
            } else if(s.startsWith("Transparency=")){
                final String s2 = s.substring(13).trim();
                transparency = Float.parseFloat(s2);
            } else if(s.startsWith("PenWidth=")){
                final String s2 = s.substring(9).trim();
                penWidth = Integer.parseInt(s2);
            } else if(s.startsWith("Background=")){
                final String s2 = s.substring(11).trim();
                backGround = s2;
            } else if(s.startsWith("RingColor=")){
                final String s2 = s.substring(10).trim();
                circleRing = s2;
            } else if(s.startsWith("EdgeColorMaternal=")){
                final String[] sa = s.substring(18).trim().split("\\;");
                edgeColorAsMum[0] = Float.parseFloat(sa[0].trim());
                edgeColorAsMum[1] = Float.parseFloat(sa[1].trim());
                edgeColorAsMum[2] = Float.parseFloat(sa[2].trim());
            } else if(s.startsWith("EdgeColorPaternal=")){
                final String[] sa = s.substring(18).trim().split("\\;");
                edgeColorAsDad[0] = Float.parseFloat(sa[0].trim());
                edgeColorAsDad[1] = Float.parseFloat(sa[1].trim());
                edgeColorAsDad[2] = Float.parseFloat(sa[2].trim());
            } else if(s.startsWith("VertexColorFinalPool=")){
                final String[] sa = s.substring(21).trim().split("\\;");
                vertexColorFinalPool[0] = Float.parseFloat(sa[0].trim());
                vertexColorFinalPool[1] = Float.parseFloat(sa[1].trim());
                vertexColorFinalPool[2] = Float.parseFloat(sa[2].trim());
            } else if(s.startsWith("OutputFile=")){
                dotOut = s.substring(11).trim();
            } else if(s.startsWith("PoolFile=")){
                poolFile = s.substring(9).trim();
            } else if(s.equalsIgnoreCase("VertexSizeChildrenBased")){
                vertexSizeChildrenBased = true;
            } else if(s.equalsIgnoreCase("VerboseLabel")){
                verboseLabel = true;
            } else if(s.startsWith("OutputMode=")){
                final String mode = s.substring(11).trim();
                switch(mode){
                    case "onlywithchildren":
                        removeAllNonChildren = true;
                        addAllFinalPoolIndividuals = false;
                        break;
                    case "onlywithchildrenandfinalpool":
                        removeAllNonChildren = true;
                        addAllFinalPoolIndividuals = true;
                        break;
                    case "onlyaccepted":
                        removeAllNonChildren = false;
                        addAllFinalPoolIndividuals = false;
                        break;
                    case "onlyacceptedandfinalpool":
                        removeAllNonChildren = false;
                        addAllFinalPoolIndividuals = true;
                        break;
                    default:
                        throw new Exception("Illegal output mode " + mode);
                }
            } else {
                throw new Exception("Unknown input option " + s);
            }
        }
    }
}
