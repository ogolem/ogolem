/**
Copyright (c) 2011-2014, J. M. Dieterich
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

      This product includes software developed at the Universities of
      Kiel, Goettingen (Germany) and Princeton University (USA) by its 
      contributors: J. M. Dieterich and B. Hartke.

    * Neither the name of the University of Kiel, the University of Goettingen,
      Princeton University nor the names of its contributors may be used to
      endorse or promote products derived from this software without specific
      prior written permission.

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
package contrib.bobyqa;

/**
 * An abstract implementation of the BOBYQA method interface. Handles the normalization
 * business.
 * @author Johannes Dieterich
 * @version 2014-03-26
 */
public abstract class AbstractBOBYQAMethod implements BOBYQAMethod{
    
    
    @Override
    public abstract double computeObjectiveValue(double[] point);
    
    @Override
    public abstract boolean doesNormalize();
    
    @Override
    public void normalizeFromBounds(final double[][] bounds, final double[] point, final double[] normalized){
        
        assert(bounds != null);
        assert(point != null);
        assert(normalized != null);
        assert(bounds.length == 2);
        assert(bounds[0].length == bounds[1].length);
        assert(bounds[0].length == point.length);
        assert(bounds[0].length == normalized.length);
        for(int i = 0; i < point.length; i++){
            final double diff = bounds[1][i] - bounds[0][i];
            normalized[i] = (point[i]-bounds[0][i])/diff;
            assert(normalized[i] <= 1.0 && normalized[i] >= 0.0) : "normalized " + i + "\t" + normalized[i] + "\t"
                    + bounds[0][i] + "\t" + bounds[1][i];
        }
    }
    
    @Override
    public void denormalizeToBounds(final double[][] bounds, final double[] normalized, final double[] point){
        
        assert(bounds != null);
        assert(point != null);
        assert(normalized != null);
        assert(bounds.length == 2);
        assert(bounds[0].length == bounds[1].length);
        assert(bounds[0].length == point.length);
        assert(bounds[0].length == normalized.length);
        for(int i = 0; i < point.length; i++){
            final double diff = bounds[1][i] - bounds[0][i];
            point[i] = normalized[i]*diff+bounds[0][i];
            assert(point[i] <= bounds[1][i] && point[i] >= bounds[0][i]);
        }
    }
}
