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
package org.ogolem.spectral;

import java.io.Serializable;

/**
 * A peak.
 * @author Johannes Dieterich
 * @version 2014-11-29
 */
public class Peak implements Comparable<Peak>, Serializable, Cloneable {
    
    private static final long serialVersionUID = (long) 20131201;
    private final double location;
    private double intensity;

    public Peak(final double loc, final double inte) {
        this.location = loc;
        this.intensity = inte;
    }

    public Peak(final Peak orig) {
        this.intensity = orig.intensity;
        this.location = orig.location;
    }

    @Override
    public Peak clone() {
        return new Peak(this);
    }

    public double getLocation() {
        return location;
    }

    public double getIntensity() {
        return intensity;
    }

    public void setIntensity(final double intensity) {
        assert (intensity == intensity);
        assert (intensity >= 0.0);
        this.intensity = intensity;
    }

    @Override
    public int compareTo(final Peak p) {

        final int BEFORE = -1;
        final int EQUAL = 0;
        final int AFTER = 1;

        // optimization
        if (this == p) {
            return EQUAL;
        }

        if (p.getLocation() > this.getLocation()) {
            return BEFORE;
        } else if (p.getLocation() < this.getLocation()) {
            return AFTER;
        } else if (p.getLocation() == this.getLocation()) {
            return EQUAL;
        } else {
            throw new RuntimeException("ERROR: WTF situation in peak!");
        }
    }
}
