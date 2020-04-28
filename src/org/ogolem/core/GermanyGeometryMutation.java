/**
Copyright (c) 2014, J. M. Dieterich
              2020, J. M. Dieterich and B. Hartke
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

import org.ogolem.generic.GenericMutation;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;

/**
 * A geometry specific mutation genotype style.
 * @author Johannes Dieterich
 * @version 2020-04-25
 */
public class GermanyGeometryMutation implements GenericMutation<Molecule,Geometry> {
    
    public static final double DEFAULTLOWCOM = 0.5;
    public static final double DEFAULTHIGHCOM = 1.5;
    public static final double DEFAULTWIDTHCOM = 0.1;
    public static final double DEFAULTLOWEULER = 0.0;
    public static final double DEFAULTHIGHEULER = 1.0;
    public static final double DEFAULTWIDTHEULER = 0.1;
    
    private static final long serialVersionUID = (long) 20200425;
    private final Lottery random = Lottery.getInstance();
    private final int whichMolMut;
    private final double lowCOM;
    private final double highCOM;
    private final double widthCOM;
    private final double lowEuler;
    private final double highEuler;
    private final double widthEuler;
    
    GermanyGeometryMutation(final int whichMolMut, final double lowCOM, final double highCOM,
            final double widthCOM, final double lowEuler, final double highEuler,
            final double widthEuler){
        assert(highCOM > lowCOM);
        assert(highEuler > lowEuler);
        assert(widthEuler > 0);
        assert(widthCOM > 0);
        assert(whichMolMut == 0 || whichMolMut == 1 || whichMolMut == 2);
        this.highCOM = highCOM;
        this.highEuler = highEuler;
        this.lowCOM = lowCOM;
        this.lowEuler = lowEuler;
        this.widthCOM = widthCOM;
        this.widthEuler = widthEuler;
        this.whichMolMut = whichMolMut;
    }
    
    GermanyGeometryMutation(final GermanyGeometryMutation orig){
        this.highCOM = orig.highCOM;
        this.highEuler = orig.highEuler;
        this.lowCOM = orig.lowCOM;
        this.lowEuler = orig.lowEuler;
        this.widthCOM = orig.widthCOM;
        this.widthEuler = orig.widthEuler;
        this.whichMolMut = orig.whichMolMut;
    }

    @Override
    public GermanyGeometryMutation clone() {
        return new GermanyGeometryMutation(this);
    }

    @Override
    public String getMyID() {
        return "GERMANY\ngeometry genotype mutation\n\t COM Gaussian: "
                + lowCOM + " " + highCOM + " " + widthCOM + "\n\t Euler Gaussian: "
                + lowEuler + " " + highEuler + " " + widthEuler
                + "\n\t Molecule mutation style: " + whichMolMut;
    }

    @Override
    public Geometry mutate(final Geometry orig) {
        
        final Geometry mutated = new Geometry(orig);
        
        final int whichMol = random.nextInt(mutated.getNumberOfIndieParticles());
        final int atsPerMol = mutated.getMoleculeAtPosition(whichMol).getNumberOfAtoms();

        /*
         * which coordinate should be mutated?
         * 0: x
         * 1: y
         * 2: z
         * 3: Euler1 (only for anything larger than a single atom)
         * 4: Euler2
         * 5: Euler3
         */
        final int whichCoord = (atsPerMol == 1) ? random.nextInt(3) : random.nextInt(6);

        final double mutFactorCOM = RandomUtils.gaussDouble(lowCOM, highCOM, widthCOM);
        final double mutFactorEuler = RandomUtils.gaussDouble(lowEuler, highEuler, widthEuler);

        // the actual mutation of a single coordinate
        mutated.mutateACoord(whichMol, whichCoord, mutFactorCOM, mutFactorEuler, whichMolMut);

        // environment mutation
        if (mutated.containsEnvironment()) {
            mutated.mutateEnvironment();
        }

        return mutated;
    }    
}
