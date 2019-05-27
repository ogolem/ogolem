/**
Copyright (c) 2010-2012, J. M. Dieterich
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
package org.ogolem.rmi;

import java.util.List;
import org.ogolem.switches.ColorPalette;
import org.ogolem.switches.Switch;
import org.ogolem.switches.SwitchesConfig;
import org.ogolem.switches.SwitchesGlobOpt;
import org.ogolem.switches.Taboos;

/**
 * A switch global optimization task.
 * @author Johannes Dieterich
 * @version 2012-06-24
 */
final class SwitchGlobOptTask implements Task<Switch>{

    private static final long serialVersionUID = (long) 20101114;

    private final SwitchesConfig config;
    private final ColorPalette palette;
    private final Taboos taboos;
    private final List<Switch> parents;
    private final int id;
    private final long[] family;

    SwitchGlobOptTask(final int futureID, final List<Switch> parentSwitches,
            final SwitchesConfig configuration, final ColorPalette cp,
            final Taboos tab){
        this.config = configuration;
        this.parents = parentSwitches;
        this.id = futureID;
        this.palette = cp;
        this.taboos = tab;
        this.family = new long[]{parents.get(0).getID(),parents.get(1).getID(),futureID};
    }

    @Override
    public Result<Switch> executeTask(int onClient) {

        //XXX we might need to take measures so that the GC is not killing palette and taboos, remains to be seen

        final SwitchesGlobOpt globopt = new SwitchesGlobOpt(config);
        final Switch child = globopt.doTheGlobOpt(id, parents.get(0), parents.get(1));

        if(child == null){
            // can be a case due to the taboos
            return new Result<>(child,false, onClient, family);
        } else{
            return new Result<>(child,true, onClient, family);
        }
    }

    @Override
    public Result<Switch> getDummyAnswer(int onClient) {
        return new Result<>(parents.get(0),false, onClient, family);
    }
}
