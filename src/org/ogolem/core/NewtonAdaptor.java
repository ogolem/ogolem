/**
Copyright (c) 2014-2015, J. M. Dieterich
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

import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.GenericFitnessBackend;
import org.ogolem.generic.GenericLocOpt;
import org.ogolem.generic.stats.GenericDetailStatistics;

/**
 * An adaptor to fit the old Newton interface into the new generic world.
 * @author Johannes Dieterich
 * @version 2015-01-04
 */
public class NewtonAdaptor implements GenericLocOpt<Molecule,Geometry> {
    
    private static final long serialVersionUID = (long) 20140328;
    private final Newton newton;
    private int count = 0;

    NewtonAdaptor(final Newton newton){
        this.newton = newton;
    }
    
    NewtonAdaptor(final NewtonAdaptor orig){
        this.newton = orig.newton.clone();
    }
    
    @Override
    public NewtonAdaptor clone() {
        return new NewtonAdaptor(this);
    }

    @Override
    public String getMyID() {
        return "NEWTON ADAPTOR:\n" + newton.myIDandMethod();
    }

    @Override
    public Geometry fitness(final Geometry individual, final boolean forceOneEval) {
        
        if(forceOneEval){
            
            GenericDetailStatistics.incrementFitnessEvals();
            
            final CartesianFullBackend back = newton.getBackend();
            if(back == null){
                // no support on the backend side of things :-(
                throw new RuntimeException("Backend does not support single energy evaluation. Sorry.");
            }
            
            final long id = individual.getID();
            final CartesianCoordinates cartes = individual.getCartesiansWithEnvironment();
            final double[] eparts = new double[cartes.getNoOfAtoms()];
            
            count++;
            
            final double e = back.energyCalculation(id, count, cartes.getAll1DCartes(), cartes.getAllAtomTypes(),
                    cartes.getAllAtomNumbers(), cartes.getAllAtomsPerMol(), eparts, cartes.getNoOfAtoms(),
                    cartes.getAllCharges(), cartes.getAllSpins(), individual.getBondInfo());
            individual.setFitness(e);
            
            return individual;
        }
        
        GenericDetailStatistics.incrementLocalOpts();
        
        return newton.localOptimization(individual);
    }
    
    @Override
    public GenericBackend<Molecule,Geometry> getBackend(){
        return new FullyCartesianCoordinates(newton.getBackend());
    }
    
    @Override
    public GenericFitnessBackend<Molecule,Geometry> getFitnessBackend(){
        return new FullyCartesianCoordinates(newton.getBackend());
    }
    
    public Newton getNewton(){
        return newton;
    }
}
