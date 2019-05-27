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
package org.ogolem.rmi;

import org.ogolem.core.FixedValues;
import org.ogolem.generic.GenericFitnessFunction;
import org.ogolem.generic.IndividualReader;
import org.ogolem.generic.Optimizable;

/**
 * A generic seeding task.
 * @author Johannes Dieterich
 * @version 2014-04-02
 */
public class GenericSeedTask <E, T extends Optimizable<E>> implements Task<T>{
    
    private static final long serialVersionUID = (long) 20140402;
    private final IndividualReader<T> reader;
    private final long newID;
    private final T start;
    private final String file;
    private final GenericFitnessFunction<E,T> fitness;

    public GenericSeedTask(final IndividualReader<T> reader, final long newID,
            final T start, final GenericFitnessFunction<E,T> fitness, final String file){
        this.fitness = fitness;
        this.reader = reader;
        this.start = start;
        this.newID = newID;
        this.file = file;
    }
    
    @Override
    public Result<T> executeTask(final int onClient) {
        
        try{
            reader.populateIndividualFromFile(start, file);
        } catch(Exception e){
            System.err.println("Failure to initialize individual from file " + file);
            e.printStackTrace(System.err);
        }
        
        start.setID(newID);
        final T opt = fitness.fitness(start, false);
        opt.setID(newID);
        
        return new Result<>(opt,true,onClient);
    }

    @SuppressWarnings("unchecked")
    @Override
    public Result<T> getDummyAnswer(final int onClient) {
        
        final T clone = (T) start.clone();
        clone.setFitness(FixedValues.NONCONVERGEDENERGY);
        clone.setID(newID);
        
        return new Result<>(clone,false,onClient);
    }    
}
