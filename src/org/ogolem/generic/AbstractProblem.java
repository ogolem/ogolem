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
package org.ogolem.generic;

import org.ogolem.core.FixedValues;

/**
 * An abstract simple implementation of optimizable.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public abstract class AbstractProblem<E> implements Optimizable<E> {
    
    private static final long serialVersionUID = (long) 20200429;
    protected long myID = -1;
    protected long fatherID = -1;
    protected long motherID = -1;
    protected double fitness = FixedValues.NONCONVERGEDENERGY;
    
    public AbstractProblem(){
    }
    
    public AbstractProblem(final AbstractProblem<E> orig){
        this.fatherID = orig.fatherID;
        this.motherID = orig.motherID;
        this.myID = orig.myID;
        this.fitness = orig.fitness;
    }
    
    @Override
    public abstract AbstractProblem<E> copy();

    @Override
    public long getID(){
        return this.myID;
    }

    @Override
    public void setID(final long id){
        this.myID = id;
    }

    @Override
    public long getFatherID(){
        return this.fatherID;
    }

    @Override
    public void setFatherID(final long fatherID){
        this.fatherID = fatherID;
    }

    @Override
    public long getMotherID(){
        return this.motherID;
    }

    @Override
    public void setMotherID(final long motherID){
        this.motherID = motherID;
    }

    @Override
    public double getFitness(){
        return this.fitness;
    }

    @Override
    public void setFitness(final double fitness){
        this.fitness = fitness;
    }
}
