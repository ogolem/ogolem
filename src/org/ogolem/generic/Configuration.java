/**
Copyright (c) 2013-2014, J. M. Dieterich
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
package org.ogolem.generic;

import java.io.Serializable;
import java.util.List;
import org.ogolem.generic.generichistory.GenericHistoryConfig;
import org.ogolem.generic.genericpool.GenericPoolConfig;
import org.ogolem.generic.genericpool.NicheComputer;
import org.ogolem.generic.threading.TaskFactory;

/**
 * Generic interface to configuration.
 * @author Johannes Dieterich
 * @version 2015-11-16
 */
public interface Configuration<E,T extends Optimizable<E>> extends Serializable {
    
    public List<String> getFormattedConfig();
    
    public String seedFolder();
    
    public String intermediatePoolFile();
    
    public GenericFitnessFunction<E,T> getFitnessFunction() throws Exception;
    
    public GenericHistoryConfig getGenericHistoryConfig() throws Exception;
    
    public GenericPoolConfig<E,T> getGenericPoolConfig() throws Exception;
    
    public int getNumberOfGlobalSteps() throws Exception;
    
    public GenericInitializer<E,T> getInitializer() throws Exception;
    
    public <V extends GenericInitializer<E,T>> TaskFactory<E,T,V> getInitFactory() throws Exception;
    
    public GenericGlobalOptimization<E,T> getGlobalOptimization() throws Exception;
    
    public <V extends GenericGlobalOptimization<E,T>> TaskFactory<E,T,V> getGlobalFactory() throws Exception;
    
    public T getExample() throws Exception;
    
    public IndividualReader<T> getReader() throws Exception;
    
    public IndividualWriter<T> getWriter(final String folder) throws Exception;
    
    public boolean wantsDetailedStats() throws Exception;
    
    /**
     * Get, if used, the niche computer to be used.
     * @return a valid niche computer if niching should be done, null otherwise
     * @throws Exception if some internal error occurs
     */
    public NicheComputer<E,T> getNicheComputer() throws Exception;
}
