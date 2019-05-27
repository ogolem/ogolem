/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015-2016, J. M. Dieterich and B. Hartke
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

import java.rmi.Remote;
import java.rmi.RemoteException;
import java.util.List;
import org.ogolem.helpers.Tuple;
import org.ogolem.helpers.Tuple3D;

/**
 * Defines all client-server communication using RMI.
 * @author Johannes Dieterich
 * @version 2016-01-25
 * @param <T> typically an optimizable
 */
public interface RMICommunication<T> extends Remote {
    
    /**
     * Register a client with the server.
     * @param sKey
     * @return A string and an ID for the client.
     * @throws RemoteException
     */
    public Tuple<String,Integer> registerWithMaster(String sKey) throws RemoteException;

    /**
     * Get a task. If it doesn't contain anything,
     * we are done and that's it.
     * @param clientID the clients id
     * @return A task and the information needed to carry it out.
     * @throws RemoteException
     */
    public Task<T> getWorkTask(int clientID) throws RemoteException;

    /**
     * Hands over a result to the server and gets in turn the current state
     * encoded as a short.
     * @param res a result of a step
     * @return the state
     * @throws RemoteException
     */
    public RMICodes.JOBSTATE returnResult(Result<T> res) throws RemoteException;

    /**
     * Get the current state encoded as a short.
     * @param id the id of the client
     * @return next task to do
     * @throws RemoteException
     */
    public RMICodes.JOBSTATE whatNext(int id) throws RemoteException;

    /**
     * Returns whether everything is OK on the server side.
     * @param id the id of the client
     * @return true if everything is OK, false otherwise.
     * @throws RemoteException
     */
    public boolean isServerAliveAndGood(int id) throws RemoteException;

    /**
     * Returns whether all tasks are complete. This is kind of redundant with
     * the above return codes but provides for an explicit answer.
     * @return true if we are done, false otherwise.
     * @throws RemoteException
     */
    public boolean isEverythingDone() throws RemoteException;
    
    /**
     * Synchronize a proxy pool with the main pool.
     * @param myPool the local list of optimizables
     * @param id the id of the client
     * @param maxStructsBack the maximal number of structures this wants back
     * @param lastStart the start ID of the last chunk before merging
     * @return the main pool after synchronizing this one and the short id what next.
     * @throws RemoteException 
     */
    public List<T> synchronizePool(final List<T> myPool, final int id, final  int maxStructsBack, final long lastStart) throws RemoteException;
    
    public List<Task<T>> getInitialTasks(final int maxTasks, final int id) throws RemoteException;
    
    public Tuple3D<RMICodes.JOBSTATE, Long, Integer> getGlobOptChunk(final int maxTasks, final int id) throws RemoteException;
    
    public void stopClientList()throws RemoteException;
}
