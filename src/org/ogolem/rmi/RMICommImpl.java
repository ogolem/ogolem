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

import java.rmi.RemoteException;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.helpers.Tuple;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.rmi.RMICodes.JOBSTATE;

/**
 * Implements our RMI communication interface.
 * @author Johannes Dieterich
 * @version 2016-01-25
 */
public class RMICommImpl<T> implements RMICommunication<T>{

    // the UID for serialization purposes, very important!
    private static final long serialVersionUID = (long) 20141101;

    private static final boolean DEBUG = false;

    private final TaskQueue<T> queue;
    private final ClientList<T> list;

    private final boolean anyProblem;
    private final int noProxies;
    private final String keySuffix;


    /**
     * Constructs the implementation.
     * @throws RemoteException
     */
    RMICommImpl(final TaskQueue<T> taskQueue, final boolean anyProblems,
            final long clientTimeoutMillis, final long clientJobTimeoutMillis,
            final int noProxies)
            throws RemoteException {
        super();

        this.queue = taskQueue;
        this.list = new ClientList<>(clientTimeoutMillis, clientJobTimeoutMillis, this);
        this.anyProblem = anyProblems;
        this.noProxies = noProxies;
        
        String tmpSuffix = System.getenv("OGO_RMIKEY");
        if(tmpSuffix == null){
            tmpSuffix = "Super secret password: Nakatomi Socrates";
        }
        this.keySuffix = tmpSuffix;
    }

    @Override
    public Tuple<String,Integer> registerWithMaster(final String key) throws RemoteException{

        if(anyProblem){
            throw new RemoteException("Configuration problem.");
        }
        
        // we use the key mainly to make sure that we do not mess up with any other program
        if(key.equalsIgnoreCase("Client speaking, I am here. " + keySuffix)){
            return new Tuple<>("Master speaking, everything fine. " + keySuffix, list.newClient());
        } else if(key.startsWith("Client speaking, I am here. ") && !key.equalsIgnoreCase("Client speaking, I am here. " + keySuffix)){
            return new Tuple<>("Master speaking, you are contacting the wrong server. Contacted server: " + keySuffix,-1);
        } else {
            return new Tuple<>("Master speaking, absolutely wrong key.",-1);
        }
    }


    /**
     * Get a compute task to work on.
     * @return the next task to work on
     * @throws RemoteException
     */
    @Override
    public Task<T> getWorkTask(final int clientID) throws RemoteException{
        Task<T> t = queue.getNextTask();
        list.idGotAJob(t, clientID);
        return t;
    }

    /**
     * Return a result from a calculation to the server.
     * @param result
     * @return Returns encoded what the overall state is.
     * @throws java.rmi.RemoteException
     */
    @Override
    public RMICodes.JOBSTATE returnResult(final Result<T> result) throws RemoteException{

        if(DEBUG) System.out.println("DEBUG: Client " + result.getClientID() + " about to return result.");

        if(list.allowedToReturn(result.getClientID())){
            if(DEBUG) System.out.println("DEBUG: Client " + result.getClientID() + " allowed to return.");
            list.idReturnedResult(result.getClientID());
            final RMICodes.JOBSTATE next = queue.submitResult(result);
            if(next == JOBSTATE.FINISH) list.removedOneClient(result.getClientID());
            return next;
        } else{
            System.err.println("WARNING: Client " + result.getClientID() + " was not allowed to return result.");
            return whatNext(result.getClientID());
        }
    }

    void returnResultForced(final Result<T> result){
        queue.submitResult(result);
    }

    @Override
    public RMICodes.JOBSTATE whatNext(final int id) throws RemoteException{
        if(queue.queueClosed()){
            list.removedOneClient(id);
            return JOBSTATE.FINISH;
        } else if(queue.queueWaiting()){
            list.clientPinging(id);
            return JOBSTATE.WAITING;
        }

        return JOBSTATE.CONTINUE;
    }

    /**
     * Returns whether everything is OK on the server side.
     * @param id the id of the client
     * @return true if everything is OK, false otherwise.
     * @throws RemoteException
     */
    @Override
    public boolean isServerAliveAndGood(final int id) throws RemoteException{
        list.clientPinging(id);
        return !anyProblem;
    }

    /**
     * Returns whether all tasks have been fulfilled.
     * @return true if the job is finished and the clients are allowed to exit
     * @throws RemoteException
     */
    @Override
    public boolean isEverythingDone() throws RemoteException{
        if(DEBUG){System.out.println("DEBUG: Everything done? " + queue.queueClosed() + " " + list.zeroClientsWaiting());}
        return (queue.queueClosed() && list.zeroClientsWaiting());
    }

    @Override
    public void stopClientList(){
        list.done();
    }
    
    @Override
    public List<T> synchronizePool(final List<T> myPool, final int whichID, final int maxStructsBack, final long lastStart) throws RemoteException {
        
        if(DEBUG) System.out.println("DEBUG: Proxy " + whichID + " about to merge pool.");

        if(list.allowedToReturn(whichID)){
            if(DEBUG) System.out.println("DEBUG: Proxy " + whichID + " allowed to merge.");
            final List<T> merged = queue.mergePools(list.howManyJobsForID(whichID), myPool, maxStructsBack, lastStart);
            list.idReturnedResult(whichID); // important that this is AFTER the call to howManyJobs()
            
            return merged;
        } else {
            System.err.println("WARNING: Proxy " + whichID + " was not allowed to merge pool.");
            return queue.mergePools(0, new ArrayList<>(0), maxStructsBack, lastStart);
        }
    }
    
    @Override
    public synchronized List<Task<T>> getInitialTasks(final int maxTasks,final int id) throws RemoteException {
        
        if(queue.queueWaiting()){
            list.clientPinging(id);
            return new ArrayList<>();
        }
        
        final List<Task<T>> tasks = queue.nextInitTasks(maxTasks, noProxies);
        list.idGotMultipleJobs(tasks.size(), new DummyTask<>(), id);
        
        return tasks;
    }
    
    @Override
    public synchronized Tuple3D<RMICodes.JOBSTATE, Long, Integer> getGlobOptChunk(final int maxTasks, final int id) throws RemoteException {
        
        if(DEBUG){System.out.print("DEBUG: State of comm impl: queue closed? " + queue.queueClosed() + " queue waiting? " + queue.queueWaiting());}
        
        if(queue.queueClosed()){
            list.removedOneClient(id);
            return new Tuple3D<>(JOBSTATE.FINISH,0l,0);
        } else if(queue.queueWaiting()){
            list.clientPinging(id);
            return new Tuple3D<>(JOBSTATE.WAITING,0l,0);
        }

        // now actually get a chunk...
        final Tuple3D<Boolean,Long,Integer> tup = queue.nextGlobOptChunk(maxTasks);
        if(DEBUG){System.out.println("DEBUG: Chunk from queue: " + tup.getObject1() + " " + tup.getObject2() + " " + tup.getObject3());}
        
        final int numTasks = tup.getObject3();
        if(numTasks > 0){
            list.idGotMultipleJobs(tup.getObject3(), new DummyTask<>(), id);
        }
        
        if(!tup.getObject1()){
            // just not done yet?
            if(tup.getObject2() == 0l){
                list.clientPinging(id);
                return new Tuple3D<>(JOBSTATE.WAITING,0l,0);
            } else {
                // done!
                list.removedOneClient(id);
                return new Tuple3D<>(JOBSTATE.FINISH,0l,0);
            }
        }
        
        return new Tuple3D<>(JOBSTATE.CONTINUE,tup.getObject2(),tup.getObject3());
    }
}