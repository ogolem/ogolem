/**
Copyright (c) 2011-2014, J. M. Dieterich
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

import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.ConcurrentHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Timer;
import java.util.TimerTask;
import org.ogolem.helpers.Tuple3D;

/**
 * A list of clients, their state and the work they are supposed to do.
 * @author Johannes Dieterich
 * @version 2016-04-05
 */
final class ClientList<T>{

    private static final boolean DEBUG = false;
    
    private final Timer timer;
    private final long timeoutJob;
    private final long timeoutContact;
    private final RMICommImpl<T> rmi;
    private final Map<Integer,Tuple3D<Integer,Task<T>,Long>> list;
    private final Map<Integer,Long> recentlyUpdated;
    private final List<Integer> finishedClients;
    private int clientCounter = -1;

    ClientList(long clientTimeoutMillis, long clientJobTimeoutMillis, RMICommImpl<T> rmiImpl){
        this.timer = new Timer();
        timer.schedule(new ClientCleaner(), clientTimeoutMillis, //initial delay
        clientTimeoutMillis); //subsequent rate
        this.timeoutContact = clientTimeoutMillis;
        this.timeoutJob = clientJobTimeoutMillis;
        this.rmi = rmiImpl;
        this.list = new ConcurrentHashMap<>(32,0.75f,2);
                //Collections.synchronizedMap(new HashMap<Integer,Tupel<Task<T>,Long>>());
        this.recentlyUpdated = new ConcurrentHashMap<>(32,0.75f,2);
        this.finishedClients = Collections.synchronizedList(new ArrayList<Integer>());
    }

    synchronized int newClient(){
        // must be synchronized: uses mutable field
        clientCounter++;
        final long time = System.currentTimeMillis();
        recentlyUpdated.put(clientCounter, time);
        return clientCounter;
    }

    void idGotAJob(final Task<T> t, final int id){
        final long time = System.currentTimeMillis();
        if(t != null){
            // add to hashmap
            final Tuple3D<Integer,Task<T>,Long> tup = new Tuple3D<>(1,t,time);
            list.put(id, tup);
        }
        recentlyUpdated.put(id, time);
    }
    
    void idGotMultipleJobs(final int noTasks, final Task<T> t, final int id){
        final long time = System.currentTimeMillis();
        if(t != null){
            // add to hashmap
            final Tuple3D<Integer,Task<T>,Long> tup = new Tuple3D<>(noTasks,t,time);
            list.put(id, tup);
        }
        recentlyUpdated.put(id, time);
    }
    
    int howManyJobsForID(final int id){
        
        final Tuple3D<Integer,Task<T>,Long> tup = list.get(id);
        if(tup == null){
            System.err.println("ERROR: ID " + id + " did not get any jobs!");
            return 0;
        }
        
        return tup.getObject1();
    }

    boolean allowedToReturn(final int id){
        return list.containsKey(id);
    }

    void idReturnedResult(final int id){
        // delete from the one hashmap
        list.remove(id);
        // update in the other
        final long time = System.currentTimeMillis();
        recentlyUpdated.put(id, time);
    }

    void removedOneClient(final int id){
        recentlyUpdated.remove(id);
        list.remove(id);
        finishedClients.add(id);
    }

    void clientPinging(final int id){
        if(finishedClients.contains(id)){return;}
        final long time = System.currentTimeMillis();
        recentlyUpdated.put(id, time);
    }

    boolean zeroClientsWaiting(){
        if(DEBUG){
            String s = "DEBUG: Current recentlyUpdated list contains: ";
            final Set<Entry<Integer,Long>> entries = recentlyUpdated.entrySet();
            for(final Entry<Integer,Long> entry : entries){
                s += " // id " + entry.getKey() + " time " + entry.getValue();
            }
            System.out.println(s);
        }
        return recentlyUpdated.isEmpty();
    }

    void done(){
        timer.cancel();
    }

    private class ClientCleaner extends TimerTask{

        @Override
        public void run(){

            final long currTime = System.currentTimeMillis();

            for(final Map.Entry<Integer,Tuple3D<Integer,Task<T>,Long>> entry : list.entrySet()){
                final Tuple3D<Integer,Task<T>,Long> t = entry.getValue();
                final int id = entry.getKey();
                final long startTime = t.getObject3();
                if(currTime - startTime >= timeoutJob){
                    // too old: clean up
                    list.remove(id);
                    System.err.println("INFO: Removed id " + id + " with timestamp " + startTime + " at " + currTime + ".");
                    // submit dummy
                    final int noResults = t.getObject1();
                    final Result<T> dummy = t.getObject2().getDummyAnswer(id);
                    try{
                        for(int assoc = 0; assoc < noResults; assoc++){
                            rmi.returnResultForced(dummy);
                        }
                    } catch(Exception e){
                        System.err.println("ERROR: Couldn't return dummy result to master. This is a serious problem! " + e.toString());
                        e.printStackTrace(System.err);
                    }
                }
            }

            for(final Map.Entry<Integer,Long> t : recentlyUpdated.entrySet()){
                final int id = t.getKey();
                final long updateTime = t.getValue();
                if(currTime - updateTime >= timeoutContact){
                    // too old: clean up
                    recentlyUpdated.remove(id);
                    System.err.println("INFO: Removed id " + id + " from the list of frequently updated clients.");
                }
            }
        }
    }
}
