/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2011, J. M. Dieterich
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

import java.util.Random;
import java.util.Timer;
import java.util.TimerTask;

/**
 * Checks whether the server is still alive.
 * @author Johannes Dieterich
 * @version 2015-08-04
 */
final class AliveCheck{

    private final Timer timer;
    private final RMICommunication<?> comm;
    private final int myID;
    
    private boolean markedDone = false;

    /**
     * Starts the alive check as a timer.
     * @param communication The remote object.
     * @param sleepTime The sleep time (randomized for the initial delay).
     * @param id The ID of this client.
     */
    AliveCheck(final RMICommunication<?> communication, long sleepTime, int id){
        this.comm = communication;
        this.timer = new Timer();
        this.myID = id;
        final Random r = new Random();
        final long initialDelay = (long) r.nextDouble() * sleepTime;
        timer.schedule(new AliveChecker(), initialDelay, //initial delay
        sleepTime); //subsequent rate
    }

    /**
     * Stops the alive check.
     */
    public void done(){
        markedDone = true;
        timer.cancel();
    }
    
    public boolean wasMarkedDone(){
        return markedDone;
    }

    private class AliveChecker extends TimerTask{

        @Override
        public void run(){
            try {
                final boolean bAliveAndGood = comm.isServerAliveAndGood(myID);
                if(!bAliveAndGood){
                    System.err.println("ERROR: Something is wrong with the server. Aborting client. My ID is " + myID);
                    System.exit(999);
                }
            } catch (Exception e) {
                System.err.println("RMI exception: " + e.getMessage() + " in AliveCheck. Aborting client. My ID is " + myID);
                e.printStackTrace(System.err);
                System.exit(245);
            }
        }
    }
}
