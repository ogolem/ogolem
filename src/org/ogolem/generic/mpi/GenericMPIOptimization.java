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
package org.ogolem.generic.mpi;

import org.ogolem.generic.Optimizable;
import mpi.MPI;
import mpi.Status;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.rmi.Job;
import org.ogolem.rmi.Result;
import org.ogolem.rmi.Task;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A generic class for MPI based parallelization of our optimization problems.
 * Please note: you must call MPI_INIT from OUTSIDE this and decide already,
 * if the process is master or slave!
 * @author Johannes Dieterich
 * @version 2014-04-02
 */
public class GenericMPIOptimization<E, T extends Optimizable<E>> {
    
    private static final Logger log = LoggerFactory.getLogger(GenericMPIOptimization.class);
    private static final int NEXTTASK = 0;
    private static final int EXITDONE = 1;
    private static final int WAITFOR = 2;
    private static final int KICKOFF = 42;
    private static final long TIMEOUT = 5000; // 5 seconds
    
    private GenericMPIOptimization(){}
    
    @SuppressWarnings("unchecked")
    public static <T> void runAsMaster(final Job<T> job) throws Exception {
        
        log.debug("Entering generic MPI globopt as master. BRACE YOURSELF BIG TIMES!");
        
        final int noProcs = MPI.COMM_WORLD.Size();
        if(noProcs <= 1){
            throw new RuntimeException("Trying to actually work all by myself is not in my nature. Bye.");
        }
        
        /*
         * initial broadcast
         */
        final char[] initMessage = "Hello from MPI master, all is well!".toCharArray();
        MPI.COMM_WORLD.Bcast(initMessage, 0, 35, MPI.CHAR, KICKOFF);
        
        /*
         * fill up all slaves once
         */
        int taskCounter = 1;
        for(int proc = 1; proc < noProcs; proc++){
            final Task<T> task = job.nextTask();
            final String outFile = "task" + taskCounter + ".dat";
            OutputPrimitives.writeObjToBinFile(outFile, task);
            final char[] message = new char[50];
            final char[] outFM = outFile.toCharArray();
            System.arraycopy(outFM, 0, message, 0, outFM.length);
            MPI.COMM_WORLD.Send(message, 0, 50, MPI.CHAR, proc, 0);
            taskCounter++;
        }
        
        while(!job.jobFinished()){
            
            // receive something
            final char[] answer = new char[50];
            final Status rcvStat = MPI.COMM_WORLD.Recv(answer, 0, 50, MPI.CHAR, MPI.ANY_SOURCE, MPI.ANY_TAG);
            final String resPath = String.copyValueOf(answer).trim();
            
            final Result<T> result = (Result<T>) InputPrimitives.readBinInput(resPath);
            
            // delete the result file
            ManipulationPrimitives.remove(resPath);
            
            // submit it
            job.submitResult(result);
            
            // and give something out
            final Task<T> task = job.nextTask();
            if(task == null){
                final char[] message = new char[50];
                MPI.COMM_WORLD.Send(message, 0, 50, MPI.CHAR, rcvStat.source, WAITFOR);
            } else {
                final String taskPath = "task" + taskCounter + ".dat";
                OutputPrimitives.writeObjToBinFile(taskPath, task);
                final char[] message = new char[50];
                final char[] outFM = taskPath.toCharArray();
                System.arraycopy(outFM, 0, message, 0, outFM.length);
            
                MPI.COMM_WORLD.Send(message, 0, 50, MPI.CHAR, rcvStat.source, NEXTTASK);
                taskCounter++;
            }
        }
        
        
        /*
         * tell everybody to go kill themselves
         */
        final char[] finalMessage = new char[50];
        for(int proc = 1; proc < noProcs; proc++){
            MPI.COMM_WORLD.Send(finalMessage, 0, 50, MPI.CHAR, proc, EXITDONE);
        }
    }
    
    @SuppressWarnings("unchecked")
    public static <X,Y extends Optimizable<X>> void runAsSlave(final long timeout) throws Exception {
        
        log.debug("Entering generic MPI globopt as slave. BRACE YOURSELF BIG TIMES!");
        
        final int myRank = MPI.COMM_WORLD.Rank();
        
        // try to receive the initial broadcast
        final char[] initMessage = new char[35];
        MPI.COMM_WORLD.Bcast(initMessage, 0, 35, MPI.CHAR, KICKOFF);
        
        if(!initMessage.toString().trim().equalsIgnoreCase("Hello from MPI master, all is well!")){
            throw new RuntimeException("Initial message from master was not what " + myRank + "expected. Exiting. Received message: "
                + initMessage);
        }
                
        // do stuff as long as time permits
        final long startTime = System.currentTimeMillis();
        int taskCounter = 0;
        BigLoop: while((System.currentTimeMillis()-startTime) < timeout ){
            
            
            // always the same idea: ask the master for a new task
            final char[] message = new char[50];
            final Status status = MPI.COMM_WORLD.Recv(message, 0, 50, MPI.CHAR, 0, MPI.ANY_TAG);
            
            final int messTag = status.tag;
            if (messTag == EXITDONE || (messTag != WAITFOR && messTag != NEXTTASK)) {
                log.debug("Master tells me to quit: doing so! Tag was " + messTag);
                break;
            } else if(messTag == WAITFOR) {
                Thread.sleep(TIMEOUT);
                continue BigLoop;
            }
            
            // there is a task: execute it
            log.debug("There is a new task!");
            taskCounter++;

            final Task<Y> task = (Task<Y>) InputPrimitives.readBinInput(message.toString().trim());
            final Result<Y> result = task.executeTask(myRank);
            
            // delete task
            ManipulationPrimitives.remove(message.toString().trim());

            final String outputFile = "result" + taskCounter + ".bin";
            OutputPrimitives.writeObjToBinFile(outputFile, result);

            // report back
            final char[] filePath = outputFile.toCharArray();
            final char[] answer = new char[50];
            assert (filePath.length <= answer.length);
            System.arraycopy(filePath, 0, answer, 0, filePath.length);
            MPI.COMM_WORLD.Send(answer, 0, 77, MPI.CHAR, 0, myRank);
        }
    }
}
