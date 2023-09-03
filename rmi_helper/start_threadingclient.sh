#!/usr/bin/env bash
set -e
set -u

server_ip=$1

java -Djava.rmi.server.codebase=file:ogolem-snapshot.jar \
-Djava.security.policy=client.policy \
--add-modules=jdk.incubator.vector --enable-preview \
-jar ogolem-snapshot.jar \
--rmithreader \
--jobtype geom \
--inputfile ar100.ogo \
--server $server_ip \
--serverregistryport 1099 \
--threads 4 \
--sleeptime 10 \
--maxchunk 2000 \
--mergeinds -1 \
--commname RMICommunication
