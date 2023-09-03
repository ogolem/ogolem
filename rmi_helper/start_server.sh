#!/usr/bin/env bash
set -e
set -u

ip_addr=$1

java -Djava.security.main -Djava.rmi.server.codebase=file:../build/libs/ogolem-snapshot.jar \
-Djava.security.policy=polfile.txt -Djava.rmi.server.hostname=$ip_addr \
--add-modules=jdk.incubator.vector --enable-preview \
-jar ../build/libs/ogolem-snapshot.jar \
--rmiserver \
--jobtype geom \
--inputfile ar100.ogo \
--timeout 120 \
--timeoutjob 120 \
--myserverport 2500 \
--myregistryport 1099 \
--noproxies 1 \
--commname RMICommunication
