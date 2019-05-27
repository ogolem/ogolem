#!/bin/sh

set -u
set -e
#set -x

tmprel=tmpRelease

mkdir $tmprel
cd $tmprel

jar -xf ../dist/ogolem.jar
for i in ../dist/lib/*
do
   jar -xf $i
done
rm -r META-INF
if [ -e AndroidManifest.xml ]; then
	rm AndroidManifest.xml
fi
if [ -e library.properties ]; then 
	rm library.properties
fi

# remove numal for licensing reasons
rm -rf numal

$JAVA_HOME/bin/jar -cfm ../ogolem.jar ../MANIFEST_REL.MF *

cd ..
rm -rf $tmprel

budate=`date '+%Y%m%d'`
$JAVA_HOME/bin/java -jar jarjar.jar process jarjar.rules ogolem.jar ogolem-snapshot.jar
rm ogolem.jar
