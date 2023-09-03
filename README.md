# OGOLEM

The Java source of the OGOLEM framework for global optimization.

## To build OGOLEM
Prerequisites:
* Java19+ JDK
* gradle7+

Invoke `gradle build` and a shadowjar will be created in `build/libs/ogolem-snapshot.jar`.

## To run OGOLEM
Prerequisites:
* Java19+ JRE
* the OGOLEM shadowjar - either from building or from downloading a Jenkins artifact

Please see the manual for input options.

Currently, one MUST explicitly enable incubator modules used by OGOLEM:
```
java --add-modules=jdk.incubator.vector --enable-preview -jar ogolem-snapshot.jar ....
```
