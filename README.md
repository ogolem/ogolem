# OGOLEM

The Java source of the OGOLEM framework for global optimization.

## To build OGOLEM
Prerequisites:
* Java25+ JDK
* gradle 9 or later

Invoke `gradle build` and a shadowjar will be created in `build/libs/ogolem-snapshot.jar`.

## To run OGOLEM
Prerequisites:
* Java25+ JRE
* the OGOLEM shadowjar - either from building or from downloading an artifact

Please see the manual for input options.

Currently, one MUST explicitly enable incubator modules used by OGOLEM:
```
java --add-modules=jdk.incubator.vector --enable-preview -jar ogolem-snapshot.jar ....
```
