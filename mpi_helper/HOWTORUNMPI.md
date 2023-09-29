## How to run MPI jobs on with OGOLEM
### Preface
Message Passing Interface (MPI) is a de-facto standard for message passing in broad use on HPC cluster installations. Unfortunately, it comes with a multitude of limitations (vendor specific libraries, struct sizes, no easy way to recover from processes failing or adding processes at runtime, ...). We hence STRONGLY recommend using our RMI implementation which is completely platform agnostic and recovers from client failures.

If you however must use MPI, be aware that these instructions come without support and will likely require adjusting based on your MPI environment (another one of MPI's most endearing "features"). 

As of July 2021, we have been able to confirm this working with OpenMPI 4.x on macOS Big Sur, different Linux distributions, and with SLURM support enabled on multi-node Linux. We have NOT been able to get a recent oneAPI MPI backend working.

### Build instructions
You must compile the `ogompi.c` wrapper in this directory into a shared library linked against your MPI library. This is NOT a portable process thanks to MPI and must be done for every MPI runtime and version you want to use.

E.g., on macOS Big Sur with an Open MPI installed in `$HOME/software/openmpi_install`:
```
~/software/openmpi_install/bin/mpicc -dynamiclib -o libogompi.dylib ogompi.c 
```
NOTE: Above needs to be adjusted for the used environment. E.g., if on Linux, replace `-dynamiclib` with `-fPIC -shared`. Also ensure that the shared library follows correct naming convention for your OS, e.g., `libogompi.so` for FreeBSD, Linux or `libogompi.dylib` for macOS.

If one does not want to use the `mpicc` wrapper (not recommended), an equivalent on macOS Big Sur is:
```
clang -dynamiclib -o libogompi.dylib -I$HOME/software/openmpi_install/include -L$HOME/software/openmpi_install/lib -lmpi  ogompi.c
```
NOTE: Above needs to be adjusted for the used environment just as the wrapped command and additionally if, e.g., using the GNU toolchain, replace `clang` with `gcc`. E.g., for a different MPI implementation and/or install location (adjust `-I...` and `-L...`).

Place the library in a location to be provided at runtime to the `java` command. Here, we placed it in the installed Open MPI `lib` directory. If this is not possible, place in a different directory and ensure BOTH the MPI library directory as well as the directory containing the wrapper library are passed to `-Djava.library.path` below.

A JDK21 is required with support for the foreign linker preview feature.

The actual command will be highly specific to your environment (thanks to MPI) but we will highlight some common requirements using an example:

```
~/software/openmpi_install/bin/mpirun -np 4 ~/software/jdk16/jdk-16.jdk/Contents/Home/bin/java --add-modules=jdk.incubator.vector --enable-preview --enable-native-access=ALL-UNNAMED -Djava.library.path=~/software/openmpi_install/lib -jar ogolem-snapshot.jar --core ar38_sweden.ogo
```

This assumes the `ogolem-snapshot.jar` in a directory together with the input file (here: `ar38_sweden.ogo`). `--core` selects the MPI backend.

You must enable the `foreign` preview using `--enable-preview` and grant permissions to it using `--enable-native-access=ALL-UNNAMED`. Also, the vectorAPI incubator must be enabled as usual.

The location of the previously compiled `libogompi.dylib` is provided using `-Djava.library.path=~/software/openmpi_install/lib` and MUST match the `mpirun` used. Here, we use `mpirun` with 4  processes.
