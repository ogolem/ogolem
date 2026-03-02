# Global optimization of clusters

For the global optimization of cluster structures on a shared memory basis, invoke

    java -jar ogolem.jar --shmem input.ogo 2

where `input.ogo` is your input file, needing to have `.ogo` as a suffix, followed by the number of threads you want to use, 2 in this case.

If you want to make use of the MPI frontend, please see notes in `mpi_helper/HOWTORUNMPI.md` on how to link <span class="smallcaps">ogolem</span> to a MPI implementation. In general, using the shared memory parallelization will be easier and is preferred for multi-CPU, single-node jobs.

See Chapter. <a href="#rmi" data-reference-type="ref" data-reference="rmi">8</a> for a JVM-based substitute for MPI with some more advanced distribution of computing tasks and stability. However, this feature requires some setup! Please contact the authors for help with this.

## General input structure

Every input file needs to start with the line

    ###OGOLEM###

for <span class="smallcaps">ogolem</span> to recognize that this is supposed to be a valid input file. Now you can add some input keywords and the building block input. This includes the actual cluster building blocks (atoms or molecules) as well as (if applicable) spins and charges.

In case of explicit degree of freedom searches or constraints, more input blocks are necessary.

## Building block input

In general, any combination of any building blocks is possible in principle. The only thing you should have for the easiest case is the number of building blocks and xyz-files of all your units.

The input starts with an opening tag followed by the number of building blocks in the next line, 5 in this example

    <GEOMETRY>
    NumberOfParticles=5

The building blocks (named molecules) are now put there in a similar form after each other, for example

    <MOLECULE>
    MoleculePath=water.xyz
    </MOLECULE>
    <MOLECULE>
    MoleculePath=ethanol.xyz
    </MOLECULE>
    <MOLECULE>
    MoleculePath=citronacid.xyz
    </MOLECULE>
    <MOLECULE>
    MoleculePath=quinine.xyz
    </MOLECULE>
    <MOLECULE>
    MoleculePath=water.xyz
    </MOLECULE>

The number of molecules and the `NumberOfParticles` needs to agree, otherwise an exception will be raised.

If the molecule consists of a single atom, another input option might be easier

    <MOLECULE>
    Na;0.0;0.0;0.0
    </MOLECULE>

where the 0.0s are internal xyz values.

In case of the molecule having internal degrees of freedom that should be explicitly optimized, the molecule path must point to a z-Matrix and the first line in the definition is the flexy keyword

    <MOLECULE>
    flexy
    MoleculePath=methanol.zmat
    </MOLECULE>

It should be noted that the zmat-format chosen is of form

    6

    C
    O  1.4   1
    H  0.95  2  109.471 1
    H  1.089 1  109.471 2  120.0 3
    H  1.089 1  109.471 2  240.0 3
    H  1.089 1  109.471 2  0.0   3

A bond matrix for the building block may be specified by the user through the `BondMatrix=` keyword. It must only be used following the `MoleculePath=` keyword in both flexible and non-flexible input. This path is relative to the location of the `.ogo` file and must point to a file containing the bond matrix definition, either as a complete matrix or as a bond list. This input option must be used if the molecule does not only contain single bonds or if the automatic bond detection (see output) fails to correctly assign bonds. Examples:

    <MOLECULE>
    flexy
    MoleculePath=water.zmat
    BondMatrix=waterbonds.mat
    </MOLECULE>

where `waterbonds.mat` (name must not end with `.list` and is otherwise not important) contains

    1 1 1
    1 1 0
    1 0 1

where 1 denotes a single bond, 0 denotes no bond. The bonds on the diagonal are arbitrary. Allowed bonds are: 0 for no bond, 1 for a single, 2 for a double, and 3 for a triple bond. We will add more bond types as backends make use of them.

Another example usage:

    <MOLECULE>
    MoleculePath=water.zmat
    BondMatrix=waterbonds.list
    </MOLECULE>

where `waterbonds.list` (files containing bond lists must end in `.list`) contains

    0 1 1
    0 2 1

do denote a single bond between atoms 0 and 1 and 0 and 2 in the molecule.

Obviously the list-based input is less verbose/redundant and should be used in general. Please note once more that the counting starts with zero! Both list-based and matrix-based input may be used in both the flexible and non-flexible context.

The geometry definition is closed using a line with a

    </GEOMETRY>

tag. Please note that the numbering of molecules for all subsequent inputs starts from 0.

This is followed by, if applicable, input concerning spins and charges. Charges are defined as follows

    <CHARGES>
    MOLECULENO;ATOMNO;CHARGE
    </CHARGES>

for example

    <CHARGES>
    0;0;1
    1;3;-2
    </CHARGES>

Above example would put a single positive charge at the first atom of the first molecule and two negative ones at the fourth atom of the second molecule. It is important to note here that, since arrays in the Java language start with 0 and not 1 as in some prehistoric languages, our counting also starts at 0.

A similar input structure exists for the spins in the system, namely

    <SPINS>
    MOLECULENO;ATOMNO;SPIN
    </SPINS>

for example

    <SPINS>
    0;2;1
    1;3;-1
    </SPINS>

where the above example would put a single spin up at atom 0 of molecule 2 and a spin down at atom 3 of molecule 1.

Should there be no spins or charges in the studied system, no spins or charges block is necessary.

Note for larger systems build from only a few building blocks: one can add inside the `<MOLECULE>` block the `MoleculeRepetitions=` keyword and set it to the number of species of this type. Please note two things: 1) the `NumberOfParticles=` setting must still match the total number of species and 2) the spins/charges/... are automatically reproduced for all species of one type. An example:

    GEOMETRY>
    NumberOfParticles=24
    <MOLECULE>
    MoleculeRepetitions=23
    MoleculePath=water.xyz
    </MOLECULE>
    <MOLECULE>
    MoleculeRepetitions=1
    MoleculePath=water.xyz
    </MOLECULE>
    </GEOMETRY>
    <CHARGES>
    0;0;0
    1;0;-1
    </CHARGES>

would put 23 neutral water molecules and 1 negatively charged one into the cluster. The first charge line is not needed.

## Constraints and explicit degrees of freedom

Constraints can be specified based on Cartesian coordinates. Although <span class="smallcaps">ogolem</span> internally allows for arbitrary constraints to be specified, the backend carrying out the local optimization might not like that. In this case, a warning should appear.

In general, the constraints are specified similar to the spins and charges, through tags of form

    <CONSTRAINTS>
    MOLECULENO;ATOMNO;COORDINATE
    </CONSTRAINTS>

where the cartesian coordinate is chosen through an integer value (0 is x, 1 is y and 2 is z). Additionally, special input exists to constrain a whole molecule or all x, y or z coordinates of a molecule. Some example are

    <CONSTRAINTS>
    0;all
    1;allx
    2;ally
    3;allz
    4;1;0
    4;2;1
    5;0;2
    </CONSTRAINTS>

The first one constraining all cartesian coordinates of molecule 0, the second only the x coordinates of molecule 1, the third the y coordinates of molecule 2 and the fourth the z coordinates of molecule 3. Then the x coordinate of atom 1 in molecule 4 is constrained, followed by y coordinate of atom 2 in molecule 4 and the z coordinate of atom 0 in molecule 5. Any logical combination of the above described input possibilities is allowed and supported by <span class="smallcaps">ogolem</span>.

Constraints are currently only supported in the following interfaces:

- DFTB+  
  in this case the constraint is actually a restraint and the
  convergence threshold needs to be increased since the gradient along
  the restraint is included in the calculation

- lotriff

- MOPAC

- Orca

Again, if your preferred level of theory does not support constraints in <span class="smallcaps">ogolem</span> and you would like to use that feature, contact us (preferably with reference input for that particular program) or send us patches.

Explicit degrees of freedom are a rather untested feature, therefore one unfortunaltely needs to expect some rough edges. Also it should be noted that you cannot put constraints and explicit degrees of freedom on the same molecule.

As said before, you need to enable flexibility on the molecule with explicit degrees of freedom in the building block input and need to speify a z-Matrix instead of cartesian coordinates. When this has been ensured, the specification of explicit degrees of freedom is relatively straightforward. Again, it includes definitions of type

    <DOF>
    MOLECULENO;ATOMNO;COORDINATE
    </DOF>

where the coordinate is a choice for the coordinate in the z-Matrix (i.e., distance, angle, dihedral). For example

    <DOF>
    0;1;dihedral
    2;3;distance
    4;2;angle
    </DOF>

where the first line would make the dihedral of atom 1 in molecule 1, the second would make the bond distance of atom 3 in molecule 2 and the third would make the angle of atom 2 in molecule 4 an explicit DoF. Alternatively, one may use the syntax

    <DOF>
    0;all

to enable explicit degrees of freedom for all DoFs of molecule 0.

Please let us know about any experiences with explicit DoFs, may they be good or bad.

## Input keywords

There are a lot of possible tunables in the <span class="smallcaps">ogolem</span> framework. Some of them are more important than others and should actually be adjusted to the studied system.

- `DebugLevel=` the debug level: controls additional output of some code
  parts. Please set to 1 or higher <span class="smallcaps">first</span>
  if you ever have unexpected failures or problems. Defaults to 0.

- `BlowBondDetect=`  
  a collision is detected if any of the non-bonded atoms in the
  structure are getting closer to each other than this blow factor times
  their added radii. Defaults to 1.2, which is a rather conservative
  choice.

- `BlowFacDissoc=`  
  the blow factor for the dissociation detection. Defaults to 3.0 which
  should be decreased for LJ systems (e.g. 2.0 seems to be a reasonable
  choice). Notice that this setting has changed since the last version!

- `BlowInitialBonds=`  
  the blow factor for the initial bond detection. Since these bonds are
  the basis for both the post sanity check and the CD, the factor
  normally should be equal to `BlowBondDetect=`. Also defaults to 1.2.

- `CellSize=`  
  how big the initial cell should be for the initial structure
  generation. Defaults to a `9;9;9` bohrs which is too small for bigger
  systems.

- `InitialFillAlgo=`  
  how the intial individuals should be generated. Defaults to
  `packrandomly`, which is a randomized packing init w/o considering
  molecular size. If that should be yielding too densly packed (or
  spherical) structures, change to `randomwithoutdd`, which is a
  randomized init w/o dissociation detection or `randomwithdd`, which is
  a randomized init w/ DD.

- `MaxIterLocOpt=`  
  how many iterations the local optimization is allowed to carry out.
  Defaults to 2000.

- `NumberOfGlobIterations=`  
  how many global optimization steps should be carried out. Defaults to
  1000 which is just sufficient for the most trivial cases. This value
  must be increased for production jobs!

- `PoolSize=`  
  one of the most important keywords. Specifies how big the genetic pool
  is allowed to be. Defaults to 100 individuals.

- `AcceptableFitness=`  
  <span id="page1:acceptableFitness"
  label="page1:acceptableFitness"></span> what fitness is acceptable and
  once reached makes the threading version of
  <span class="smallcaps">ogolem</span> finish. Set to negative infinity
  per default and should only be set for benchmarking purposes and only
  if you are know your way around the
  <span class="smallcaps">ogolem</span> program code extremely well.

- `CollisionDetection=`  
  which collision detection to use. Defaults to `simplepairwise`, which
  is a pairwise collision detection engine. Keep.

- `CrossoverPossibility=`  
  chances for crossover in the global optimization. Defaults to 1.0
  (100%).

- `DissociationDetection=`  
  which dissociation detection to use. Defaults to `dfs`, which is a
  recursive quadratic scaling algorithm. Keep.

- `DoGeometryNiching=`  
  if a niching algorithm should be used for this run. Defaults to false.

- `GeneticRecordBufferSize=`  
  number of steps till a summary of the genetic records are flushed to
  the output file. Defaults to 1000.

- `GeneticRecordsToSerial=`  
  number of steps till the genetic records are written to disk. Defaults
  to 1000.

- `GeometriesToSerial=`  
  number of steps till the pool is written to disk. Defaults to 99,
  which is too small for big systems with a lot of steps and causes an
  1/O bottleneck and should be increased then.

- `SerializeAfterNewBest=`  
  if the pool shall be serialized whenever a new best individual has
  been added. Defaults to true.

- `GeometryChoice=`  
  how the parent individuals should be selected. Three options exist:

  - `random:`  
    both parent individuals are randomly selected. With the additional
    option `stepat=X.X` (with $`0\leq`$`X.X`$`\leq 1`$) this is turned
    into a step-function selector for the first X.X percent of the
    sorted pool. With `stepat1=X.X,stepat2=Y.Y`, you can assign
    different step functions to parent1 and parent2. When no steps are
    given/desired, this corresponds to `X.X`=`Y.Y`=1

  - `fitnessrankbased:`  
    one or both parents are selected based on a Gaussian over their
    ranked fitness. Valid, comma separated options: `bothfitness`:
    choose both parents fitness based (default: no), `gausswidth=X.X`
    sets the width of the Gaussian to X.X (default: 0.5). If you want to
    select parent1 and parent2 with Gaussians of different width, simply
    use `gausswidth1=X.X` and `gausswidth2=Y.Y`. Our default selector
    with a Gaussian width of 0.05.  
    **THIS SETTING WILL VERY LIKELY CHANGE IN THE VERY NEAR FUTURE!**

  - `fitnessvaluebased:`  
    one or both parents are selected based on their fitness value.
    Valid, comma separated options: `bothfitness`: choose both parents
    fitness based (default: no), `usegausshape`: use a Gaussian shaped
    random number for the selection (default: linear), `tailprob=X.X`
    sets the tail probability to X.X (default: 0.1). The tail
    probability is defined as the probability for the last element in
    the pool to be chosen with the head probability (i.e., the
    probability for the fittest individual) being 1.0.

- `GlobOptTries=`  
  how many tries the global optimization algorithm should carry out at
  any step that cause either collision and/or dissociation before
  proceding with the next parents. Defaults to 200 which is suitable for
  bigger systems with more degrees of freedom.

- `GrowCell=`  
  whether the cell of the starting structure is allowed to grow.
  Defaults to `false`. If you are having problems with the randomized
  structure preparation, specify a bigger cell and/or set this to
  `true`.

- `IntLocOptMaxStep=`  
  the maximum step length for the internal BFGS. Defaults to 1E0 which
  seems to be reasonable.

- `IndividualsPerNicheAtMax=`  
  how many individuals are maximally allowed per niche, if niching is
  enabled. Defaults to 50.

- `MaxBondStretch=`  
  how much a bond is allowed to be streched as a factor. Defaults to 1.3
  which seems to be reasonable.

- `MolecularCDInInit=`  
  whether there should be a molecular collision detection in the
  geometry initialization. Defaults to true.

- `MolecularCoordinateMutation=`  
  how the coordinates should be mutated. Defaults to `some`, which
  always mutates a couple of them at the same time. `one` just mutates a
  single one at any time.

- `MutationPossibility=`  
  chances for mutation of an individual during the global optimization.
  Defaults to 0.05 (5%).

- `NichingAddsToStats=`  
  how many individuals shall be added to niches (if enabled) before
  statistics are printed to command line. Defaults to 50.

- `PostSanityCD=`  
  if post-local optimization sanity checks are enabled, should the
  resulting structure be checked for collisions? Defaults to false.

- `PostSanityDD=`  
  if post-local optimization sanity checks are enabled, should the
  resulting structure be checked for dissociation? Defaults to false.

- `PostSanityCheck=`  
  whether the structure should be checked for sanity after the local
  optimization. Defaults to `true` which is reasonable for most cases.

- `PrintGeomsBeforeLocOpt=`  
  whether the geometries should be printed to the command line after the
  global optimization operators, before the local optimization. Useful
  for setup purposes. Defaults to false.

- `RatioExplDoFInit=`  
  what fraction of the explicit degrees of freedom should be randomized
  in the initialization. Applied per molecule with explicit DoFs.
  Defaults to 1.0, which might be too big for big molecules with a lot
  of explicit DoFs. If <span class="smallcaps">ogolem</span> seems to be
  stalled in the init, this might be a reason.

- `SeedingPath=`  
  if the initial pool should be seeded with some structures set this to
  the folder containing these structures in the xyz-file format.
  <span class="smallcaps">ogolem</span> will try to read in all files in
  this folder so there must only be xyz files present.

- `ThreshLocOptCoord=`  
  which threshold is acceptable for the coordinates in the local
  optimization. Defaults to 1E-8 which should be kept, if one is
  uncertain.

- `ThreshLocOptGradient=`  
  which threshold is acceptable for the gradient in the local
  optimization. Defaults to 1E-8 which should be kept, if one is
  uncertain.

Some special keywords for restarts are

- `restart=`  
  Whether this run is supposed to be restarted or not. Defaults to
  `false`. If set to `true`, some other keywords need to be set as well

- `InitialFill=`  
  Defaults to `true` and defines whether we should make an initial fill.
  Will be automagically set by the restart keyword.

- `InitialLocOpt=`  
  Defaults to `true` and defines whether the read in pool should be
  locally optimized in a first step.

- `GlobalOptimization=`  
  Defaults to `true` and defines whether there should be any global
  optimization steps.

- `NamingOffset=`  
  which ID to start from. Increase to the number of steps of your
  previous run.

However, we have found it superior to use (part of) the results of a previous run as seeds for a subsequent run. Hence, these options may vanish in the future.

Within the block of keywords, one can comment lines using either `#` or `//`.

## Local optimization and external backends

<span class="smallcaps">ogolem</span> supports a variety of program packages through interfaces. We hope that this satisfies most users. If there should be a (preferably free) program that you want to see supported and/or the methods listed here to do not satisfy you, feel free to write a mail to developers@ogolem.org .

In general, most methods are easily called through input of form

    LocOptAlgo=PACKAGE:METHOD

In general, there are thoughts to change this input for something more flexible so future versions might behave different. Also it should be noted that of course <span class="smallcaps">ogolem</span> might not support any revision of the programs mentioned if their input changes. If you should hit such a case and are not working with some ancient version (e.g. there is for sure no PM6 in MOPAC7), please contact us and wait for a fix.

### Interface to ADFSuite

An experimental interface to the ADFSuite of programs exists. It is chosen by setting `LocOptAlgo=adf:JOBKEY`, where `JOBKEY` can either be a build-in job identifier, such as `GO` or `DFTB-GO`, or the file name of an adf stub file, e.g., `mylocopt.adf`. In the later case, this file name must end in `.adf` and be present in the working directory.

A stub file can be created using, e.g., `adfjobs` and tuned to the particular computational needs of the system under study.

**IMPORTANT:** <span class="smallcaps">ogolem</span> will **NOT** adjust anything except for the geometry coordinates in the job stub file. Therefore, a charged system or a system with spin does require creation of a fitting job stub file. Likewise, convergence thresholds will not be propagated from within <span class="smallcaps">ogolem</span> to ADF and must be set in the job stub file!

Please note, you must have all the usual ADF/SCM environment variables properly setup, e.g., `SCMLICENSE` and `ADFBIN`. This can be done by sourcing the correct setup file for your platform. Obviously, a valid license for the codes to be run is required.

### Interface to AMBER

The interface to AMBER is called through `LocOptAlgo=amber`, followed by the optional choice `:implicitwater`, if implicit solvation is wanted or needed.

The `sander` executable needs to be in the path and named such.

Additionally, there needs to be a fitting amber parameter file named `input.prmtop`, where `input` is the prefix of your `input.ogo` file.

The AMBER interface in its current state is not suitable for molecule only optimizations. This does currently not have any effect on <span class="smallcaps">ogolem</span>.

### Interface to CP2K

The interface to CP2K is called through `LocOptAlgo=cp2k:`, followed by a sequence of keywords to specify the level of theory. Supported levels of theory at the moment are

- `B3LYP`

- `BLYP`

- `BP`

- `HCTH120`

- `OLYP`

- `PADE`

- `PBE`

- `PBE0`

- `TPSS`

The `cp2k.sopt` executable needs to be in the path and named such.

Additionally, an `input-cp2k.aux` file is required, where `input` is the prefix of your `ogo`-file. This auxiliary file needs to start with a line

    ###OGOLEMAUX###

followed by a range of lines specifying in three whitespace-separated columns which atom should use which basis set and which potential. These potentials and basis sets need to be present in the directory.

### Interface to CPMD

The interface to CPMD is called through `LocOptAlgo=cpmd:`, followed by a sequence of keywords to specify the level of theory. Supported levels of theory at the moment are

- `SONLY`

- `LDA`

- `BONLY`

- `BP`

- `BLYP`

- `XLYP`

- `GGA`

- `PBE`

- `PBES`

- `REVPBE`

- `HCTH`

- `OPTX`

- `OLYP`

- `TPSS`

- `PBE0`

- `B1LYP`

- `B3LYP`

- `X3LYP`

- `CUSTOM`

The `cpmd.x` executable needs to be in the path and named such.

Additionally, an `input-cpmd.aux` file is required, where `input` is the prefix of your `ogo`-file. This auxiliary file needs to start with

    ###OGOLEMAUX###
    70
    15;15;15

Where the second line is the cutoff in Rydberg and the third line is the cell size to be used (if line is empty it will be estimated for each geometry, based on the assumption of a clustered system!). This is followed by a range of lines specifying in three whitespace-separated columns which atom should use which pseudopotential and which lmax. The pseudopotentials need to be present in the directory. Please note that no extra lines are allowed.

### Interface to DFTB+

The interface to DFTB+ is called through `LocOptAlgo=dftb+:`, followed by a sequence of keywords to specify some needed parameters for the local optimization employed by DFTB+

- `steepest`  
  use a steepest descent algorithm

- `conjugate`  
  use a conjugate gradient algorithm

- `steepest,scc`  
  use a steepest descent algorithm with SCC approximation

- `conjugate,scc`  
  use a conjugate gradient algorithm with SCC approximation

The DFTB+ executable needs to be in your path and named `dftb+`.

Additionally, you need to make sure to have all Slater-Koster files for all interactions in your system in the directory your `ogo`-input file is. This is needed for <span class="smallcaps">ogolem</span> using the SK-files for the DFTB+ input.

### Interface to GROMACS

The interface to GROMACS is called through `LocOptAlgo=gromacs`. It is followed by an optional choice of the QM backend to be used. Choices are

- empty/whitespace  
  no QM backend

- `:mopac`  
  use MOPAC7 as a backend. Please note that MOPAC7 only provides AM1 and
  PM3.

- `:orca`  
  use Orca as a backend.

- `:gaussian`  
  use Gaussian as a backend.

- `:gamess-uk`  
  use Gamess-UK as a backend.

In all cases, the GROMACS binaries need to be recompiled according to the GROMACS documentation depending on the chosen backend.

The path to the `grompp` and `mdrun` excecutables needs to be exported through the environment variables `OGO_GROMACSMPP` and `OGO_GROMACSCMD`, e.g. in bash syntax

    export OGO_GROMACSMPP=$HOME/software/gromacs/grompp
    export OGO_GROMACSCMD=$HOME/software/gromacs/mdrun

Additionally, the files

- `gromacs-topology.top`,

- `gromacs-parameters.mdp` and

- `gromacs-index.ndx`

need to present and correct for the system under study. These files are also responsible for all GROMACS related configuration, e.g. constraints and QM/MM partitioning. The file

    INPUTNAME-gromacs.aux

needs to contain in the first line the identifier `###OGOLEMAUX###` followed by the first half of a `.g96` file of the system (residue and atom-IDs).

If one uses Orca as a QM-backend for Gromacs, one also needs to export the path to the orca excecutables and the basename of the optimization job. This is e.g. in bash syntax accomplished through

    export ORCA_PATH=$HOME/software/orca
    export BASENAME=topol

where the `BASENAME` is always `topol`! Additionally, the file `gromacs-orca.aux` needs to contain the input keywords, as described in the GROMACS manual.

### Interface to lotriff

Lotriff is an in-house development by Bernd Hartke. If you are in desperate need of this FF engine, contact him via mail. If you obtained the `lotriff` executable you can use `LocOptAlgo=lotriff`.

The path to the lotriff executable needs to be exported in the `OGO_LOTRIFFCMD` environment variable, e.g. in bash syntax

    export OGO_LOTRIFFCMD=$HOME/software/lotriff/lotriff

If this variable should be not set, <span class="smallcaps">ogolem</span> will assume the binary to be in your path.

Additionally, a correct `lotriff.in` file and a correct `hetero.def` file need to be present in your directory, the constrained building blocks are supposed to be specified first in the building block input.

A mapping of atom types in an auxiliary file `input-lotriff.aux`, where input is the prefix of your input file name, needs to be present. First line of this file needs to be

    ###OGOLEMAUX###

followed by a column of force field types as present in your lotriff input.

The lotriff interface in its current state is not suitable for molecule only optimizations. This does currently not have any effect on <span class="smallcaps">ogolem</span>.

### Interface to MNDO

The interface to MNDO is called through `LocOptAlgo=mndo:`, followed by a sequence of keywords to specify the method. Builtin support exists at the moment for the following levels of theory

- `mndod`

- `om1`

- `om2`

- `om3`

- `pm3`

- `am1`

- `mndoc`

- `mndo`

- `mindo3`

- `cndo2`

- `scc-dftb`

- `scc-dftbj`  
  Jørgensen DFTB

all of the above can make use of implicit solvation using COSMO when appending `:cosmo`.

Since the MNDO program is rather unique when it comes to input and output, it was ultimately decided to call not MNDO directly but instead to always call a shell-script which just wraps that. This obviously just works on UNIX machines but MNDO is also just available for them to our knowledge.

The shell script, named `mndo.sh`, needs to be in the path and include something alike to the following

    #!/bin/sh
    $1 > mndo > $2

### Interface to MOLPRO

The interface to MOLPRO is called through `LocOptAlgo=molpro:`, followed by a sequence of keywords to specify the method. Builtin support exists at the moment for the following levels of theory

- `hf/vdz`

- `b86/vdz`

- `mp2/avtz`

- `custom`  
  Requiring you to provide an input stub `YOURINPUT-molpro.aux`
  containing all relevant input and instead of the geometry the key
  `OGOCOORDS`.

- `custom,nolocopt`  
  Just like above, an auxiliary file is required, with the `nolocopt`
  flag signaling that no local optimization takes place, interesting for
  chained local optimizations.

Additionally, the location of the MOLPRO executable needs to be exported in the environment variable `OGO_MOLPROCMD`, e.g. in bash syntax

    export OGO_MOLPROCMD=$HOME/software/molpro2010.1/bin/molpro

If this variable should be not set, <span class="smallcaps">ogolem</span> will assume the binary to be in your path.

Since MOLPRO is a closed-source program it is only supported on the platforms the official MOLPRO release supports.

### Interface to MOPAC

The interface to MOPAC is called through `LocOptAlgo=mopac:`, followed by a sequence of keywords to specify the method. Builtin support exists at the moment for the following levels of theory

- `mndo`

- `am1`

- `pm3`

- `pm5`

- `pm6`

- `mndo,cosmo(ethanol)`

- `am1,cosmo(ethanol)`

- `pm3,cosmo(ethanol)`

- `pm5,cosmo(ethanol)`

- `pm6,cosmo(ethanol)`

- `mndo,cosmo(water)`

- `am1,cosmo(water)`

- `pm3,cosmo(water)`

- `pm5,cosmo(water)`

- `pm6,cosmo(water)`

Additionally, the location of the MOPAC executable needs to be exported in the environment variable `OGO_MOPACCMD`, e.g. in bash syntax

    export OGO_MOPACCMD=$HOME/software/mopac2010/mopac

If this variable should be not set, <span class="smallcaps">ogolem</span> will assume the binary to be in your path.

Since MOPAC is a closed-source program it is only supported on the platforms the official MOPAC release supports.

### Interface to NAMD

The interface to NAMD is called through `LocOptAlgo=namd`.

The `namd2` executable needs to be in the path and named such.

Additionally, there needs to be a fitting amber parameter file (yes, this is not an error!) named `input.prmtop`, where `input` is the prefix of your `input.ogo` file.

The NAMD interface in its current state is not suitable for molecule only optimizations. This does currently not have any effect on <span class="smallcaps">ogolem</span>.

### Interface to OpenBabel

The interface to OpenBabel is called through `LocOptAlgo=openbabel:`, followed by a sequence of keywords to specify which force field to use. Currently supported are

- `ghemical`  
  the usage of this FF is highly encouraged since it is finnish
  technology ;-)

- `mmff94s`

- `mmff94`

- `uff`

The `obminimize` executable needs to be in the path and named so.

As a side note, the atom type recognition is done by OpenBabel automagically. This might be error prone and it is recommended to be checked in advance to a full global optimization.

### Interface to Orca

The interface to Orca is called through `LocOptAlgo=orca:`, followed by a sequence of keywords to specify the method. Builtin support exists at the moment for the following levels of theory

- `mndo`

- `am1`

- `pm3`

- `b3lyp/vdz`

- `bp86/svp`

- `bp86/tzvp`

- `ubp86/tzvp`  
  forces unrestricted Kohn-Sham

- `b3lyp/6-31+g**`

- `mp2/aug-cc-pVDZ`

- `ccsd(t)/tzvpp`  
  only use this if you have way too much computing time

- `custom`  
  Requiring you to provide an input stub `YOURINPUT-orca.aux` containing
  all relevant input before the geometry block. Please note, that this
  input needs to contain the constraints, if wanted.

- `custom,nolocopt`  
  Just like above, an auxiliary file is required, with the `nolopt` flag
  signaling that no local optimization takes place, interesting for
  chained local optimizations.

Additionally, the location of the orca executable needs to be exported in the environment variable `OGO_ORCACMD`, e.g. in bash syntax

    export OGO_ORCACMD=$HOME/software/orca_v2.8.0/orca

If this variable should be not set, <span class="smallcaps">ogolem</span> will assume the binary to be in your path.

Since Orca is a free, albeit closed-source program, it is only supported on the platforms the official Orca release supports (at the moment: Linux, Solaris, Windows).

### Interface to Tinker

The interface to Tinker is called through `LocOptAlgo=tinker:`, followed by a sequence of keywords to specify some needed parameters for the local optimization employed by Tinker

- `minimize`  
  use Tinkers minimize program for local optimization

- `newton`  
  use Tinkers newton program for local optimization

- `optimize`  
  use Tinkers optimize program for local optimization

- `custom`  
  use a custom command and options for the local optimization

For the earlier three, the Tinker executables need to be in your path and named `minimize`, `newton` and `optimize`. For the last option, one needs to export the system variables `OGO_TINKERCMD` and `OGO_TINKEROPT`, e.g. in bash syntax

    export OGO_TINKERCMD=$HOME/software/tinker_patched/optimize
    export OGO_TINKEROPT="0.01"

Additionally, you need to have a file called `input-tinker.aux`, where `input` is the prefix of your input file in the same directory as the `ogo`-file. This file needs to start with a line

    ###OGOLEMAUX###

followed by a line with the name of the tinker parameter file to use (e.g. `mm3.prm`), which also needs to be in the directory. Then the force field types for all atoms in your systems (same order as in the building block definition). The first column being some arbirtray identifier, followed by whitespace and the force field type ID as in the parameter file.

Occasionally, the setup described above does not work but results in error messages like this:

    WARNING: Error in Tinker: optimize locopt. Returning non-optimized geometry. 
    org.ogolem.core.ConvergenceException: Tinker has a problem (local optimization).

If this happens, first make sure that Tinker in stand-alone mode (w/o <span class="smallcaps">ogolem</span>) can really do the intended calculations, using the input files automatically generated by <span class="smallcaps">ogolem</span>. If Tinker does work in stand-alone mode, the above error may be related to system I/O issues on certain Linux systems. As a possible fix, you are then encouraged to try this workaround:

- replace the direct invocation of the desired Tinker routine (e.g.,
  `optimize`), with a script of this contents:

      #!/bin/sh
      base=`basename $1 .xyz`
      out=$base".out"
      /path/to/your/tinker/executables/optimize $@ > $out

  For example, you can name this script “optimize” and put it in a   directory that occurs early in your `PATH`-variable (e.g.,   `$HOME/bin`), before the normal Tinker executables (to be on the safe   side, remove the Tinker executables from your `PATH`). Alternatively,   you can also let `OGO_TINKERCMD` point to this script.

- additionally, in the $`\ast`$.ogo input file, call the Tinker
  interface with

      LocOptAlgo=tinker:optimize,fileout

This (or obvious variations of this recipe for the other Tinker routines) may get rid of the above error.

The Tinker interface in its current state is not suitable for molecule only optimizations. This does currently not have any effect on <span class="smallcaps">ogolem</span>.

### Interface to xtb

The interface to xtb is called through `LocOptAlgo=xtb:`, followed by a sequence of optional keywords adjusting its behavior.

- `method=`  
  Set the xtb method to be used. Options are: `gfn-ff`, `gfn0-xtb`,
  `gfn1-xtb`, `gfn2-xtb`. Default is GFN2.

- `optlevel=`  
  Set the convergence level of local optimizations in xtb. Options are:
  `crude`, `sloppy`, `loose`, `normal`, `tight`, `verytigth`. Default is
  `normal`.

- `dontserialize`  
  Do not set the environment to force xtb to only use one thread.

- `xcontrol=`  
  Path to an xtb control file to be used for the local optimizations.
  Path is relative to <span class="smallcaps">ogolem</span> input file
  location. Default: no additional input file is used.

### Interfacing to **any** other program (package)

Although <span class="smallcaps">ogolem</span> provides a numerous set of interfaces to programs implementing methods ranging from force fields to wavefunction methods, we cannot implement an interface to every code that we find interesting. Additionally, if one of the interfaced codes changes input or output structure in a non-compatible manner, we will need to chase that.

But fear not as there is an easy solution to interface **any** other code to <span class="smallcaps">ogolem</span> with a little bit of simple scripting. This interface can be requested by setting `LocOptAlgo=generic`. Additionally, the environment variables `OGO_GENERALCMD` will need to be set to the command to be called for the actual local optimization. Additionally, if `OGO_GENERALOPTS` are set, these options will be appended to the called command.

If this interface is chosen, <span class="smallcaps">ogolem</span> will create individual directories for each local optimization task. In such directory, a file `input.xyz` within the initial guess coordinates in xyz-format will be created. Subsequently, the `OGO_GENERALCMD` will be executed within the directory (optinally with the options appended) which is supposed to read in the coordinates, run the local optimization and return the results in an xyz-formatted file `output.xyz`. In this file, the second line is supposed to contain the fitness of the optimized structure in kJ/mol as the second token.

Happy <span class="smallcaps">ogolem</span> scripting!

### Chaining multiple local optimizations

Any number of local optimizations can be chained. This is easily achieved by specifying `LocOptAlgo= chained:` followed by a -separated (pipe-separated) list of the local optimization options (inluding their possible options) to be chained.

### Build-in local optimizations

A variety of local optimization algorithms are available build-in in <span class="smallcaps">ogolem</span>. They all have in common that they require a backend specification as described in the following sections.

Please note that not all of the local optimization algorithms now described are classically seen as local optimization but some of them are used in other groups as global optimizations engines. Despite the semantics, we are convinced that coupling mutliple global optimization engines can be benefical. No bigger offense should be taken by naming them local optimizations.

A word on syntax: unless otherwise noted, options to local optimizations are semicolon separated. Classic choices for the local optimizations include:

- `chained:XXX|YYY` which chains the local optimizations `XXX` and
  `YYY`. Any number of local optimizations can be chained and are *pipe*
  \| separated.

- `none:backend=XXX` which uses the backend `XXX` for single-point only
  evaluations.

- `lbfgs:` to use a L-BFGS optimization through RISO’s version of the
  original Nocedal implementation. Recommended choice. Options are:

  - `backend=XXX` uses backend `XXX` for this optimization. Mandatory
    option.

  - `maxiter=XXX` sets the maximal number of iterations for this local
    optimization to `XXX`. Default: global default. Optional input.

  - `lineiter=XXX` sets the number of line search iterations to `XXX`.
    Default: 50. Optional input.

  - `gtol=X.X` sets the gradient tolerance to `X.X`. Default: 0.9.
    Optional input.

  - `nocorrs=X` sets the number of corrections to `X`. Default: 7.
    Optional input.

  - `absconvthresh=X.XX` sets the absolute convergence threshold to
    `X.XX`. Default: global default. Optional input.

  - `tryresets=X` sets the number of resets of the L-BFGS state to `X`.
    Default: 3. Optional input.

- `apachecg:` to use Apache’s conjugate gradient local optimization,
  options include:

  - `backend=XXX` uses backend `XXX` for this optimization. Mandatory
    option.

  - `maxiter=XXX` sets the maximal number of iterations for this local
    optimization to `XXX`. Default: global default. Optional input.

  - `usefletcher=true` enables Fletcher-Reeves updates over
    Polak-Ribiere updates. Default: Polak-Ribiere. Optional input.

  - `absconvthresh=X.XX` sets the absolute convergence threshold to
    `X.XX`. Default: global default. Optional input.

- `fire:` to use the Fast Inertial Relaxation Engine (FIRE). Options
  are:

  - `backend=XXX` uses backend `XXX` for this optimization. Mandatory
    option. Note that FIRE benefits greatly from backends that support
    history functionality.

  - `maxiter=XXX` sets the maximal number of iterations for this local
    optimization to `XXX`. Default: global default. Optional input.

  - `fmax=X.X` sets the convergence criterion to `X.X`. Default: global
    default. Optional input.

  - `maxmove=X.X` sets the maximum movement in one iteration to `X.X`.
    Default: 0.2. Optional input.

  - `dt=X.X` sets the initial time step size to `X.X`. Default: 0.1.
    Optional input.

  - `dtmax=X.X` sets the maximum time step size to `X.X`. Default: 1.0.
    Optional input.

  - `nmin=X` sets the number of steps without step increase to `X`.
    Default: 5. Optional input.

  - `finc=X.X` sets the factor for time step increase to `X.X`. Default:
    1.1. Optional input.

  - `fdec=X.X` sets the factor for time step decrease to `X.X`. Default:
    0.5. Optional input.

  - `astart=X.X` sets the starting acceleration mixing factor to `X.X`.
    Default: 0.1. Optional input.

  - `fa=X.X` sets the factor for acceleration mixing increase to `X.X`.
    Default: 0.99. Optional input.

  - `a=X.X` sets the initial acceleration mixing factor to `X.X`.
    Default: 0.1. Optional input.

  - `tryresets=X` sets the maximum number of resets if too many trials
    fail to improve the result to `X`. Default: 2. Optional Input.

  - `maxtrials=X` sets the number of trials without improvement before a
    reset is initiated to `X`. Default: 20. Optional input.

  - `resettostable=false` turns off resets to stable states. Default:
    true. Optional input.

  - `resettobestpoint=true` turns on resets to remembered points in
    history. Requires a backend that supports history functionality.
    Default: false. Optional input.

- `flemin:` to use NUMAL’s flemin optimization. Options:

  - `backend=XXX` uses backend `XXX` for this optimization. Mandatory
    option.

  - `maxiter=XXX` sets the maximal number of iterations for this local
    optimization to `XXX`. Default: global default. Optional input.

  - `convthresh=X.XX` sets the convergence threshold in the fitness to
    `X.XX`. Default: global default. Optional input.

  - `convthreshgrad=X.XX` sets the convergence threshold in the gradient
    to `X.XX`. Default: global default. Optional input.

- `rnk1min:` to use NUMAL’s rnk1min optimization. Options:

  - `backend=XXX` uses backend `XXX` for this optimization. Mandatory
    option.

  - `maxiter=XXX` sets the maximal number of iterations for this local
    optimization to `XXX`. Default: global default. Optional input.

  - `convthresh=X.XX` sets the convergence threshold in the fitness to
    `X.XX`. Default: global default. Optional input.

  - `maxstep=X.X` sets the maximal step to `X.XX`. Default: 0.1.
    Optional input.

  - `convthreshgrad=X.XX` sets the convergence threshold in the gradient
    to `X.XX`. Default: global default. Optional input.

<!-- -->

- 

Please note that the `flemin`, `praxis`, `rnk1min` algorithms may only be used if you are in possession of Lau’s Java Numal book(Lau 2003) (and therefore the accompanying source code and jar-file). In this case, you can add the numal.jar to the classpath prior to <span class="smallcaps">ogolem</span> excecution and make use of the optimization algorithms implemented in Numal.

### Local heat pulses

One of the global optimization algorithms included in <span class="smallcaps">ogolem</span> are local heat pulses, as originally developed by Arnulf Moebius. The version in <span class="smallcaps">ogolem</span> is fairly customized in comparison to the original, and was described and tested in detail (Dieterich and Hartke 2017). Local heat pulses can be used either as a geometry mutation (cf. section <a href="#sec:geomglobopt" data-reference-type="ref" data-reference="sec:geomglobopt">2.7</a>), or as pre-stage to the actual local optimization. In the former case, it will come into effect only with the mutation probability (and possibly as one of several mutations); in the latter case it will influence every newly generated geometry. Both choices lead to the same algorithmic machinery used in the end, hence the options are the same (described below). There are quite a few options because the originally proposed heat pulses are a full-blown global optimization algorithm in themselves, borrowing features from MCM/basin-hopping and simulated annealing. The latter features lead to annealing-style options.

You can specify local head pulses as local optimization pre-stage by inserting `localheat:` with its options, followed by the actual local optimization (may be any of the ones available in <span class="smallcaps">ogolem</span>) after a “;”, in this way: `localheat:OPTIONS;lbfgs:backend=xyz:lennardjones`. The options are keyword=value pairs, from the following list:

- `choosemode=` the mode for the local heat pulses. We support the
  following choices: (default `pick5`)

  - `pick5`: move 5 randomly choosen atoms/COMs

  - `upto5`: move 0 to 5 atoms/COMs

  - `10percent`: move 10% of the atoms/COMs

  - `upto10percent`: move 0% to 10% of the atoms/COMs

  - `onsphere`: move atoms/COMs within a 3D-Gaussian centered at a point
    on a spherical surface enclosing the cluster

  - `insphere`: move atoms/COMs within a 3D-Gaussian centered somewhere
    inside the whole cluster

  - `incenter`: move atoms/COMs within a 3D-Gaussian around the cluster
    center

- `doCDDD=` if we should check the geometry’s sanity after heat pulse
  perturbation, using CD/DD. Default: false.

- `movemode=` `coms` for moving the COMs of the building blocks, or
  `cartesians` for moving atoms. Default: `coms`

- `iters=` the number of iterative repetitions of the local heat pulses.
  Default: 10.

- `eqIters=` the number of equilibration iterations before cooling the
  movement amplitudes. Of course, we should normally have `iters` $`>`$
  `eqIters`. Default: 3.

- `startAmpl=` what the start amplitude for atom/COM perturbation should
  be. Default: 2Å

- `startEuler=` the initial amplitude for Euler angle perturbation.
  Default: 0.5.

- `scaleFac=` how we scale amplitudes, temperature and Euler angle
  amplitudes in annealing-style cooling. Default: 0.9.

- `useMetropolis=` if we use the Metropolis criterion (also allowing
  energy-upward moves) instead of the simplistic downward-only
  criterion. Default: true (use Metropolis).

- `temperature=` when we use a Metropolis criterion, what the start
  temperature (in Kelvin) should be. Default: 50.

- `sigmaX=`,`sigmaY=`,`sigmaZ=` the x,y,z-width of the 3D Gaussian in
  the modes `onsphere`, `insphere` and `incenter`. Default: 0.1.

## Internal and external backends

As discussed above, all internal local optimizations require a backend definition. This backend definition can be either directly inlined or pre-defined and tagged. Please note that we prefer the pre-defined version discussed here over direct inlining.

In order to define a backend a `CLUSTERBACKEND` block must be in the input. Options include:

- `BackendTag=XYZ` to name this backend defintion `XYZ`. Tag must be
  unique. Mandatory option.

- `Backend=XXX` to provide the fully qualified backend definition.
  Mandatory option.

As an example:

    <CLUSTERBACKEND>
    BackendTag=lj
    Backend=xyz:lennardjones
    </CLUSTERBACKEND>

would be used to assign a Cartesian-based LJ force-field the tag `lj`. Subsequently, the tag `lj` could be used to refer to this backend in the local optimization definitions.

Actual backend definitions will be covered in the next sections.

#### Internal LJ force field

Very simple, based on standard parameters for *homogenous* LJ clusters of noble gases (excluding radon). Again, this does *not* include any mixing and can therefore only be used for homogenous clusters. Specify `lennardjones` as the backend.

#### Internal mixed LJ force field

Also simple, based on standard parameters for *heterogenous* LJ clusters of noble gases. Through the use of Lorentz-Berthelot mixing rules, any mixed cluster of noble gases (excluding radon) can be studied. Specify `mixedlj` as the backend.

#### scaTTM3F

The <span class="smallcaps">ogolem</span> contains a Scala-based implementation of the TTM3F force field by Xantheas *et al.*. Despite being significantly faster than the reference implementation, this force field is equivalent in behaviour to said original. To use it, simply specify `scattm3f:X,Y.Y`, where `X` is the number of water molecules your system contains (our implementation also pre-allocates some scratch arrays) and `Y.Y` is an energy cutoff in $`E_h`$. This cutoff is of importance, as (like the reference implementation) the induced dipols may converge to unphysical solutions which can be very low in energy. It is recommended to scale the cutoff with system size.

Please also note that the implementation assumes your cluster to only contain water molecules with internal atom order O/H/H.

#### (EVB-)QMDFF

Thanks Stefan Grimme, we include his QMD force field natively in <span class="smallcaps">ogolem</span>. Additionally, we support the EVB-QMDFF extension by Bernd Hartke and Stefan Grimme. Before running <span class="smallcaps">ogolem</span>, a force field parameter must be created file with the official tools from Stefan Grimme. Subsequently, QMDFF can be configured as `qmdff:` followed by a , (comma) separated list of options:

- `params=XXX` sets the path to the parameter file to `XXX`. Default:
  `solvent`.

- `addterm=XXX` add an additional term to QMDFF. Currently only
  supported: `electrostatic:PATHTOCHARGEFILE`. Default: no additional
  term.

To configure EVB-QMDFF, specify `evb-qmdff:` followed by any of the following options:

- `evb=`

- `qmdffs=`

Please note that we are in the testing phase for both EVB-QMDFF and QMDFF for structure global optimization and cannot make any guarantees as to the physical sensibility of obtained results!

#### Adaptive choices

Adaptive choices are in more detail described in section <a href="#parameterfit" data-reference-type="ref" data-reference="parameterfit">[parameterfit]</a> *Global parametrization of potentials*.

## Global optimization algorithms

All global optimization algorithms are chosen through the `GlobOptAlgo=` keyword. The grammar is as follows: `cluster`$`\{`$`xover(GEOMXOVERALGO:OPTS)mutation(GEOMMUTATIONALGO:OPTS)`$`\}`$`molecules`$`\{`$` xover(MOLXOVERALGO:OPTS)mutation(MOLMUTATIONALGO:OPTS)`$`\}`$. The cluster algorithm specification is required, the molecule specification (for explicit degrees of freedom) is optional. Algorithmic choices must end with a colon, options are optional and comma-separated.

Valid choices for `GEOMXOVERALGO` are:

- `germany:` a very basic one-point genotype crossover. Options:

  - `gausswidth=X.X` sets the width of the Gaussian deciding for the
    cutting spot to X.X, default 0.3.

- `portugal:` a $`N`$-point genotype crossover.

  - `nocuts=X` sets the number of cuts to `X`. Default: 1.

- `sweden:` our $`N`$-species phenotype crossover.

  - `cutstyle=` specifies the cutting style, valid choices 0, 1 or 2. 0:
    cutting plane is always through the total COM of the cluster. 1:
    Gauss-distributed cutting plane. 2: normal-distributed cutting
    plane. Default: 2.

- `lapland:` a phenotype crossover with implicit exchange.

  - `cutstyle=`, valid choices 0, 1 or 2. Default: 2.

- `iceland:` a merging phenotype operator.

  - `cutstyle=X`, valid choices 0, 1 or 2. Default: 2.

  - `colldetect=X` Choses a collision detection engine. See global
    option for valid choices. Default: global setting.

  - `blowcoll=X.X` specifies the collision blow factor to be `X.X`.
    Default: global setting.

  - `blowdiss=X.X` specifies the dissociation blow factor to be `X.X`.
    Default: global setting.

  - `inittrust=X.X` specifies the initial trust radius for the merge.
    Default: 0.1.

  - `stoptrust=X.X` specifies the stopping trust radius for the merge.
    Default: 1E-5.

  - `mergeiter=X` specifies the number of iteration for the merge.
    Default: 100.

  - `nointerpoints=X` specifies the number of interpolation points for
    the merging. Default: 5.

  - `optbounds=` specifies the optimization bounds for the merging as a
    $`/`$-separated array. Order: lower translation bound, higher
    translation bound, lower rotation bound, higher rotation bound.
    Defaults: `0.0/15.0/0.0/2*pi`

- `aland:` a fitness-based phenotype operator with optional merging.

  - `cutstyle=X`, valid choices 0, 1 or 2. Default: 2.

  - `mergingpheno`, enables merging. Default: off.

  - `colldetect=X` Choses a collision detection engine. See global
    option for valid choices. Default: global setting.

  - `blowcoll=X.X` specifies the collision blow factor to be `X.X`.
    Default: global setting.

  - `blowdiss=X.X` specifies the dissociation blow factor to be `X.X`.
    Default: global setting.

  - `inittrust=X.X` specifies the initial trust radius for the merge.
    Default: 0.1.

  - `stoptrust=X.X` specifies the stopping trust radius for the merge.
    Default: 1E-5.

  - `mergeiter=X` specifies the number of iteration for the merge.
    Default: 100.

  - `nointerpoints=X` specifies the number of interpolation points for
    the merging. Default: 5.

  - `optbounds=` specifies the optimization bounds for the merging as a
    $`/`$-separated array. Order: lower translation bound, higher
    translation bound, lower rotation bound, higher rotation bound.
    Defaults: `0.0/15.0/0.0/2*pi`

- `snaefellsjoekull:`

  - `cutstyle=X` specifies a cut style. Allowed: 0 (full random), 1
    (Gauss-distributed), 2 (inverted Gauss-distributed). Default: 0.

  - `gausswidth=X.X` the standard deviation if a Gauss-cut is chosen.
    Default: 1.0.

  - `adjustradius=true/false` if the radius of the spheres should be
    adjusted. Default: false.

- `noxover:` disables crossovers in the global optimization.

- `mutationasxover:GEOMMUTATIONALGO:OPTS` uses the mutation specified by
  `GEOMMUTATIONALGO` as a crossover.

- `multiple:XX%GEOMXOVER1:OPTS|YY%GEOMXOVER2:OPTS...` couples multiple
  crossover operators together. Each definition starts with the percent
  probabilty for this operator to be used. Different operators are
  separated by pipes \|. The percentages must add up to 100%.

- `chained:XX%GEOMXOVER1:OPTS|YY%GEOMXOVER2:OPTS...` chains multiple
  crossover operators together. Each definition starts with the percent
  probabilty for this operator to be used in the chain. Different
  operators are separated by pipes \|. Each percentage is independent.

We do recommend to choose `sweden` if one is not certain what the best crossover is.

Valid choices for `GEOMMUTATIONALGO` are:

- `germany:` for a basic genotype mutation of the COMs and rotations for
  the molecules.

  - `lowcom=X.X` the minimal translation for the COMs. Default: 0.5
    bohr.

  - `highcom=X.X` the maximal translation for the COMs. Default: 1.5
    bohr.

  - `widthCOM=X.X` the width of the Gaussian centered between the two.
    Default: 0.1.

  - `loweuler=X.X` the minimal rotation of the COMs. Default: 0.0 (in
    the specific Euler angles definition.

  - `higheuler=X.X` the maximal rotation of the COM. Default: 1.0 (in
    the specific Euler angles definition.

  - `widtheuler=X.X` the width of the Gaussian centered between the two
    Euler borders. Default: 0.1.

- `montecarlo:` for a Monte-Carlo based mutation in the atomic space.

  - `mode=` to adjust the move mode. Valid choices: `all`, `some`,
    `one`.

  - `maxmove=X` specifies the maximum move. Default: 0.2 (in bohr).

- `finland:` for a one-parent phenotype crossover.

  - `cutstyle=X`, valid choices 0, 1, 2. Default: 2.

  - `blowcoll=X.X` to set the blow factor for the collision detection.
    Default: global setting.

  - `blowdiss=X.X` to set the blow factor for the dissociation
    detection. Default: global setting.

- `norway:` for a one-parent packing operator.

  - `colldetect=` to set the collision detection algorithm. Default:
    global setting.

  - `dissdetect=` to set the dissociation detection algorithm. Default:
    global setting.

  - `blowcoll=X.X` to set the blow factor for the collision detection.
    Default: global setting.

  - `blowdiss=X.X` to set the blow factor for the dissociation
    detection. Default: global setting.

- `xchangemut:` for an explicit exchange mutation.

  - `mode=` the exchange mode. Valid choices: `single`, `multiple`.
    Default: `single`.

  - `gausswidth=X.X` if multiple is chosen, the width of the Gaussian
    used to determine how many should be exchanged. Default: 0.4.

- `graphbasedmut:` the original GDM doing a graph-based analysis of the
  cluster to determine the least connected building block and move it to
  a better position.

  - 

- `localheat:` local heat pulses can be selected as local
  pre-optimization or as mutation, cf. section
  <a href="#sec:localheat" data-reference-type="ref"
  data-reference="sec:localheat">2.5.19</a>

To be continued...

### HACTAR

HACTAR is...... **not yet ready for prime time. Please do not use!**

Input syntax: `hactar`$`\{`$`HACTAROPTIONS`$`\}`$`cluster`$`\{`$`xover()mutation()`$`\}`$. The cluster specifications for both `xover()` and `mutation()` can (and should) contain more than one algorithmic definition. These definitions follow the `multiple:` syntax explained above, without the leading `multiple:` keyword.

Currently no molecular HACTAR exists. `HACTAROPTION` are a series of optional, comma-separated options. The following input keywords exists:

- `mode=` Wich moving mode to use. Possible modes are `all`, `some`,
  `one`. Default: `some`.

- `dometropolis` To enable the Metropolis criterion. Default: accept
  steps that lower the energy, no Metropolis criterion.

- `beta=X.X` If the Metropolis criterion is enabled, which beta factor
  to use. Default: 1.0.

- `maxtoreset=X` After how many unsuccessful operations a reset to the
  initial algorithmic configuration should take place. Default: 100.

- `stepstopring=X` After how many steps we print the current algorithmic
  configuration to standard out. Default: 1000.

- `maxpercstep=X.X` What is the maximum percentage the configuration is
  allowed to be moved in one step. Default: 1%.

- `dopercentages` evaluate the improvements as a percentage, not as the
  absolute.

## Diversity checkers

In order to maintain a sufficient diversity in the genetic pool, we use so called diversity checkers. Their quality (and performance) are crucial for the success of your global optimization.

The default is a fitness-based diversity criterion with a threshold of $`1\cdot10^{-6}\:E_h`$. This threshold may be too tight for your optimization and is related to the convergence threshold of the local optimization. To change the default, the `DiversityCheck=` keyword must be used. Options are:

- `fitnessbased:X.X`  
  which is the simple fitness-based diversity checker with a threshold
  of `X.X`.

- `percfitnessbased:X.X`  
  which is the simple fitness diversity checker based on the percentage
  difference w.r.t. the individual with the smaller fitness. `X.X` is
  the percentage in the interval \[0,100\].

- `hundtoverlap:CONFIG`  
  to use Hundt’s overlap check as a diversity criterion using the
  configuration options in `CONFIG`, semicolon separated.

- `fitnesshundtoverlap:X.X,CONFIG`  
  to use a fitness-prescreened version (with threshold `X.X`) of Hundt’s
  overlap check as a diversity criterion using the configuration options
  in `CONFIG`, semicolon separated.

## Running the job

After calling, e.g.,

    java -jar ogolem.jar --shmem MYJOB.ogo 2

there will, potentially, appear (a lot of files) in the directory the job was started from. Additionally, depending upon the number of threads allowed, a java process will start eating up CPU time. Concerning the latter: it is not recommended to use all CPUs on a workstation which is supposed to be still used for anything else.

<span class="smallcaps">ogolem</span> will clean up all automagically generated input files and/or directories for other programs in case of a sucessfull call to that program. In case that something unexpected (e.g. failing convergence) happens, the files will be left for inspection.

Interesting datafiles are the folder named like your input file, where structures and output (the latter one is to be inspected to see whether all input keywords are correctly read in) are collected. Additionally, there is a log file containing timestamps and informations whenever a new individual took the first rank in the pool.

Once the run is complete, there will be a lot of files of form `rankNUMBERgeometryNUMBER.xyz` appearing in that folder. This is the final pool of structures, where the first number is the rank in the pool (best candidate has rank 0) and the second is the ID of the geometry to give a measure when it was found in the run.

<span class="smallcaps">ogolem</span> might throw warnings or errors. These will be written to System.out and System.err and are of interest.

## Cluster analysis

<span class="smallcaps">ogolem</span> provides some, admittingly very trivial, means to analyze the outcome of a global optimization of clusters and, perhaps more important, to check a running job using the written restart files. For this, you call the `ogolem.jar` with the

    --clusters

option.

The analysis program needs a binary pool, either the final `pool.bin` or a temporary `IntermediateClusterPool.bin`. This is specified using the `-i` flag followed by whitespace, e.g.

    java -jar ogolem.jar --clusters -i IntermediateClusterPool.bin

Some other options are possible

- `-binstructs`  
  bin structures by number of structures instead of energies (which is
  the default)

- `-dissblow`  
  the blow factor for the dissociation detection, defaults to 1.8.

- `-getstructs`  
  create a directory named `structs` and write the structures of this
  pool there in the xyz-file format. Good option for intermediate
  checks.

- `-ljstrains`  
  calculate strain energies for LJ clusters. Defaults to off.

- `-nocomdiffs`  
  disable analysis of COM differences (on by default).

- `-nodissdetect`  
  disable dissociation detection (on by default).

- `-noenergies`  
  disable scanning of energies (on by default).

- `-noinertias`  
  disable moments of inertia analysis (on by default).

- `-noofbins`  
  the number of bins for structure binning. Defaults to 100.

- `-threshdiff`  
  threshold for moments of inertia, defining when they are considered to
  be different. Defaults to 2.5E-1.

- `-threshequal`  
  threshold for moments of inertia, defining when they are considered to
  be same. Defaults to 1E-1.

We want to encourage you here to contribute to the <span class="smallcaps">ogolem</span> by proposing cluster analysis techniques that are of interest to you!

