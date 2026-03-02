# Global parametrization of potentials, pseudopotentials, ...

This part of <span class="smallcaps">ogolem</span> started as a way to globally optimize potentials in a system-specific fashion to later use them for cluster structure optimization. It has since significantly evolved from that basis towards optimizing potentials, pseudopotentials, and test benchmark functions. However, some traces of its ancestry are still present, i.e., in the names of keywords.

**IMPORTANT:   This part of the manual is not complete yet and does, e.g., not cover all possible adaptivable choices! It also is outdated as the input syntax has been adjusted.**

The part for the global optimization of potentials in the <span class="smallcaps">ogolem</span> was primarily designed with the idea in mind to allow for a system-specific parametrization of potentials against a reference for later usage in the geometry structure optimization part. Therefore, it is strongly coupled to that part.

Of importance is the definition of the fitness in this part. In the easiest case (no weights, no penalty functions, exact computation), it is ``` math \begin{equation} F=\sum_{i=0}^{N}abs\left( E_{i,fit}-E_{i,ref} \right) \end{equation} ``` where $`N`$ is the number of reference points with their reference energies $`E_{i,ref}`$ and the fitted energies $`E_{i,fit}`$. Therefore, an exact agreement between reference and fitted energies causes the fitness to vanish, making 0.0 the optimal one.

<span id="parameterfit" label="parameterfit"></span>

## General input structure

Due to the tight interaction between this part and the geometry optimization part, all direct input for the parameter optimization resides in the `.ogo` file. Therefore, there *needs* to be a geometry definition. This is useful if one wants to use the feature of automatical generated reference structures and is otherwise only a small nuissance.

All direct input is encapsulated into `<ADAPTIVE>` tags, e.g.,

    ###OGOLEM###
    <GEOMETRY>
    NumberOfParticles=2
    <MOLECULE>
    Ar;0.0;0.0;0.0
    </MOLECULE>
    <MOLECULE>
    Ar;0.0;0.0;0.0
    </MOLECULE>
    </GEOMETRY>

    <ADAPTIVE>
    AcceptableFitness=1E-8
    <REFERENCE>
    <PATH>
    arar_2.xyz
    </PATH>
    <ENERGY>
    1E-5
    </ENERGY>
    </REFERENCE>
    <REFERENCE>
    <PATH>
    arar_22.xyz
    </PATH>
    <ENERGY>
    1.2E-5
    </ENERGY>
    </REFERENCE>
    </ADAPTIVE>

If one does not use the by default disabled feature of automatic reference point generation, the geometry definition does not need to be identical to the reference points. Also, in general the reference points may (depending on the adaptivable chosen) have different compositions and number of atoms.

Additionally, there exists an optional file specific to the parametrization part. This file defines lower and upper boundaries to each parameter and will be used for this purpose if present, otherwise default borders will be used which might be very inefficient for your specific purpose.

The file should be named `yourinputname-borders.aux`, where `yourinputname` is the prefix of the used `.ogo` file, and the general format is

    ###OGOLEMAUX###
    -3.0        2.0
    0.1         0.2
    1E3         1E5
    -1E3        -1E1

where for each individual parameter a lower boundary is assigned on the left and a upper one on the right side. These bounds typically (exceptions are mentioned in the corresponding paragraph of the adaptivable choice) are in atomic units.

It is important to note that to-date the length of the specified array of borders needs to agree with the overall number of parameters for the fit and the order needs to agree with the internal representation. Depending on the actual interface chosen, it may therefore by reasonable to start <span class="smallcaps">ogolem</span> with default borders to figure these things out.

The automatic calculation, how many parameters and which parameters are needed for a given set of reference points may be too conservative and generate too many parameters. This of course highers the problem dimensionality and makes the genetic operations and local optimizations be less efficient. If this is the case for you, you may provide a file named `yourinputname-stub.aux`, where `yourinputname` is the prefix of the used `.ogo` file. This file will then be used as a reference stub for the global optimization. Please be aware that any errors here can cause really nasty behaviour of the optimization. In general, the format of this file should follow the later described parameter output format. The numerical values of the parameters given can be any floating point number and do not have any effect on the optimization. Again, this potentially very harmful option should only be used if one knows what one is doing!

## Input of reference points

The reference points are handed over through input of form

    <REFERENCE>
    <PATH>
    pathToCartesianReference.xyz
    </PATH>
    <ENERGY>
    0.1
    </ENERGY>
    <REFCHARGES>
    0;1;2
    </REFCHARGES>
    <REFSPINS>
    0;2;1
    </REFSPINS>
    <REFBONDS>
    bonds.list
    </REFBONDS>
    <REFWEIGHT>
    2.4
    </REFWEIGHT>
    <REFMAXDIFF>
    1E-5
    </REFMAXDIFF>
    </REFERENCE>

All adaptive backends support energy fits currently (this may change in the future) and for these `<PATH>` and `<ENERGY>` tags are mandatory and hold the path to the Cartesian structure and its energy (in hartree) of this reference point. The other tags are optional with both `<REFSPINS>` and `<REFCHARGES>` following the exact same input syntax as described above for the geometry structure optimization.

Please note that the `<PATH>` tags must always be located before any other tags!

The `<REFBONDS>` tags may contain the path to either a list or a full matrix of the bonds in this reference geometry. If the file suffix is `.list`, a list will be assumed, otherwise a full matrix needs to be in the file (in most cases, one wants to use the list option). If these tags are not given, <span class="smallcaps">ogolem</span> makes estimates based on a distance criterion. The list needs to be in format

    0 1
    6 5
    2 3
    4 5

where each row specifies the two atoms (once more counting starts at 0) between a bond exists. Empty lines are permitted, no order needs to be pertained. If a full matrix needs to be read in, the file should contain

    true true true
    true true false
    true false true

for (in this case) a water molecule. No empty lines are permitted.

The `<REFWEIGHT>` tag assigns a weight to this specific point, a feature that might be usable if e.g. the description of a certain region is important. The `<REFMAXDIFF>` tag defines how much the fitted energy is allowed to differ from the reference before a penalty is applied. This last tag obviously must be used with a penalty function enabled separatly.

There is growing support for other properties. Please check if your problem supports the following. If not, patches are welcome!

### Parametrizing energy order

Energy order refers to the ordering of different atomic arrangements (e.g., crystal structures). It can be input as follows:

    <REFERENCE>
    <ENERGYORDER>
    ConsiderRelativeEnergies
    First;pathtofirst.xyz;-945.64
    Second;pathtosecond.xyz;-1213.3
    </ENERGYORDER>
    </REFERENCE>

where `ConsiderRelativeEnergies` is an optional tag (must be in first line if used) to not use total energies for the parametrization but only relative ones. Followed by arrangements having a tag as the first element, a semicolon as a separator, the path (relative or absolute) to the atomic positions in either xyz or z-matrix format, another semicolon separator, and the energy (typically in Hartree, can be relative if above relative energies were requested).

### Parametrizing bulk modulus

Bulk moduli are of importance for condensed matter applications. They can be input as follows:

    <REFERENCE>
    <BULKMODULUS>
    Li(bcc):150.624
    </BULKMODULUS>
    </REFERENCE>

where `Li(bcc)` refers to the bcc phase of lithium, followed by colon and the bulk modulus (typically in atomic units).

### Parametrizing equilibrium volume

Equilibrium cell volumes are of importance for condensed matter applications. They can be input as follows:

    <REFERENCE>
    <CELLVOLUME>
    Li(bcc):42.4242
    </CELLVOLUME>
    </REFERENCE>

where `Li(bcc)` refers to the bcc phase of lithium, followed by colon and the equilibrium cell volume (typically in atomic units).

### Parametrizing nuclear forces

Nuclear forces can be input as follows:

    <REFERENCE>
    <PATH>
    ...
    </PATH>
    <FORCES>
    ForcesFile=pathtoforces.dat
    </FORCES>
    </REFERENCE>

where `ForcesFile=` points to a file (relative or absolute path) containing the forces for the current reference in a Cartesian xyz format. The forces input requires also a `<PATH>` environment with atomic positions as defined above.

### Parametrizing $`\delta`$ gauge

$`\delta`$ gauges are of importance for condensed matter applications. They can be input as follows:

    <REFERENCE>
    <PATH>
    ...
    </PATH>
    <DELTAGAUGE>
    Li(bcc):0.0
    </DELTAGAUGE>
    </REFERENCE>

where `Li(bcc)` refers to the bcc phase of lithium, followed by colon and the $`\delta`$ gauge (typically 0.0 in atomic units). Currently, the $`\delta`$ gauge requires also a `<PATH>` environment with atomic positions as defined above.

### Parametrizing electron density

Electron density (on a grid) can be input as follows:

    <REFERENCE>
    <PATH>
    ...
    </PATH>
    <DENSITY>
    DensityFile=pathtodensity.dat
    ReferenceTag=myTag
    </DENSITY>
    </REFERENCE>

where `DensityFile=` points to a file (relative or absolute path) containing the density on a grid for the current reference. First line contains the Cartesian grid dimensions (whitespace separated) and subsequent lines contain one data point on that grid. The electron density input requires also a `<PATH>` environment with atomic positions as defined above. The second line of the `<REFERENCE>` environment must define a tag for this density using `ReferenceTag=`.

### Mixing properties for the same reference

Some of the above properties can be mixed if the same reference is meant (the definition of that is somewhat backend dependent!). For this, simply add multiple property tags (e.g., `<CELLVOLUME>` `BULKMODULUS` in the same `<REFERENCE>` environment.

## Input keywords

Again, there are a lot of tunables to this part of <span class="smallcaps">ogolem</span> with some of them more important than others. It is important to note that a certain setting might be optimal for one fit but not for another.

- `AcceptableFitness=`  
  <span id="page2:acceptableFitness"
  label="page2:acceptableFitness"></span> which fitness is considered to
  be acceptable and causes the optimization to stop. Defaults to 0.0
  which is by our definition total agreement between all references and
  the fit.

- `AdaptivableChoice=`  
  what to adapt/parametrize. See details below.

- `AllRefsSameChargesAndSpins=`  
  set to true if all reference geometries have the same charge and spin.
  Use with extraordinary care only!

- `AnyParamHistory=`  
  set to true if the genetic history of parameters should be tracked.
  This implies a higher use of memory.

- `DoNiching=`  

- `DoubleMorseDMax=`  
  if the double morse adaptivable is used, the $`d_{max}`$ value.
  Deprecated! set to true if niching should be enabled.

- `NichesPerDim=`  
  if niching is enabled, the number of niches per (static) grid
  dimension.

- `MaxIndividualsPerNiche=`  
  if niching is enabled, the maximum number of individuals per niche.

- `MaxNonConvSteps=`  
  maximum number of steps w/o fitness progress in the internal local
  optimization. Defaults to 100.

- `MaxTasksToSubmit=`  
  maximum number of tasks that are submitted to the threadpool at once.
  Defaults to 1000000. Please keep this setting if you are not
  absolutely certain that you benefit from a change.

- `ParamBorderPrint=`  
  set to false if the parameter borders should be printed.

- `ParameterDetailedStats=`  
  set to true if detailed statistics for the global optimization are
  wanted.

- `ParamEnergyDiv=`  
  the minimal fitness diversity needing to exist between two individuals
  in the pool. Defaults to 1E-8 which might be (depending on e.g. the
  convergence criterion of the local optimization) too small and cause
  premature convergence of the pool.

- `ParamGlobOptIter=`  
  the number of global optimization steps. Defaults to 9900. Must be
  increased for production runs!

- `ParamLinesearchEnergyDecr=`  
  threshold for the fitness decrease in the linesearch. Defaults to
  1E-10.

- `ParameterGlobOpt=`  
  which global optimization ingredients to use. Syntax is
  `parameters{xover()mutation()}`. Choices include: `portugal:` with
  optional `nocuts=N` argument for an N-point genotype xover. and
  `germany:` for both mutation and xover. Details below.

- `ParamLocOptIter=`  
  maximum number of iterations for the parameter local optimization.
  Defaults to 5000.

- `ParameterLocOpt=`  
  which local optimization engine should be used for the local
  optimization. Choices include: apachecg (for Apache’s conjugate
  gradient implementation), lbfgs (for an L-BFGS), and bobyqa (for
  Powell’s BOBYQA).

- `ParamLocOptMaxStep=`  
  maximum stepsize in the local optimization. Defaults to 1E-5.

- `ParamSeedingFolder=`  
  specify the folder where seeds for a seeding initialization are
  located.

- `ParamSerializeAfterNewBest=`  
  set to false if the genetic pool should not be stored to disk after
  each new best individual found.

- `ParamThreshDiv=`  
  the minimal difference between the same gene of two individuals.
  Applies only if the corresponding diversity check is enabled. Defaults
  to 1E-8 which might be (depending on e.g. the convergence criterion of
  the local optimization) too small and cause premature convergence of
  the pool.

- `ParamThreshLocOptGradient=`  
  threshold for the gradient convergence. Defaults to 1E-8.

- `ParamThreshLocOptParam=`  
  threshold for the parameter convergence. Defaults to 1E-8.

- `ParamsToBorders=`  
  whether the parameters should be put back into their borders.
  Otherwise, a local optimization may cause the parameters to run out of
  the borders. Defaults to false.

- `ParamsToSerial=`  
  After how may steps a binary intermediate parameter pool is written to
  disk. Defaults to 10000;

- `PercentageForDiffs=`  
  if percentage-scaled differences are enabled, what percentage should
  be used. Default: 1 \[in percent\].

- `PopulationSize=`  
  the pool size for the parameter fit. Defaults to 100.

- `UsePercentageScaledDiffs=`  
  set to true to enable percentage-scaled differences for the maximum
  allowed difference.

- `WhichNicher=`  
  selects the niching protocol being used. First an integer selection
  for the niching algorithm, followed by a semicolon and the nicher.
  Choices are 0 for static grid niches without ratio consideration and 1
  for static grid niches with ratio consideration.

- `WhichDiversityCheck=`  
  integer choice which diversity check is to be used. Defaults to 0,
  which is a fitness based diversity check. Also possible: 2, which
  compares the individuals based on their gene values.

## Adaptivable choices

<span class="smallcaps">ogolem</span> supports a to-date only a few adaptivable choices. If there should be a (preferably free) program that you want to see supported and/or the methods listed here to do not satisfy you, feel free to write a mail to dieterich@ogolem.org . It is necessary to provide reference in- and output to the program you want to see support for and (in case that it is no freeware) at least a binary of the program/engine is required.

In general, all adaptivable choices are easily called through input of form

    AdaptivableChoice=IDENTIFIER:FLAGS

### Adaptive GUPTA potential

The adaptive GUPTA potential can be called using `adaptivegupta:DISTBLOW,ENABLECACHE`, where `DISTBLOW` is a blow factor defining when to neglect contributions. `ENABLECACHE` is a boolean defining whether or not to enable the internal caching procedures. The later should only be used when all geometries are of precisely the same system. In general, it can and should be anabled for cluster global optimization and may be enable for parameter global optimization if (and *only* if) all reference structures are of the same system. Using caching provides a significant speedup in comparision to the uncached procedures.

The used function definition is as follows: ``` math \begin{equation} E(a,b)=A(a,b)\cdot \exp{\left[-p(a,b)\left(\dfrac{r_{ab}}{r_0(a,b)}-1\right)\right]} -\sqrt{\chi(a,b)^2\cdot\exp{\left[-2\cdot q(a,b) \left(\dfrac{r_{ab}}{r_0(a,b)}-1\right)\right]}} \end{equation} ``` where $`r_{ab}`$ is the distance between atoms $`a`$ and $`b`$.

The implementation loops over all atoms in the system, independent whether they are bonded or not. Therefore, this should primarily be used for e.g. metal clusters. The output parameters are in the order $`A(a,b)`$, $`p(a,b)`$, $`r_0(a,b)`$, $`\chi(a,b)`$ and $`q(a,b)`$ for each possible pair in atomic units. The implementation does feature an analytical gradient for both coordinates and parameters for enhanced performance of the local optimizations.

### Adaptive LJ potential

### Adaptive Morse potential

The adaptive Morse potential can be called using `adaptivemorse:CLOSEBLOW,DISTBLOW`, where `CLOSEBLOW` is a blow factor defining when to use a cutoff for close cases and `DISTBLOW` defines when to neglect contributions.

The implementation loops over all atoms in the system, independent whether they are bonded or not. Therefore, this should primarily be used for e.g. metal clusters. The output parameters are in the order $`D_e`$, $`a`$ and $`r_e`$ for each possible pair in atomic units. The implementation does feature an analytical gradient for both coordinates and parameters (suboptimal implementation) for enhanced performance of the local optimizations.

Please note that currently no caching is implemented in the Morse potential. If you want to use this potential for bigger systems, please contact the author(s) whether this can be implemented.

### Benchmark: Ackley

### Benchmark: Schaffer F6

### Benchmark: Schaffer F7

### Benchmark: Lunacek

### goLPS parametrization using PROFESS

<span class="smallcaps">ogolem</span> can be used to fit local pseudopotentials following the goLPS formalism. Our current implementation requires PROFESS (open-source) for the actual OFDFT calculations with the modified pseudopotential and a set of helper scripts designed to translate from <span class="smallcaps">ogolem</span> output to PROFESS input and from PROFESS output back to <span class="smallcaps">ogolem</span> input. We will document here only the data exchange protocol, not the scripts themselves. The current recommendation is to write them based on the problem to be fitted yourself. If assistance with this step is required, please contact the author(s) of the original goLPS publication.

In the following, the properties to-be-fitted are individually discussed. Please note that this feature is experimental and we do NOT recommend fitting goLPSss for multiple atoms simulaneously.

#### Energy

If a reference contains the energy (i.e., total energy), a directory following the scheme `professenergycalc_RUNNINGNUMNER_TIMESTAMP` will be created. In it, the following files will be created:

- a file `pseudopot_ATOMNAME.recpot` for each ATOMNAME to be optimized,
  containing the current goLPS trial in PROFESS format

- a file `atompos.xyz` containing the atom positions of the reference in
  xyz format.

The following files will be copied from the directory in which <span class="smallcaps">ogolem</span> is executed

- `profess_energy_calc.sh` a script to be called (name can be overridden
  by setting the environment variable `OGO_PROFESSENERGYSCRIPT` to the
  wanted name).

Subsequently, the script will be executed, using `XXX` as the first argument. When execution of this script finishes, <span class="smallcaps">ogolem</span> expects the calculated energy (in the same units as the reference) as a simple number in a file `energy.dat`).

#### Energy order

If a reference contains the energy order of different atomic arrangements (e.g., crystal structures), for each crystal structure a directory following the scheme `professenergycalc_RUNNINGNUMNER_TIMESTAMP` will be created. In it, the following files will be created:

- a file `pseudopot_ATOMNAME.recpot` for each ATOMNAME to be optimized,
  containing the current goLPS trial in PROFESS format

- a file `atompos.xyz` containing the atom positions of the reference in
  xyz format.

The following files will be copied from the directory in which <span class="smallcaps">ogolem</span> is executed

- `profess_energy_calc.sh` a script to be called (name can be overridden
  by setting the environment variable `OGO_PROFESSENERGYSCRIPT` to the
  wanted name).

Subsequently, the script will be executed for each structure, using the defined key for this arrangement as the first argument. When execution of this script finishes, <span class="smallcaps">ogolem</span> expects the calculated energy (in the same units as the reference) as a simple number in a file `energy.dat`). It will internally rank the energies.

#### Bulk modulus

If a reference contains the bulk modulus, a directory following the scheme `professbulkcalc_RUNNINGNUMNER_TIMESTAMP` will be created. In it, the following files will be created:

- a file `pseudopot_ATOMNAME.recpot` for each ATOMNAME to be optimized,
  containing the current goLPS trial in PROFESS format.

The following files will be copied from the directory in which <span class="smallcaps">ogolem</span> is executed

- `profess_bulk_calc.sh` a script to be called (name can be overridden
  by setting the environment variable `OGO_PROFESSBULKSCRIPT` to the
  wanted name).

Subsequently, the script will be executed, using the crystal symmetry (e.g., fcc) of the reference as the first argument. When execution of this script finishes, <span class="smallcaps">ogolem</span> expects the calculated bulk modulus (in the same units as the reference) as a simple number in a file `bulk.modulus`).

#### Equilibrium volume

If a reference contains the equilibrium volume, a directory following the scheme `professcellvolcalc_RUNNINGNUMNER_TIMESTAMP` will be created. In it, the following files will be created:

- a file `pseudopot_ATOMNAME.recpot` for each ATOMNAME to be optimized,
  containing the current goLPS trial in PROFESS format.

The following files will be copied from the directory in which <span class="smallcaps">ogolem</span> is executed

- `profess_cellvol_calc.sh` a script to be called (name can be
  overridden by setting the environment variable
  `OGO_PROFESSCELLVOLSCRIPT` to the wanted name).

Subsequently, the script will be executed, using the crystal symmetry (e.g., fcc) of the reference as the first argument. When execution of this script finishes, <span class="smallcaps">ogolem</span> expects the calculated equilibrium volume (in the same units as the reference) as a simple number in a file `cell.volume`).

#### Combined energy, bulk modulus, and equilibrium volume calculations

If the same reference contains total energy, bulk modulus, and equilibrium volume, they will be evaluated simulaneously. A directory following the scheme `professcombcalc_RUNNINGNUMNER_TIMESTAMP` will be created. In it, the following files will be created:

- a file `pseudopot_ATOMNAME.recpot` for each ATOMNAME to be optimized,
  containing the current goLPS trial in PROFESS format

- a file `atompos.xyz` containing the atom positions of the reference in
  xyz format.

The following files will be copied from the directory in which <span class="smallcaps">ogolem</span> is executed

- `profess_combined_calc.sh` a script to be called (name cannot be
  overridden).

Subsequently, the script will be executed, using the crystal symmetry (e.g., fcc) of the reference as the first argument. When execution of this script finishes, <span class="smallcaps">ogolem</span> expects the calculated equilibrium volume (in the same units as the reference) as a simple number in a file `cell.volume`), the calculated bulk modulus (in the same units as the reference) as a simple number in a file `bulk.modulus`), and the calculated energy (in the same units as the reference) as a simple number in a file `energy.dat`).

#### $`\delta`$ gauge

If a reference contains the $`\delta`$ gauge, a directory following the scheme `professdeltagaugecalc_RUNNINGNUMNER_TIMESTAMP` will be created. In it, the following files will be created:

- a file `pseudopot_ATOMNAME.recpot` for each ATOMNAME to be optimized,
  containing the current goLPS trial in PROFESS format.

The following files will be copied from the directory in which <span class="smallcaps">ogolem</span> is executed

- `profess_deltagauge_calc.sh` a script to be called (name can be
  overridden by setting the environment variable
  `OGO_PROFESSDELTAGAUGESCRIPT` to the wanted name).

Subsequently, the script will be executed, using the crystal symmetry (e.g., fcc) of the reference as the first argument. When execution of this script finishes, <span class="smallcaps">ogolem</span> expects the calculated $`\delta`$ gauge (in the same units as the reference) as a simple number in a file `Delta.out`).

#### Forces

If a reference contains nuclear forces, a directory following the scheme `professforcescalc_RUNNINGNUMNER_TIMESTAMP` will be created. In it, the following files will be created:

- a file `pseudopot_ATOMNAME.recpot` for each ATOMNAME to be optimized,
  containing the current goLPS trial in PROFESS format

- a file `atompos.xyz` containing the atom positions of the reference in
  xyz format.

The following files will be copied from the directory in which <span class="smallcaps">ogolem</span> is executed

- `profess_forces_calc.sh` a script to be called (name can be overridden
  by setting the environment variable `OGO_PROFESSFORCESSCRIPT` to the
  wanted name).

Subsequently, the script will be executed, using `XXX` as the first argument. When execution of this script finishes, <span class="smallcaps">ogolem</span> expects the calculated nuclear forces (in the same units as the reference) as a simple number in a file `forces.out`). The format is as follows: first five lines are ignored, sixth line onwards contains the forces (same order as reference) in columns two (x-component) to four (z-component).

#### Electron density

If a reference contains the electron density, a directory following the scheme `professdensitycalc_RUNNINGNUMNER_TIMESTAMP` will be created. In it, the following files will be created:

- a file `pseudopot_ATOMNAME.recpot` for each ATOMNAME to be optimized,
  containing the current goLPS trial in PROFESS format

- a file `atompos.xyz` containing the atom positions of the reference in
  xyz format.

The following files will be copied from the directory in which <span class="smallcaps">ogolem</span> is executed

- `profess_density_calc.sh` a script to be called (name can be
  overridden by setting the environment variable
  `OGO_PROFESSDENSITYSCRIPT` to the wanted name).

Subsequently, the script will be executed, using the defined tag for this reference as the first argument. When execution of this script finishes, <span class="smallcaps">ogolem</span> expects the calculated energy (in the same units as the reference) as a simple number in a file `density.out`). The format is as follows (please note, this is PROFESS’s fault!): a single (!), white space separated line, containing the grid dimensions as values two, four, six, and density values (same order as reference) as values 11 and later. Negative density values will be ignored.

### Interface to Polaris

A word of prelimenary caution is needed here. This interface was written for a single and special purpose. Therefore, it will not work for general-purpose tasks by design.

The interface to Polaris is called through `AdaptivableChoice=polaris`, followed by the choices `:totalenergy` or `:waterinteraction` defining which energy holds as reference. The references need to be adjusted accordingly and must be in $`E_h`$.

Additionally, some auxiliary input and file stubs are required. Both a `MET.FFM` and `REP.FFM` file need to be present in the input directory, where all parameters that are to be fitted are replaced by the string `OGOPARAM`. A reference `donee` file needs to be present, starting from

    -----------------------
    |   Options energie   |
    -----------------------

and named `donee-ref.aux`. If there are other files to be copied for the Polaris jobs, they need to be present in the input directory and their name needs to be mentioned in a file `adaptive-polaris.aux` of form

    ###OGOLEMAUX###
    filenametocopy1
    filenametocopy2

Since Polaris requires more than just a simple xyz-file as input, we assume that for each reference point the following files exist:

- `GEOM.xyz_NUMBER` and

- `GEOM.sle_NUMBER`

where number is the the position in the references starting from zero.

It is highly recommended to use a custom border file as explained above where the parameters start with the `MET.FFM` parameters followed by the `REP.FFM` parameters. Due to the parameters not being known by <span class="smallcaps">ogolem</span> in advance, the borders and parameter values in <span class="smallcaps">ogolem</span> are in the same units as in Polaris.

Last but not least, the path to the Polaris command needs to be exported through an environment variable named `OGO_POLARISCMD`.

Again, this interface was written for a single usage scenario and is not designed to be used for other purposes.

## Global optimization algorithms

All global optimization algorithms are chosen through the `ParameterGlobOpt=` keyword. The input syntax follows the grammar as described in Section <a href="#sec:geomglobopt" data-reference-type="ref" data-reference="sec:geomglobopt">2.7</a>, it is `parameters``xover()mutation()`.

Possible crossover choices are:

- `germany:` a basic single cut, Gaussian distributed genotype
  crossover. Valid options:

  - `gausswidth=X.X` Sets the width of the Gaussian to `X.X`. Default:
    0.3.

- `portugal:` a genotype, normal-distributed crossover.

  - `nocuts=X` sets the number of genotype cuts to `X`. Default: 1.

- `noxover:` disables crossover operations.

- `multiple:` uses multiple crossover operators. Syntax identical to the
  one described in Section
  <a href="#sec:geomglobopt" data-reference-type="ref"
  data-reference="sec:geomglobopt">2.7</a>.

- `chained:` chains multiple crossover operators. Syntax identical to
  the one described in Section
  <a href="#sec:geomglobopt" data-reference-type="ref"
  data-reference="sec:geomglobopt">2.7</a>.

- `mutationasxover:` uses a mutation operator instead of a crossover
  operator for the xover step. Keyword must be followed by a valid
  mutation as described below.

After some prelimenary testing we are of the opinion that the exact crossover (when using local optimization steps) is not important. Therefore we recommend to use either `germany` or `portugal` with a low number of cuts.

In some cases (e.g. for a hypersurface with a lot of minima) it seems to be a good idea to aggressively mutate the parameters, making `france` or `spain` with a low number of cuts the methods of choice.

Possible mutation choices are:

- `germany:` a basic genotype mutation. Valid options:

  - `mode=` sets the strength of the mutation. Valid input: `all` or
    `one`. Default: `one`.

- `nomutation:` disables mutation operations.

- `multiple:` uses multiple mutation operators. Syntax identical to the
  one described in Section
  <a href="#sec:geomglobopt" data-reference-type="ref"
  data-reference="sec:geomglobopt">2.7</a>

- `chained:` chains multiple crossover operators. Syntax identical to
  the one described in Section
  <a href="#sec:geomglobopt" data-reference-type="ref"
  data-reference="sec:geomglobopt">2.7</a>.

### HACTAR

Please see Section <a href="#sec:geomhactar" data-reference-type="ref" data-reference="sec:geomhactar">2.7.1</a> for a description of the input format.

## Local optimization algorithms

All local optimization algorithms are chosen through the `ParameterLocOpt=` keyword. In general, the same choices are available as mentioned in Sec. <a href="#sec:geomlocopt" data-reference-type="ref" data-reference="sec:geomlocopt">2.5.18</a>. Please note that the same restrictions concerning the methods obtained through the NUMAL interface apply as mentioned there.

In difference to the cluster structure optimization, the `backend=` definition present in the local optimization definitions has a somewhat restricted feature set. Either, one can specify `alldefault` to default to the classic setup of the fitness function. Note however, this option may well chance in the future! More future-proof is a pre-definition of fitness functions as defined in Sec. <a href="#sec:paramfitfuncdef" data-reference-type="ref" data-reference="sec:paramfitfuncdef">3.7</a>.

## Defining fitness functions

Fitness functions can be defined in between `<PARAMFITNESSFUNCTION>` tags. Per fitness function definition, one set of tags is required. Configuration options inside the tags are:

- `FitFunctionTag=XYZ` assigns tag `XYZ` to this fitness function.
  Mandatory.

- `FitnessFunction=` specifies the fitness function to be used. At the
  moment, only `default` is allowed. Mandatory.

- `FitnessFunctionStyle=` specifies the fitness function style to be
  used. `full` choses all parameters to be exposed always. `ranged:`
  followed by range definitions (example: 1/3-4/7-8) specifies only
  these parts of the parameters to be exposed. Optional, default:
  `full`.

## Parameter output format

Both the inital and final population are written in a human-readable format. This format, as most input and output formats of <span class="smallcaps">ogolem</span> uses an XML-like formatting. The general format is as follows

    <PARAMETERSFOR>
    methodID
    </PARAMETERSFOR>
    <PARAMETERFITNESS>
    1.23E-1
    </PARAMETERFITNESS>
    <PARAMETERS>
    <ATOM>
    Key
    </ATOM>
    <VALUES>
     -0.4636119053304735
     -0.5097016790211515
     -1.5220016461819896
     -0.5629664946696742
     -3.541843106605196
     -3.5394934242499287
     0.5378639455382299
     -3.5096479279350183
     0.2614123652818187
     -2.399411743051027
     -3.51017085778294
     1.4941457394022188
     -2.539798139301017
     -4.5084460220195215
     -0.5152036256889749
     -3.511654320068149
     -4.456450012737577
     -3.52578457936217
     -2.539843561431608
     1.477287999076277
    </VALUES>
    </PARAMETERS>

where the method ID resembles the fitted potential. The parameter fitness is in hartree and a *reasonable* fitness depends on the settings, the dimensionality and form of the problem. E.g. the fitness of a high dimensional problem with a lot of reference points employing a steep penalty function accepting only a small deviation from the reference is surely higher than, e.g., a two-body LJ(6,12,6) fit.

In general, as pointed out before, the best fitness is *per definitionem* 0.0. This is, except for benchmark functions, in general not reachable.

This is followed by (depending on the potential and reference data) a number of `<PARAMETER>` tags each uniquely identified by the key in between the `<ATOM>` tags. The tag is somewhat misleading since depending on the potential, it might not be a atom but, e.g., a three-body interaction. The actual key transparently informs about such cases. Each `<PARAMETER>` section contains in between the `<VALUES>` tags all raw values for the specific key. Which parameter means what depends on the potential and is explained in detail in the corresponding section of this manual.

Although we consider this format to be rather final, we do not guarantee that it stays unchanged in later revisions of <span class="smallcaps">ogolem</span>.

## Running the job

After calling, e.g.,

    java -jar ogolem.jar --adaptive input.ogo 2

Depending on the adaptivable backend chosen, there might be a lot of files appearing in the directory and the specified amount of CPUs be used. Concerning the latter: it is not recommended to use all CPUs on a workstation which is supposed to be still used for something else.

<span class="smallcaps">ogolem</span> will clean up all automagically generated input files and/or directories for other programs in case of a sucessfull call to that program. In case that something unexpected (e.g. failing convergence) happens, the files will be left for inspection.

Interesting datafiles are the intermediate binary pool (`IntermediateParamPool.bin`) as well as the human-readable parameter initial files `locrankRANKparamID.prm` and the final files `rankRANKparamID.prm`, where `RANK` is the rank of the parameter set (0 being the global minimum) and `ID` denotes when the specified set was found.

<span class="smallcaps">ogolem</span> might throw warnings or errors. These will be written to System.out and System.err and are of interest.

## Parameter analysis

<span class="smallcaps">ogolem</span> provides some possibilities for parameter analysis. The most useful application of them to-date is probably the intermediate check of parameters and genetic progress at runtime. For this, you call the `ogolem.jar` with the

    --parameters

option.

The analysis program needs a binary pool, either the final `paramPool.bin` or the temporary `IntermediateParamPool.bin`. This is specified using the `-pool=` flag directly followed by the filename, e.g.

    java -jar ogolem.jar --parameters -pool=IntermediateParamPool.bin

Another option is possible

- `-dontgetparams`  
  Do not write the parameters out. Currently a rather useless option.

Once more, we want to encourage you here to contribute to the <span class="smallcaps">ogolem</span> by proposing some analysis techniques!

