# Meta-benchmarking <span class="smallcaps">ogolem</span> using Lionbench

## Introduction

Lionbench is a Scala-code developed by Marvin Kammler and Johannes Dieterich at the Georg-August-University Göttingen and has been embedded into the ogolem.jar. It is designed to automatically create <span class="smallcaps">ogolem</span> input, run and analyze the performance based on all combinations of input settings. It is therefore possible to e.g. determine optimal input settings for a series of smaller cluster structure optimizations to use them for larger optimizations and save considerable amounts of computing time. Also, it allows to benchmark new algorithms against a set of known problems and in comparison to other algorithms.

Currently, no set of benchmark problems exists. Any contribution in this respect is very welcome!

## Usage within <span class="smallcaps">ogolem</span>

Access to Lionbench is simply realized by calling the ogolem.jar with the `--lionbench` flag followed by the Lionbench configuration file and the number of CPUs to be used for all <span class="smallcaps">ogolem</span>runs.

Lionbench is calling ogolem through a runOGOLEM file which must be present and executable. It automatically sniffs through all subdirectories of the directory it has been started in and searches for <span class="smallcaps">ogolem</span> input files to run them.

## Input file format: Main Lionbench configuration

In general, the Lionbench configuration file contains keywords followed by whitespace and the option(s).

- `NumberOfSamples`  
  The number of independent runs on each input configuration for
  statistics. Default: 5.

- `ReadLocalOpts`  
  If the number of local optimizations should be read. Default: false.

- `ExitOnFirstNotFound`  
  If Lionbench should exit when a global minimum is not found. Default:
  true.

- `ExitOnProblem`  
  If Lionbench should exit ungracefully upon the first encountered
  problem. Default: true.

- `Blacklist`  
  A comma-separated list of subdirectories that shall be ignored.

- `PartialBlacklist`  
  A comma-separated list of to-be-ignored inputs.

- `Whitelist`  
  A comma-separated list of subdirectories that shall be used (inverts
  the `Blacklist`)

- `VerbosePerFolder`  
  If verbose output for every folder of tests is wanted. Default: true.

- `VerbosePerInput`  
  If verbose output for every test input is wanted. Default: true.

- `DeleteAllFiles`  
  If files should be deleted after successful excecution. Default: true.

- `SaveErrorErrout`  
  If error and output stream of every
  <span class="smallcaps">ogolem</span> run should be save. Default:
  true.

- `AllowedDiffToGlobMin`  
  By how much the solution found by any input configuration is allowed
  to differ from the known global minimum. The known global minimum
  value is set in the .ogo-file with the `AcceptableFitness` option
  (cf. pp.  and ). Default: 1E-6.

## Input file format: Individual <span class="smallcaps">ogolem</span> and `.lion` input files

Lionbench works based on a stub for an <span class="smallcaps">ogolem</span> input (e.g. test.ogo) and a corresponding Lionbench file (test.lion) defining how input parameters for ogolem should be manipulated. For convinience, it is possible to have one .lion file for all <span class="smallcaps">ogolem</span> inputs if they should be altered in the same way. This file must be named `generic.lion` and will be used if present.

Every <span class="smallcaps">ogolem</span> input variable to be varied with Lionbench will be replaced by the `LIONVAR` identifier in the `.ogo` file, e.g., `GlobOptAlgo=LIONVAR`. It is possible to vary multiple variables at once, Lionbench will automatically test all variations of them. The `.lion` file contains information on how those variables should be substituted and which files are needed for this <span class="smallcaps">ogolem</span> run to complete (e.g. coordinate files, auxiliary files,...). For example:

    <LIONOPTS>
    germany AND portugal:1 AND holland
    </LIONOPTS>
    <LIONFILES>
    dummy.xyz
    </LIONFILES>

where the first `LIONVAR` marked input option will be replaced with the identifiers `germany`, `portugal:1` and `holland`. If more than one input option are to be changed in the same Lionbench run, their order in both the `.ogo` and the `.lion` file must match.

In case of numerical variables, another syntax is available. A start to end with increment can be specified as `XXX...YYY...ZZZ` where `XXX` is the numerical start value, `YYY` is the end value and `ZZZ` is the increment.

Every directory with `.ogo`/`.lion` combinations must contain an executable shell script run ogolem named `runOGOLEM`. As an example:

    #!/bin/sh
    $JAVA_HOME/bin/java -jar $OGOLEMDIR/ogolem.jar --adaptive $1 $2 > errout 2>&1

with the usual call of the `ogolem.jar`, the correct switch for the job to run and taking the input parameters \$1 and \$2 for the input file name and the number of threads, respectively.

## Output format

*Please note that the current output format is rather verbose and will be subject to change in later versions of Lionbench!*

Depending on the chosen verbosity, Lionbench will report the average and standard deviation of each input combination to console. Also, individual reports for each run are possible. There is a Python-script written by Marvin Kammler available to parse and beautify the Lionbench output. Eventually, the output will be altered to provide a better out-of-the-box solution.

