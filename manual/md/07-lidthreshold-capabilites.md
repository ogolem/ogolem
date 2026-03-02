# Lid/Threshold capabilites

<span class="smallcaps">ogolem</span>features an implementation of Schön and Sibanis lid/threshold algorithm that can be used to study funnels of the energy landscape of clusters and molecules. Please read the original publications for a description of this elegant algorithm.

Our lid implementation is called with the `--lid` flag followed by the path to the .ogo file and the number of threads you intent to use for the lid search. The .ogo file must contain a geometry definition and a choice for the local optimization algorithm (please make sure to use something that also doubles as a backend!) and contains all settings for the lid run in `<LIDCONFIG>` tags.

Relevant input keywords are

- `ThreshMetaRuns=`  
  The number of meta runs for the lid algorithm (read: how often the
  total lid cycle will be run). Default: 2.

- `ThreshUseKearsleyOverlap`  
  Choose Kearsley’s aligner for structural comparison. Works only with
  single molecules. Mutually exclusive with `TheshUseHundtOverlap`.
  Default: off.

- `ThreshUseHundtOverlap`  
  Choose Hundt’s aligner for structural comparison. Mutually exclusive
  with `ThreshUseKearsleyOverlap`. Default: off.

- `ThreshKearslerRMSD=`  
  Maximal RMSD in bohr for structures to be considered equal. Default:
  0.1.

- `ThreshOverlapCOnfig=`  
  Modifies the settings of Hundt’s aligner. See description there.

- `ThreshEnergyPrescreen=`  
  The energy threshold in Hartree for structural comparison or, if a
  secondary aligner is used, for prescreening. Default: 1E-4.

- `NoVerboseSerialization`  
  Disable storage of intermediate results after every meta run. Default:
  on.

- `ThreshMoveclass=`  
  Choosing a move class for the MC steps of the lid run. Options:

  - `mcmove:` MC moves in the Cartesian coordinate space. mcmove must be
    followed by comma separated choices of move mode (0: move one atom,
    1: move random number of atoms, 2: move all atoms) and the maximum
    move in bohr.

  - `mcextmove:` MC moves of fragements in the COM/Euler space.
    mcextmove must be followed by comma separated choices of move mode
    (0: move one fragment, 1: move random number of fragements, 2: move
    all fragements), the maximum move in the COM space in bohr and the
    maximum move in the Euler space as a fraction of the coordinate
    defintion (e.g. 0.2 as 20% of the Euler coordinate).

  - `mcintmove:` MC moves in the internal coordinate space. THIS MUST BE
    USED ONLY FOR SYSTEMS CONTAINING ONE MOLECULE WITH A Z-MATRIX INPUT
    PRESENT AND PROPAGATED IN DIHEDRALS! mcintmove must be followed by
    comma separated choices of the usual mode definition and the maximum
    move in the dihedral space, e.g., 0.2 for 20% of the $`-\pi`$ to
    $`+\pi`$ defined dihedral space. This input options relies on the
    previously mentioned `DOF` definition to allow specific dihedrals
    for propagation.

- `ThresholdProvider=`  
  Defines which form the threshold increase should have. Currently, you
  must use `simple:` followed by a comma separated list of the starting
  energy, the end energy and the linear energy increment in Hartree
  (e.g. `simple:0.02,0.20,0.01` to go from 0.02 to 0.2 Hartree in 0.01
  Hartree increments).

- `ThresholdConfig=`  
  Can be used to modify the internals of the threshold algorithm. Tokens
  are semicolon separated and may contain the following keywords

  - `threshruns=` the number of propagations at a certain lid level.
    Default: 30.

  - `nomoves=` the number of MC moves in one threshold run. Default:
    50000.

  - `freqanalysis` carry out a frequency analysis to check each obtained
    structure if it is a real minimum. NOT WORKING YET!

Additionally, at least one starting structure is required. The defintion is enclosed in `<LIDMINIMA>` tags, containing one structure (represented by the path to an xyz file) per line, any number of lines is supported.

When running the lid algorithm, intermediate data will be printed out (unless disabled) with the `intermediate` suffix and the number of the meta run finished. Ultimately, the final results will be printed. A folder containing all structures found by the algorithm, a dot file representing the graph and the same graph as an <span class="smallcaps">ogolem</span>internal binary for later reuse.

