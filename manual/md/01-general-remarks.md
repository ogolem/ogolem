# General remarks

## This manual

Is alpha and might not include everything you would like to know. If this is true, contact <span class="smallcaps">ogolem</span> via webpage and mail. If you happen to find any kind of error or uncertainty, please do also contact us.

## Licensing

<span class="smallcaps">ogolem</span> is distributed under the following license

    Copyright (c) 2009-2023, J. M. Dieterich and B. Hartke
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

        * Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.

        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.

        * All advertising materials mentioning features or use of this software
          must display the following acknowledgement:

          This product includes software of the ogolem.org project developed by
          J. M. Dieterich and B. Hartke (Christian-Albrechts-University Kiel,
          Germany) and contributors.

        * Neither the name of the ogolem.org project, the University of Kiel
          nor the names of its contributors may be used to endorse or promote
          products derived from this software without specific prior written
          permission.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

in the hope to help people crack problems and enable interesting science. We do ask you to cite our papers introducing the framework and techniques you use in your publications(Dieterich and Hartke 2010, 2011, 2012, 2014; Carstensen et al. 2011; Dieterich et al. 2011; Buck et al. 2014; Forck et al. 2012; Dieterich, Clever, et al. 2012; Frank et al. 2013; Dieterich, Gerke, et al. 2012). Depending on what is used in or through <span class="smallcaps">ogolem</span> additional citations to papers describing these methods may be required.

We encourage all collaborations, external implementations, usages of <span class="smallcaps">ogolem</span> as a library for any purposes. We are open for discussions regarding uses of <span class="smallcaps">ogolem</span> under `developers@ogolem.org`.

## External libraries and code contributions

<span class="smallcaps">ogolem</span> is grateful for the following codes it uses or took inspiration from (in random oder):

- Atomdroid by J. Feldt, R. A. Mata and J. M. Dieterich (proprietary
  license)

- Lionbench by M. Kammler and J. M. Dieterich (4-clause BSD)

- L-BFGS from RISO by J. Nocedal and R. Dodier (Apache/public domain)

- BOBYQA from the Apache Software Foundation (Apache)

- NEWUOA code by J. Feldt and J. M. Dieterich, original by M. J. D.
  Powell (public domain)

- JAMA

- Apache Commons Math

- Scala

- openJDK

- Numal

- Arpack/Netlib

- SL4J

- JGraphT

- mpiJava

- phenix v0-42

- Tinker

## Financial Acknowledgements

The authors would like to acknowledge the following affiliations:

- J. M. Dieterich

  - *2008-2010:*  
    Institute for Physical Chemistry  
    Christian-Albrechts-University Kiel  
    Germany, Europe

  - *2010-2012:*  
    Institute for Physical Chemistry  
    Georg-August-University Göttingen  
    Germany, Europe

  - *2015-2018*  
    Department of Mechanical&Aerospace Engineering  
    Princeton University  
    USA

- B. Hartke

  - *2008-now:*  
    Institute for Physical Chemistry  
    Christian-Albrechts-University Kiel  
    Germany, Europe

## The <span class="smallcaps">ogolem</span>-jar

<span class="smallcaps">ogolem</span> is distributed as a single jar file. By invoking

    java -jar ogolem.jar

followed by a single mandatory flag deciding which program part to use. These flags are

       --adaptive
         for global parameter optimization
       --beswitched
         for extracting intermediate switch optimization results
       --clusteroverlap
         for checking if two cluster structures are identical
       --core
         for MPI parallelized cluster structure optimization
       --corrfunc
         to compute correlation functions from MD trajectories
       --help
         to display this overview
       --md
         for running simple molecular dynamics
       --molecules
         to optimize molecules with respect to a or multiple properties
       --parameters
         for extracting intermediate parametrization results
       --rmiserver
         for starting an OGOLEM RMI server
       --rmiproxy
         for starting an OGOLEM RMI proxy
       --rmiclient
         for starting an OGOLEM RMI client
       --rmithreader
         for starting an OGOLEM RMI threading client
       --shmem
         for thread-based cluster structure optimization
       --switches
         for thread-based molecular switch optimization
       --tests
         for testing and development purposes
       --version
         for the OGOLEM snapshot version

Please note here that you should have installed a Java runtime environment (JRE) complying at least to the JAVA21 standard. Additionally, for above command to work, the java binary should be in your path, otherwise you will need to explicitly specify the path to it.

If you have a choice, use the newest openJDK (currently openJDK21) in server mode (the latter is really important for general performance and seems to be the default for most JREs for most UNIX incarnations). Also, the author(s) want to note that they could not see any performance gain on standard hardware with proprietary JREs like JRockit or IBM JDK.

## General notes

### Units

All units within <span class="smallcaps">ogolem</span> are (if applicable) atomic units. As a secondary system, we provide SI units as a convinience. Therefore, input and output is (unless explicitly stated otherwise) in atomic units. Obviously, this only applies for the case of data having a physical meaning, otherwise input/output is system dependent and typically unit-free.

## Bugs and problems

Since <span class="smallcaps">ogolem</span> is a relatively new program that is still changing rapidly, we unfortunately can only guarantee that you will probably hit bugs and problems.

In case you hit a bug:

- Relax.

- Mail `developers@ogolem.org` describing the problem and how to
  reproduce it. Please attach input and error messages (in particular:
  stack traces).

- Wait for the bugfix which will be rolled automatically in the next
  autobuild.

