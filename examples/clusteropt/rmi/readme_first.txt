For these RMI-parallelization examples, do consult the manual first,
and optimally also the first RMI paper:
   J.M.Dieterich and B.Hartke, J. Chem. Theory Comput. XX (2016) XXXX,
   DOI: 10.1021/acs.jctc.6b00716
Otherwise, there is little chance of understanding what this is all about.

In the corresponding subdirectories, there are examples for
- simple-client mode
- and threading-client mode.
They both do global optimization of LJ42 clusters, which of course only is a
demonstration example, nothing even remotely close to a case that would need
high-performance computing in a production setting. Therefore, this example
input has been stripped down to the bare essentials.

In both cases, it should suffice
- to start an RMI server job first, via the script "run.server",
- and then to start an RMI client, via the script "run.client(1/2)".
In simple-client mode, both scripts contain only one line. In threading-client
mode, putting all this in scripts makes somewhat more sense, since the client
script first moves itself into a subdirectory. As pointed out in the manual,
this is necessary for the outputs of the server on one hand and of proxies or
threading clients on the other hand to not interfere with each other. In case
of several proxies/threading clients, each should have its own
subdirectory. Of course, setting up and changing into subdirectories is not
necessary for RMI subjobs running in different file systems.

All scripts start all jobs on "localhost", and initiate RMI communication to
target "localhost". Again, this is for demonstration purposes only, and to
make the scripts (possibly) actually executable. In real RMI applications,
most RMI jobs will run on different machines/nodes, and then "localhost"
entries have to be replaced by appropriate, case-specific IP addresses.

Upon normal execution and finishing, expect outputs as shown in the
corresponding *_output subdirectories. Note that
- in simple client mode, the server does almost all of the output, and it
should be similar to that produced by a corresponding --shmem run; the clients
normally only output one single line at the end of execution, reporting their
regular shutdown and their ID number.
- in threading client mode, each threading client produces output of its own,
similar to (but not identical to) that of the server.
