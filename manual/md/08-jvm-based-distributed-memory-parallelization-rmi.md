# JVM-based distributed-memory parallelization (RMI)

The cluster, parameter and switch optimization parts of <span class="smallcaps">ogolem</span> can also use a JVM-based parallelization via RMI (remote method invocation). While possessing some unique features of its own, its targeted use is replacing the MPI-based parallelization layer. See Ref.  for a full discussion of the background rationale for this RMI parallelization mode, of its different modi, and of its performance, including in particular demonstrations of its excellent parallel scaling. The most important benefits include the ability to dynamically add or subtract clients at will during runtime, also resulting in almost total immunity against hardware failures, and parallelization across arbitrarily heterogenous ad-hoc networks, only needing a standard JVM on each node and standard TCP/IP connections between them.

In order for this parallelization scheme to work, a standard Java distribution must be used, there need be two open TCP ports dedicated on the server and (if applicable) proxy nodes. While not recommended, one can run a server and (multiple) proxies on the same machine if different TCP ports are used.

Also, some installation-dependent flags to the JVM must be passed for security policies, CLASSPATH settings etc. These are included in the startup syntax discussions below but are not yet explained further in this version of this manual.

There are three possible modes:

simple clients:   1 server communicating directly with $`N`$ clients,

proxies:   1 server communicating with $`M`$ proxies, which in turn communicate with $`N_M`$ simple clients each,

threading clients:   1 server communicating with $`N`$ tread-parallel clients.

Essential similarities and differences between these modes (fully explained in Ref. ) include the following:

- The clients in the “simple clients” and “proxies” modes are simple in
  several ways:

  - tasks they receive for execution always consist of just one single
    global optimization step, i.e., crossover and mutation of one or two
    “parents”, and local optimization of the resulting one or two
    “children”;

  - these tasks are executed serially (no parallelization);

  - from the user, these clients do not receive any problem-specific
    input (which is provided at runtime by the server).

  In contrast,

  - both proxies and threading clients receive chunks of tasks,
    consisting of 1,2,…,$`m`$ global optimization steps;

  - threading clients execute their task chunks in thread-parallel
    manner;

  - and both proxies and threading clients need a full `ogo`-file as
    input.

- From the viewpoint of the RMI server, there is no difference in
  communication to proxies or to threaded clients. From the viewpoint of
  the simple clients, there is no difference in communication to proxies
  or directly to the server.

Threading clients can be understood as combining the roles of one proxy and all the simple clients connected to this proxy into one single job that is then executed in shared-memory/thread-based parallelization. Alternatively, one proxy and all simple clients connected to it can be understood as a threading client executed in RMI parallelization.

## Starting RMI job-hierarchies

As things are set up inside <span class="smallcaps">ogolem</span>, all three of the above modes of RMI-mode hierarchies (simple clients, proxies, threading clients) consist of several (at least two) separate <span class="smallcaps">ogolem</span> jobs that need to be started separately, and jobs “higher” in the hierarchy (servers, proxies) always need to be started before jobs “lower” in the hierarchy (proxies, clients).

### Starting the RMI server

Hence, in every mode, first a server process must be started. This can be accomplished using the following (in one single, long line):

    java -Djava.security.main -Djava.rmi.server.codebase=file:<path/to/ogolem-jar-file> 
    -Djava.security.policy=polfile.txt -Djava.rmi.server.hostname=<IP-address> 
    -jar <path/to/ogolem-jar-file> 
    --rmiserver 
    --jobtype <job type> 
    --inputfile <input.ogo> 
    --timeout <timeout> 
    --timeoutjob <timeoutjob> 
    --myserverport <server-port#> 
    --myregistryport <registry-port#> 
    --noproxies <#proxies>
    --commname <nameofregistryobject>

where

- `-Djava.security.main` enables Java’s security manager (mandatory)

- `-Djava.rmi.server.codebase=file:` tells the security manager where
  the ogolem.jar lies (mandatory)

- `-Djava.security.policy=` specifies the location of the security
  policy file (mandatory)

- `polfile.txt` is the filename of a java security policy file, with a
  content like this one (as the simplest and most permissive version):

      grant{
              permission java.security.AllPermission;
      };

  We strongly advice you to only grant a minimum of permissions!

- `-Djava.rmi.server.hostname=` sets the internally used hostname
  (optional, use if there are problems with name resolution)

- `IP-address` is a standard IP address for the machine this server job
  is running on, e.g., a fully qualified numerical IPv4 address like
  `121.232.212.202`, or a string-based one like “machine.domain.org” or
  possibly even “localhost” (depending on the IP-setup of this machine);

- `job type` is the optimization problem to be solved. Options are:
  `geom` (cluster structure optimization), `param` (parameter fitting),
  `switch` (switch design), this flag is mandatory;

- `input.ogo` is a regular <span class="smallcaps">ogolem</span> input
  file for the problem to be optimized (mandatory);

- `timeout` is a timeout in seconds after which a client will be purged
  from the server-internal list of active clients, provided there has
  been no communication with or heart beat from this client within this
  time interval. Default: 10s. Can be increased on stable networks.
  (optional);

- `timeoutjob` is a timeout in seconds after which a job will be purged.
  This should correspond to the expected time a global optimization step
  (for simple clients) or a chunk of steps (proxies, threading clients)
  needs for normal execution, hence this timeout must be increased for
  demanding problems. Default: 10s. (optional);

- `server-port#` is the IP port number that the server will use for task
  communication with the clients, typically/default 2500 (optional);

- `registry-port#` is the IP port number where the RMI registry is
  running on, typically/default 1099 (optional);

- `#proxies` is the number of proxies or threading clients that are
  initially expected to talk to this server. If there will be no
  proxies, specify -1 (optional). Note that `#proxies` does *not* limit
  (or otherwise affect) the number of actual proxies or threading
  clients that can connect to this server at runtime. `#proxies = -1`
  denotes the “simple client” mode, and `#proxies >= 0` denotes a
  “proxies” or “threading clients” mode, where the initial value of
  `#proxies` is only used to organize the initialization phase more
  sensibly.

- `commname#` is the String name of the RMI communication for RMI object
  export it must match on server, proxy, and client, typically/default
  RMICommunication (optional);

After `--rmiserver`, the following flags can occur in any order.

The server process will now setup and configure itself and listen to incoming connections on `server-port#` from clients.

### Starting RMI proxies or threading clients

Subsequently, proxy processes or threading clients can be started if wished.

The proxy invocation is as follows (again with everything in one line):

    java -Djava.rmi.server.codebase=file:<path/to/ogolem-jar-file> 
    -Djava.security.policy=client.policy 
    -jar <path/to/ogolem-jar-file> 
    --rmiproxy
    --jobtype <job type>
    --inputfile <input.ogo>
    --server <IP-address>
    --serverregistryport <registry-port#>
    --myserverport <myserver-port#>
    --registryport <myregistry-port#>
    --timeout <timeout> 
    --timeoutjob <timeoutjob> 
    --sleeptime <sleeptime>
    --maxchunk <chunked steps>
    --mergeinds <#indivs2merge>
    --commname <nameofregistryobject>

where

- `IP-address`: As above; however, on a remote machine, the user should
  make sure beforehand (with commands like `ping`, `telnet`, etc.) that
  the address provided here is actually resolved in the desired way;
  this is a mandatory flag;

- `registry-port#` is the port number of the registry on the server
  (optional). The additional `server-port#` provided at server startup
  need not be re-specified here since it is communicated to the client
  from the server.

- `myserver-port#` and `myregistry-port#` are for the “server-like role”
  this proxy plays towards the RMI clients to be associated with it, and
  hence are fully analogous to the corresponding flags in the actual
  server invocation, including their default values and optional status.
  Only if proxy and server are located on the same physical machine, the
  ports must differ between them;

- `sleeptime` has a dual role: It specifies the time delay between
  “heartbeat” signals sent to the server, indicating that the proxy is
  still “alive”. And it also is the time this proxy waits after
  contacting the server and finding that the server is busy (e.g.,
  waiting for another client);

- `chunked steps` is the number of global optimization steps packed
  together into one chunk of tasks communicated from the server to the
  proxy, and of results communicated back from the proxy to the server;
  actual chunks may be smaller during initialization or in the final
  phase.

- `#indivs2merge` indicates how many individuals of the proxy’s pool (at
  the time results of a chunk are communicated back to the server)
  should be merged with the server’s pool. Default: -1 (full pool)
  (optional). Note that the proxy pool size need not be equal to the
  server pool size, since the proxy has its own input (where a pool size
  is specified).

- `client.policy` is an optional setting to specify restrictions on the
  client JVM. Please see JVM documentation for options and syntax.

After `--rmiproxy`, the following flags can occur in any order.

Note that both server and proxy jobs do not do any actual calculations (global optimization steps) themselves, but initially they do parse their <span class="smallcaps">ogolem</span> input files. Hence, before simple or threading clients are started,

- there will be no optimization progress;

- but server and/or proxy jobs should report successful input parsing,

- or may exit prematurely in case of input parsing errors (in which case
  client jobs being started subsequently may fail to connect).

Please note that both proxies and threading clients may output log files, binary snapshots etc. Hence, they should be started in a directory other than the server process to avoid random overwrites.

The threading client invocation is as follows:

    java -Djava.rmi.server.codebase=file:<path/to/ogolem-jar-file> 
    -Djava.security.policy=client.policy 
    -jar <path/to/ogolem-jar-file> 
    --rmithreader 
    --jobtype <job type> 
    --inputfile <input.ogo> 
    --server <IP-address>
    --serverregistryport <registry-port#>
    --threads <#threads> 
    --sleeptime <sleeptime>
    --maxchunk <chunked steps>
    --mergeinds <#indiv2merge>
    --commname <nameofregistryobject>

Several of these flags have already been explained above. The remaining ones are:

- `registry-port#` is the port number of the registry on the server. The
  additional `server-port#` provided at server startup need not be
  re-specified here since it is communicated to the client from the
  server. In contrast to a proxy, this threading client does not play a
  dual RMI role, hence further ports for interactions with RMI clients
  need not be specified.

- `#threads` specifies the number of shared-memory-parallel threads this
  client should use on this machine. Note that a value of 1 (one) is
  also possible, which establishes a mixture of features from simple
  clients (serial execution) with those of threaded clients (chunking,
  separate input file).

### Starting RMI simple clients

Lastly, one can start simple client processes. As explained above, in “simple client” mode, the simple clients connect directly to the server. In “proxies” mode, there typically are several groups of simple clients, each of which connects to one proxy, and from the viewpoint of this simple client group, this proxy takes the role of a server.

As explained above, simple clients need only very little information from the user, hence their invocation is very simple, too:

    java -Djava.rmi.server.codebase=file:<path/to/ogolem-jar-file> 
    -Djava.security.policy=client.policy 
    -jar <path/to/ogolem-jar-file> 
    --rmiclient 
    --server <IP-address>
    --sleeptime <sleeptime>
    --serverregistryport <server-port#>
    --taskwaittime <taskwaittime>

with

- `taskwaittime` the time in seconds the client will wait for the server
  upon receiving a wait signal before contacting again. Default: 5s
  (optional).

and all other options have already been explained above, and again all options after `--rmiclient` can appear in any order. Naturally, `--server` is mandatory, but the other options are optional. If proxies are used, `--server` and `--serverregistryport` denote the IP address and registry port of the desired proxy (which takes the role of the server here); the “actual server” is then effectively invisible to the simple clients.

