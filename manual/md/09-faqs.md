# FAQs

## <span class="smallcaps">ogolem</span> doesn’t run and/or behaves crazy

Ensure that your environment is properly setup. In particular reduce the input to a minimal example and check the following things:

- A JRE (preferably openJDK/JRE or Sun/Oracle-JDK/JRE) of at least
  version 7 exists and works on your system.

- If an external program is called, check whether is is properly
  installed and all conditions specified in the corresponding section
  are met.

- Check whether your input is sane (in particular blow factors for
  global geometry optimizations) and makes sense.

- set DebugLevel= to 1 or higher and check the additional output for
  hints.

## I have problems running an external program

Most of the times, this is due to the fact that the program cannot be found in the \$PATH environment variable (for unixoid systems). Setup a wrapper script

    #!/bin/sh
    PROGRAM $@

in your \$HOME/bin directory named like the `PROGRAM` you want to call, mark it executable and see if it works (you might need to setup necessary environment variables and/or load modules in this script as well). Debug statements might help to make sure the wrapper is being executed:

    #!/bin/sh
    echo "Calling external program..." >> $HOME/debug_ogolem
    PROGRAM $@
    echo "Finished calling external program." >> $HOME/debug_ogolem

## I want to report a bug, how does it work and what is required?

In case your data is not confidential (or you can reduce the data set such that it is not confidential), please file a bug report to `developers@ogolem.org`, providing your input(s) and a description of the problems encountered. We will treat the data confidentially.

Once the bug has been found and fixed, you can obtain a new snapshot of the ogolem.jar from the <span class="smallcaps">ogolem</span> homepage.

## Does <span class="smallcaps">ogolem</span> loose or change features over time?

In general, we strive to remain backwards compatible with regards to both features and input. However, we will remove features that have been superseded by superior solutions, with the same rationale applying to input formats. Either of the two has so far happened only most seldomly. If a regression should be introduced, please notify the authors as we will want to fix it.