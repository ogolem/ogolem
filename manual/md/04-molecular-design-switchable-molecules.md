# Molecular design: Switchable molecules

**IMPORTANT:   This part of the manual is not complete yet!**

## Switch analysis

<span class="smallcaps">ogolem</span> provides means to check a running switch design job using the written restart files. For this, you call the `ogolem.jar` with the

    --beswitched

option.

The analysis program needs a binary pool, either the final `switchespool.bin` or the temporary `snapswitchpool.bin`. This is specified using the `-i` flag followed by whitespace, e.g.

    java -jar ogolem.jar --beswitched -i snapswitchpool.bin

Some other options are possible

- `-collblow`  
  the blow factor for the collision detection (needed for sidechain
  glueing), defaults to 1.4.

- `-getswitches`  
  create a directory named `switches` and write the switched of this
  pool there in the human-readable format.

