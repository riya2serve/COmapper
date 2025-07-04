#!/bin/bash

# the home directory of git repository
dir_home=/Users/riyarampalli/hacks/COmapper

# the number of available threads
n_parallel=4

# The name of input file (fq.gz)
input=sim_lyrata.fq.gz

# Location of minimap reference of TAIR10.
#ref=/path/to/TAIR10.mmi
ref=A_thaliana.mmi

# Location of masksam.awk. 
awkloc=resources/masksam_Q14.awk

export dir_home
export n_parallel
export input
export ref
export awkloc
