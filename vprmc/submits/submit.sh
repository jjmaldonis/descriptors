#!/bin/bash

#$ -N hrmc

#$ -q all.q
##$ -pe orte 1

#$ -l h_rt=INFINITY
#$ -cwd
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -V

# Write the date and time to outfile.
echo "Submitted on:"
date '+%s'
echo $@

echo "Got $NSLOTS processors."

./rmc $JOB_ID src/temp.txt ../testmodels/Zr50Cu35Al15_start.xyz

echo "Finished on:"
date '+%s'
