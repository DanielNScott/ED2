#!/bin/sh

#--------------------------------------------------#
# A simple script for installing via slurm queue.  #
#--------------------------------------------------#
#SBATCH -t 15              # Time required
#SBATCH --mem=2000         # Memory required
#SBATCH -p moorcroft_6100  # Queue Name
#SBATCH -J ed_install      # Job Name
#--------------------------------------------------#

./install.sh 1>qinst.out 2>qinst.err
