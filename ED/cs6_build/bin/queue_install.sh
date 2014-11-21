#!/bin/sh

#SBATCH -n 1
#SBATCH -t 0
#SBATCH -p moorcroft_6100
#SBATCH --mem=5000
#SBATCH -J Install_r85-DS-Iso-E

install.sh 1>qinst.out 2>qinst.err
