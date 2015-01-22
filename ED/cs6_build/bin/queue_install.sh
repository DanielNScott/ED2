#!/bin/sh

#SBATCH -n 1
#SBATCH -t 180
#SBATCH -p moorcroft_6100
#SBATCH --mem=0
#SBATCH -J Install_r85-SMP

install.sh 1>qinst.out 2>qinst.err
