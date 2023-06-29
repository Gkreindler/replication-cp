#!/bin/bash
#SBATCH -o log_het_boot.out
#SBATCH -e log_het_boot.err
#SBATCH -p bigmem
#SBATCH -N 1
#SBATCH -n 43
#SBATCH -t 1-00:00
#SBATCH --mem=450G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gabriel.kreindler@gmail.com

export PATH=$PATH:/n/holylabs/LABS/kreindler_lab/Lab/software/julia-1.9.1/bin

julia runsims-het-boot.jl
