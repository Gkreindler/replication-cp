#!/bin/bash
#SBATCH -o log_robustness_proptowage.out
#SBATCH -e log_robustness_proptowage.err
#SBATCH -p bigmem
#SBATCH -N 1
#SBATCH -n 44
#SBATCH -t 0-24:00
#SBATCH --mem=300G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gabriel.kreindler@gmail.com

export PATH=$PATH:/n/holylabs/LABS/kreindler_lab/Lab/software/julia-1.9.1/bin

julia estimation_robust5_propto_wage.jl

