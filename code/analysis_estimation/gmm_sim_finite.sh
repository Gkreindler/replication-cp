#!/bin/bash
#SBATCH -o log_sim_finite.out
#SBATCH -e log_sim_finite.err
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 44
#SBATCH -t 0-08:00
#SBATCH --mem=160G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gabriel.kreindler@gmail.com

export PATH=$PATH:/n/holylabs/LABS/kreindler_lab/Lab/software/julia-1.9.1/bin

julia numerical_identification/estimation_simulated_data_fin.jl

