#!/bin/bash
#SBATCH -o log_delta.out
#SBATCH -e log_delta.err
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 44
#SBATCH -t 2-00:00
#SBATCH --mem=180G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gabriel.kreindler@gmail.com

module load Julia/1.7.1-linux-x86_64
julia estimation_delta_d99.jl

julia estimation_delta_estimation.jl

julia estimation_delta_other.jl
