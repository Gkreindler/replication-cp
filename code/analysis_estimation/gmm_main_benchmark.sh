#!/bin/bash
#SBATCH -o main_benchmark.out
#SBATCH -e main_benchmark.err
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 44
#SBATCH -t 1-00:00
#SBATCH --mem=180G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gabriel.kreindler@gmail.com

module load Julia/1.7.1-linux-x86_64
julia estimation_main_benchmark.jl
