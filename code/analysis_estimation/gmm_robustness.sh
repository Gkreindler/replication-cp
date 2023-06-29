#!/bin/bash
#SBATCH -o log_robustness.out
#SBATCH -e log_robustness.err
#SBATCH -p bigmem
#SBATCH -N 1
#SBATCH -n 44
#SBATCH -t 3-00:00
#SBATCH --mem=360G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gabriel.kreindler@gmail.com

module load Julia/1.7.1-linux-x86_64
julia estimation_robust1_static.jl

julia estimation_robust2_asym_switch.jl

julia estimation_robust3_FE.jl

julia estimation_robust4_half_attention.jl

julia estimation_robust5_propto_wage.jl

julia estimation_robust6_single_arrival_time.jl
