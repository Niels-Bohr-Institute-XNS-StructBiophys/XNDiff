#!/bin/bash -i
#
#SBATCH --job-name=SSS_180x360_P020_012_ST_225_050_N50_OSL_SPST_123
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --partition=SB3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

export OMP_NUM_THREADS=2

/home/martins/projects/XNDiff/XNDiff -l -openmp +o SSS_180x360_P020_012_ST_225_050_N50_OSL_SPST_123 -X +cif betaSSS_abc.cif cif_core.dic -sym symop_betaSSS_abc.txt -cis 0.052586372247525714 0.004510639736581632 0.42054301527185284 0.0 0.0 0.0 +av 0 0.0 90.0 0.5 0.0 360.0 1.0 -sin -init_n_rand -10112 -conc def -distr +par 2.0 1.2 310.0 350.0 333.0 +f stackcpp 200.0 40.0 150.0 30.0 0.0 0.0 0.0 20.0 0.0 20.0 22.5 0.50 0 2 0 0.001 0.45 450 50 0 0 1 0 0 1 6 5 -silent

cd out/
tar cfvz SSS_180x360_P020_012_ST_225_050_N50_OSL_SPST_123.tar.gz SSS_180x360_P020_012_ST_225_050_N50_OSL_SPST_123* --remove-files
cd ../
