#!/bin/bash -i
#
#SBATCH --job-name=SSS_180x360_P028_040_ST_360_400_N50_OSL_HUGE
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --partition=SB3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

export OMP_NUM_THREADS=2

/home/martins/projects/XNDiff/XNDiff -l -openmp +o SSS_180x360_P028_040_ST_360_400_N50_OSL_HUGE -X +cif betaSSS_abc.cif cif_core.dic -sym symop_betaSSS_abc.txt +av 0 0.0 90.0 0.5 0.0 360.0 1.0 -sin -init_n_rand -10112 -conc def -distr +par 2.8 4.0 310.0 350.0 333.0 +f stackcpp 400.0 80.0 300.0 60.0 2.522 0.382 0.0 20.0 0.0 20.0 36.0 4.0 0 1 0 0.001 0.45 450 50 0 0 1 0 0 1 7 1 -silent -z 10 -conv 0

cd out/
tar cfvz SSS_180x360_P028_040_ST_360_400_N50_OSL_HUGE.tar.gz SSS_180x360_P028_040_ST_360_400_N50_OSL_HUGE* --remove-files
cd ../
