#!/bin/bash -i
#
#SBATCH --job-name=SSS_180x360_P020_040_STMOD3_nis_totd_N50_OSL_ASSYM_5p0BC_native_v2
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --partition=SB3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

export OMP_NUM_THREADS=2

/home/martins/projects/XNDiff/XNDiff -l -openmp +o SSS_180x360_P020_040_STMOD3_nis_totd_N50_OSL_ASSYM_5p0BC_native_v2 -X +cif betaSSS_abc.cif cif_core.dic -sym symop_betaSSS_abc.txt -nis 0.364901 0.149400 0.365224 0.090911 0.029456 0.000109 +av 0 0.0 90.0 0.5 0.0 360.0 1.0 -sin -init_n_rand -10112 -conc def -distr +par 2.0 4.0 310.0 350.0 333.0 +f stackcpp 400.0 100.0 120.0 40.0 0.0 0.0 0.0 50 0.0 20.0 0.0 0.0 0 2 3 0.001 0.45 450 50 0 0 1 0 0 1 6 5 -silent

cd out/
tar cfvz SSS_180x360_P020_040_STMOD3_nis_totd_N50_OSL_ASSYM_5p0BC_native_v2.tar.gz --exclude=*.tar.gz SSS_180x360_P020_040_STMOD3_nis_totd_N50_OSL_ASSYM_5p0BC_native_v2* --remove-files
cd ../
