#!/bin/bash -i
#
# use -i flag here to use aliases like XNDiff
#
# screen -d -m ./20190308_betaSSS_abc_R_sigma_R_disl_dosl_STMOD4_scan.sh > 20190308_betaSSS_abc_R_sigma_R_disl_dosl_STMOD4_scan.log



# more MathematicaOut/SSS_0p0BC_1to5dil_OSL_005/SSS_180x360_P020_012_ST_360_400_N50_OSL.log
# c1 = 0.052586372247525714
# c2 = 0.004510639736581632
# c3 = 0.42054301527185284
# c4 = 0.35619147143577523
# c5 = 0.11911260635744349
# c6 = 0.047055894950821076



# 3*45+2*32=199 !< min(D)=213 Å
# 4*45+2*32=244 > min(D)=213 Å !!!
# -> use 1 to 4 platelet thicknesses with new stack mode 4, which pushes D[k] so much that the particles do not overlap

# R=21.3:0.2:21.7, dR=1.0,1.5 dosl=disl=8:4:24 with dtot<=32

ddisl=0.4
disl0=0.8
imax=4

ddosl=0.4
dosl0=0.8
jmax=4

dR0=1.0
ddR=1.4142
kmax=1

D0=21.3
dD=0.2
lmax=2


XNDiff='/home/martins/projects/XNDiff/XNDiff'

echo 0*$ddisl+$disl0 | bc
echo 0*$ddosl+$dosl0 | bc

echo $imax*$ddisl+$disl0 | bc
echo $jmax*$ddosl+$dosl0 | bc


count=0

for (( i = 0 ; i<=$imax; i++ ))
do

for (( j = 0 ; j<=$jmax; j++ ))
do

for (( k = 0 ; k<=$kmax; k++ ))
do

for (( l = 0 ; l<=$lmax; l++ ))
do


disl=$(echo $i*$ddisl+$disl0 | bc)
dosl=$(echo $j*$ddosl+$dosl0 | bc)

dtot=$(echo $disl+$dosl | bc)

# < 32.0 only
upper_check=$(echo "$dtot >  3.2" | bc)
lower_check=$(echo "$dtot <  0.8" | bc)


if [ $upper_check -eq 1 -o $lower_check -eq 1 ]
then
continue
fi


dR=$(echo $ddR^$k*$dR0 | bc)
D=$(echo $l*$dD+$D0 | bc)

(( count++ ))

disl_str=$(printf "%.1lf" $disl)
dosl_str=$(printf "%.1lf" $dosl)
D_str=$(printf "%.1lf" $D)
dR_str=$(printf "%01.2lf" $dR)
echo "disl=$disl_str dosl=$dosl_str D=$D_str dR=$dR_str"

disl_str2=$(printf "%04.1lf" $disl)
dosl_str2=$(printf "%04.1lf" $dosl)
disl_str2=$(echo $disl_str2 | sed 's#\.##g')
dosl_str2=$(echo $dosl_str2 | sed 's#\.##g')
D_str2=$(echo $D_str | sed 's/\.//g')
dR_str2=$(echo $dR_str | sed 's/\.//g')


jobname="SSS_180x360_P${disl_str2}_${dosl_str2}_ST_${D_str2}_${dR_str2}_N50_OSL_SPST_1234_STMOD4"
echo $jobname


output="$jobname.sh"

echo "#"'!'"/bin/bash -i" > $output
echo "#" >> $output
echo "#SBATCH --job-name=$jobname" >> $output
echo "#SBATCH --nodes=1" >> $output
echo "#SBATCH --time=48:00:00" >> $output
echo "#SBATCH --partition=SB4" >> $output
echo "#SBATCH --ntasks=1" >> $output
echo "#SBATCH --cpus-per-task=2" >> $output
echo "" >> $output
echo "export OMP_NUM_THREADS=2" >> $output
echo "" >> $output
echo "$XNDiff -l -openmp +o $jobname -X +cif betaSSS_abc.cif cif_core.dic -sym symop_betaSSS_abc.txt -cis 0.052586372247525714 0.004510639736581632 0.42054301527185284 0.35619147143577523 0.0 0.0 +av 0 0.0 90.0 0.5 0.0 360.0 1.0 -sin -init_n_rand -10112 -conc def -distr +par $disl_str $dosl_str 310.0 350.0 333.0 +f stackcpp 200.0 40.0 150.0 30.0 0.0 0.0 0.0 20.0 0.0 20.0 $D_str $dR_str 0 2 4 0.001 0.45 450 50 0 0 1 0 0 1 6 5 -silent" >> $output
echo "" >> $output
echo "cd out/" >> $output
echo "tar cfvz $jobname.tar.gz $jobname* --remove-files" >> $output
echo "cd ../" >> $output

chmod +x $output

echo "sbatch $output"
sbatch $output
#sleep 0.5

done

done

done

done

echo $count






