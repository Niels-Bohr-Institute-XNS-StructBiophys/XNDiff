#!/bin/bash -i
#
# use -i flag here to use aliases like XNDiff
#
# screen -d -m ./20190320_betaSSS_abc_disl_dosl_scan_LARGER.sh > 20190320_betaSSS_abc_disl_dosl_scan_LARGER.log


ddisl=0.4
disl0=0.4
imax=35

ddosl=0.4
dosl0=0.4
jmax=35

#dR0=1
#ddR=2.0
#kmax=2


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

#for (( k = 0 ; k<=$kmax; k++ ))
#do



disl=$(echo $i*$ddisl+$disl0 | bc)
dosl=$(echo $j*$ddosl+$dosl0 | bc)

dtot=$(echo $disl+$dosl | bc)

# 7.2-10.0
upper_check=$(echo "$dtot > 10.0" | bc)
lower_check=$(echo "$dtot <  7.2" | bc)


if [ $upper_check -eq 1 -o $lower_check -eq 1 ]
then
continue
fi

# cut old ones with disl<=6.4, dosl<=6.4 and dtot<=8.4
#octave:17> 19*5+2*6+5+4+3+2+1
#ans =  122

upper_check=$(echo "$dtot <  8.4" | bc)
lower_check=$(echo "$disl <  6.8" | bc)
lower_check2=$(echo "$dosl < 6.8" | bc)

if [[ $upper_check -eq 1 && $lower_check -eq 1 && $lower_check2 -eq 1 ]]
then
continue
fi

#dR=$(echo $ddR^$k*$dR0 | bc)

(( count++ ))

disl_str=$(printf "%.1lf" $disl)
dosl_str=$(printf "%.1lf" $dosl)
#dR_str=$(printf "%.2lf" $dR)
#echo "disl=$disl_str dosl=$dosl_str dR=$dR_str"
echo "disl=$disl_str dosl=$dosl_str"

disl_str2=$(printf "%04.1lf" $disl)
dosl_str2=$(printf "%04.1lf" $dosl)
disl_str2=$(echo $disl_str2 | sed 's#\.##g')
dosl_str2=$(echo $dosl_str2 | sed 's#\.##g')
#dR_str2=$(echo $dR_str | sed 's#\.##g')


jobname="SSS_180x360_P${disl_str2}_${dosl_str2}_ST_360_400_N50_OSL_LARGER"


output="$jobname.sh"

echo "#"'!'"/bin/bash -i" > $output
echo "#" >> $output
echo "#SBATCH --job-name=$jobname" >> $output
echo "#SBATCH --nodes=1" >> $output
echo "#SBATCH --time=48:00:00" >> $output
echo "#SBATCH --partition=SB3" >> $output
echo "#SBATCH --ntasks=1" >> $output
echo "#SBATCH --cpus-per-task=2" >> $output
echo "" >> $output
echo "export OMP_NUM_THREADS=2" >> $output
echo "" >> $output
echo "$XNDiff -l -openmp +o $jobname -X +cif betaSSS_abc.cif cif_core.dic -sym symop_betaSSS_abc.txt +av 0 0.0 90.0 0.5 0.0 360.0 1.0 -sin -init_n_rand -10112 -conc def -distr +par $disl_str $dosl_str 310.0 350.0 333.0 +f stackcpp 282.8 57.6 212.1 42.4 2.522 0.382 0.0 20.0 0.0 20.0 36.0 4.0 0 1 0 0.001 0.45 450 50 0 0 1 0 0 1 6 1 -silent -z 10 -conv 0" >> $output
echo "" >> $output
echo "cd out/" >> $output
echo "tar cfvz $jobname.tar.gz $jobname* --remove-files" >> $output
echo "cd ../" >> $output

chmod +x $output

echo "sbatch $output"
sbatch $output

#sleep 0.5

#done

done

done

echo $count




