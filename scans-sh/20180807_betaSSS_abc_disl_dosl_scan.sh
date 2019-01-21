#!/bin/bash -i
#
# use -i flag here to use aliases like XNDiff
#
# screen -m -d ./betaSSS_abc_disl_dosl_scan.sh > betaSSS_abc_disl_dosl_scan.log


ddisl=0.2
disl0=0.4
imax=6

ddosl=0.2
dosl0=0.4
jmax=8

echo 0*$ddisl+$disl0 | bc
echo 0*$ddosl+$dosl0 | bc

echo $imax*$ddisl+$disl0 | bc
echo $jmax*$ddosl+$dosl0 | bc


count=0

for (( i = 0 ; i<=$imax; i++ ))
do

for (( j = 0 ; j<=$jmax; j++ ))
do

disl=$(echo $i*$ddisl+$disl0 | bc)
dosl=$(echo $j*$ddosl+$dosl0 | bc)

#dtot=$(echo $disl+$dosl | bc)

(( count++ ))

disl_str=$(printf "%.1lf" $disl)
dosl_str=$(printf "%.1lf" $dosl)

echo "disl=$disl_str dosl=$dosl_str"

disl_str2=$(echo $disl_str | sed 's#\.##g')
dosl_str2=$(echo $dosl_str | sed 's#\.##g')

jobname="SSS_180x360_P${disl_str2}_${dosl_str2}_ST_360_400_N100"
echo "run $jobname"

#XNDiff -l -openmp +o $jobname +cif betaSSS_abc.cif cif_core.dic -sym symop_betaSSS_abc.txt +av 0 0.0 90.0 0.5 0.0 360.0 1.0 -sin -init_n_rand -10112 -conc def -distr +par $disl_str $dosl_str 310.0 350.0 333.0 -0.4939 1.1995 6.3664 0.0 0.0 0.0 +f stackcpp 200.0 40.0 150.0 30.0 2.522 0.382 0.0 20.0 0.0 20.0 36.0 4.00 0 1 0 0.001 0.6 600 100 0 0 1 0 0 1 7 2 -silent -z 20

#sleep 0.5

cd out/
tar cfvz $jobname.tar.gz $jobname* --remove-files
#tar xzvf $jobname.tar.gz
cd ../

#sleep 0.5

done

done

echo $count




