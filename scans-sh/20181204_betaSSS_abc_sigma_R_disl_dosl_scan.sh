#!/bin/bash -i
#
# use -i flag here to use aliases like XNDiff
#
# screen -d -m ./20181204_betaSSS_abc_sigma_R_disl_dosl_scan.sh > 20181204_betaSSS_abc_sigma_R_disl_dosl_scan.log
#
# note that "./" before the *.sh file is mandatory to run the shell script !!!


ddisl=0.2
disl0=0.4
imax=8

ddosl=0.2
dosl0=0.4
jmax=8

dR0=1
ddR=2.0
kmax=2


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

disl=$(echo $i*$ddisl+$disl0 | bc)
dosl=$(echo $j*$ddosl+$dosl0 | bc)

dtot=$(echo $disl+$dosl | bc)

upper_check=$(echo "$dtot > 3.0" | bc)
lower_check=$(echo "$dtot < 1.0" | bc)

if [ $upper_check -eq 1 -o $lower_check -eq 1 ]
then
continue
fi

dR=$(echo $ddR^$k*$dR0 | bc)

(( count++ ))

disl_str=$(printf "%.1lf" $disl)
dosl_str=$(printf "%.1lf" $dosl)
dR_str=$(printf "%.2lf" $dR)

echo "disl=$disl_str dosl=$dosl_str dR=$dR_str"

disl_str2=$(echo $disl_str | sed 's#\.##g')
dosl_str2=$(echo $dosl_str | sed 's#\.##g')
dR_str2=$(echo $dR_str | sed 's#\.##g')

jobname="SSS_180x360_P${disl_str2}_${dosl_str2}_ST_230_${dR_str2}_N100"
echo "run $jobname"

# cis: 0.0507 0.00696 0.452 0.317 0.0846 0.0885 -> skip all except 3,4 use only nsp=5, nst=5
XNDiff -l -openmp +o $jobname -X +cif betaSSS_abc.cif cif_core.dic -sym symop_betaSSS_abc.txt -cis 0.0 0.0 0.452 0.317 0.0 +av 0 0.0 90.0 0.5 0.0 360.0 1.0 -sin -init_n_rand -10112 -conc def -distr +par $disl_str $dosl_str 310.0 350.0 333.0 +f stackcpp 200.0 40.0 150.0 30.0 0.0 0.0 0.0 20.0 0.0 20.0 23.0 $dR_str 0 2 0 0.001 0.45 450 100 0 0 1 0 0 1 5 5 -silent -z 50

sleep 0.5

cd out/
tar cfvz $jobname.tar.gz $jobname* --remove-files
cd ../

sleep 0.5

done

done

done

echo $count




