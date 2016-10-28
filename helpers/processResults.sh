#!/bin/bash

#phases of the schur solve process
phase=(Init Schur Solution)

for match in "${phase[@]}"
do

    printf "Results for ----${match}---- phase:\n"
    printf "Id Nodes Time Ipopt_iter Schur_iter Residual\n"
    
    id=0
    
    for f in `ls $1/job_*`
    do
        
        #Match NP N NS in job_NP_N_NS.out filename
        NP=`echo $f | sed -e 's/.*job_\([0-9]*\)_\([0-9]*\)_\([0-9]*\).*out/\1/g'`
        
        T=`cat $f | grep "${match} "  | sed -e "s/${match}[ ]*: \([0-9].[0-9]*e+[0-9]*\) sec/\1/g" |\
            awk '{ total += $1; count++ } END { print total/count }'`
        
        #how many times the schur solve has been performed
        SCHUR_ITER=`cat $f | grep "${match} "  | sed -e "s/${match}[ ]*: \([0-9].[0-9]*e+[0-9]*\) sec/\1/g" | wc -l`
    
        #how many IPOPT iteration were performed
        ITER=`cat $f | grep "Number of Iterations" | sed -e "s/^[^0-9]*\([0-9]*\)/\1/g"`
        
        #Get compute time
        #perl: match non-digit characters \D+, forget everything \K, match decimal digit with decimal dot
        #cat $f | grep "Total CPU secs in IPOPT" | perl -lne 'print $& if /\D+\K\d+.\d+/'
        
        id=$((id + 1))
        
        printf "${id} ${NP} ${T} ${ITER} ${SCHUR_ITER}\n"
    done
    printf "\n"

done

#Id Nodes Max_iter Time Residual
#1 1 1 1800.586205 4.29742e-09
#2 2 28 457.139046 0.0006036
#3 4 38 124.963360 0.000707
#4 8 50 33.557306 0.000797479
