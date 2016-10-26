#!/bin/bash

for f in `ls jobs/const/job_*`
do

#Match NP N NS in job_NP_N_NS.out filename
echo $f | sed -e 's/job_\([0-9]*\)_\([0-9]*\)_\([0-9]*\).*out/\1/g'

#Get compute time
#perl: match non-digit characters \D+, forget everything \K, match decimal digit with decimal dot
#cat $f | grep "Total CPU secs in IPOPT" | perl -lne 'print $& if /\D+\K\d+.\d+/'

#----------------------------------------------------------------------

#cat $f | grep "Init "  |sed -e 's/Init      : \([0-9].[0-9]*e+[0-9]*\) sec/\1/g'
cat $f | grep "Schur "  |sed -e 's/Schur     : \([0-9].[0-9]*e+[0-9]*\) sec/\1/g'
#cat $f | grep "Solution "  |sed -e 's/Solution  : \([0-9].[0-9]*e+[0-9]*\) sec/\1/g'
done
