#!/bin/bash

N=300  # dimension of the grid
NS=128 # number of scenarios
SOLVER=schur #linear solver
SUFFIX="" #e.g. "_weak" suffix to the name of the output job
WEAK=0

##############################################################################
# This scipt sets proper options in ipopt.opt required by SchurSolver.
# Particularly the dimension of the grid N and the number of stochastic scenarios NS.
# These parameter has to match with runtime parameters of the solve_problems

MKLROOT=/apps/intel/17.0.0/mkl
SCHUR_DIR=/home/kardos/PowerGrid

#location of the executable
PREFIX=/home/kardos/optimalControl/build/Ipopt/examples/ScalableProblems

# check command line arguments
#if [ $# -ne 4 ]; then
#    echo 'Usage: $./run.sh NP N NS SOLVER'
#    echo 'SOLVER is pardiso or schur'
#    exit
#fi

for comp_nodes in 1 2 4 8; do

#master process is not computing, only distributes work
tot_nodes=$((comp_nodes + 1))

#adjust number of scenarios for weak scaling test
if [ $WEAK -eq 1 ]
then
    ns=$((comp_nodes * NS)) #constant NS per node for weak scaling
else
    ns=$NS
fi

#run always single process for pardiso
if [ "$SOLVER" = "pardiso" ]
then
    tot_nodes=1
fi

#set ipopt.opt options N NS and solver
if [ -e ${PREFIX}/ipopt.opt ]; then
    sed -i "s/^problem_dimension\ [0-9]*[[:space:]]*$/problem_dimension\ $N/g" ${PREFIX}/ipopt.opt
    sed -i "s/^problem_scenarios\ [0-9]*[[:space:]]*$/problem_scenarios\ $ns/g" ${PREFIX}/ipopt.opt
    sed -i "s/^linear_solver\ [a-z]*[[:space:]]*$/linear_solver\ $SOLVER/g" ${PREFIX}/ipopt.opt
else
    echo "problem_dimension $N" > ${PREFIX}/ipopt.opt
    echo "problem_scenarios $NS" >> ${PREFIX}/ipopt.opt
    echo "linear_solver $SOLVER" >> ${PREFIX}/ipopt.opt
fi

sbatch <<-_EOF
#!/bin/bash
#SBATCH --job-name=OC_${comp_nodes}_${N}_${NS}
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=${tot_nodes}
#SBATCH --time=03:30:00
#SBATCH --output=job_${comp_nodes}_${N}_${NS}${SUFFIX}.out

export OMP_NUM_THREADS=8
export PARDISOLICMESSAGE=1

LD_LIBRARY_PATH=${SCHUR_DIR}/lib:${MKLROOT}/lib/intel64:${LD_LIBRARY_PATH} \
mpirun -np ${tot_nodes} -N 1 ${PREFIX}/solve_problem $N $ns

_EOF
done

