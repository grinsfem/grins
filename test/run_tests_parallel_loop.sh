#!/bin/sh

# This script will run `make check` with various settings of LIBMESH_RUN.
# In particular, LIBMESH_RUN will start at mpiexec -np 1 and the number of
# processors will increase. The max value can be set in the envrionment using
# GRINS_TESTING_MAX_N_PROCS or it can be the first argument supplied to
# this script. The argument is preferred if both are specified. It is an
# error if neither is specified. Note that the -j argument to check is halved
# each time and starts at the max.
#
# An optional second argument is allowed. If second argument is present,
# then all tests are run each time. If you only want to run a subset of tests,
# then supply the argument TESTS="<space separated list of tests>" (exactly
# as you supply to normal make check). For example, if we only wanted to run the
# following two tests for each pass, supply the following as the second argument:
# TESTS="regression/thermally_driven_2d_flow.sh regression/reacting_low_mach_antioch_statmech_constant.sh"
#
# The output goes GRINS_TESTING_OUTPUT_FILE if it's set in the environment
# or to parallel_loop_check_tests.log if it's not. A handy command for parsing
# the output is:
# grep -e LIBMESH_RUN -e FAIL <file> | grep -v XFAIL

set -e

export GRINS_TESTING_OUTPUT_VALUE=${GRINS_TESTING_OUTPUT_FILE:-parallel_loop_check_tests.log}

echo $GRINS_TESTING_OUTPUT_VALUE

export GRINS_TESTING_MAX_N_PROCS_VALUE=$GRINS_TESTING_MAX_N_PROCS

if [ -n $1 ]; then
   export GRINS_TESTING_MAX_N_PROCS_VALUE=$1
fi

if [ -n $2 ]; then
   export GRINS_TESTING_TESTS_SUBSET=$2
fi

for i in $(seq 1 ${GRINS_TESTING_MAX_N_PROCS_VALUE:?undefined}); \
   do export LIBMESH_RUN="mpiexec -np $i"; \
   j=$((${GRINS_TESTING_MAX_N_PROCS_VALUE:?undefined}/$i)); \
   echo LIBMESH_RUN is $LIBMESH_RUN | tee -a ${GRINS_TESTING_OUTPUT_VALUE}; \
   make -j $j $GRINS_TESTING_TESTS_SUBSET check 2>&1 | tee -a ${GRINS_TESTING_OUTPUT_VALUE}; \
done
