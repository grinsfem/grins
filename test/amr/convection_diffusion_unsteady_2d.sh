#!/bin/bash

set -e

INPUT="${GRINS_TEST_SRCDIR_DIR}/amr/input_files/convection_diffusion_unsteady_2d_amr.in"

TESTDATA="./convection_diffusion_unsteady_2d_amr"
MESHDATA=$TESTDATA

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT

# Untar the files into a local directory in the build tree
LOCALTESTDIR="./convection_diffusion_unsteady_2d_amr_test_tmp"

mkdir $LOCALTESTDIR
tar -xzf ${GRINS_TEST_SRCDIR_DIR}/amr/gold_data/convection_diffusion_unsteady_2d/files.tar.gz --directory $LOCALTESTDIR

# Now give prefixes for the gold data
GOLDDATA=$LOCALTESTDIR/convection_diffusion_unsteady_2d_amr
GOLDMESH=$GOLDDATA

# Now run the test part to make sure we're getting the correct thing
${GRINS_TEST_DIR}/generic_amr_testing_app \
                 input=$INPUT \
                 vars='u' \
                 norms='L2' \
                 tol='1.0e-10' \
                 test_data_prefix=$TESTDATA \
                 mesh_data_prefix=$MESHDATA \
                 gold_data_prefix=$GOLDDATA \
                 gold_mesh_prefix=$GOLDMESH \
                 n_steps=25

# Now remove the test turd
rm -rf $LOCALTESTDIR
rm $TESTDATA*.xdr $MESHDATA*.xda
