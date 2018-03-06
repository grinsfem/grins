#!/bin/sh
Dir=ConstraintRunYO2_0$1YO3_0$2Mdot_0$3 
mkdir $Dir
cd ./$Dir
echo "[RunVars]" >> ozone_constraint.in
echo "      Tmax = '700'" >> ozone_constraint.in
echo "      YO2_Val = '$1'" >> ozone_constraint.in
echo "      YO3_Val = '$2'" >> ozone_constraint.in
echo "      M_dot_Value = '$3'" >> ozone_constraint.in 
echo "[]" >> ozone_constraint.in
cat ../ozone_constraint.base >> ozone_constraint.in




echo "#!/bin/sh" >> Run_Constraint.sh
echo "GRINS_RUN=\${GRINS_RUN:-\$LIBMESH_RUN}" >> Run_Constraint.sh
echo "DEFAULT_SOLVER_OPTIONS=\"-ksp_type gmeres -pc_type bjacobi -sub_pc_type lu -sub_pc_factor_shift_type nonzero\"" >> Run_Constraint.sh
echo "GRINS_SOLVER_OPTIONS=\${GRINS_SOLVER_OPTIONS:-\$LIBMESH_OPTIONS:\$DEFAULT_SOLVER_OPTIONS}" >> Run_Constraint.sh
echo "\$GRINS_RUN /zoidberg1/data/shared/klbudzin/WorkingGrins/OptBuild/bin/grins /zoidberg1/data/shared/klbudzin/Flames/OzoneOD/$Dir/ozone_constraint.in \$GRINS_SOLVER_OPTIONS" >> Run_Constraint.sh
chmod u+x Run_Constraint.sh


