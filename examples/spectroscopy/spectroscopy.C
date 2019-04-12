//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

// GRINS
#include "grins/runner.h"
#include "grins/string_utils.h"
#include "grins/absorption_coeff.h"
#include "grins/simulation.h"
#include "grins/multiphysics_sys.h"
#include "grins/composite_qoi.h"
#include "grins/spectroscopic_absorption.h"

#if GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

int main(int argc, char* argv[])
{
#if GRINS_HAVE_ANTIOCH
  GRINS::Runner runner(argc,argv);
  runner.init();

  // Parse the wavenumber range to plot over from the input file
  const GetPot & input = runner.get_input_file();

  libMesh::Real nu_min  = input("SpectroscopyExample/wavenumber_range",-1.0,0);
  libMesh::Real nu_max  = input("SpectroscopyExample/wavenumber_range",-1.0,1);
  libMesh::Real nu_step = input("SpectroscopyExample/wavenumber_range",-1.0,2);

  if ( (nu_min < 0) || (nu_max < 0) || (nu_step < 0) )
    libmesh_error_msg("ERROR: please specify a wavenumber range for the SpectroscopyExample as 'nu_min nu_max nu_step'");

  std::string filename = input("SpectroscopyExample/output_prefix","");
  if (filename == "")
    libmesh_error_msg("ERROR: please specify an output file prefix for the SpectroscopyExample");

  runner.run();

  GRINS::Simulation & sim = runner.get_simulation();
  GRINS::MultiphysicsSystem * system = sim.get_multiphysics_system();

  GRINS::CompositeQoI * comp_qoi = libMesh::cast_ptr<GRINS::CompositeQoI*>(system->get_qoi());
  GRINS::SpectroscopicAbsorption & qoi = libMesh::cast_ref<GRINS::SpectroscopicAbsorption &>(comp_qoi->get_qoi(0));
  GRINS::AbsorptionCoeff<GRINS::AntiochChemistry> & abs_coeff = libMesh::cast_ref<GRINS::AbsorptionCoeff<GRINS::AntiochChemistry> &>(qoi.get_function());

  libMesh::QoISet qs;
  qs.add_index(0);

  std::ofstream output;
  if (system->get_mesh().comm().rank() == 0)
    {
      output.open(filename+".dat",std::ofstream::app);
    }

  for (libMesh::Real nu = nu_min; nu <= nu_max; nu += nu_step)
    {
      abs_coeff.set_wavenumber(nu);

      system->assemble_qoi(qs);
      libMesh::Real qoi = sim.get_qoi_value(0);

      if (system->get_mesh().comm().rank() == 0)
        output <<std::fixed <<std::setprecision(8) <<nu <<"," <<std::setprecision(16) <<qoi <<std::endl;
    }

  if (system->get_mesh().comm().rank() == 0)
    output.close();
#else
  libmesh_error_msg("ERROR: GRINS must be built with Antioch to use the Spectroscopy example. Please reconfigure your build to include the Antioch library.");
#endif
  return 0;
}

