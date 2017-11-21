//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

// This class
#include "grins/qoi_factory.h"

// GRINS
#include "grins/string_utils.h"
#include "grins/physics_naming.h"
#include "grins/qoi_names.h"
#include "grins/average_nusselt_number.h"
#include "grins/vorticity.h"
#include "grins/parsed_boundary_qoi.h"
#include "grins/parsed_interior_qoi.h"
#include "grins/weighted_flux_qoi.h"
#include "grins/integrated_function.h"
#include "grins/hitran.h"
#include "grins/absorption_coeff.h"
#include "grins/spectroscopic_absorption.h"
#include "grins/chemistry_builder.h"

#if GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

#if GRINS_HAVE_CANTERA
#include "grins/cantera_mixture.h"
#endif

// libMesh
#include "libmesh/parsed_function.h"

namespace GRINS
{
  QoIFactory::QoIFactory()
  {
    return;
  }

  QoIFactory::~QoIFactory()
  {
    return;
  }

  SharedPtr<CompositeQoI> QoIFactory::build(const GetPot& input)
  {
    std::string qoi_list = input("QoI/enabled_qois", "none" );

    std::vector<std::string> qoi_names;

    if( qoi_list != std::string("none") )
      {
        StringUtilities::split_string( qoi_list, std::string(" "), qoi_names );
      }

    SharedPtr<CompositeQoI> qois( new CompositeQoI );

    if( !qoi_names.empty() )
      {
        for( std::vector<std::string>::const_iterator name = qoi_names.begin();
             name != qoi_names.end(); ++name )
          {
            this->add_qoi( input, *name, qois );

            this->check_qoi_physics_consistency( input, *name );
          }

        if( input( "screen-options/echo_qoi", false ) )
          {
            this->echo_qoi_list( qois );
          }
      }

    return qois;
  }

  void QoIFactory::add_qoi( const GetPot& input, const std::string& qoi_name, SharedPtr<CompositeQoI>& qois )
  {
    QoIBase* qoi = NULL;

    if( qoi_name == avg_nusselt )
      {
        qoi = new AverageNusseltNumber( avg_nusselt );
      }

    else if( qoi_name == parsed_boundary )
      {
        qoi =  new ParsedBoundaryQoI( parsed_boundary );
      }

    else if( qoi_name == parsed_interior )
      {
        qoi =  new ParsedInteriorQoI( parsed_interior );
      }

    else if( qoi_name == vorticity )
      {
        qoi =  new Vorticity( vorticity );
      }

    else if( qoi_name == weighted_flux )
      {
        qoi =  new WeightedFluxQoI( weighted_flux );
      }

    else if( qoi_name == integrated_function )
      {
        std::string function;
        if (input.have_variable("QoI/IntegratedFunction/function"))
          function = input("QoI/IntegratedFunction/function", "");
        else
          libmesh_error_msg("ERROR: Could not find function to integrate");

        unsigned int p_level = input("QoI/IntegratedFunction/quadrature_level", 2);

        SharedPtr<libMesh::FunctionBase<libMesh::Real> > f = new libMesh::ParsedFunction<libMesh::Real>(function);

        qoi =  new IntegratedFunction<libMesh::FunctionBase<libMesh::Real> >(input,p_level,f,"IntegratedFunction",integrated_function);
      }

    else if ( qoi_name == spectroscopic_absorption )
      {
        std::string material;
        this->get_var_value<std::string>(input,material,"QoI/SpectroscopicAbsorption/material","NoMaterial!");

        std::string hitran_data;
        this->get_var_value<std::string>(input,hitran_data,"QoI/SpectroscopicAbsorption/hitran_data_file","");

        std::string hitran_partition;
        this->get_var_value<std::string>(input,hitran_partition,"QoI/SpectroscopicAbsorption/hitran_partition_function_file","");

        libMesh::Real T_min,T_max,T_step;
        std::string partition_temp_var = "QoI/SpectroscopicAbsorption/partition_temperatures";
        if (input.have_variable(partition_temp_var))
          {
            T_min = input(partition_temp_var, 0.0, 0);
            T_max = input(partition_temp_var, 0.0, 1);
            T_step = input(partition_temp_var, 0.0, 2);
          }
        else
          libmesh_error_msg("ERROR: Could not find tenmperature range specification for partition functions: "+partition_temp_var+" 'T_min T_max T_step'");

        SharedPtr<HITRAN> hitran( new HITRAN(hitran_data,hitran_partition,T_min,T_max,T_step) );

        std::string species;
        this->get_var_value<std::string>(input,species,"QoI/SpectroscopicAbsorption/species_of_interest","");

        SharedPtr<FEMFunctionAndDerivativeBase<libMesh::Real> > absorb;

        libMesh::Real thermo_pressure = -1.0;
        bool calc_thermo_pressure = input("QoI/SpectroscopicAbsorption/calc_thermo_pressure", false );
        if (!calc_thermo_pressure)
          thermo_pressure = input("Materials/"+material+"/ThermodynamicPressure/value", 0.0 );

        libMesh::Real nu_desired;
        this->get_var_value<libMesh::Real>(input,nu_desired,"QoI/SpectroscopicAbsorption/desired_wavenumber",0.0);

        libMesh::Real nu_min;
        this->get_var_value<libMesh::Real>(input,nu_min,"QoI/SpectroscopicAbsorption/min_wavenumber",0.0);

        libMesh::Real nu_max;
        this->get_var_value<libMesh::Real>(input,nu_max,"QoI/SpectroscopicAbsorption/max_wavenumber",0.0);

        ChemistryBuilder chem_builder;

#if GRINS_HAVE_ANTIOCH
        libMesh::UniquePtr<AntiochChemistry> chem_ptr;
        chem_builder.build_chemistry(input,material,chem_ptr);
        SharedPtr<AntiochChemistry> chem(chem_ptr.release());

        absorb = new AbsorptionCoeff<AntiochChemistry>(chem,hitran,nu_min,nu_max,nu_desired,species,thermo_pressure);
#elif GRINS_HAVE_CANTERA
        libMesh::UniquePtr<CanteraMixture> chem_ptr;
        chem_builder.build_chemistry(input,material,chem_ptr);
        SharedPtr<CanteraMixture> chem(chem_ptr.release());

        absorb = new AbsorptionCoeff<CanteraMixture>(chem,hitran,nu_min,nu_max,nu_desired,species,thermo_pressure);
#else
        libmesh_error_msg("ERROR: GRINS must be built with either Antioch or Cantera to use the SpectroscopicAbsorption QoI");
#endif

        qoi = new SpectroscopicAbsorption(input,qoi_name,absorb);
      }

    else
      {
        libMesh::err << "Error: Invalid QoI name " << qoi_name << std::endl;
        libmesh_error();
      }

    libmesh_assert(qoi);

    qois->add_qoi( *qoi );

    return;
  }

  void QoIFactory::check_qoi_physics_consistency( const GetPot& input,
                                                  const std::string& qoi_name )
  {
    int num_physics =  input.vector_variable_size("Physics/enabled_physics");

    // This should be checked other places, but let's be double sure.
    libmesh_assert(num_physics > 0);

    std::set<std::string> requested_physics;
    std::set<std::string> required_physics;

    // Build Physics name set
    for( int i = 0; i < num_physics; i++ )
      {
        requested_physics.insert( input("Physics/enabled_physics", "NULL", i ) );
      }

    /* If it's Nusselt, we'd better have HeatTransfer or LowMachNavierStokes.
       HeatTransfer implicitly requires fluids, so no need to check for those. `*/
    if( qoi_name == avg_nusselt )
      {
        required_physics.insert(PhysicsNaming::heat_transfer());
        required_physics.insert(PhysicsNaming::low_mach_navier_stokes());
        this->consistency_helper( requested_physics, required_physics, qoi_name );
      }

    return;
  }

  void QoIFactory::echo_qoi_list( SharedPtr<CompositeQoI>& qois )
  {
    /*! \todo Generalize to multiple QoI case when CompositeQoI is implemented in libMesh */
    std::cout << "==========================================================" << std::endl
              << "List of Enabled QoIs:" << std::endl;

    for( unsigned int q = 0; q < qois->n_qois(); q++ )
      {
        std::cout << qois->get_qoi(q).name() << std::endl;
      }

    std::cout <<  "==========================================================" << std::endl;

    return;
  }

  void QoIFactory::consistency_helper( const std::set<std::string>& requested_physics,
                                       const std::set<std::string>& required_physics,
                                       const std::string& qoi_name )
  {
    bool physics_found = false;
    for( std::set<std::string>::const_iterator name = required_physics.begin();
         name != required_physics.end();
         name++ )
      {
        if( requested_physics.find( (*name) ) != requested_physics.end() )
          physics_found = true;
      }

    if( !physics_found )
      this->consistency_error_msg( qoi_name, required_physics );

    return;
  }

  void QoIFactory::consistency_error_msg( const std::string& qoi_name,
                                          const std::set<std::string>& required_physics )
  {
    libMesh::err << "================================================================" << std::endl
                 << "QoI " << qoi_name << std::endl
                 << "requires one of the following physics which were not found:" <<std::endl;

    for( std::set<std::string>::const_iterator name = required_physics.begin();
         name != required_physics.end();
         name++ )
      {
        libMesh::err << *name << std::endl;
      }

    libMesh::err << "================================================================" << std::endl;

    libmesh_error();
  }

  template<typename T>
  void QoIFactory::get_var_value( const GetPot & input, T & value, std::string input_var, T default_value )
  {
    if (input.have_variable(input_var))
      value = input(input_var, default_value);
    else
      libmesh_error_msg("ERROR: Could not find required input parameter: "+input_var);
  }

  template void QoIFactory::get_var_value<std::string>( const GetPot &, std::string &, std::string, std::string);
  template void QoIFactory::get_var_value<libMesh::Real>( const GetPot &, libMesh::Real &, std::string, libMesh::Real);

} //namespace GRINS
