//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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
#include "grins/spectroscopic_transmission.h"
#include "grins/spectroscopic_absorption.h"
#include "grins/chemistry_builder.h"
#include "grins/constant_laser_intensity_profile.h"
#include "grins/collimated_gaussian_laser_intensity_profile.h"
#include "grins/laser_absorption.h"

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

  std::shared_ptr<CompositeQoI> QoIFactory::build(const GetPot& input)
  {
    std::string qoi_list = input("QoI/enabled_qois", "none" );

    std::vector<std::string> qoi_names;

    if( qoi_list != std::string("none") )
      {
        StringUtilities::split_string( qoi_list, std::string(" "), qoi_names );
      }

    std::shared_ptr<CompositeQoI> qois( new CompositeQoI );

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

  void QoIFactory::add_qoi( const GetPot& input, const std::string& qoi_name, std::shared_ptr<CompositeQoI>& qois )
  {
    QoIBase* qoi = NULL;

    bool do_final_add = true;

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

        std::shared_ptr<libMesh::FunctionBase<libMesh::Real> >
          f( new libMesh::ParsedFunction<libMesh::Real>(function) );

        std::shared_ptr<RayfireMesh> rayfire( new RayfireMesh(input,"IntegratedFunction") );

        qoi =  new IntegratedFunction<libMesh::FunctionBase<libMesh::Real> >(p_level,f,rayfire,qoi_name);
      }

    else if ( qoi_name == spectroscopic_transmission )
      {
        do_final_add = this->create_spectroscopic_qoi(input,qoi_name,"SpectroscopicTransmission",&qoi,qois);
      }

    else if ( qoi_name == spectroscopic_absorption )
      {
        do_final_add = this->create_spectroscopic_qoi(input,qoi_name,"SpectroscopicAbsorption",&qoi,qois);
      }

    else if ( qoi_name == laser_absorption )
      {
        std::string qoi_string = "LaserAbsorption";

        std::shared_ptr<FEMFunctionAndDerivativeBase<libMesh::Real>> absorb = this->create_absorption_coeff(input,qoi_string);

        unsigned int dim = 2;
        if (input.have_variable("QoI/"+qoi_string+"/Rayfire/phi"))
          dim = 3;

        libMesh::Point top_origin;
        top_origin(0) = input("QoI/"+qoi_string+"/top_origin", 0.0, 0);
        top_origin(1) = input("QoI/"+qoi_string+"/top_origin", 0.0, 1);
        if (dim == 3)
          top_origin(2) = input("QoI/"+qoi_string+"/top_origin", 0.0, 2);

        libMesh::Point centerline_origin;
        centerline_origin(0) = input("QoI/"+qoi_string+"/centerline_origin", 0.0, 0);
        centerline_origin(1) = input("QoI/"+qoi_string+"/centerline_origin", 0.0, 1);
        if (dim == 3)
          centerline_origin(2) = input("QoI/"+qoi_string+"/centerline_origin", 0.0, 2);

        libMesh::Point bottom_origin;
        bottom_origin(0) = input("QoI/"+qoi_string+"/bottom_origin", 0.0, 0);
        bottom_origin(1) = input("QoI/"+qoi_string+"/bottom_origin", 0.0, 1);
        if (dim == 3)
          bottom_origin(2) = input("QoI/"+qoi_string+"/bottom_origin", 0.0, 2);

        libMesh::Real theta = input("QoI/"+qoi_string+"/Rayfire/theta", -7.0);
        libMesh::Real phi = -1.0;
        if (dim == 3)
          this->get_var_value<libMesh::Real>(input,phi,"QoI/"+qoi_string+"/Rayfire/phi",-1.0);

        unsigned int n_qp = input("QoI/"+qoi_string+"/n_quadrature_points", 1); // default to single rayfire
        
        std::shared_ptr<LaserIntensityProfileBase> intensity_profile;
        std::string profile_name;
        this->get_var_value<std::string>(input,profile_name,"QoI/"+qoi_string+"/intensity_profile","");
        if (profile_name == "constant")
          {
            libMesh::Real I0;
            this->get_var_value<libMesh::Real>(input,I0,"QoI/"+qoi_string+"/I0",0.0);
            intensity_profile.reset( new ConstantLaserIntensityProfile(I0) );
          }
        else if (profile_name == "collimated gaussian")
          {
            libMesh::Real w;
            this->get_var_value<libMesh::Real>(input,w,"QoI/"+qoi_string+"/w",0.0);
            intensity_profile.reset( new CollimatedGaussianLaserIntensityProfile(w) );
          }
        else
          {
            libmesh_error_msg("ERROR: Please specify either 'constant' or 'collimated gaussian' intensity profile");
          }

        if (dim == 3)
          qoi = new LaserAbsorption(absorb,top_origin,centerline_origin,bottom_origin,theta,phi,n_qp,intensity_profile,qoi_name);
        else
          qoi = new LaserAbsorption(absorb,top_origin,centerline_origin,bottom_origin,theta,n_qp,intensity_profile,qoi_name);
      }

    else
      {
        libMesh::err << "Error: Invalid QoI name " << qoi_name << std::endl;
        libmesh_error();
      }

    libmesh_assert(qoi);

    if (do_final_add)
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

  void QoIFactory::echo_qoi_list( std::shared_ptr<CompositeQoI>& qois )
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

  bool QoIFactory::create_spectroscopic_qoi(const GetPot & input, const std::string & qoi_name, const std::string & qoi_string, QoIBase ** qoi, std::shared_ptr<CompositeQoI> & qois)
  {
    bool do_final_add = true;

    std::shared_ptr<FEMFunctionAndDerivativeBase<libMesh::Real>> absorb = this->create_absorption_coeff(input,qoi_string);

    std::shared_ptr<RayfireMesh> rayfire( new RayfireMesh(input,qoi_string) );

    // This is the wavenumber of the "laser" or a range of wavenumbers to scan over
    std::string desired_wavenumber_var = "QoI/"+qoi_string+"/desired_wavenumber";
    unsigned int num_wavenumbers = input.vector_variable_size(desired_wavenumber_var);
    if (num_wavenumbers == 1)
      {
        // Calculating QoI at a single wavenumber

        // Output the QoI value normally
        bool output_as_csv = false;
        
        if (qoi_name == spectroscopic_absorption)
          *qoi = new SpectroscopicAbsorption(absorb,rayfire,qoi_name,output_as_csv);
        else // spectroscopic_transmission
          *qoi = new SpectroscopicTransmission(absorb,rayfire,qoi_name,output_as_csv);

      }
    else if (num_wavenumbers == 3)
      {
        // Calculating QoI within a range of wavenumbers (i.e. absorption plot)

        // Output wavenumber,absorption (ideally to a file) for easy plotting
        bool output_as_csv = true;

        // We don't want to duplicate the last qoi
        do_final_add = false;
        
        libMesh::Real min_nu = input(desired_wavenumber_var,0.0,0);
        libMesh::Real max_nu = input(desired_wavenumber_var,0.0,1);
        libMesh::Real step_nu = input(desired_wavenumber_var,0.0,2);
        
        // sanity check the given range
        if ( (min_nu>max_nu) || (step_nu>(max_nu-min_nu)) )
          {
            std::stringstream ss;
            ss <<"Invalid specification of desired_wavenumber range:" <<std::endl
               <<"nu_min: " <<min_nu <<std::endl
               <<"nu_max: " <<max_nu <<std::endl
               <<"nu_step: " <<step_nu <<std::endl;
            libmesh_error_msg(ss.str());
          }
        
        for (libMesh::Real nu = min_nu; nu <= max_nu; nu += step_nu)
          {            
            libMesh::Real nu_desired = nu;

            AbsorptionCoeffBase * coeff = libMesh::cast_ptr<AbsorptionCoeffBase *>(absorb->clone().release());
            coeff->set_wavenumber(nu_desired);
            absorb.reset( coeff );

            if (qoi_name == spectroscopic_absorption)
              *qoi = new SpectroscopicAbsorption(absorb,rayfire,qoi_name,output_as_csv);
            else // spectroscopic_transmission
              *qoi = new SpectroscopicTransmission(absorb,rayfire,qoi_name,output_as_csv);

            qois->add_qoi( **qoi );
          }
      }

    return do_final_add;
  }

  std::shared_ptr<FEMFunctionAndDerivativeBase<libMesh::Real>> QoIFactory::create_absorption_coeff(const GetPot & input, const std::string & qoi_string)
  {
    std::string material;
    this->get_var_value<std::string>(input,material,"QoI/"+qoi_string+"/material","NoMaterial!");

    std::string hitran_data;
    this->get_var_value<std::string>(input,hitran_data,"QoI/"+qoi_string+"/hitran_data_file","");

    std::string hitran_partition;
    this->get_var_value<std::string>(input,hitran_partition,"QoI/"+qoi_string+"/hitran_partition_function_file","");

    libMesh::Real T_min,T_max,T_step;
    std::string partition_temp_var = "QoI/"+qoi_string+"/partition_temperatures";
    if (input.have_variable(partition_temp_var))
      {
        T_min = input(partition_temp_var, 0.0, 0);
        T_max = input(partition_temp_var, 0.0, 1);
        T_step = input(partition_temp_var, 0.0, 2);
      }
    else
      libmesh_error_msg("ERROR: Could not find temperature range specification for partition functions: "+partition_temp_var+" 'T_min T_max T_step'");

    std::shared_ptr<HITRAN> hitran( new HITRAN(hitran_data,hitran_partition,T_min,T_max,T_step) );

    std::string species;
    this->get_var_value<std::string>(input,species,"QoI/"+qoi_string+"/species_of_interest","");

    libMesh::Real thermo_pressure = -1.0;
    bool calc_thermo_pressure = input("QoI/"+qoi_string+"/calc_thermo_pressure", false );
    if (!calc_thermo_pressure)
      thermo_pressure = input("Materials/"+material+"/ThermodynamicPressure/value", 0.0 );


    // These options are for the linecenter wavenumbers included in the calculation of kv
    libMesh::Real nu_data_min,nu_data_max;
    this->get_var_value<libMesh::Real>(input,nu_data_min,"QoI/"+qoi_string+"/min_wavenumber",0.0);
    this->get_var_value<libMesh::Real>(input,nu_data_max,"QoI/"+qoi_string+"/max_wavenumber",0.0);

    ChemistryBuilder chem_builder;
    std::shared_ptr<FEMFunctionAndDerivativeBase<libMesh::Real> > absorb;

#if GRINS_HAVE_ANTIOCH
    std::unique_ptr<AntiochChemistry> chem_ptr;
    chem_builder.build_chemistry(input,material,chem_ptr);
    std::shared_ptr<AntiochChemistry> chem(chem_ptr.release());
#elif GRINS_HAVE_CANTERA
    std::unique_ptr<CanteraMixture> chem_ptr;
    chem_builder.build_chemistry(input,material,chem_ptr);
    std::shared_ptr<CanteraMixture> chem(chem_ptr.release());
#else
    libmesh_error_msg("ERROR: GRINS must be built with either Antioch or Cantera to use the LaserAbsorption QoI");
#endif

    libMesh::Real nu_desired;
    this->get_var_value<libMesh::Real>(input,nu_desired,"QoI/"+qoi_string+"/desired_wavenumber",0.0);

#if GRINS_HAVE_ANTIOCH
    absorb.reset( new AbsorptionCoeff<AntiochChemistry>(chem,hitran,nu_data_min,nu_data_max,nu_desired,species,thermo_pressure) );
#elif GRINS_HAVE_CANTERA
    absorb.reset( new AbsorptionCoeff<CanteraMixture>(chem,hitran,nu_data_min,nu_data_max,nu_desired,species,thermo_pressure) );
#else
    libmesh_error_msg("ERROR: GRINS must be built with either Antioch or Cantera to use the LaserAbsorption QoI");
#endif
    
    return absorb;
  }

  template void QoIFactory::get_var_value<std::string>( const GetPot &, std::string &, std::string, std::string);
  template void QoIFactory::get_var_value<libMesh::Real>( const GetPot &, libMesh::Real &, std::string, libMesh::Real);

} //namespace GRINS
