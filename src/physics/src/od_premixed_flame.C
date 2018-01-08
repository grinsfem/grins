// This class
#include "grins/od_premixed_flame.h"

#include "grins_config.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/cantera_mixture.h"
#include "grins/grins_enums.h"
#include "grins/antioch_mixture.h"
#include "grins/materials_parsing.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"
#include "grins/generic_ic_handler.h"
#include "grins/postprocessed_quantities.h"

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  template<class Mixture, class Evaluator>

  ODPremixedFlame<Mixture,Evaluator>::ODPremixedFlame(const std::string& physics_name,
								  const GetPot& input,
								  std::unique_ptr<Mixture> & gas_mix)
    : Physics(physics_name,input),
      _temp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>(VariablesParsing::temp_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _species_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>(VariablesParsing::species_mass_frac_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _mass_flux_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<SingleVariable>(VariableParsing::single_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _n_species(_species_vars.n_species()),
      _fixed_density( input("Physics/"+PhysicsNaming::od_premixed_flame+"/fixed_density", false ) ),
      _fixed_rho_value(0.0),
      _gas_mixture(gas_mix.release()),
      _rho_index(0),
      _k_index(0),           
      _cp_index(0)
  {
    this->set_parameter                    
      (_fixed_rho_value, input,
       "Physics/"+PhysicsNaming::od_premixed_flame+"fixed_rho_value", 0.0 );

    this->read_input_options(input);
                                                                      
    this->check_var_subdomain_consistency(_mass_flux_vars);
    this->check_var_subdomain_consistency(_temp_vars);
    this->check_var_subdomain_consistency(_species_vars);

   
    this->_ic_handler = new GenericICHandler( physics_name, input );
  }

  //Reading and setting the Thermodynamic Pressure
  void ODPremixedFlame::read_input_options( const GetPot& input )
  {
    // Read thermodynamic pressure info
    MaterialsParsing::read_property( input,
				     "ThermodynamicPressure",
                                     PhysicsNaming::od_premixed_flame(),                     
                                     (*this),
                                     _p0 );
  }

  
 
  void ODPremixedFlame::set_time_evolving_vars( libMesh::FEMSystem* system )
  {    
    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
        system->time_evolving( _species_vars.species(i), 1 );
      }
        system->time_evolving(_temp_vars.T(), 1);
  }


  void ODPremixedFlame::init_context( AssemblyContext& context )   
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_species_vars.species(0))->get_JxW();
    context.get_element_fe(_species_vars.species(0))->get_phi();
    context.get_element_fe(_species_vars.species(0))->get_dphi();
    context.get_element_fe(_species_vars.species(0))->get_xyz();

    context.get_element_fe(_mass_flux_vars.var())->get_JxW();
    context.get_element_fe(_mass_flux_vars.var())->get_phi();
    context.get_element_fe(_mass_flux_vars.var())->get_dphi();
    context.get_element_fe(_mass_flux_vars.var())->get_xyz();
    
    context.get_element_fe(_temp_vars.T())->get_JxW();
    context.get_element_fe(_temp_vars.T())->get_phi();
    context.get_element_fe(_temp_vars.T())->get_dphi();
    context.get_element_fe(_temp_vars.T())->get_xyz();


    // We also need the side shape functions, etc.

    context.get_side_fe(_mass_flux_vars.var())->get_JxW();
    context.get_side_fe(_mass_flux_vars.var())->get_phi();
    context.get_side_fe(_mass_flux_vars.var())->get_dphi();
    context.get_side_fe(_mass_flux_vars.var())->get_xyz();

    context.get_side_fe(this->_temp_vars.T())->get_JxW();
    context.get_side_fe(this->_temp_vars.T())->get_phi();
    context.get_side_fe(this->_temp_vars.T())->get_dphi();
    context.get_side_fe(this->_temp_vars.T())->get_xyz();
  }		   


 
  template<typename Mixture, typename Evaluator>
  void ODPremixedFlame<Mixture,Evaluator>::register_postprocessing_vars( const GetPot& input,
									       PostProcessedQuantities<libMesh::Real>& postprocessing )
  {
    std::string section = "physics/"+PhysicsNaming::od_premixed_flame()+"/output_vars";
    
    if( input.have_variable(section) )
      {
	unsigned int n_vars = input.vector_variable_size(section);

	for( unsigned int v=0; v <n_vars; v++ )
	  {
	    std::string name = input(section,"DIE!",v);
	    
	    if( name == std::string("rho") )
	      {
		this->_rho_index = postprocessing.register_quantity( name );
	      }
	    else if( name == std::string("k") ) 
	      {
		this->_k_index = postprocessing.register_quantity( name );
	      }
	    else if( name == std::string("cp") )
	      {
		this->_cp_index = post.processing.register_quantity( name );
	      }
	    else if( name == (std::string("u")) )
	      {
		this->_u_index = postprocessing.register_quantity( name );
	      }
		
	    //time for species specific values

	      else if(name == std::string("mole_fractions") )
	      {
		this ->_mole_fractions_index.resize(this->n_species());
	    
		for(unsigned int s=0; s < this->n_species(); s++)
		  {
		    this->mole_fractions_index[s] = postprocessing.register_quantity("X_"+this->_gas_mixture->species_name(s) );
		  }
	      }
	    else if(name == std::string("h_s") )
	      {
		this->_h_s_index.resize(this->n_species());
		
		for(unsigned int s=0; s < this->n_species(); s++)
		  {
		    this->_h_s_index[s] = postprocessing.register_quantity( "h_"+this->_gas_mixture->species_name(s) );
		  }
	      }
	    else if(name == std::string("omega_dot") )
	      {
		this->_omega_dot_index.resize(this->n_species());
		
		for(unsigned int s=0; s < this->n_species(); s++)
		  {
		    this->_omega_dot_index[s] = postprocessing.register_quantity( "omega_dot_"+this->_gas_mixture->species_name(s) );
		  }
	      }
	    else if( name == std::string("D_s") )
	      {
		this->_Ds_index.resize(this->n_species());
		
		for(unsigned int s=0;s < this->n_species(); s++)
		  {
		    this->_Ds_index[s] = postprocessing.register_quantity( "D_"+this->_gas_mixture->species_name(s) );
		  }
	      }
	      /* else if( name == std::string("cp_s") )
	      {
		this->_cp_s_index.resize(this->n_species());
		
		for(unsigned int s=0;s < this->n_species(); s++)
		  {
		    this ->_cps_index[s] = postprocessing.register_quantity( "cp_"+this->_gas_mixture->species_name(s) );
		  }
		  }*/
	    else
	      {
		std::cerr << "Error: Invalid output_vars value for ODPremixedFlame " << std::endl
                          << "       Found " << name << std::endl
                          << "       Acceptable values are: rho" << std::endl
                          << "                              D_s" << std::endl
                          << "                              k" << std::endl
			  << "                              h_s" << std::endl
                          << "                              cp" << std::endl
                          << "                              mole_fractions" << std::endl
                          << "                              omega_dot" << std::endl
			  << "                              u" << std::endl;
                libmesh_error();
	      }
	  }
      }
    
    return;
  }    //end post processing vars
  


  template<typename Mixture, typename Evaluator>
  void ODPremixedFlame<Mixture,Evaluator>::element_time_derivative
( bool compute_jacobian, AssemblyContext &  context )
  {
    if( compute_jacobian )
      libmesh_not_implemented();
    
    //Convenience
    const VariableIndex s0_var = this->_species_vars.species(0);
    
    // The number of local degrees of freedom in each variable                                  //Variables are Species, Mass flux, and Temperature
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();
    const unsigned int n_M_dofs = context.get_dof_indices(this->_mass_flux_vars.var()).size();
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();
    
    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real>& JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();
    
    //the species shape function at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> > & s_phi = context.get_element_fe(s0_var)->get_phi();
    
    //the species shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::Gradient> > & s_dphi = context.get_element_fe(s0_var)->get_dphi();
    
    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();
    
    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_dphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();    
    
    
    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T()); // R_{T}
    
    //species set in residual calculations since the need of a for loop
    
    unsigned int n_qpoints = context.get_element_qrule().n_points();
    for(unsigned int qp = 0;qp != n_qpoints;qp++)
      {
	libMesh::Real T, M_dot;
	libMesh::Gradient Dx_T;
       
	libmesh::Real R, M, k, cp, rho, p0, mu;
	
	T = context.interior_value(this->_temp_vars_T(),qp);
	Dx_T = context.interior_gradient(this->_temp_vars.T(), qp);
	
	M_dot = context.interior_value(this->_mass_flux_vars.var(), qp);
	
	p0 = this->get_p0();
	
	
	std::vector<libMesh::Real> mass_fractions, h, D, omega_dot;
	std::vector<libMesh::Gradient> Dx_mass_fractions;
	
	mass_fractions.resize(this->_n_species);
	Dx_mass_fractions.resize(this->_n_species);
	h.resize(this->_n_species);
	
	for (unsigned int s = 0; s < this->_n_species; s++)
	  {
	    mass_fractions[s] = std::max( context.interior_value(this->_species_vars.species(s),qp),0.0);
	    Dx_mass_fractions[s] = context.interior_gradient(this->_species_vars.species(s),qp);
	    h[s] = gas_evaluator.h_s( T, s );
	  }
	M = gas_evaluator.M_mix( mass_fractions );
	
	R = gas_evaluator.R_mix( mass_fractions );
	
	rho = this->rho( T, p0, R );
	
	cp = gas_evaluator.cp( T, p0, mass_fractions);
	
	D.resize(this->_n_species);
	
	gas_evaluator.mu_and_k_and_D( T, rho, cp, mass_fractions, 
				      mu, k, D );
	
	omega_dot.resize(this->_n_species);
	
	gas_evaluator.omega_dot( T, rho, mass_fractions, omega_dot );
	
	libMesh::Real jac = JxW[qp];
	
	//Energy equation Residual
	libMesh::Real chem_term = 0.0;
	for ( unsigned int s=0; s< this->_n_species; s++ )
	  {
	    chem_term +=h[s]*(this->_gas_mixture->M(s))*omega_dot[s];
	  }
	for (unsigned int i=0;i != n_T_dofs; i++ )
	  {
	    FT(i) += ( ( -cp*M_dot *Dx_T - chem_term )*T_phi[i][qp]
		       -k*Dx_T*T_dphi[i][qp])*jac;
	    
	  }
	
	//Species equation Residuals
	for (unsigned int s = 0; s < _n_species; s++)
	  {
	    libMesh::DenseSubVector<libMesh::Number> &FS =
	      contest.get_elem_residual(this->_species_vars.species(s)); //R_{s}
	    
	    const libMesh::Real term1 = -M_dot*Dx_mass_fractions[s] + omega_dot[s]*(this->_gas_mixture->M(s));
	    const libMesh::gradient term2 = -rho*D[s]*Dx_mass_fractions[s];
	    
	    for (unsigned int i =0;i != n_s_dofs;i++)
	      {
		FS(i) += ( term1 * s_phi[i][qp] + term2 * s_dphi[i][qp] )*jac;
	      }
	  }
      } // end of quadrature loop
    return;
  }            // end element time derivative



  template<typename Mixture, typename Evaluator>
  void ODPremixedFlame<Mixture,Evaluator>::mass_residual
  (bool compute_jacobian, AssemblyContext & context )
  {
    const VariableIndex s0_var = this->_species_vars.species(0);
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();
    const unsigned int n_M_dofs = context.get_dof_indices(this->_mass_flux_vars.var()).size();
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();
  
    
    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real>& JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();
    
    //the species shape function at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> > & s_phi = context.get_element_fe(s0_var)->get_phi();
 
    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();
 
    //The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F_T = context.get_elem_residual(this->_temp_vars.T());

    libMesh::DenseSubVector<libMesh::Number> &F_M = context.get_elem_residual(this->_mass_flux_vars.var());

    //Get the number of quadrature points
    unsigned int n_qpoints = context.get_element_qrule().n_points();
    
    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
	libMesh::Real T_dot;
	context.interior_rate(this->_temp_vars.T(), qp, T_dot);
	
        libMesh::Real T = context.interior_value(this->_temp_vars.T(), qp);

	std::vector<libMesh::Real> mass_fractions(this->n_species());
        for(unsigned int s=0; s < this->_n_species; s++ )
          mass_fractions[s] = context.interior_value(this->_species_vars.species(s), qp);

	Evaluator gas_evaluator(*(this-> _gas_mixture));
	const libMesh::Real R_mix = gas_evaluator.R_mix(mass_fractions);
        const libMesh::Real p0 = this->get_p0();
        const libMesh::Real rho = this->rho(T, p0, R_mix);
        const libMesh::Real cp = gas_evaluator.cp(T,p0,mass_fractions);
        const libMesh::Real M = gas_evaluator.M_mix(mass_fractions);

	libMesh::Real jac = JxW[qp];
	
	// Species residual
	for(unsigned int s=0; s < this->n_species(); s++)
	  {
	    libMesh::DenseSubVector<libMesh::Number> & F_s =
	      context.get_elem_residual(this->_species_vars.species(s));

	    libMesh::Real mass_fractions_dot;
	    context.interior_rate(this->_species_vars.species(s));

	    for (unsigned int i = 0; i != n_s_dofs; ++i)
	      {
	      F_s(i) -= rho*mass_fractions_dot*s_phi[i][qp]*jac;
	      }
	  }
	
	//Energy Residual
	
	for (unsigned int i = 0; i!= n_T_dofs)
	  {
	    F_T(i) = rho*cp*T_dot*T_phi[i][qp]*jac;
	  }
	
	if( compute_jacobian ) 
	  libmesh_not_implemented();


      } // end Quadrature loop
  }                   //end Mass Residual



  template<typename Mixture, typename Evaluator>
    void ODPremixedFlame<Mixture,Evaluator>::element_constraint
    ( bool compute_jacobian,
      AssemblyContext & context )
    {
      if(compute_jacobian)
	libmesh_not_implemented();
      
      const std::vector<libMesh::Real> &JxW =
	context.get_element_fe(this->_temp_vars.T())->get_JxW();
      
      // The Mass Flux shape functions at interior quadrature points.
      const std::vector<std::vector<libMesh::Real> >& M_phi =
	context.get_element_fe(this->_Mass_flux_vars.var())->get_phi();

    

      const unsigned int n_M_dofs = context.get_dof_indices(this->_mass_flux_vars.var()).size();

      libMesh::DenseSubVector<libMesh::Number> &Fm = context.get_elem_residual(this->_Mass_flux_vars.var()); // R_{M}

      for(unsigned int qp=0; qp!=n_qpoints; qp++)
	{
	  libmesh::Gradient Dx_M_dot = context.interior_gradient(this->_mass_flux_vars.var(), qp);
	  libMesh::Real jac = JxW[qp];
	  for(unsigned int i=0;i != n_M_dofs; i++)
	    {
	      Fm(i) += Dx_M_dot*M_phi[i][qp]*jac;
	    }
	}
    }   //end Element Constraint

      





  template<typename Mixture, typename Evaluator>
  void ODPremixedFlame<Mixture,Evaluator>::compute_postprocessed_quantity( unsigned int quantity_index, 
										       const Assembly Context& context,
										       const libMesh::Point& point,
										       libMesh::Real & value )
  {
    Evaluator gas_evaluator( *(this->_gas_mixture) );

    if( quantity_index == this->_rho_index )
      {
	std::vector<libMesh::Real> Y( this->_n_species );
	libMesh Real T = this->T(point,context);
	libMesh::Real p0 = this->get_p0();
	this->mass_fractions( point, context, Y );
	
	value = this->rho(T,p0, gas_evaluator.R_mix(Y) );
      }
    else if( quantity_index == this->_k_index )
      {
	std::vector<libMesh::Real> Y(this->_n_species );
	
	libMesh::Real T = this->T(point,context);
	this->mass_fractions( point,context, Y);
	libMesh::Real p0 = this->get_p0();
	
	libMesh::Real cp = gas_evaluator.cp( T, p0, Y );
	
	libMesh::Real rho = this->rho(T, p0, gas_evaluator.R_mix(Y) );

	libMesh::Real mu,k;
	
	gas_evaluator.mu_and_k_and_D( T, rho, cp, Y, mu, k, D );

	value = k;
	return;
      }
    else if( quantity_index == this->_cp_index )
      {
	std::vector<libMesh::Real> Y( this->_n_species );
	libMesh::Real T = this->T(point,context);
	this->mass_fractions( point, context, Y);
	libMesh::Real p0 = this0>get_p0();
	
	value = gas_evaluator.cp( T, p0, Y );
      }
    else if ( quantity_index == this->_u_index )
      {
	libMesh::Real M_dot = this->M_dot(point,context);

	std::vector<libMesh::Real> Y( this->_n_species );
	libMesh Real T = this->T(point,context);
	libMesh::Real p0 = this->get_p0();
	this->mass_fractions( point, context, Y );
	rho = this->rho(T,p0,Y );

	value = M_dot/rho;
	
      }
    //now onto the species dependent stuff
    
    else 
      {
	if( !this->_mole_fractions_index.empty() )
	  {
	    libmesh_assert_equal_to( _mole_fractions_index.size(), this->n_species() );
	    
	    for( unsigned int s = 0; s < this->n_species(); s++ )
	      {
		if( quantity_index == this->_mole_fractions_index[s] )
		  {
		    std::vector<libMesh::Real> Y( this->_n_species );
		    this->mass_fractions( point, context, Y );
		    
		    libMesh::Real M = gas_evaluator.M_mix(Y);
		    
		    value = gas_evaluator.X( s, M, Y[s]);
		    return;
		  }
	      }
	  }
	
	if( !this->_h_s_index.empty() )
	  {
	    libmesh_assert_equal_to( _h_s_index.size(), this->n_species() );
	    
	    for( unsigned int s=0; s < this->n_species(); s++)
	      {
		if( quantity_index == this->_h_s_index[s] )
		  {
		    libMesh::Real T = this->T(point,context);

		    value = gas_evaluator.h_s( T, s );
		    return;
		  }
	      }
	  }
	
	if( !this->_omega_dot_index.empty() )
	  {
	    libmesh_assert_equal_to( _omega_dot_index.size(), this->n_species() );

	    for(unsigned int s=0; s < this->n_species(); s++)
	      {
		if( quantity_index == this->_omega_dot_index[s] )
		  {
		    std::vector<libMesh::Real> Y( this->n_species() );
		    this->mass_fractions( point, context, Y );
		    
		    libMesh::Real T = this->T(point,context);
		    
		    libMesh::Real p0 = this->get_p0();
		    
		    libMesh::Real rho = this->rho(T,p0, gas_evaluator.R_mix(Y) );
		    
		    std::vector<libMesh::Real> omega_dot( this->n_species() );
		    gas_evaluator.omega_dot( T, rho, Y, omega_dot );
		    
		    value = omega_dot[s]
		      return;
		  }
	      }
	  }

	if( !this->_Ds_index.empty() )
	  {
	    libmesh_assert_equal_to( _Ds_index.size(), this->n_species() );

	    for( unsigned int s = 0; s < this->n_species(); s++ )
	      {
		if(qunatity_index == this->_Ds_index[s] )
		  {
		    std::vector<libMesh::Real> Y( this->_n_species );
		    
		    libMesh::Real T = this->T(point,context);
		    this->mass_fractions(point,context, Y );
		    libMesh::Real p0 = this->get_p0();

                    libMesh::Real cp = gas_evaluator.cp( T, p0, Y );

                    libMesh::Real rho = this->rho( T, p0, gas_evaluator.R_mix(Y) );

                    libMesh::Real mu, k;
                    std::vector<libMesh::Real> D( this->_n_species );

                    gas_evaluator.mu_and_k_and_D( T, rho, cp, Y, mu, k, D );

                    value = D[s];
                    return;
		  }
	      }
	  }
      }//if/else quantity_index
    
    return;
  }         //end Postproc


}//end Namespace Grins
	       
	  
	
    
	 
     
     

     
