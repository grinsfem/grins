#ifndef GRINS_OD_PREMIXED_FLAME_H
#define GRINS_OD_PREMIXED_FLAME_H

// Grins items, can see not needing multi-component-vector variable but will leave in until i've written more of the code.
#include "grins_config.h"
#include "grins/grins_enums.h"
#include "grins/physics.h"
#include "grins/assembly_context.h"


#include "grins/single_variable.h"
#include "grins/materials_parsing.h"

namespace GRINS
{
  class ODPremixedFlame : public Physics
  {
  public:
    ODPremixedFlame(const PhysicsName& physics_name,
		    const GetPot & input,
		    std::unique_ptr<Mixture> & gas_mix);
    
    virtual ~ODPremixedFlame(){};




    //! Sets  variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    unsigned int n_species() const;

    libMesh::Real T( const libMesh::Point& p, const AssemblyContext& c ) const;

    libMesh::Real M_dot( const libMesh::Point& p, const AssemblyContext& c ) const;

    void mass_fractions( const libMesh::Point& p, const AssemblyContext& c,
                         std::vector<libMesh::Real>& mass_fracs ) const;

    libMesh::Real rho( libMesh::Real T, libMesh::Real p0, libMesh::Real R_mix) const;

    libMesh::Real get_p0() const;



    //

    
    //Register postprocessing variables for OD Premixed Flame 
    virtual void register_postprocessing_vars( const GetPot& input,
					       PostProcessedQuantities<libMesh::Real>& postprocessing);

    // Context Initializations
    virual void init_context( AssemblyContext& context );

    //Time dependent part(s)          // The Energy and species, and Mass equations will be included here, F(u,v)
    virtual void element_time_derivative( bool compute_jacobian,
					  AssemblyContext & context);
  

    //Mass matrix part(s)             //This would by my LHS of the species and energy, and mass residuals  M(u,v)
    virtual void mass_residual( bool compute_jacobian,
				AssemblyContext & context );

    virtual void element_constraint(bool compute_jacobian,
				    AssemblyContext & context );

    virtual void compute_postprocessed_quantity( unsigned int quantity_index, 
						 const AssemblyContext& context,
						 const libMesh::Point& point,
						 libMesh::Real& value );
    //Ditched the Caching of the time derivatives

   

    









  protected:
    //Variables 
    PrimitiveTempFEVariables& _temp_vars;
    SingleVariable& _mass_flux_vars;
    SpeciesMassFractionsVariable& _species_vars;


    //! Number of species
    unsigned int _n_species;
    

    bool _fixed_density;

    libMesh::Real _fixed_rho_value;






    //! Index from registering this quantity
    unsigned int _rho_index;               

    //! Index from registering this quantity                        //This is the coefficient of thermal conductivity, or lambda in my case
    unsigned int _k_index;

    //!Index from registering this quantity
    unsigned int _cp_index;
    
    //!Index from registering this quantity
    libMesh::Number _p0;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _h_s_index;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _omega_dot_index;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _Ds_index;

    /* //! Index from registering this quantity. Each species will have it's own index.
       std::vector<unsigned int> _cp_s_index;*/                

 
  private:

    void read_input_options( const GetPot& input );

    ODPremixedFlame();

  }; //Class ODPremixedFlame


  inline
    unsigned int ODPremixedFlame::n_species() const
    { return _n_species }


  inline
    libMesh::Real ODPremixedFlame::T( const libMesh::Point& p,
                                                        const AssemblyContext& c ) const
  { return c.point_value(_temp_vars.T(),p); }


  inline
    libMesh::Real M_dot( const libMesh::Point& p,
			 const AssemblyContext& c ) const
  { return c.point_value(_mass_flux_vars.M(),p); }

  
  inline
    void ODPremixedFlame::mass_fractions( const libMesh::Point& p,
						const AssemblyContext& c,
						std::vector<libMesh::Real>& mass_fracs ) const
  {
    libmesh_assert_equal_to(mass_fracs.size(), this->_n_species);
    
    for( unsigned int var = 0; var < this->_n_species; var++ )
      {
        mass_fracs[var] = c.point_value(_species_vars.species(var),p);
      }
  }

  
  inline
    libMesh::Real ODPremixedFlame::rho( libMesh::Real T,
					      libMesh::Real p0,
					      libMesh::Real R_mix) const
  {
    libMesh::Real value = 0;
    if( this->_fixed_density )
      value = this->_fixed_rho_value;
    else
      value = p0/(R_mix*T);

    return value;
  }


  inline
    libMesh::Real get_p0() const
  {return _p0;}



} // namespace Grins


#endif // OD_Premixed_Flame
