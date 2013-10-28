//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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

#ifndef GRINS_REACTING_LOW_MACH_NAVIER_STOKES_BC_HANDLING_H
#define GRINS_REACTING_LOW_MACH_NAVIER_STOKES_BC_HANDLING_H

// shared_ptr
#include "boost/tr1/memory.hpp"

// GRINS
#include "grins/low_mach_navier_stokes_bc_handling.h"
#include "grins/catalytic_wall.h"

namespace GRINS
{
  template<typename Chemistry>
  class ReactingLowMachNavierStokesBCHandling : public LowMachNavierStokesBCHandling
  {
  public:

    ReactingLowMachNavierStokesBCHandling( const std::string& physics_name, const GetPot& input,
					   const Chemistry& chem_mixture );

    virtual ~ReactingLowMachNavierStokesBCHandling();

    virtual int string_to_int( const std::string& bc_type_in ) const;

    virtual void init_bc_data( const libMesh::FEMSystem& system );
    
    virtual void init_bc_types( const GRINS::BoundaryID bc_id, 
			        const std::string& bc_id_string, 
			        const int bc_type, 
				const std::string& bc_vars, 
				const std::string& bc_value, 
			        const GetPot& input );

    virtual void user_init_dirichlet_bcs( libMesh::FEMSystem* system,
					  libMesh::DofMap& dof_map,
					  GRINS::BoundaryID bc_id,
					  GRINS::BCType bc_type ) const;

    virtual void init_dirichlet_bcs( libMesh::FEMSystem* system ) const;

    virtual void user_apply_neumann_bcs( AssemblyContext& context,
					 const GRINS::CachedValues& cache,
					 const bool request_jacobian,
					 const GRINS::BoundaryID bc_id,
					 const GRINS::BCType bc_type ) const;

    void set_species_bc_type( GRINS::BoundaryID bc_id, int bc_type );

    void set_species_bc_values( GRINS::BoundaryID bc_id, const std::vector<libMesh::Real>& species_values );

    libMesh::Real get_species_bc_value( GRINS::BoundaryID bc_id, unsigned int species ) const;

  protected:

    void build_catalycities( const GetPot& input,
                             const std::string& reactant,
                             const std::string& bc_id_string,
                             const BoundaryID bc_id,
                             boost::scoped_ptr<CatalycityBase>& gamma_r,
                             boost::scoped_ptr<CatalycityBase>& gamma_p );

     // We need a another container to stash dirichlet values for the speccies
    std::map< GRINS::BoundaryID, std::vector<libMesh::Real> > _species_bc_values;

    // We also need another map container
    std::vector<std::pair<BoundaryID,BCType> > _species_bc_map;

    unsigned int _n_species;
    std::vector<std::string> _species_var_names;
    std::vector<GRINS::VariableIndex> _species_vars;

    //! Temporarily cache the pointers to the CatalyticWall functors
    /*! We need to be able to init them before they get inserted into the BCHandling
        system, so we cached the pointers here and then clean up when we're done with
        them so they're only stored in one place. */
    std::map<BoundaryID,std::vector<std::pair<unsigned int,std::tr1::shared_ptr<CatalyticWall<Chemistry> > > > > _catalytic_walls;

    //! Cache which species are catalytic
    /*! \todo Currently restricted to one reaction per species */
    std::map<BoundaryID,std::set<unsigned int> > _catalytic_species;

    const Chemistry& _chemistry;

  private:

    ReactingLowMachNavierStokesBCHandling();

    // Needs to start larger than the LMNS_BC_TYPES end
    enum RLMNS_BC_TYPES{ ZERO_SPECIES_FLUX=20, 
			 PRESCRIBED_SPECIES,
                         PRESCRIBED_MOLE_FRACTIONS,
			 CATALYTIC_WALL,
			 GENERAL_SPECIES };

  };
}

#endif // GRINS_REACTING_LOW_MACH_NAVIER_STOKES_BC_HANDLING_H
