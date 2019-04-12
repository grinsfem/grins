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

// This class
#include "grins/gas_solid_catalytic_wall.h"

// GRINS
#include "grins/physics.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{
  template<typename Chemistry>
  GasSolidCatalyticWall<Chemistry>::GasSolidCatalyticWall(std::shared_ptr<Chemistry>& chem,
                                                          std::shared_ptr<CatalycityBase>& gamma,
                                                          const std::vector<VariableIndex>& species_vars,
                                                          VariableIndex T_var,
                                                          libMesh::Real p0,
                                                          unsigned int reactant_gas_species_idx,
                                                          unsigned int reactant_solid_species_idx,
                                                          unsigned int product_species_idx)
  : CatalyticWallBase<Chemistry>(chem,gamma,species_vars,T_var,p0,reactant_gas_species_idx),
    _reactant_gas_species_idx(reactant_gas_species_idx),
    _reactant_gas_var_idx(species_vars[reactant_gas_species_idx]),
    _reactant_solid_species_idx(reactant_solid_species_idx),
    _product_species_idx(product_species_idx),
    _product_var_idx(species_vars[product_species_idx])
  {}

  template<typename Chemistry>
  GasSolidCatalyticWall<Chemistry>::GasSolidCatalyticWall( const Chemistry& chem_mixture,
                                                           CatalycityBase& gamma,
                                                           const unsigned int reactant_gas_species_idx,
                                                           const unsigned int reactant_solid_species_idx,
                                                           const unsigned int product_species_idx )
    : CatalyticWallBase<Chemistry>(chem_mixture,gamma,reactant_gas_species_idx),
    _reactant_gas_species_idx(reactant_gas_species_idx),
    _reactant_solid_species_idx(reactant_solid_species_idx),
    _product_species_idx(product_species_idx)
  {
    libmesh_deprecated();
  }

  template<typename Chemistry>
  void GasSolidCatalyticWall<Chemistry>::init( const libMesh::FEMSystem& system )
  {
    libmesh_do_once(libmesh_deprecated());

    const std::string r_var_name = std::string("w_"+this->_chemistry.species_name( this->_reactant_gas_species_idx ) );

    const std::string p_var_name = std::string("w_"+this->_chemistry.species_name( this->_product_species_idx ) );

    libmesh_assert( system.has_variable( r_var_name ) );
    libmesh_assert( system.has_variable( p_var_name ) );

    this->_reactant_gas_var_idx = system.variable_number( r_var_name );

    this->_product_var_idx = system.variable_number( p_var_name );

    return;
  }

  template<typename Chemistry>
  bool GasSolidCatalyticWall<Chemistry>::eval_flux( bool compute_jacobian,
                                                    AssemblyContext& context,
                                                    libMesh::Real sign,
                                                    bool is_axisymmetric )
  {
    libMesh::FEGenericBase<libMesh::Real>* side_fe = NULL;
    context.get_side_fe( _reactant_gas_var_idx, side_fe );

    // The number of local degrees of freedom in each variable.
    const unsigned int n_var_dofs = context.get_dof_indices(_reactant_gas_var_idx).size();

    libmesh_assert_equal_to( n_var_dofs, context.get_dof_indices(_product_var_idx).size() );

    // Element Jacobian * quadrature weight for side integration.
    const std::vector<libMesh::Real> &JxW_side = side_fe->get_JxW();

    // The var shape functions at side quadrature points.
    const std::vector<std::vector<libMesh::Real> >& var_phi_side = side_fe->get_phi();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& var_qpoint = side_fe->get_xyz();

    // reactant residual
    libMesh::DenseSubVector<libMesh::Number> &F_r_var = context.get_elem_residual(_reactant_gas_var_idx);

    // product residual
    libMesh::DenseSubVector<libMesh::Number> &F_p_var = context.get_elem_residual(_product_var_idx);

    unsigned int n_qpoints = context.get_side_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real jac = JxW_side[qp];

        if(is_axisymmetric)
          {
            const libMesh::Number r = var_qpoint[qp](0);
            jac *= r;
          }

        std::vector<libMesh::Real> mass_fractions(this->_chem_ptr->n_species());
        for( unsigned int s = 0; s < this->_chem_ptr->n_species(); s++ )
          mass_fractions[s] = context.side_value(this->_species_vars[s], qp);

        libMesh::Real T =  context.side_value(this->_T_var, qp);
        libMesh::Real R_mix = this->_chem_ptr->R_mix(mass_fractions);
        libMesh::Real rho = this->rho( T, this->_p0, R_mix );

        libMesh::Real Y_r = mass_fractions[this->_reactant_gas_species_idx];

        const libMesh::Real r_value = this->compute_reactant_gas_mass_flux(rho, Y_r, T);

        const libMesh::Real p_value = this->compute_product_mass_flux(rho, Y_r, T);

        for (unsigned int i=0; i != n_var_dofs; i++)
          {
            F_r_var(i) += sign*r_value*var_phi_side[i][qp]*jac;

            F_p_var(i) += sign*p_value*var_phi_side[i][qp]*jac;

            if( compute_jacobian )
              {
                libmesh_not_implemented();
              }
          }
      }

    // We're not computing the Jacobian yet
    return false;
  }

  template<typename Chemistry>
  void GasSolidCatalyticWall<Chemistry>::apply_fluxes( AssemblyContext& context,
                                                       const CachedValues& cache,
                                                       const bool request_jacobian )
  {
    libmesh_do_once(libmesh_deprecated());

    libMesh::FEGenericBase<libMesh::Real>* side_fe = NULL;
    context.get_side_fe( _reactant_gas_var_idx, side_fe );

    // The number of local degrees of freedom in each variable.
    const unsigned int n_var_dofs = context.get_dof_indices(_reactant_gas_var_idx).size();

    libmesh_assert_equal_to( n_var_dofs, context.get_dof_indices(_product_var_idx).size() );

    // Element Jacobian * quadrature weight for side integration.
    const std::vector<libMesh::Real> &JxW_side = side_fe->get_JxW();

    // The var shape functions at side quadrature points.
    const std::vector<std::vector<libMesh::Real> >& var_phi_side = side_fe->get_phi();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& var_qpoint = side_fe->get_xyz();

    // reactant residual
    libMesh::DenseSubVector<libMesh::Number> &F_r_var = context.get_elem_residual(_reactant_gas_var_idx);

    // product residual
    libMesh::DenseSubVector<libMesh::Number> &F_p_var = context.get_elem_residual(_product_var_idx);

    unsigned int n_qpoints = context.get_side_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real jac = JxW_side[qp];

        if(Physics::is_axisymmetric())
          {
            const libMesh::Number r = var_qpoint[qp](0);
            jac *= r;
          }

        const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];

        const libMesh::Real Y_r = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp][this->_reactant_gas_species_idx];

        const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

        const libMesh::Real r_value = this->compute_reactant_gas_mass_flux(rho, Y_r, T);

        const libMesh::Real p_value = this->compute_product_mass_flux(rho, Y_r, T);

        for (unsigned int i=0; i != n_var_dofs; i++)
          {
            F_r_var(i) += r_value*var_phi_side[i][qp]*jac;

            F_p_var(i) += p_value*var_phi_side[i][qp]*jac;

            if( request_jacobian )
              {
                libmesh_not_implemented();
              }
          }
      }

    return;
  }

} // end namespace GRINS
