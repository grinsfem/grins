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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_POSTPROCESSED_QUANTITIES_H
#define GRINS_POSTPROCESSED_QUANTITIES_H

//libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/equation_systems.h"

//GRINS
#include "grins/multiphysics_sys.h"

namespace GRINS
{
  template<class NumericType>
  class PostProcessedQuantities : public libMesh::FEMFunctionBase<NumericType>
  {
  public:

    PostProcessedQuantities( const GetPot& input );
    virtual ~PostProcessedQuantities();

    /* Methods to override from FEMFunctionBase needed for libMesh-based evaluations */
    virtual void init_context( const libMesh::FEMContext & context);

    virtual AutoPtr<libMesh::FEMFunctionBase<NumericType> > clone() const
    { return AutoPtr<libMesh::FEMFunctionBase<NumericType> >( new PostProcessedQuantities(*this) ); }

    virtual NumericType operator()( const libMesh::FEMContext& context, 
				    const libMesh::Point& p,
				    const Real time = 0. );

    virtual void operator()( const libMesh::FEMContext& context, 
			     const libMesh::Point& p,
			     const Real time,
			     libMesh::DenseVector<NumericType>& output );

    virtual NumericType component( const libMesh::FEMContext& context, 
				   unsigned int i,
				   const libMesh::Point& p,
				   Real time=0. );

    /* Methods for GRINS usage below */
    virtual void initialize( MultiphysicsSystem& system,
			     libMesh::EquationSystems& equation_systems );

    virtual void update_quantities( libMesh::EquationSystems& equation_systems );

    unsigned int n_quantities() const
    {return _quantities.size();}

  protected:

    virtual void build_name_map();

    virtual void init_quantities( const MultiphysicsSystem& multiphysics_system,
				  libMesh::System& output_system,
				  const unsigned int component );

    virtual NumericType compute_quantities( const unsigned int component ) const;

    std::vector<unsigned int> _quantities;
    std::map<std::string, unsigned int> _quantity_name_map;
    std::map<VariableIndex,unsigned int> _quantity_var_map;
    
    MultiphysicsSystem* _multiphysics_sys;
    std::tr1::shared_ptr<libMesh::FEMContext> _multiphysics_context;

    CachedValues _cache;

    libMesh::Point _prev_point;

    //! Place to cache species names for species-dependent quantities.
    std::vector<std::string> _species_names;

    //! Cache the map between species-components variable indices and species number
    std::map<VariableIndex, unsigned int> _species_var_map;

  private:

    PostProcessedQuantities();

    /* Perfect gas ==> LowMachNavierStokes
     * Species/Mixture ==> ReactingLowMachNavierStokes */
    enum QuantityList{ PERFECT_GAS_DENSITY = 0,
                       MIXTURE_DENSITY,
                       PERFECT_GAS_VISCOSITY,
                       SPECIES_VISCOSITY,
                       MIXTURE_VISCOSITY,
                       PERFECT_GAS_THERMAL_CONDUCTIVITY,
                       SPECIES_THERMAL_CONDUCTIVITY,
                       MIXTURE_THERMAL_CONDUCTIVITY,
                       PERFECT_GAS_SPECIFIC_HEAT_P,
                       SPECIES_SPECIFIC_HEAT_P,
                       MIXTURE_SPECIFIC_HEAT_P,
                       PERFECT_GAS_SPECIFIC_HEAT_V,
                       SPECIES_SPECIFIC_HEAT_V,
                       MIXTURE_SPECIFIC_HEAT_V,
                       MOLE_FRACTIONS,
                       SPECIES_ENTHALPY,
                       OMEGA_DOT };

  };

} // namespace GRINS

#endif //GRINS_POSTPROCESSED_QUANTITIES_H
