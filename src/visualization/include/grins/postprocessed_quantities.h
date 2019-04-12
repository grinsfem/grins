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

    virtual std::unique_ptr<libMesh::FEMFunctionBase<NumericType> >
    clone() const
    {
      return std::unique_ptr<libMesh::FEMFunctionBase<NumericType> >
        ( new PostProcessedQuantities(*this) );
    }

    virtual NumericType operator()( const libMesh::FEMContext& context,
                                    const libMesh::Point& p,
                                    const libMesh::Real time = 0. );

    virtual void operator()( const libMesh::FEMContext& context,
                             const libMesh::Point& p,
                             const libMesh::Real time,
                             libMesh::DenseVector<NumericType>& output );

    virtual NumericType component( const libMesh::FEMContext& context,
                                   unsigned int i,
                                   const libMesh::Point& p,
                                   libMesh::Real time=0. );

    /* Methods for GRINS usage below */

    //! Register quantity to be postprocessed
    /*!
      This method returns an index that will correspond to the name provided.
      Note that this assumes all quantities are scalar. Thus, if there's
      a gradient quantity, each component will need to be separately registered.
    */
    unsigned int register_quantity( std::string name );

    virtual void initialize( MultiphysicsSystem& system,
                             libMesh::EquationSystems& equation_systems );

    virtual void update_quantities( libMesh::EquationSystems& equation_systems );

  protected:

    std::map<std::string, unsigned int> _quantity_name_index_map;
    std::map<VariableIndex, unsigned int> _quantity_index_var_map;

    MultiphysicsSystem* _multiphysics_sys;
    std::shared_ptr<AssemblyContext> _multiphysics_context;

  private:

    PostProcessedQuantities();

  };

} // namespace GRINS

#endif //GRINS_POSTPROCESSED_QUANTITIES_H
