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
#include "grins/integrated_function.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"
#include "grins/materials_parsing.h"
#include "grins/rayfire_mesh.h"
#include "grins/fem_function_and_derivative_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/fe_type.h"
#include "libmesh/function_base.h"

namespace GRINS
{
  template<typename Function>
  IntegratedFunction<Function>::IntegratedFunction(unsigned int p_level, const std::shared_ptr<Function> & f, const std::shared_ptr<RayfireMesh> & rayfire, const std::string & qoi_name) :
    QoIBase(qoi_name),
    _p_level(p_level),
    _f(f),
    _rayfire(rayfire)
  {}

  template<typename Function>
  QoIBase* IntegratedFunction<Function>::clone() const
  {
    return new IntegratedFunction<Function>( *this );
  }

  template<typename Function>
  void IntegratedFunction<Function>::init
  (const GetPot & /*input*/,
   const MultiphysicsSystem & system,
   unsigned int /*qoi_num*/ )
  {
    _multiphysics_system = &( const_cast<MultiphysicsSystem &>(system) );
    _rayfire->init(system.get_mesh());
  }

  template<typename Function>
  void IntegratedFunction<Function>::reinit(MultiphysicsSystem & system)
  {
    _rayfire->reinit(system.get_mesh());
  }

  template<typename Function>
  void IntegratedFunction<Function>::element_qoi( AssemblyContext & context,
                                                  const unsigned int qoi_index )
  {
    const libMesh::Elem & original_elem = context.get_elem();
    const libMesh::Elem * rayfire_elem = _rayfire->map_to_rayfire_elem(original_elem.id());

    // rayfire_elem will be NULL if the main_elem
    // is not in the rayfire
    if (rayfire_elem)
      {
        // create and init the quadrature base on the rayfire elem
        libMesh::QGauss qbase(rayfire_elem->dim(),libMesh::Order(_p_level));
        qbase.init(rayfire_elem->type(),libMesh::Order(_p_level));

        // need the QP coordinates and JxW
        std::unique_ptr< libMesh::FEGenericBase<libMesh::Real> > fe = libMesh::FEGenericBase<libMesh::Real>::build(rayfire_elem->dim(), libMesh::FEType(libMesh::FIRST, libMesh::LAGRANGE) );

        fe->attach_quadrature_rule( &qbase );
        const std::vector<libMesh::Real> & JxW = fe->get_JxW();
        const std::vector<libMesh::Point> & xyz = fe->get_xyz();

        fe->reinit(rayfire_elem);

        const unsigned int n_qpoints = fe->n_quadrature_points();

        libMesh::Number & qoi = context.get_qois()[qoi_index];

        for (unsigned int qp = 0; qp != n_qpoints; ++qp)
          qoi += this->qoi_value((*_f),context,xyz[qp])*JxW[qp];

      }
  }

  template<typename Function>
  void IntegratedFunction<Function>::element_qoi_derivative( AssemblyContext & context,
                                                             const unsigned int qoi_index )
  {
    const libMesh::Elem& original_elem = context.get_elem();
    const libMesh::Elem* rayfire_elem = _rayfire->map_to_rayfire_elem(original_elem.id());

    // rayfire_elem will be NULL if the main_elem
    // is not in the rayfire
    if (rayfire_elem)
      {
        // create and init the quadrature base on the rayfire elem
        libMesh::QGauss qbase(rayfire_elem->dim(),libMesh::Order(_p_level));
        qbase.init(rayfire_elem->type(),libMesh::Order(_p_level));

        // need the QP coordinates and JxW
        std::unique_ptr< libMesh::FEGenericBase<libMesh::Real> > fe = libMesh::FEGenericBase<libMesh::Real>::build(rayfire_elem->dim(), libMesh::FEType(libMesh::FIRST, libMesh::LAGRANGE) );

        fe->attach_quadrature_rule( &qbase );
        const std::vector<libMesh::Real> & JxW = fe->get_JxW();
        const std::vector<libMesh::Point> & xyz = fe->get_xyz();

        fe->reinit(rayfire_elem);

        const unsigned int n_qpoints = fe->n_quadrature_points();

        // FIXME the inverse_map() is ugly AF
        for (unsigned int qp = 0; qp != n_qpoints; qp++)
          this->qoi_derivative((*_f),context,xyz[qp],JxW[qp],qoi_index);

      }
  }

  // speciaizations of the qoi_value() function
  template<>
  libMesh::Real IntegratedFunction<FEMFunctionAndDerivativeBase<libMesh::Real> >::qoi_value(FEMFunctionAndDerivativeBase<libMesh::Real> & f, AssemblyContext & context, const libMesh::Point & xyz)
  {
    return f(context,xyz);
  }

  template<>
  libMesh::Real IntegratedFunction<libMesh::FunctionBase<libMesh::Real> >::qoi_value(libMesh::FunctionBase<libMesh::Real> & f, AssemblyContext & /*context*/, const libMesh::Point & xyz)
  {
    return f(xyz);
  }

  // speciaizations of the qoi_derivative() function
  template<>
  void IntegratedFunction<FEMFunctionAndDerivativeBase<libMesh::Real> >::qoi_derivative( FEMFunctionAndDerivativeBase<libMesh::Real> & f, AssemblyContext & context,
                                                                                         const libMesh::Point & qp_xyz, const libMesh::Real JxW, const unsigned int qoi_index)
  {
    f.derivatives(context,qp_xyz,JxW,qoi_index);
  }

  template<>
  void IntegratedFunction<libMesh::FunctionBase<libMesh::Real> >::qoi_derivative( libMesh::FunctionBase<libMesh::Real> & /*f*/, AssemblyContext & /*context*/,
                                                                                  const libMesh::Point & /*qp_xyz*/, const libMesh::Real /*JxW*/, const unsigned int /*qoi_index*/)
  {
    // derivatives are always zero for FunctionBase
  }

  template<>
  void IntegratedFunction<libMesh::FunctionBase<libMesh::Real>>::register_active_vars
  ( std::set<unsigned int> & /*element_vars*/,
    std::set<unsigned int> & /*side_vars*/ )
  {}

  template<>
  void IntegratedFunction<FEMFunctionAndDerivativeBase<libMesh::Real>>::register_active_vars
  ( std::set<unsigned int> & element_vars,
    std::set<unsigned int> & side_vars )
  {
    _f->register_active_vars(element_vars,side_vars);
  }

  template class IntegratedFunction<libMesh::FunctionBase<libMesh::Real> >;
  template class IntegratedFunction<FEMFunctionAndDerivativeBase<libMesh::Real> >;

} //namespace GRINS
