//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/parsed_interior_qoi.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/fe_base.h"
#include "libmesh/parsed_fem_function.h"

namespace GRINS
{
  ParsedInteriorQoI::ParsedInteriorQoI( const std::string& qoi_name )
    : QoIBase(qoi_name) {}

  ParsedInteriorQoI::ParsedInteriorQoI( const ParsedInteriorQoI& original )
    : QoIBase(original)
  {
    if (original.qoi_functional.get())
      {
        this->qoi_functional = original.qoi_functional->clone();
        this->move_parameter
          (*libMesh::libmesh_cast_ptr<libMesh::ParsedFEMFunction<libMesh::Number>*>
             (original.qoi_functional.get()),
           *libMesh::libmesh_cast_ptr<libMesh::ParsedFEMFunction<libMesh::Number>*>
             (this->qoi_functional.get()));
      }
  }

  ParsedInteriorQoI::~ParsedInteriorQoI() {}

  QoIBase* ParsedInteriorQoI::clone() const
  {
    return new ParsedInteriorQoI( *this );
  }

  void ParsedInteriorQoI::init
    (const GetPot& input,
     const MultiphysicsSystem& system,
     unsigned int /*qoi_num*/ )
  {
    libMesh::ParsedFEMFunction<libMesh::Number> *qf
      (new libMesh::ParsedFEMFunction<libMesh::Number>
         (system, ""));
    this->qoi_functional.reset(qf);

    this->set_parameter(*qf, input,
                        "QoI/ParsedInterior/qoi_functional", "DIE!");
  }

  void ParsedInteriorQoI::init_context( AssemblyContext& context )
  {
    libMesh::FEBase* element_fe;
    context.get_element_fe<libMesh::Real>(0, element_fe);
    element_fe->get_JxW();
    element_fe->get_xyz();

    qoi_functional->init_context(context);
  }

  void ParsedInteriorQoI::element_qoi( AssemblyContext& context,
                                       const unsigned int qoi_index )
  {
    libMesh::FEBase* element_fe;
    context.get_element_fe<libMesh::Real>(0, element_fe);
    const std::vector<libMesh::Real> &JxW = element_fe->get_JxW();

    const std::vector<libMesh::Point>& x_qp = element_fe->get_xyz();

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    /*! \todo Need to generalize this to the multiple QoI case */
    libMesh::Number& qoi = context.get_qois()[qoi_index];

    for( unsigned int qp = 0; qp != n_qpoints; qp++ )
      {
        const libMesh::Number func_val =
          (*qoi_functional)(context, x_qp[qp], context.get_time());

        qoi += func_val * JxW[qp];
      }
  }

  void ParsedInteriorQoI::element_qoi_derivative( AssemblyContext& context,
                                                  const unsigned int qoi_index )
  {
    libMesh::FEBase* element_fe;
    context.get_element_fe<libMesh::Real>(0, element_fe);
    const std::vector<libMesh::Real> &JxW = element_fe->get_JxW();

    const std::vector<libMesh::Point>& x_qp = element_fe->get_xyz();

    // Local DOF count and quadrature point count
    const unsigned int n_u_dofs = context.get_dof_indices().size();

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // Local solution vector - non-const version for finite
    // differenting purposes
    libMesh::DenseVector<libMesh::Number>& elem_solution =
      const_cast<libMesh::DenseVector<libMesh::Number>&>
        (context.get_elem_solution());

    /*! \todo Need to generalize this to the multiple QoI case */
    libMesh::DenseVector<libMesh::Number> &Qu =
      context.get_qoi_derivatives()[qoi_index];

    for( unsigned int qp = 0; qp != n_qpoints; qp++ )
      {
        // Central finite differencing to approximate derivatives.
        // FIXME - we should hook the FParserAD stuff into
        // ParsedFEMFunction

        for( unsigned int i = 0; i != n_u_dofs; ++i )
          {
            libMesh::Number &current_solution = elem_solution(i);
            const libMesh::Number original_solution = current_solution;

            current_solution = original_solution + libMesh::TOLERANCE;

            const libMesh::Number plus_val =
              (*qoi_functional)(context, x_qp[qp], context.get_time());

            current_solution = original_solution - libMesh::TOLERANCE;

            const libMesh::Number minus_val =
              (*qoi_functional)(context, x_qp[qp], context.get_time());

            Qu(i) += (plus_val - minus_val) *
                     (0.5 / libMesh::TOLERANCE) * JxW[qp];

            // Don't forget to restore the correct solution...
            current_solution = original_solution;
          }
      }
  }

} //namespace GRINS
