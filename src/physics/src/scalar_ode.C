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
#include "grins/scalar_ode.h"

// GRINS
#include "grins/generic_ic_handler.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/boundary_info.h"
#include "libmesh/parsed_fem_function.h"

namespace GRINS
{

  ScalarODE::ScalarODE( const std::string& physics_name, const GetPot& input )
    : Physics(physics_name, input),
      _epsilon(1e-6),
      _input(input),
      _var(GRINSPrivate::VariableWarehouse::get_variable_subclass<ScalarVariable>(VariablesParsing::scalar_variable_name(input,physics_name,VariablesParsing::PHYSICS)))
  {
    this->read_input_options(input);

    this->_ic_handler = new GenericICHandler( physics_name, input );
  }

  void ScalarODE::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    system->time_evolving(this->scalar_ode_var(),1);

    // FIXME: this doesn't fit here at all, but it's the only time
    // we've clearly got a System to grab hold of with all it's
    // variables initialized.

    {
      auto tdf = std::make_unique<libMesh::ParsedFEMFunction<libMesh::Number>>(*system, "");

      this->set_parameter(*tdf, _input,
                          "Physics/"+PhysicsNaming::scalar_ode()+"/time_deriv",
                          std::string("0"));

      if (tdf->expression() == "0")
        std::cout << "Warning! Zero time_deriv function specified!" << std::endl;

      this->time_deriv_function = std::move(tdf);
    }

    {
      auto mrf = std::make_unique<libMesh::ParsedFEMFunction<libMesh::Number>>(*system, "");

      if (mrf->expression() == "0")
        std::cout << "Warning! Zero mass_residual function specified!" << std::endl;

      this->set_parameter(*mrf, _input,
                          "Physics/"+PhysicsNaming::scalar_ode()+"/mass_residual",
                          std::string("0"));

      this->mass_residual_function = std::move(mrf);
    }

    {
      auto cf = std::make_unique<libMesh::ParsedFEMFunction<libMesh::Number>>(*system, "");

      this->set_parameter(*cf, _input,
                          "Physics/"+PhysicsNaming::scalar_ode()+"/constraint",
                          std::string("0"));

      this->constraint_function = std::move(cf);
    }
  }


  void ScalarODE::read_input_options( const GetPot& input )
  {
    this->_epsilon = input("Physics/"+PhysicsNaming::scalar_ode()+"/epsilon", 1e-6);
  }


  void ScalarODE::init_context( AssemblyContext& context )
  {
    mass_residual_function->init_context(context);
    time_deriv_function->init_context(context);
    constraint_function->init_context(context);

    context.get_side_fe(this->scalar_ode_var())->get_nothing();
  }


  void ScalarODE::nonlocal_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    libMesh::DenseSubMatrix<libMesh::Number> &Kss =
      context.get_elem_jacobian(this->scalar_ode_var(), this->scalar_ode_var()); // R_{s},{s}

    libMesh::DenseSubVector<libMesh::Number> &Fs =
      context.get_elem_residual(this->scalar_ode_var()); // R_{s}

    const libMesh::Number time_deriv =
      (*time_deriv_function)(context, libMesh::Point(0),
                             context.get_time());

    Fs(0) += time_deriv;

    if (compute_jacobian)
      {
        // FIXME: we should replace this hacky FDM with a hook to the
        // AD fparser stuff
        libMesh::DenseSubVector<libMesh::Number> &Us =
          const_cast<libMesh::DenseSubVector<libMesh::Number>&>
          (context.get_elem_solution(this->scalar_ode_var())); // U_{s}

        const libMesh::Number s = Us(0);
        Us(0) = s + this->_epsilon;
        libMesh::Number time_deriv_jacobian =
          (*time_deriv_function)(context, libMesh::Point(0),
                                 context.get_time());

        Us(0) = s - this->_epsilon;
        time_deriv_jacobian -=
          (*time_deriv_function)(context, libMesh::Point(0),
                                 context.get_time());

        Us(0) = s;
        time_deriv_jacobian /= (2*this->_epsilon);

        Kss(0,0) += time_deriv_jacobian *
          context.get_elem_solution_derivative();
      }
  }


  void ScalarODE::nonlocal_mass_residual
  (bool compute_jacobian, AssemblyContext & context )
  {
    libMesh::DenseSubMatrix<libMesh::Number> &Kss =
      context.get_elem_jacobian(this->scalar_ode_var(), this->scalar_ode_var()); // R_{s},{s}

    libMesh::DenseSubVector<libMesh::Number> &Fs =
      context.get_elem_residual(this->scalar_ode_var()); // R_{s}

    const libMesh::Number mass_res =
      (*mass_residual_function)(context, libMesh::Point(0),
                                context.get_time());

    Fs(0) -= mass_res;

    if (compute_jacobian)
      {
        // FIXME: we should replace this hacky FDM with a hook to the
        // AD fparser stuff
        libMesh::DenseSubVector<libMesh::Number> &Us =
          const_cast<libMesh::DenseSubVector<libMesh::Number>&>
          (context.get_elem_solution_rate(this->scalar_ode_var())); // U_{s}

        const libMesh::Number s = Us(0);
        Us(0) = s + this->_epsilon;
        libMesh::Number mass_residual_jacobian =
          (*mass_residual_function)(context, libMesh::Point(0),
                                    context.get_time());

        Us(0) = s - this->_epsilon;
        mass_residual_jacobian -=
          (*mass_residual_function)(context, libMesh::Point(0),
                                    context.get_time());

        Us(0) = s;
        mass_residual_jacobian /= (2*this->_epsilon);

        Kss(0,0) -= mass_residual_jacobian *
          context.get_elem_solution_rate_derivative();
      }
  }


  void ScalarODE::nonlocal_constraint
  ( bool compute_jacobian, AssemblyContext & context )
  {
    libMesh::DenseSubMatrix<libMesh::Number> &Kss =
      context.get_elem_jacobian(this->scalar_ode_var(), this->scalar_ode_var()); // R_{s},{s}

    libMesh::DenseSubVector<libMesh::Number> &Fs =
      context.get_elem_residual(this->scalar_ode_var()); // R_{s}

    const libMesh::Number constraint =
      (*constraint_function)(context, libMesh::Point(0),
                             context.get_time());

    Fs(0) += constraint;

    if (compute_jacobian)
      {
        // FIXME: we should replace this hacky FDM with a hook to the
        // AD fparser stuff
        libMesh::DenseSubVector<libMesh::Number> &Us =
          const_cast<libMesh::DenseSubVector<libMesh::Number>&>
          (context.get_elem_solution(this->scalar_ode_var())); // U_{s}

        const libMesh::Number s = Us(0);
        Us(0) = s + this->_epsilon;
        libMesh::Number constraint_jacobian =
          (*constraint_function)(context, libMesh::Point(0),
                                 context.get_time());

        Us(0) = s - this->_epsilon;
        constraint_jacobian -=
          (*constraint_function)(context, libMesh::Point(0),
                                 context.get_time());

        Us(0) = s;
        constraint_jacobian /= (2*this->_epsilon);

        Kss(0,0) += constraint_jacobian *
          context.get_elem_solution_derivative();
      }
  }

} // namespace GRINS
