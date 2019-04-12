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


#include "grins_config.h"

#include <iostream>

// GRINS
#include "grins/runner.h"
#include "grins/mesh_builder.h"
#include "grins/multiphysics_sys.h"
#include "grins/dirichlet_bc_factory_function_old_style_base.h"
#include "grins/var_typedefs.h"

//libMesh
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh_function.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/gmv_io.h"
#include "libmesh/exact_solution.h"
#include "libmesh/zero_function.h"
#include "libmesh/composite_function.h"
#include "libmesh/zero_function.h"
#include "libmesh/enum_xdr_mode.h"

namespace GRINS
{
  class MultiphysicsSystem;
}

void test_error_norm( libMesh::ExactSolution& exact_sol,
                      const std::string& system_name,
                      const std::string& var,
                      const std::string& norm,
                      const double tol,
                      int& return_flag );

namespace GRINSTesting
{

  class TurbBoundFuncBase : public libMesh::FunctionBase<libMesh::Number>
  {
  public:
    TurbBoundFuncBase (libMesh::MeshFunction* turbulent_bc_values)
      : _turbulent_bc_values(turbulent_bc_values)
    { this->_initialized = true; }

    virtual libMesh::Number operator() (const libMesh::Point&, const libMesh::Real = 0)
    { libmesh_not_implemented(); }

    virtual libMesh::Number compute_final_value( const libMesh::DenseVector<libMesh::Number>& u_nu_values ) =0;

    virtual void operator() (const libMesh::Point& p,
                             const libMesh::Real t,
                             libMesh::DenseVector<libMesh::Number>& output)
    {
      output.resize(1);
      output.zero();

      // Since the turbulent_bc_values object has a solution from a 1-d problem, we have to zero out the y coordinate of p
      libMesh::Point p_copy(p);
      // Also, the 1-d solution provided is on the domain [0, 1] on the x axis and we need to map this to the corresponding point on the y axis
      p_copy(0) = p_copy(1);
      p_copy(1)= 0.0;
      // Also, the 1-d solution provided is actually a symmetry solution, so we have to make the following map
      // x_GRINS < 0.5 => x_meshfunction = 2*x_GRINS , x_GRINS >= 0.5 => x_GRINS = 1 - x_GRINS, x_meshfunction = 2*x_GRINS
      if(p_copy(0) > 0.5)
        {
          p_copy(0) = 1 - p_copy(0);
        }
      p_copy(0) = 2*p_copy(0);

      libMesh::DenseVector<libMesh::Number> u_nu_values;
      _turbulent_bc_values->operator()(p_copy, t, u_nu_values);

      output(0) = this->compute_final_value(u_nu_values);
    }

  protected:

    libMesh::MeshFunction* _turbulent_bc_values;

  };

  // Class to construct the Dirichlet boundary object and operator for the inlet u velocity and nu profiles
  class TurbulentBdyFunctionU : public TurbBoundFuncBase
  {
  public:
    TurbulentBdyFunctionU (libMesh::MeshFunction* turbulent_bc_values)
      : TurbBoundFuncBase(turbulent_bc_values)
    {}

    virtual libMesh::Number compute_final_value( const libMesh::DenseVector<libMesh::Number>& u_nu_values )
    { return  u_nu_values(0)/21.995539; }

    virtual std::unique_ptr<libMesh::FunctionBase<libMesh::Number> > clone() const
    { return std::unique_ptr<libMesh::FunctionBase<libMesh::Number> > (new TurbulentBdyFunctionU(_turbulent_bc_values)); }
  };

  // Class to construct the Dirichlet boundary object and operator for the inlet u velocity and nu profiles
  class TurbulentBdyFunctionNu : public TurbBoundFuncBase
  {
  public:
    TurbulentBdyFunctionNu (libMesh::MeshFunction* turbulent_bc_values)
      : TurbBoundFuncBase(turbulent_bc_values)
    {}

    virtual libMesh::Number compute_final_value( const libMesh::DenseVector<libMesh::Number>& u_nu_values )
    { return  u_nu_values(1)/(2.0*21.995539); }

    virtual std::unique_ptr<libMesh::FunctionBase<libMesh::Number> > clone() const
    { return std::unique_ptr<libMesh::FunctionBase<libMesh::Number> > (new TurbulentBdyFunctionNu(_turbulent_bc_values)); }

  };

  class SATurbBCFactoryBase : public GRINS::DirichletBCFactoryFunctionBase<libMesh::FunctionBase<libMesh::Number> >
  {
  public:
    SATurbBCFactoryBase( const std::string& bc_type_name )
      : DirichletBCFactoryFunctionBase<libMesh::FunctionBase<libMesh::Number> >(bc_type_name)
    {}

    static void set_turb_bc_values( libMesh::MeshFunction* turbulent_bc_values )
    { _turbulent_bc_values = turbulent_bc_values; }

  protected:

    static libMesh::MeshFunction* _turbulent_bc_values;
  };

  libMesh::MeshFunction* SATurbBCFactoryBase::_turbulent_bc_values = NULL;

  class SATurbUBCFactory : public SATurbBCFactoryBase
  {
  public:

    SATurbUBCFactory( const std::string& bc_type_name )
      : SATurbBCFactoryBase(bc_type_name)
    {}

  protected:

    virtual std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
    build_func( const GetPot& /*input*/,
                GRINS::MultiphysicsSystem& system,
                std::vector<std::string>& var_names,
                const std::string& /*section*/ )
    {
      libmesh_assert_equal_to(var_names.size(), 2);
      libmesh_assert_equal_to(var_names[0], std::string("u"));
      libmesh_assert_equal_to(var_names[1], std::string("v"));

      std::unique_ptr<libMesh::CompositeFunction<libMesh::Number> >
        composite_func( new libMesh::CompositeFunction<libMesh::Number> );

      std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
        bound_func(new TurbulentBdyFunctionU(_turbulent_bc_values) );

      std::vector<GRINS::VariableIndex> vars(1);
      vars[0] = system.variable_number(var_names[0]);
      composite_func->attach_subfunction(*bound_func, vars);

      vars[0] = system.variable_number(var_names[1]);
      composite_func->attach_subfunction(libMesh::ZeroFunction<libMesh::Number>(), vars);

      return std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >( composite_func.release() );
    }
  };

  class SATurbNuBCFactory : public SATurbBCFactoryBase
  {
  public:

    SATurbNuBCFactory( const std::string& bc_type_name )
      : SATurbBCFactoryBase(bc_type_name)
    {}

  protected:

    virtual std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
    build_func( const GetPot& /*input*/,
                GRINS::MultiphysicsSystem& system,
                std::vector<std::string>& var_names,
                const std::string& /*section*/ )
    {
      libmesh_assert_equal_to(var_names.size(), 1);
      libmesh_assert_equal_to(var_names[0], std::string("nu"));

      std::unique_ptr<libMesh::CompositeFunction<libMesh::Number> >
        composite_func( new libMesh::CompositeFunction<libMesh::Number> );

      std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
        bound_func( new TurbulentBdyFunctionNu(_turbulent_bc_values) );

      std::vector<GRINS::VariableIndex> vars(1);
      vars[0] = system.variable_number(var_names[0]);
      composite_func->attach_subfunction(*bound_func, vars);

      return std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >( composite_func.release() );
    }
  };

} // end namespace GRINSTesting

int main(int argc, char* argv[])
{
  // Factories needed for run
  GRINSTesting::SATurbUBCFactory grins_factory_testing_turb_u_bc("testing_turb_u");
  GRINSTesting::SATurbNuBCFactory grins_factory_testing_turb_nu_bc("testing_turb_nu");

  // Need only
  GRINS::Runner grins(argc,argv);

  const GetPot & command_line = grins.get_command_line();

  const GetPot & inputfile = grins.get_input_file();

  // Don't flag our command-line-specific variables as UFOs later
  inputfile.have_variable("mesh-1d");
  inputfile.have_variable("data-1d");
  inputfile.have_variable("soln-data");
  inputfile.have_variable("vars");
  inputfile.have_variable("norms");
  inputfile.have_variable("tol");

  const libMesh::LibMeshInit & libmesh_init = grins.get_libmesh_init();

  // Build a 1-d turbulent_bc_system to get the bc data from files
  libMesh::SerialMesh mesh(libmesh_init.comm());

  std::string oned_mesh = command_line("mesh-1d", "DIE!");
  mesh.read(oned_mesh);

  // And an EquationSystems to run on it
  libMesh::EquationSystems equation_systems (mesh);

  std::string oned_data = command_line("data-1d", "DIE!");
  equation_systems.read(oned_data, libMesh::XdrMODE::READ,
                        libMesh::EquationSystems::READ_HEADER |
                        libMesh::EquationSystems::READ_DATA |
                        libMesh::EquationSystems::READ_ADDITIONAL_DATA);

  // Get a reference to the system so that we can call update() on it
  libMesh::System & turbulent_bc_system = equation_systems.get_system<libMesh::System>("flow");

  // We need to call update to put system in a consistent state
  // with the solution that was read in
  turbulent_bc_system.update();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Prepare a global solution and a MeshFunction of the Turbulent system
  std::unique_ptr<libMesh::MeshFunction> turbulent_bc_values;

  std::unique_ptr<libMesh::NumericVector<libMesh::Number> > turbulent_bc_soln = libMesh::NumericVector<libMesh::Number>::build(equation_systems.comm());

  std::vector<libMesh::Number> flow_soln;

  turbulent_bc_system.update_global_solution(flow_soln);

  turbulent_bc_soln->init(turbulent_bc_system.solution->size(), true, libMesh::SERIAL);

  (*turbulent_bc_soln) = flow_soln;

  std::vector<unsigned int>turbulent_bc_system_variables;
  turbulent_bc_system_variables.push_back(0);
  turbulent_bc_system_variables.push_back(1);

  turbulent_bc_values = std::unique_ptr<libMesh::MeshFunction>
    (new libMesh::MeshFunction(equation_systems,
                               *turbulent_bc_soln,
                               turbulent_bc_system.get_dof_map(),
                               turbulent_bc_system_variables ));

  turbulent_bc_values->init();

  GRINSTesting::SATurbBCFactoryBase::set_turb_bc_values( turbulent_bc_values.get() );

  // Initialize
  grins.init();

  // Solve
  grins.run();

  GRINS::Simulation & sim = grins.get_simulation();

  // Get equation systems to create ExactSolution object
  std::shared_ptr<libMesh::EquationSystems> es = sim.get_equation_system();

  // Create Exact solution object and attach exact solution quantities
  //libMesh::ExactSolution exact_sol(*es);

  // Create Exact solution object and attach exact solution quantities
  libMesh::ExactSolution exact_sol(*es);

  libMesh::EquationSystems es_ref( es->get_mesh() );

  // Filename of file where comparison solution is stashed
  std::string solution_file = command_line("soln-data", "DIE!");
  es_ref.read( solution_file, libMesh::XdrMODE::DECODE,
               libMesh::EquationSystems::READ_HEADER |
               libMesh::EquationSystems::READ_DATA |
               libMesh::EquationSystems::READ_ADDITIONAL_DATA);

  exact_sol.attach_reference_solution( &es_ref );

  // Now grab the variables for which we want to compare
  unsigned int n_vars = command_line.vector_variable_size("vars");
  std::vector<std::string> vars(n_vars);
  for( unsigned int v = 0; v < n_vars; v++ )
    {
      vars[v] = command_line("vars", "DIE!", v);
    }

  // Now grab the norms to compute for each variable error
  unsigned int n_norms = command_line.vector_variable_size("norms");
  std::vector<std::string> norms(n_norms);
  for( unsigned int n = 0; n < n_norms; n++ )
    {
      norms[n] = command_line("norms", "DIE!", n);
      if( norms[n] != std::string("L2") &&
          norms[n] != std::string("H1") )
        {
          std::cerr << "ERROR: Invalid norm input " << norms[n] << std::endl
                    << "       Valid values are: L2" << std::endl
                    << "                         H1" << std::endl;
        }
    }

  const std::string& system_name = sim.get_multiphysics_system_name();

  // Now compute error for each variable
  for( unsigned int v = 0; v < n_vars; v++ )
    {
      exact_sol.compute_error(system_name, vars[v]);
    }

  int return_flag = 0;

  double tol = command_line("tol", 1.0e-10);

  // Now test error for each variable, for each norm
  for( unsigned int v = 0; v < n_vars; v++ )
    {
      for( unsigned int n = 0; n < n_norms; n++ )
        {
          test_error_norm( exact_sol, system_name, vars[v], norms[n], tol, return_flag );
        }
    }

  return return_flag;
}

void test_error_norm( libMesh::ExactSolution& exact_sol,
                      const std::string& system_name,
                      const std::string& var,
                      const std::string& norm,
                      const double tol,
                      int& return_flag )
{
  // We don't set return_flag unless we are setting it 1
  // since this function gets called multiple times and we don't
  // want to overwrite a previous "fail" (return_flag = 1) with
  // a "pass" (return_flag = 0)

  double error = 0.0;

  std::cout << "==========================================================" << std::endl
            << "Checking variable " << var << " using error norm " << norm << " with tol " << tol << "...";

  if( norm == std::string("L2") )
    {
      error = exact_sol.l2_error(system_name, var);
    }
  else if( norm == std::string("H1") )
    {
      error = exact_sol.h1_error(system_name, var);
    }
  else
    {
      std::cerr << "ERROR: Invalid norm " << norm << std::endl;
      exit(1);
    }

  if( error > tol )
    {
      return_flag = 1;

      std::cerr << "Tolerance exceeded for generic regression test!" << std::endl
                << "tolerance     = " << tol << std::endl
                << "norm of error = " << error << std::endl
                << "norm type     = " << norm << std::endl
                << "var           = " << var << std::endl;
    }
  else
    {
      std::cout << "PASSED!" << std::endl
                << "==========================================================" << std::endl;
    }

  return;
}
