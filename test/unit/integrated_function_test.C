//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

#ifdef GRINS_HAVE_CPPUNIT

#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include "test_comm.h"
#include "grins_test_paths.h"
#include "regression_helper.h"

// GRINS
#include "grins/math_constants.h"
#include "grins/grins_enums.h"
#include "grins/integrated_function.h"
#include "grins/composite_qoi.h"
#include "grins/simulation_builder.h"
#include "grins/simulation.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/parsed_function.h"
#include "libmesh/qoi_set.h"
#include "libmesh/steady_solver.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/elem.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class IntegratedFunctionTest : public CppUnit::TestCase,
                                 public RegressionHelper
  {
  public:
    CPPUNIT_TEST_SUITE( IntegratedFunctionTest );

    CPPUNIT_TEST( test_exact_answer );
    CPPUNIT_TEST( test_convergence );
    CPPUNIT_TEST( qoi_from_input_file );
    CPPUNIT_TEST( reinit_through_system );

    CPPUNIT_TEST_SUITE_END();

  public:

    void tearDown()
    {
      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

    //! Functions that can be integrated exactly using Gauss Quadrature
    void test_exact_answer()
    {
      const std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/integrated_function_quad9.in";
      this->init_sim(filename);

      libMesh::Point origin(0.0,0.0);

      // angles for rayfire
      std::vector<libMesh::Real> theta;
      theta.push_back(0.0);
      theta.push_back(0.25);

      for (unsigned int t=0; t<theta.size(); t++)
        {
          libMesh::Real L = 3.0/std::cos(theta[t]);
          libMesh::Real sintheta = std::sin(theta[t]);
          libMesh::Real costheta = std::cos(theta[t]);

          std::vector<std::string> functions;
          std::vector<libMesh::Real> calc_answers;

          // constant
          functions.push_back("1");
          calc_answers.push_back(L);

          // linear
          functions.push_back("x");
          calc_answers.push_back(0.5/costheta*3.0*3.0);
          functions.push_back("y");
          calc_answers.push_back( sintheta/2.0*L*L );

          // polynomial
          functions.push_back("x*(y^2)");
          calc_answers.push_back(0.25*std::pow(L,4)*std::pow(sintheta,2)*costheta);
          functions.push_back("(4/3)*(x^3)+10");
          calc_answers.push_back( (4.0/12.0)*pow(costheta,3)*std::pow(L,4)+10*L );

          GRINS::MultiphysicsSystem* system = _sim->get_multiphysics_system();

          CPPUNIT_ASSERT_EQUAL(functions.size(), calc_answers.size());

          // and now the CompositeQoI to do the evalutaion
          GRINS::CompositeQoI comp_qoi;

          for (unsigned int i=0; i<functions.size(); i++)
            {
              GRINS::RayfireMesh * rayfire = new GRINS::RayfireMesh(origin,theta[t]);
              GRINS::IntegratedFunction<libMesh::FunctionBase<libMesh::Real> > integ_func((unsigned int)2,new libMesh::ParsedFunction<libMesh::Real>(functions[i]),rayfire,"integrated_function");
              comp_qoi.add_qoi(integ_func);
            }

          comp_qoi.init(*_input,*system);
          system->attach_qoi(&comp_qoi);
          system->assemble_qoi();

          for (unsigned int i=0; i<calc_answers.size(); i++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL( calc_answers[i], _sim->get_qoi_value(i),libMesh::TOLERANCE  );

        } // for t
    }

    //! These functions cannot be integrated exactly wih Gauss Quadrature,
    //! and so we look for quartic convergence to the analytical solution
    /*!
      A quartic convergence rate is identified by graphing log(error) vs. log(h),
      where h is the length of the longest rayfire elem.
      Linear regression is then used to fit a line to that data.
      For quartic convergence, this line should have a slope of 4.
      A 2% tolerance for the slope value is allowed in order to keep
      the test runtime from becoming excessive.
    */
    void test_convergence()
    {
      const std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/integrated_function_quad9.in";
      this->init_sim(filename);

      libMesh::Point origin(0.0,0.0);

      // angles for rayfire
      std::vector<libMesh::Real> theta;
      theta.push_back(0.0);
      theta.push_back(0.25);

      for (unsigned int t=0; t<theta.size(); t++)
        {
          libMesh::Real costheta = std::cos(theta[t]);
          libMesh::Real sintheta = std::sin(theta[t]);
          libMesh::Real L = 3.0/costheta;

          std::vector<std::string> functions; // parsed functions
          std::vector<libMesh::Real> calc_answers; // analytical solutions
          std::vector<libMesh::Real> tolerance; // absolute error tolerance
          std::vector<bool> converged;

          // trig
          functions.push_back("sin(x)+cos(y)");

          if (theta[t] == 0.0)
            calc_answers.push_back( L-std::cos(L)+1.0 );
          else
            calc_answers.push_back( ( (std::sin(L*sintheta))/sintheta ) - ( (std::cos(L*costheta)-1.0)/costheta ) );

          tolerance.push_back(1e-9);
          converged.push_back(false);

          // exponential
          functions.push_back("exp(x)");
          calc_answers.push_back( (std::exp(3.0) - 1.0) / costheta);
          tolerance.push_back(1e-8);
          converged.push_back(false);

          std::vector<std::vector<libMesh::Real> > errors(functions.size());
          std::vector<std::vector<libMesh::Real> > h_vals(functions.size());

          GRINS::MultiphysicsSystem* system = _sim->get_multiphysics_system();
          GRINS::SharedPtr<libMesh::EquationSystems> es = _sim->get_equation_system();

          CPPUNIT_ASSERT_EQUAL(functions.size(), calc_answers.size());

          // and now the CompositeQoI to do the evalutaion
          GRINS::CompositeQoI comp_qoi;

          // all functions use the same mesh, so they have identical rayfires.
          // instead of trying to access each of them directly, we create a reference rayfire that
          // can be used for checking the number of rayfire elements and
          // finding the longest rayfire element for calculating convergence
          GRINS::SharedPtr<GRINS::RayfireMesh> ref_rayfire(new GRINS::RayfireMesh(origin,theta[t]));
          ref_rayfire->init( system->get_mesh() );

          for (unsigned int i=0; i<functions.size(); i++)
            {
              GRINS::RayfireMesh * rayfire = new GRINS::RayfireMesh(origin,theta[t]);
              GRINS::IntegratedFunction<libMesh::FunctionBase<libMesh::Real> > integ_func((unsigned int)3,new libMesh::ParsedFunction<libMesh::Real>(functions[i]),rayfire,"integrated_function");
              comp_qoi.add_qoi(integ_func);
            }

          comp_qoi.init(*_input,*system);
          system->attach_qoi(&comp_qoi);

          libMesh::MeshRefinement mr(system->get_mesh());

          unsigned int num_converged = 0;
          unsigned int iter = 0;

          // limit iterations to prevent excessive runtime
          unsigned int max_iter = 7;

          // calculate qoi's, check error
          // refine until all qoi's converge within their tolerance
          do
            {
              std::vector<libMesh::dof_id_type> elems_in_rayfire;
              ref_rayfire->elem_ids_in_rayfire(elems_in_rayfire);

              system->assemble_qoi();

              // get the length of the longest rayfire elem
              libMesh::Real h = -1.0;
              for (unsigned int i=0; i<elems_in_rayfire.size(); i++)
                {
                  libMesh::Real l = (ref_rayfire->map_to_rayfire_elem(elems_in_rayfire[i]))->length(0,1);
                  if (l>h)
                    h=l;
                }

              CPPUNIT_ASSERT(h > 0.0);

              // check for convergence
              for (unsigned int i=0; i<calc_answers.size(); i++)
                {
                  if (!converged[i])
                    {
                      libMesh::Real err = std::abs( _sim->get_qoi_value(i) - calc_answers[i] );
                      if (err < tolerance[i])
                        {
                          num_converged++;
                          converged[i] = true;
                          h_vals[i].push_back(std::log10(h));
                          errors[i].push_back(std::log10(err));
                        }
                      else
                        {
                          h_vals[i].push_back(std::log10(h));
                          errors[i].push_back(std::log10(err));
                        }
                    }
                }

              // if all functions have not yet converged, refine along the rayfire
              if (num_converged < calc_answers.size())
                {
                  if (iter++ > max_iter+1 )
                    libmesh_error_msg("Exceeded maximum iterations");

                  for (unsigned int i=0; i<elems_in_rayfire.size(); i++)
                    system->get_mesh().elem(elems_in_rayfire[i])->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);

                  mr.refine_elements();

                  // ensure all elems marked for refinement were actually refined
                  for (unsigned int i=0; i<elems_in_rayfire.size(); i++)
                    CPPUNIT_ASSERT( !( system->get_mesh().elem(elems_in_rayfire[i])->active() ) );

                  // need to manually reinit the reference rayfire
                  ref_rayfire->reinit(system->get_mesh());

                  // and reinit the CompositeQoI (which will reinit all internal rayfires)
                  libMesh::DifferentiableQoI* diff_qoi = system->get_qoi();
                  GRINS::CompositeQoI* qoi = libMesh::cast_ptr<GRINS::CompositeQoI*>(diff_qoi);
                  qoi->reinit(*system);

                  es->reinit();
                }

            } while(num_converged < calc_answers.size());

          // verify that all functions had quartic convergence within 2%
          for (unsigned int i=0; i<functions.size(); i++)
            this->check_convergence_rate(h_vals[i],errors[i],0,errors[i].size()-1,4.0,0.08);

        } //for t
    }

    //! Tests that an IntegratedFunction can be initialized and computed
    //! directly from the input file
    void qoi_from_input_file()
    {
      const std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/integrated_function_qoi_quad9.in";
      this->init_sim(filename);

      libMesh::Real theta = 0.25;
      libMesh::Real L = 10.0/std::cos(theta);
      libMesh::Real costheta = std::cos(theta);

      // function: "(4/3)*(x^3)+10"

      libMesh::Real calc_answer = (4.0/12.0)*pow(costheta,3)*std::pow(L,4)+10*L;

      _sim->run();

      CPPUNIT_ASSERT_DOUBLES_EQUAL( calc_answer, _sim->get_qoi_value(0),libMesh::TOLERANCE  );
    }

    //! Here, we use CompositeQoI::reinit() to reinitialize the rayfire post-refinement
    void reinit_through_system()
    {
      const std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/integrated_function_quad9.in";
      this->init_sim(filename);

      libMesh::Point origin(0.0,0.0);
      libMesh::Real theta = 0.0;

      GRINS::MultiphysicsSystem* system = _sim->get_multiphysics_system();
      GRINS::SharedPtr<libMesh::EquationSystems> es = _sim->get_equation_system();

      libMesh::MeshRefinement mr(system->get_mesh());

      // build rayfire
      GRINS::RayfireMesh * rayfire = new GRINS::RayfireMesh(origin,theta);

      // and now the CompositeQoI to do the evalutaion
      GRINS::CompositeQoI comp_qoi;

      libMesh::Real costheta = std::cos(theta);

      // simple function
      std::string function = "x";
      libMesh::Real calc_answer = 0.5/costheta*3.0*3.0;

      GRINS::IntegratedFunction<libMesh::FunctionBase<libMesh::Real> > integ_func((unsigned int)3,new libMesh::ParsedFunction<libMesh::Real>(function),rayfire,"integrated_function");
      comp_qoi.add_qoi(integ_func);

      comp_qoi.init(*_input,*system);
      system->attach_qoi(&comp_qoi);

      // evaluate qoi
      system->assemble_qoi();

      // ensure we get the correct answer
      CPPUNIT_ASSERT_DOUBLES_EQUAL( calc_answer, _sim->get_qoi_value(0),libMesh::TOLERANCE  );

      // get the IntegratedFunction object, since comp_qoi would have already been cloned
      libMesh::DifferentiableQoI* diff_qoi = system->get_qoi();
      GRINS::CompositeQoI* qoi = libMesh::cast_ptr<GRINS::CompositeQoI*>(diff_qoi);
      GRINS::IntegratedFunction<libMesh::FunctionBase<libMesh::Real> > * integrated_func = libMesh::cast_ptr<GRINS::IntegratedFunction<libMesh::FunctionBase<libMesh::Real> > * >( &(qoi->get_qoi(0)) );
      
      // make sure we have exactly 3 elements along the rayfire
      std::vector<libMesh::dof_id_type> elems_in_rayfire;
      integrated_func->get_rayfire().elem_ids_in_rayfire(elems_in_rayfire);
      unsigned int num_rayfire_elems = elems_in_rayfire.size();
      CPPUNIT_ASSERT_EQUAL((unsigned int)3,num_rayfire_elems);

      for (unsigned int i=0; i<elems_in_rayfire.size(); i++)
        system->get_mesh().elem(elems_in_rayfire[i])->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);

      mr.refine_elements();

      // this should trigger RayfireMesh::reinit()
      es->reinit();

      // recalculate the qoi to make sure
      // we still get the same answer
      system->assemble_qoi();

      // after the refinement, we should now have 6 rayfire elements
      std::vector<libMesh::dof_id_type> refined_elems_in_rayfire;
      integrated_func->get_rayfire().elem_ids_in_rayfire(refined_elems_in_rayfire);
      unsigned int num_refined_rayfire_elems = refined_elems_in_rayfire.size();
      CPPUNIT_ASSERT_EQUAL((unsigned int)6,num_refined_rayfire_elems);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( calc_answer, _sim->get_qoi_value(0),libMesh::TOLERANCE );
    }


  private:
    GRINS::SharedPtr<GRINS::Simulation> _sim;
    GRINS::SharedPtr<GetPot> _input;

    //! Initialize the GetPot and Simulation class objects
    void init_sim(const std::string& filename)
    {
      _input.reset(new GetPot(filename));

      const char* const argv = "unit_driver";
      GetPot empty_command_line( (const int)1,&argv );
      GRINS::SimulationBuilder sim_builder;

      _sim = new GRINS::Simulation(*_input,
                                   empty_command_line,
                                   sim_builder,
                                   *TestCommWorld );
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( IntegratedFunctionTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
