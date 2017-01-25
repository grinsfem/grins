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

// GRINS
#include "grins/math_constants.h"
#include "grins/grins_enums.h"
#include "grins/simulation_builder.h"
#include "grins/simulation.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/parsed_function.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class SpectroscopicAbsorptionTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( SpectroscopicAbsorptionTest );
#if GRINS_HAVE_ANTIOCH
    CPPUNIT_TEST( single_elem_mesh );
    CPPUNIT_TEST( multi_elem_mesh );
#endif
    CPPUNIT_TEST_SUITE_END();

  public:

    void tearDown()
    {
      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

    //! Single QUAD4 elem, uniform T,P,Y
    void single_elem_mesh()
    {
      const std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/spectroscopic_absorption_qoi.in";
      libMesh::Real calc_answer = 0.520403290868787;

      this->run_test(filename,calc_answer);
    }

    //! 10x10 mesh, uniform T,P,Y
    void multi_elem_mesh()
    {
      const std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/spectroscopic_absorption_qoi_fine.in";
      libMesh::Real calc_answer = 0.520403290868787; // same physical conditions as single_elem_mesh, so answer should not change

      this->run_test(filename,calc_answer);
    }

  private:
    GRINS::SharedPtr<GRINS::Simulation> _sim;
    GRINS::SharedPtr<GetPot> _input;

    //! Run the test on a given input file and calculated answer
    void run_test(const std::string filename, libMesh::Real calc_answer)
    {
      this->init_sim(filename);

      _sim->run();

      CPPUNIT_ASSERT_DOUBLES_EQUAL( calc_answer, _sim->get_qoi_value(0),libMesh::TOLERANCE  );
    }

    //! Initialize the GetPot and Simulation class objects
    void init_sim(const std::string & filename)
    {
      _input.reset(new GetPot(filename));

      const char * const argv = "unit_driver";
      GetPot empty_command_line( (const int)1,&argv );
      GRINS::SimulationBuilder sim_builder;

      _sim = new GRINS::Simulation(*_input,
                                    empty_command_line,
                                    sim_builder,
                                   *TestCommWorld );
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( SpectroscopicAbsorptionTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
