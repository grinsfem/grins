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

#ifdef GRINS_HAVE_CPPUNIT

#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// Testing headers
#include "test_comm.h"
#include "grins_test_paths.h"
#include "system_helper.h"

// GRINS
#include "grins/grins_enums.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/multicomponent_variable.h"
#include "grins/single_variable.h"
#include "grins/variable_warehouse.h"
#include "grins/variable_builder.h"
#include "grins/variables_parsing.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class VariablesTest : public CppUnit::TestCase,
                        public SystemHelper
  {
  public:
    CPPUNIT_TEST_SUITE( VariablesTest );

    CPPUNIT_TEST( test_variable_builder );
    CPPUNIT_TEST( test_var_constraint );
    CPPUNIT_TEST( test_variable_arbitrary_names );
    CPPUNIT_TEST( test_variable_species_from_cantera );
    CPPUNIT_TEST( test_variable_species_from_antioch );

    CPPUNIT_TEST_SUITE_END();

  public:

    void tearDown()
    {
      this->reset_all();
    }

    void test_variable_builder()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_2d.in";
      this->setup_multiphysics_system(filename);

      GRINS::VariableBuilder::build_variables((*_input),(*_system));

      this->test_all_variables( GRINS::VariablesParsing::velocity_section(),
                                GRINS::VariablesParsing::temperature_section(),
                                GRINS::VariablesParsing::pressure_section(),
                                GRINS::VariablesParsing::single_var_section() );

      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

    void test_var_constraint()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_3d.in";
      this->setup_multiphysics_system(filename);

      GRINS::VariableBuilder::build_variables((*_input),(*_system));

      const GRINS::FEVariablesBase& vel_vars =
        GRINS::GRINSPrivate::VariableWarehouse::get_variable(GRINS::VariablesParsing::velocity_section());

      GRINS::FEVariablesBase& press_vars =
        GRINS::GRINSPrivate::VariableWarehouse::get_variable(GRINS::VariablesParsing::pressure_section());
      press_vars.set_is_constraint_var(true);

      CPPUNIT_ASSERT(!vel_vars.is_constraint_var());
      CPPUNIT_ASSERT(press_vars.is_constraint_var());

      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

    void test_variable_arbitrary_names()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_arbitrary_names.in";
      this->setup_multiphysics_system(filename);

      GRINS::VariableBuilder::build_variables((*_input),(*_system));

      this->test_all_variables( "MySpeed",
                                "TestingIsSoHot",
                                "SoMuchPressure",
                                "ForeverAlone" );

      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

    void test_variable_species_from_cantera()
    {
#ifdef GRIN_HAVE_CANTERA
      std::string filename =
        std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_species_cantera.in";

      this->setup_multiphysics_system(filename);

      GRINS::VariableBuilder::build_variables((*_input),(*_system));

      const GRINS::FEVariablesBase& species_vars =
          GRINS::GRINSPrivate::VariableWarehouse::get_variable("SpeciesMassFractions");

        const std::vector<std::string>& var_names = species_vars.active_var_names();
        CPPUNIT_ASSERT_EQUAL( std::string("Y_O"),  var_names[0] );
        CPPUNIT_ASSERT_EQUAL( std::string("Y_O2"), var_names[1] );
        CPPUNIT_ASSERT_EQUAL( std::string("Y_O3"), var_names[2] );
#endif // GRIN_HAVE_CANTERA
    }

    void test_variable_species_from_antioch()
    {
#ifdef GRIN_HAVE_ANTIOCH
      std::string filename =
        std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_species_antioch.in";

      this->setup_multiphysics_system(filename);

      GRINS::VariableBuilder::build_variables((*_input),(*_system));

      // We're also testing the arbitrary naming of species variables with this test
      const GRINS::FEVariablesBase& species_vars =
          GRINS::GRINSPrivate::VariableWarehouse::get_variable("Barf");

        const std::vector<std::string>& var_names = species_vars.active_var_names();
        CPPUNIT_ASSERT_EQUAL( std::string("Y_O"),  var_names[0] );
        CPPUNIT_ASSERT_EQUAL( std::string("Y_O2"), var_names[1] );
        CPPUNIT_ASSERT_EQUAL( std::string("Y_O3"), var_names[2] );
#endif // GRIN_HAVE_ANTIOCH
    }

  private:

    void test_all_variables( const std::string& velocity_name,
                             const std::string& temp_name,
                             const std::string& press_name,
                             const std::string& single_var_name )
    {
      // There should be 7 variables generated from that input file
      CPPUNIT_ASSERT_EQUAL((unsigned int)5,_system->n_vars());

      // Check Velocity variables
      {
        const GRINS::FEVariablesBase& vel_vars =
          GRINS::GRINSPrivate::VariableWarehouse::get_variable(velocity_name);

        const std::vector<std::string>& var_names = vel_vars.active_var_names();
        this->test_vel_var_names_2d(var_names);

        // Verify the FE part
        this->test_vel_fe_2d(*_system);
      }

      // Check Temperature variables
      {
        const GRINS::FEVariablesBase& temp_vars =
          GRINS::GRINSPrivate::VariableWarehouse::get_variable(temp_name);

        const std::vector<std::string>& var_names = temp_vars.active_var_names();
        this->test_temp_var_names(var_names);

        // Verify the FE part
        this->test_temp_fe(*_system);
      }

      // Check Pressure variable
      {
        const GRINS::FEVariablesBase& press_vars =
          GRINS::GRINSPrivate::VariableWarehouse::get_variable(press_name);

        const std::vector<std::string>& var_names = press_vars.active_var_names();
        this->test_press_var_names(var_names);

        // Verify the FE part
        this->test_press_fe(*_system);
      }

      // Check Single variable
      {
        const GRINS::FEVariablesBase& single_var =
          GRINS::GRINSPrivate::VariableWarehouse::get_variable(single_var_name);

        const std::vector<std::string>& var_names = single_var.active_var_names();
        CPPUNIT_ASSERT_EQUAL(1,(int)var_names.size());
        CPPUNIT_ASSERT_EQUAL(std::string("u"),var_names[0]);

        // Verify the FE part
        libMesh::Order order = _system->variable_type("u").order;
        CPPUNIT_ASSERT_EQUAL(GRINSEnums::LAGRANGE,_system->variable_type("u").family);
        CPPUNIT_ASSERT_EQUAL(GRINSEnums::FIRST,order);
      }
    }

    void test_vel_var_names_2d( const std::vector<std::string>& var_names )
    {
      // For 2-D, we should only have 2 components
      CPPUNIT_ASSERT_EQUAL(2,(int)var_names.size());
      CPPUNIT_ASSERT_EQUAL(std::string("Ux"),var_names[0]);
      CPPUNIT_ASSERT_EQUAL(std::string("Uy"),var_names[1]);
    }

    void test_vel_var_names_3d( const std::vector<std::string>& var_names )
    {
      // For 3-D, we should have 3 components
      CPPUNIT_ASSERT_EQUAL(3,(int)var_names.size());
      CPPUNIT_ASSERT_EQUAL(std::string("Ux"),var_names[0]);
      CPPUNIT_ASSERT_EQUAL(std::string("Uy"),var_names[1]);
      CPPUNIT_ASSERT_EQUAL(std::string("Uz"),var_names[2]);
    }

    void test_temp_var_names( const std::vector<std::string>& var_names )
    {
      CPPUNIT_ASSERT_EQUAL(1,(int)var_names.size());
      CPPUNIT_ASSERT_EQUAL(std::string("T"),var_names[0]);
    }

    void test_species_var_names( const std::vector<std::string>& var_names )
    {
      CPPUNIT_ASSERT_EQUAL(2,(int)var_names.size());
      CPPUNIT_ASSERT_EQUAL(std::string("Y_N2"),var_names[0]);
      CPPUNIT_ASSERT_EQUAL(std::string("Y_N"),var_names[1]);
    }

    void test_press_var_names( const std::vector<std::string>& var_names )
    {
      CPPUNIT_ASSERT_EQUAL(1,(int)var_names.size());
      CPPUNIT_ASSERT_EQUAL(std::string("p"),var_names[0]);
    }

    void test_species_var_names_no_order( const std::vector<std::string>& var_names )
    {
      // For this one, we can't guarantee the order, so we check to
      // make sure both the species are there.
      CPPUNIT_ASSERT_EQUAL(2,(int)var_names.size());
      CPPUNIT_ASSERT( std::find( var_names.begin(), var_names.end(),"Y_N2")
                      != var_names.end() );
      CPPUNIT_ASSERT( std::find( var_names.begin(), var_names.end(),"Y_N")
                      != var_names.end() );
    }

    void test_vel_fe_2d( const libMesh::System& system )
    {
      libMesh::Order order = system.variable_type("Ux").order;
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::LAGRANGE,system.variable_type("Ux").family);
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::FIRST,order);

      order = system.variable_type("Uy").order;
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::LAGRANGE,system.variable_type("Uy").family);
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::FIRST,order);
    }

    void test_vel_fe_3d( const libMesh::System& system )
    {
      this->test_vel_fe_2d(system);
      libMesh::Order order = system.variable_type("Uz").order;
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::LAGRANGE,system.variable_type("Uz").family);
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::FIRST,order);
    }

    void test_temp_fe( const libMesh::System& system )
    {
      libMesh::Order order = system.variable_type("T").order;
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::LAGRANGE,system.variable_type("T").family);
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::FIRST,order);
    }

    void test_press_fe( const libMesh::System& system )
    {
      libMesh::Order order = system.variable_type("p").order;
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::LAGRANGE,system.variable_type("p").family);
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::FIRST,order);
    }


    void test_species_fe( const libMesh::System& system )
    {
      libMesh::Order order = system.variable_type("Y_N2").order;
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::LAGRANGE,system.variable_type("Y_N2").family);
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::SECOND,order);

      order = system.variable_type("Y_N").order;
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::LAGRANGE,system.variable_type("Y_N").family);
      CPPUNIT_ASSERT_EQUAL(GRINSEnums::SECOND,order);
    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( VariablesTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
