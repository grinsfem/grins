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

#include "grins_config.h"

#ifdef GRINS_HAVE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

// Testing headers
#include "test_comm.h"
#include "grins_test_paths.h"
#include "system_helper.h"

// GRINS
#include "grins/grins_enums.h"
#include "grins/velocity_fe_variables.h"
#include "grins/primitive_temp_fe_variables.h"
#include "grins/species_mass_fracs_fe_variables.h"
#include "grins/pressure_fe_variable.h"
#include "grins/variable_warehouse.h"

namespace GRINSTesting
{
  class VariablesTest : public CppUnit::TestCase,
                        public SystemHelper
  {
  public:
    CPPUNIT_TEST_SUITE( VariablesTest );

    CPPUNIT_TEST( test_velocity_2d );
    CPPUNIT_TEST( test_velocity_3d );
    CPPUNIT_TEST( test_temp );
    CPPUNIT_TEST( test_species_mass_fracs );
    CPPUNIT_TEST( test_pressure );
    CPPUNIT_TEST( test_var_constraint_and_warehouse );

    CPPUNIT_TEST_SUITE_END();

  public:

    void tearDown()
    {
      this->reset_all();
    }

    void test_velocity_2d()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_2d.in";
      this->setup_multiphysics_system(filename);

      // This will add the variables to the system
      GRINS::VelocityFEVariables vel_vars(*_input,"PhysicsNameIsDUMMYForThisTest");
      vel_vars.init(_system);
      CPPUNIT_ASSERT_EQUAL((unsigned int)2,_system->n_vars());

      const std::vector<std::string>& var_names = vel_vars.active_var_names();
      this->test_vel_var_names_2d(var_names);

      // Verify the FE part
      this->test_vel_fe_2d(*_system);
    }

    void test_velocity_3d()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_3d.in";
      this->setup_multiphysics_system(filename);

      // This will add the variables to the system
      GRINS::VelocityFEVariables vel_vars(*_input,"PhysicsNameIsDUMMYForThisTest");
      vel_vars.init(_system);
      CPPUNIT_ASSERT_EQUAL((unsigned int)3,_system->n_vars());

      const std::vector<std::string>& var_names = vel_vars.active_var_names();
      this->test_vel_var_names_3d(var_names);

      // Verify the FE part
      this->test_vel_fe_3d(*_system);
    }

    void test_temp()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_2d.in";
      this->setup_multiphysics_system(filename);

      // This will add the variables to the system
      GRINS::PrimitiveTempFEVariables temp_vars(*_input,"PhysicsNameIsDUMMYForThisTest");
      temp_vars.init(_system);
      CPPUNIT_ASSERT_EQUAL((unsigned int)1,_system->n_vars());

      const std::vector<std::string>& var_names = temp_vars.active_var_names();
      this->test_temp_var_names(var_names);

      // Verify the FE part
      this->test_temp_fe(*_system);
    }

    void test_species_mass_fracs()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_2d.in";
      this->setup_multiphysics_system(filename);

      // This will add the variables to the system
      GRINS::SpeciesMassFractionsFEVariables species_vars(*_input,"TestSpeciesMassFractionsVariables");
      species_vars.init(_system);
      CPPUNIT_ASSERT_EQUAL((unsigned int)2,_system->n_vars());

      const std::vector<std::string>& var_names = species_vars.active_var_names();
      this->test_species_var_names(var_names);

      // Verify the FE part
      this->test_species_fe(*_system);
    }

    void test_var_constraint_and_warehouse()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_3d.in";
      this->setup_multiphysics_system(filename);

      GRINS::VelocityFEVariables vel_vars(*_input,"PhysicsNameIsDUMMYForThisTest");
      GRINS::PressureFEVariable press_vars(*_input,"PhysicsNameIsDUMMYForThisTest",true /*is_constraint_variable*/ );

      GRINS::GRINSPrivate::VariableWarehouse::register_variable("Velocity", vel_vars);
      GRINS::GRINSPrivate::VariableWarehouse::register_variable("Pressure", press_vars);

      // This should not throw an error
      GRINS::GRINSPrivate::VariableWarehouse::check_and_register_variable("Velocity", vel_vars);
      GRINS::GRINSPrivate::VariableWarehouse::check_and_register_variable("Pressure", press_vars);

      vel_vars.init(_system);
      press_vars.init(_system);

      {
        const GRINS::FEVariablesBase& vel_base =
          GRINS::GRINSPrivate::VariableWarehouse::get_variable("Velocity");
        CPPUNIT_ASSERT(!vel_base.is_constraint_var());

        const GRINS::FEVariablesBase& press_base =
          GRINS::GRINSPrivate::VariableWarehouse::get_variable("Pressure");
        CPPUNIT_ASSERT(press_base.is_constraint_var());
      }

      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

    void test_pressure()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_2d.in";
      this->setup_multiphysics_system(filename);

      // This will add the variables to the system
      GRINS::PressureFEVariable press_vars(*_input,"PhysicsNameIsDUMMYForThisTest");
      press_vars.init(_system);
      CPPUNIT_ASSERT_EQUAL((unsigned int)1,_system->n_vars());

      const std::vector<std::string>& var_names = press_vars.active_var_names();
      this->test_press_var_names(var_names);

      // Verify the FE part
      this->test_press_fe(*_system);
    }

  private:

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
