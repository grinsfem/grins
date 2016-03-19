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
      {
        GRINS::VelocityFEVariables vel_vars(*_input,"PhysicsNameIsDUMMYForThisTest");
        vel_vars.init(_system);
        CPPUNIT_ASSERT_EQUAL((unsigned int)2,_system->n_vars());

        const std::vector<std::string>& var_names = vel_vars.active_var_names();
        this->test_vel_var_names_2d(var_names);

        // Verify the FE part
        this->test_vel_fe_2d(*_system);
      }

      // Now we should be able to also use a basic VelocityVariables class
      // and get the right var names out once the variables are added to the
      // system
      {
        GRINS::VelocityVariables vel_vars(*_input);
        vel_vars.init_vars(_system);

        const std::vector<std::string>& var_names = vel_vars.active_var_names();
        this->test_vel_var_names_2d(var_names);
      }
    }

    void test_velocity_3d()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_3d.in";
      this->setup_multiphysics_system(filename);

      // This will add the variables to the system
      {
        GRINS::VelocityFEVariables vel_vars(*_input,"PhysicsNameIsDUMMYForThisTest");
        vel_vars.init(_system);
        CPPUNIT_ASSERT_EQUAL((unsigned int)3,_system->n_vars());

        const std::vector<std::string>& var_names = vel_vars.active_var_names();
        this->test_vel_var_names_3d(var_names);

        // Verify the FE part
        this->test_vel_fe_3d(*_system);
      }

      // Now we should be able to also use a basic VelocityVariables class
      // and get the right var names out once the variables are added to the
      // system
      {
        GRINS::VelocityVariables vel_vars(*_input);
        vel_vars.init_vars(_system);

        const std::vector<std::string>& var_names = vel_vars.active_var_names();
        this->test_vel_var_names_3d(var_names);
      }
    }

    void test_temp()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/variables_2d.in";
      this->setup_multiphysics_system(filename);

      // This will add the variables to the system
      {
        GRINS::PrimitiveTempFEVariables temp_vars(*_input,"PhysicsNameIsDUMMYForThisTest");
        temp_vars.init(_system);
        CPPUNIT_ASSERT_EQUAL((unsigned int)1,_system->n_vars());

        const std::vector<std::string>& var_names = temp_vars.active_var_names();
        this->test_temp_var_names(var_names);

        // Verify the FE part
        this->test_temp_fe(*_system);
      }

      // Now we should be able to also use a basic PrimitiveTempVariables class
      // and get the right var names out once the variables are added to the
      // system
      {
        GRINS::PrimitiveTempVariables temp_vars(*_input);
        temp_vars.init_vars(_system);

        const std::vector<std::string>& var_names = temp_vars.active_var_names();
        this->test_temp_var_names(var_names);
      }
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

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( VariablesTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
