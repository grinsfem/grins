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
#include "grins/default_bc_builder.h"
#include "grins/boundary_condition_names.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class DefaultBCBuilderTest : public CppUnit::TestCase,
                               public SystemHelper,
                               public GRINS::DefaultBCBuilder // So we can test proctected methods
  {
  public:
    CPPUNIT_TEST_SUITE( DefaultBCBuilderTest );

    CPPUNIT_TEST( test_parse_and_build_bc_id_map );
    CPPUNIT_TEST( test_verify_bc_ids_with_mesh );
    CPPUNIT_TEST( test_parse_periodic_master_slave_ids );
    CPPUNIT_TEST( test_parse_periodic_offset );

    CPPUNIT_TEST_SUITE_END();

  public:

    void tearDown()
    {
      this->reset_all();
    }

    void test_parse_and_build_bc_id_map()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/default_bc_builder.in";
      GetPot input(filename);

      std::map<std::string,std::set<GRINS::BoundaryID> > bc_id_map;
      this->parse_and_build_bc_id_map(input,bc_id_map);

      // Make sure we have the right sections
      CPPUNIT_ASSERT_EQUAL(3,(int)bc_id_map.size());
      CPPUNIT_ASSERT( bc_id_map.find("Hot") != bc_id_map.end() );
      CPPUNIT_ASSERT( bc_id_map.find("Together") != bc_id_map.end() );
      CPPUNIT_ASSERT( bc_id_map.find("Cold") != bc_id_map.end() );

      // Make sure we have the right number and values of the bc ids
      {
        std::set<GRINS::BoundaryID> bc_ids = bc_id_map["Hot"];
        CPPUNIT_ASSERT_EQUAL(1,(int)bc_ids.size());
        CPPUNIT_ASSERT(bc_ids.find(0) != bc_ids.end());
      }

      // Make sure we have the right number and values of the bc ids
      {
        std::set<GRINS::BoundaryID> bc_ids = bc_id_map["Together"];
        CPPUNIT_ASSERT_EQUAL(2,(int)bc_ids.size());
        CPPUNIT_ASSERT(bc_ids.find(1) != bc_ids.end());
        CPPUNIT_ASSERT(bc_ids.find(2) != bc_ids.end());
      }

      // Make sure we have the right number and values of the bc ids
      {
        std::set<GRINS::BoundaryID> bc_ids = bc_id_map["Cold"];
        CPPUNIT_ASSERT_EQUAL(1,(int)bc_ids.size());
        CPPUNIT_ASSERT(bc_ids.find(3) != bc_ids.end());
      }
    }

    void test_verify_bc_ids_with_mesh()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/default_bc_builder.in";
      this->setup_multiphysics_system(filename);

      std::map<std::string,std::set<GRINS::BoundaryID> > bc_id_map;
      this->parse_and_build_bc_id_map(*_input,bc_id_map);

      // This shouldn't error
      this->verify_bc_ids_with_mesh( *_system, bc_id_map );
    }



    void test_parse_periodic_master_slave_ids()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/default_bc_builder.in";
      this->setup_multiphysics_system(filename);

      libMesh::boundary_id_type invalid_bid =
        std::numeric_limits<libMesh::boundary_id_type>::max();

      libMesh::boundary_id_type master_id = invalid_bid;
      libMesh::boundary_id_type slave_id = invalid_bid;

      std::string section = GRINS::BoundaryConditionNames::bc_section()+"/Together";

      this->parse_periodic_master_slave_ids(*_input,section,master_id,slave_id);

      CPPUNIT_ASSERT_EQUAL(1,(int)master_id);
      CPPUNIT_ASSERT_EQUAL(2,(int)slave_id);
    }

    void test_parse_periodic_offset()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/default_bc_builder.in";
      this->setup_multiphysics_system(filename);

      std::string section = GRINS::BoundaryConditionNames::bc_section()+"/Together";

      libMesh::RealVectorValue offset =
        this->parse_periodic_offset(*_input,section);

      libMesh::Real tol = std::numeric_limits<libMesh::Real>::epsilon()*10;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.21,offset(0),tol);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0,offset(1),tol);
    }

  private:

    void test_for_var_name( const std::vector<std::string>& var_names,
                            const std::string& var_to_find )
    {
      CPPUNIT_ASSERT( std::find( var_names.begin(), var_names.end(), var_to_find ) != var_names.end() );
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( DefaultBCBuilderTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
