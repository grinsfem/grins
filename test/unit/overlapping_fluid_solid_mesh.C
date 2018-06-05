//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

#include "test_comm.h"
#include "grins_test_paths.h"
#include "system_helper.h"

// GRINS
#include "grins/variable_builder.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/variable_warehouse.h"
#include "grins/variable_builder.h"
#include "grins/variables_parsing.h"
#include "grins/overlapping_fluid_solid_map.h"

// libMesh
#include "libmesh/steady_solver.h"
#include "libmesh/elem.h"

namespace GRINSTesting
{

  class OverlappingFluidSolidMeshTest : public CppUnit::TestCase,
                                        SystemHelper
  {
  public:
    CPPUNIT_TEST_SUITE( OverlappingFluidSolidMeshTest );

    CPPUNIT_TEST( one_overlapping_element_test );

    CPPUNIT_TEST_SUITE_END();

  public:

    void tearDown()
    {
      this->reset_all();
    }

    void one_overlapping_element_test()
    {
      std::string input_file = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/one_elem_overlap.in";
      this->setup_multiphysics_system(input_file);

      // This isn't used internally in the OverlappingFluidSolidMap, but the FEMContext
      // therein wants a time solver to be defined.
      libMesh::SteadySolver* time_solver = new libMesh::SteadySolver( (*_system) );
      _system->time_solver = libMesh::UniquePtr<libMesh::TimeSolver>(time_solver);

      std::set<libMesh::subdomain_id_type> solid_ids;
      solid_ids.insert(1);

      std::set<libMesh::subdomain_id_type> fluid_ids;
      fluid_ids.insert(2);

      GRINS::VariableBuilder::build_variables((*_input),(*_system));

      const GRINS::DisplacementVariable & disp_vars =
        GRINS::GRINSPrivate::VariableWarehouse::get_variable_subclass<GRINS::DisplacementVariable>
        (GRINS::VariablesParsing::displacement_section());

      _es->init();

      libMesh::UniquePtr<libMesh::PointLocatorBase> point_locator = _mesh->sub_point_locator();

      GRINS::OverlappingFluidSolidMap mesh_overlap( (*_system), (*point_locator), solid_ids, fluid_ids, disp_vars );

      /*
        Overlapping element in the interior: id = 0, node 0 = (0.5,0.5), length is 1 on each edge
        Right interior element: id = 1, node 0 = (1,1), width = 1,  height is 2
        Left interior element: id = 2, node 0 = (0,1), width = 1,  height is 2

        (-1,1) ------------------------ (1,1)
               |          |          |
               |          |          |
               |      ---------      |
               |      |       |      |
               |      |       |      |
               |      ---------      |
               |          |          |
               |          |          |
       (-1,-1) ------------------------ (-1,1)

       */
      for( libMesh::MeshBase::element_iterator e = (*_mesh).active_local_elements_begin();
           e != (*_mesh).active_local_elements_end();
           ++e )
        {
          // Convenience
          const libMesh::Elem * elem = *e;
          elem->print_info();
        }


      // There should be 1 "solid" element mapping to other fluid two elements
      CPPUNIT_ASSERT_EQUAL(1,(int)mesh_overlap.solid_map().size());
      CPPUNIT_ASSERT_EQUAL(2,(int)(mesh_overlap.solid_map().find(0)->second).size());

      // These are a little silly. Morally equivalent to map.find() != map.end()
      CPPUNIT_ASSERT_EQUAL(1,(int)(mesh_overlap.solid_map().find(0)->second).find(1)->first);
      CPPUNIT_ASSERT_EQUAL(2,(int)(mesh_overlap.solid_map().find(0)->second).find(2)->first);

      // Check for the correct number and values of solid quadrature point indices
      // associated with each of the fluid elements
      CPPUNIT_ASSERT_EQUAL(6,(int)((mesh_overlap.solid_map().find(0)->second).find(1)->second).size());
      CPPUNIT_ASSERT_EQUAL(3,(int)((mesh_overlap.solid_map().find(0)->second).find(2)->second).size());

      std::vector<unsigned int> elem1_qps(6);
      elem1_qps[0] = 0; elem1_qps[1] = 1; elem1_qps[2] = 3; elem1_qps[3] = 4; elem1_qps[4] = 6; elem1_qps[5] = 7;

      std::vector<unsigned int> elem2_qps(3);
      elem2_qps[0] = 2; elem2_qps[1] = 5; elem2_qps[2] = 8;

      for( unsigned int i = 0; i < 6; i++ )
        CPPUNIT_ASSERT_EQUAL(elem1_qps[i],((mesh_overlap.solid_map().find(0)->second).find(1)->second)[i]);

      for( unsigned int i = 0; i < 3; i++ )
        CPPUNIT_ASSERT_EQUAL(elem2_qps[i],((mesh_overlap.solid_map().find(0)->second).find(2)->second)[i]);

      // There should be 2 "fluid" elements, each mapping to the solid element
      CPPUNIT_ASSERT_EQUAL(2,(int)mesh_overlap.fluid_map().size());
      CPPUNIT_ASSERT_EQUAL(1,(int)(mesh_overlap.fluid_map().find(1)->second).size());
      CPPUNIT_ASSERT_EQUAL(1,(int)(mesh_overlap.fluid_map().find(2)->second).size());
      CPPUNIT_ASSERT_EQUAL(0,(int)(mesh_overlap.fluid_map().find(1)->second).find(0)->first);
      CPPUNIT_ASSERT_EQUAL(0,(int)(mesh_overlap.fluid_map().find(2)->second).find(0)->first);
      CPPUNIT_ASSERT_EQUAL(6,(int)((mesh_overlap.fluid_map().find(1)->second).find(0)->second).size());
      CPPUNIT_ASSERT_EQUAL(3,(int)((mesh_overlap.fluid_map().find(2)->second).find(0)->second).size());

      for( unsigned int i = 0; i < 6; i++ )
        CPPUNIT_ASSERT_EQUAL(elem1_qps[i],((mesh_overlap.fluid_map().find(1)->second).find(0)->second)[i]);

      for( unsigned int i = 0; i < 3; i++ )
        CPPUNIT_ASSERT_EQUAL(elem2_qps[i],((mesh_overlap.fluid_map().find(2)->second).find(0)->second)[i]);

      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

  private:


  };

  CPPUNIT_TEST_SUITE_REGISTRATION( OverlappingFluidSolidMeshTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
