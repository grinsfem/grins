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

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include "test_comm.h"
#include "grins_test_paths.h"
#include "system_helper.h"

// GRINS
#include "grins/mesh_builder.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/variable_warehouse.h"
#include "grins/variable_builder.h"
#include "grins/variables_parsing.h"
#include "grins/overlapping_fluid_solid_map.h"

// libMesh
#include "libmesh/steady_solver.h"
#include "libmesh/elem.h"
#include "libmesh/face_quad9.h"
#include "libmesh/face_tri6.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/serial_mesh.h"

#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{

  class OverlappingFluidSolidMeshTest : public CppUnit::TestCase,
                                        SystemHelper
  {
  public:
    CPPUNIT_TEST_SUITE( OverlappingFluidSolidMeshTest );

    CPPUNIT_TEST( one_overlapping_element_test );
    CPPUNIT_TEST( quad_on_quad_overlapping_test );
    CPPUNIT_TEST( tri_on_quad_overlapping_test );
    CPPUNIT_TEST( quad_on_tri_overlapping_test );
    CPPUNIT_TEST( tri_on_tri_overlapping_test );

    CPPUNIT_TEST_SUITE_END();

  public:

    void tearDown()
    {
      this->reset_all();
    }

    //---------------------------------------------QUAD-ON-QUAD-TEST(.exo)------------------------------------------------- 

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

      // There should be 1 "solid" element mapping to other two fluid elements
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

    //---------------------------------------------QUAD-ON-QUAD-TEST-------------------------------------------------

    void quad_on_quad_overlapping_test()
    {
      std::string input_file = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/overlap_variables.in";
      
      std::unique_ptr<GetPot> _input;
      _input.reset( new GetPot(input_file) );

      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

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

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.5,0.5),0 );
      mesh->add_point( libMesh::Point(-0.5,0.5),1 );
      mesh->add_point( libMesh::Point(-0.5,-0.5),2 );
      mesh->add_point( libMesh::Point(0.5,-0.5),3 );
      mesh->add_point( libMesh::Point(0.0,0.5),4 );
      mesh->add_point( libMesh::Point(-0.5,0.0),5 );
      mesh->add_point( libMesh::Point(0.0,-0.5),6 );
      mesh->add_point( libMesh::Point(0.5,0.0),7 );
      mesh->add_point( libMesh::Point(0.0,0.0),8 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Quad9 );
      elem->set_id(0);
      elem->subdomain_id() = 1;
      for (unsigned int n=0; n<9; n++)
        elem->set_node(n) = mesh->node_ptr(n);
      
      mesh->add_point( libMesh::Point(1.0,1.0),9 );
      mesh->add_point( libMesh::Point(0.0,1.0),10 );
      mesh->add_point( libMesh::Point(0.0,-1.0),11 );
      mesh->add_point( libMesh::Point(1.0,-1.0),12 );
      mesh->add_point( libMesh::Point(0.5,1.0),13 );
      mesh->add_point( libMesh::Point(0.0,0.0),14 );
      mesh->add_point( libMesh::Point(0.5,-1.0),15 );
      mesh->add_point( libMesh::Point(1.0,0.0),16 );
      mesh->add_point( libMesh::Point(0.5,0.0),17 );

      elem = mesh->add_elem( new libMesh::Quad9 );
      elem->set_id(1);
      elem->subdomain_id() = 2;
      elem->set_node(0) = mesh->node_ptr(9);
      elem->set_node(1) = mesh->node_ptr(10);
      elem->set_node(2) = mesh->node_ptr(11);
      elem->set_node(3) = mesh->node_ptr(12);
      elem->set_node(4) = mesh->node_ptr(13);
      elem->set_node(5) = mesh->node_ptr(14);
      elem->set_node(6) = mesh->node_ptr(15);
      elem->set_node(7) = mesh->node_ptr(16);
      elem->set_node(8) = mesh->node_ptr(17);
      
      //mesh->add_point( libMesh::Point(0.0,1.0),10 );
      mesh->add_point( libMesh::Point(-1.0,1.0),18 );
      mesh->add_point( libMesh::Point(-1.0,-1.0),19 );
      //mesh->add_point( libMesh::Point(0.0,-1.0),11 );
      mesh->add_point( libMesh::Point(-0.5,1.0),20 );
      mesh->add_point( libMesh::Point(-1.0,0.0),21 );
      mesh->add_point( libMesh::Point(-0.5,-1.0),22 );
      //mesh->add_point( libMesh::Point(0.0,0.0),14 );
      mesh->add_point( libMesh::Point(-0.5,0.0),23 );

      elem = mesh->add_elem( new libMesh::Quad9 );
      elem->set_id(2);
      elem->subdomain_id() = 2;
      elem->set_node(0) = mesh->node_ptr(10);
      elem->set_node(1) = mesh->node_ptr(18);
      elem->set_node(2) = mesh->node_ptr(19);
      elem->set_node(3) = mesh->node_ptr(11);
      elem->set_node(4) = mesh->node_ptr(20);
      elem->set_node(5) = mesh->node_ptr(21);
      elem->set_node(6) = mesh->node_ptr(22);
      elem->set_node(7) = mesh->node_ptr(14);
      elem->set_node(8) = mesh->node_ptr(23);
      
      mesh->prepare_for_use();

      std::unique_ptr<libMesh::EquationSystems> _es;
      _es.reset( new libMesh::EquationSystems(*mesh) );

      GRINS::MultiphysicsSystem* _system;
      _system = &_es->add_system<GRINS::MultiphysicsSystem>( "GRINS-TEST" );
      _system->read_input_options( (*_input) );

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
      
      libMesh::UniquePtr<libMesh::PointLocatorBase> point_locator = mesh->sub_point_locator();

      GRINS::OverlappingFluidSolidMap mesh_overlap( (*_system), (*point_locator), solid_ids, fluid_ids, disp_vars );

      CPPUNIT_ASSERT_EQUAL(1,(int)mesh_overlap.solid_map().size());
      CPPUNIT_ASSERT_EQUAL(2,(int)(mesh_overlap.solid_map().find(0)->second).size());

      CPPUNIT_ASSERT_EQUAL(1,(int)(mesh_overlap.solid_map().find(0)->second).find(1)->first);
      CPPUNIT_ASSERT_EQUAL(2,(int)(mesh_overlap.solid_map().find(0)->second).find(2)->first);

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

    //---------------------------------------------TRI-ON-QUAD-TEST-------------------------------------------------

    void tri_on_quad_overlapping_test()
    {
      std::string input_file = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/overlap_variables.in";
      
      std::unique_ptr<GetPot> _input;
      _input.reset( new GetPot(input_file) );

      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      /*
        Overlapping element in the interior (TRI6): id = 0, node 0 = (0.0,0.5)
        Right interior element (QUAD9): id = 1, node 0 = (1,1), width = 1,  height is 2
        Left interior element (QUAD9): id = 2, node 0 = (0,1), width = 1,  height is 2
     
        (-1,1) ----------------------- (1,1)
               |          |          |
               |          |          |
               |         / \         |
               |        /   \        |
               |       /     \       |
               |      ---------      |
               |          |          |
               |          |          |
       (-1,-1) ----------------------- (-1,1)

      */

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.5),0 );
      mesh->add_point( libMesh::Point(-0.5,-0.5),1 );
      mesh->add_point( libMesh::Point(0.5,-0.5),2 );
      mesh->add_point( libMesh::Point(-0.25,0.0),3 );
      mesh->add_point( libMesh::Point(0.0,-0.5),4 );
      mesh->add_point( libMesh::Point(0.25,0.0),5 );
      
      libMesh::Elem* elem = mesh->add_elem( new libMesh::Tri6 );
      elem->set_id(0);
      elem->subdomain_id() = 1;
      for (unsigned int n=0; n<6; n++)
        elem->set_node(n) = mesh->node_ptr(n);
      
      mesh->add_point( libMesh::Point(1.0,1.0),6 );
      mesh->add_point( libMesh::Point(0.0,1.0),7 );
      mesh->add_point( libMesh::Point(0.0,-1.0),8 );
      mesh->add_point( libMesh::Point(1.0,-1.0),9 );
      mesh->add_point( libMesh::Point(0.5,1.0),10 );
      mesh->add_point( libMesh::Point(0.0,0.0),11 );
      mesh->add_point( libMesh::Point(0.5,-1.0),12 );
      mesh->add_point( libMesh::Point(1.0,0.0),13 );
      mesh->add_point( libMesh::Point(0.5,0.0),14 );

      elem = mesh->add_elem( new libMesh::Quad9 );
      elem->set_id(1);
      elem->subdomain_id() = 2;
      elem->set_node(0) = mesh->node_ptr(6);
      elem->set_node(1) = mesh->node_ptr(7);
      elem->set_node(2) = mesh->node_ptr(8);
      elem->set_node(3) = mesh->node_ptr(9);
      elem->set_node(4) = mesh->node_ptr(10);
      elem->set_node(5) = mesh->node_ptr(11);
      elem->set_node(6) = mesh->node_ptr(12);
      elem->set_node(7) = mesh->node_ptr(13);
      elem->set_node(8) = mesh->node_ptr(14);
 
      //mesh->add_point( libMesh::Point(0.0,1.0),7 );
      mesh->add_point( libMesh::Point(-1.0,1.0),15 );
      mesh->add_point( libMesh::Point(-1.0,-1.0),16 );
      //mesh->add_point( libMesh::Point(0.0,-1.0),8 );
      mesh->add_point( libMesh::Point(-0.5,1.0),17 );
      mesh->add_point( libMesh::Point(-1.0,0.0),18 );
      mesh->add_point( libMesh::Point(-0.5,-1.0),19 );
      //mesh->add_point( libMesh::Point(0.0,0.0),11 );
      mesh->add_point( libMesh::Point(-0.5,0.0),20 );

      elem = mesh->add_elem( new libMesh::Quad9 );
      elem->set_id(2);
      elem->subdomain_id() = 2;
      elem->set_node(0) = mesh->node_ptr(7);
      elem->set_node(1) = mesh->node_ptr(15);
      elem->set_node(2) = mesh->node_ptr(16);
      elem->set_node(3) = mesh->node_ptr(8);
      elem->set_node(4) = mesh->node_ptr(17);
      elem->set_node(5) = mesh->node_ptr(18);
      elem->set_node(6) = mesh->node_ptr(19);
      elem->set_node(7) = mesh->node_ptr(11);
      elem->set_node(8) = mesh->node_ptr(20);
      
      mesh->prepare_for_use();
    
      std::unique_ptr<libMesh::EquationSystems> _es;
      _es.reset( new libMesh::EquationSystems(*mesh) );

      GRINS::MultiphysicsSystem* _system;
      _system = &_es->add_system<GRINS::MultiphysicsSystem>( "GRINS-TEST" );
      _system->read_input_options( (*_input) );

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

      libMesh::UniquePtr<libMesh::PointLocatorBase> point_locator = mesh->sub_point_locator();

      GRINS::OverlappingFluidSolidMap mesh_overlap( (*_system), (*point_locator), solid_ids, fluid_ids, disp_vars );

      CPPUNIT_ASSERT_EQUAL(1,(int)mesh_overlap.solid_map().size());
      CPPUNIT_ASSERT_EQUAL(2,(int)(mesh_overlap.solid_map().find(0)->second).size());

      CPPUNIT_ASSERT_EQUAL(1,(int)(mesh_overlap.solid_map().find(0)->second).find(1)->first);
      CPPUNIT_ASSERT_EQUAL(2,(int)(mesh_overlap.solid_map().find(0)->second).find(2)->first);

      CPPUNIT_ASSERT_EQUAL(5,(int)((mesh_overlap.solid_map().find(0)->second).find(1)->second).size());
      CPPUNIT_ASSERT_EQUAL(2,(int)((mesh_overlap.solid_map().find(0)->second).find(2)->second).size());

      std::vector<unsigned int> elem1_qps(5);
      elem1_qps[0] = 0; elem1_qps[1] = 1; elem1_qps[2] = 3; elem1_qps[3] = 4; elem1_qps[4] = 5; 

      std::vector<unsigned int> elem2_qps(2);
      elem2_qps[0] = 2; elem2_qps[1] = 6; 

      for( unsigned int i = 0; i < 5; i++ )
        CPPUNIT_ASSERT_EQUAL(elem1_qps[i],((mesh_overlap.solid_map().find(0)->second).find(1)->second)[i]);

      for( unsigned int i = 0; i < 2; i++ )
        CPPUNIT_ASSERT_EQUAL(elem2_qps[i],((mesh_overlap.solid_map().find(0)->second).find(2)->second)[i]);


      CPPUNIT_ASSERT_EQUAL(2,(int)mesh_overlap.fluid_map().size());
      CPPUNIT_ASSERT_EQUAL(1,(int)(mesh_overlap.fluid_map().find(1)->second).size());
      CPPUNIT_ASSERT_EQUAL(1,(int)(mesh_overlap.fluid_map().find(2)->second).size());
      CPPUNIT_ASSERT_EQUAL(0,(int)(mesh_overlap.fluid_map().find(1)->second).find(0)->first);
      CPPUNIT_ASSERT_EQUAL(0,(int)(mesh_overlap.fluid_map().find(2)->second).find(0)->first);
      CPPUNIT_ASSERT_EQUAL(5,(int)((mesh_overlap.fluid_map().find(1)->second).find(0)->second).size());
      CPPUNIT_ASSERT_EQUAL(2,(int)((mesh_overlap.fluid_map().find(2)->second).find(0)->second).size());

      for( unsigned int i = 0; i < 5; i++ )
        CPPUNIT_ASSERT_EQUAL(elem1_qps[i],((mesh_overlap.fluid_map().find(1)->second).find(0)->second)[i]);

      for( unsigned int i = 0; i < 2; i++ )
        CPPUNIT_ASSERT_EQUAL(elem2_qps[i],((mesh_overlap.fluid_map().find(2)->second).find(0)->second)[i]);

      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

    //---------------------------------------------QUAD-ON-TRI-TEST-------------------------------------------------
 
    void quad_on_tri_overlapping_test()
    {
      std::string input_file = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/overlap_variables.in";
      
      std::unique_ptr<GetPot> _input;
      _input.reset( new GetPot(input_file) );

      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );
      
      /*
        Overlapping element in the interior(QUAD9): id = 0, node 0 = (0.5,0.5), length is 1 on each edge
        Right interior element(TRI6): id = 1, node 0 = (1,1) 
        Left interior element(TRI6): id = 2, node 0 = (1,1)
        
        (-1,1) ---------------------- (1,1)
               |                   /|
               |                  / |
               |                 /  |
               |                /   |
               |               /    |
               |       _______/     |
               |      |       |     |
               |      |       |     |
               |      |       |     |
               |      |_______|     |
               |     /              |
               |    /               |
               |   /                |
               |  /                 |
               | /                  |
               |/                   |
       (-1,-1) ---------------------- (-1,1)
      
      */
      
      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.5,0.5),0 );
      mesh->add_point( libMesh::Point(-0.5,0.5),1 );
      mesh->add_point( libMesh::Point(-0.5,-0.5),2 );
      mesh->add_point( libMesh::Point(0.5,-0.5),3 );
      mesh->add_point( libMesh::Point(0.0,0.5),4 );
      mesh->add_point( libMesh::Point(-0.5,0.0),5 );
      mesh->add_point( libMesh::Point(0.0,-0.5),6 );
      mesh->add_point( libMesh::Point(0.5,0.0),7 );
      mesh->add_point( libMesh::Point(0.0,0.0),8 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Quad9 );
      elem->set_id(0);
      elem->subdomain_id() = 1;
      for (unsigned int n=0; n<9; n++)
        elem->set_node(n) = mesh->node_ptr(n);
      
      mesh->add_point( libMesh::Point(1.0,1.0),9 );
      mesh->add_point( libMesh::Point(-1.0,-1.0),10 );
      mesh->add_point( libMesh::Point(1.0,-1.0),11 );
      mesh->add_point( libMesh::Point(0.0,0.0),12 );
      mesh->add_point( libMesh::Point(0.0,-1.0),13 );
      mesh->add_point( libMesh::Point(1.0,0.0),14 );
      
      elem = mesh->add_elem( new libMesh::Tri6 );
      elem->set_id(1);
      elem->subdomain_id() = 2;
      elem->set_node(0) = mesh->node_ptr(9);
      elem->set_node(1) = mesh->node_ptr(10);
      elem->set_node(2) = mesh->node_ptr(11);
      elem->set_node(3) = mesh->node_ptr(12);
      elem->set_node(4) = mesh->node_ptr(13);
      elem->set_node(5) = mesh->node_ptr(14);
            
      //mesh->add_point( libMesh::Point(1.0,1.0),9 );
      mesh->add_point( libMesh::Point(-1.0,1.0),15 );
      //mesh->add_point( libMesh::Point(-1.0,-1.0),10 );
      mesh->add_point( libMesh::Point(0.0,1.0),16 );
      mesh->add_point( libMesh::Point(-1.0,0.0),17 );
      //mesh->add_point( libMesh::Point(0.0,0.0),12 );
      
      elem = mesh->add_elem( new libMesh::Tri6 );
      elem->set_id(2);
      elem->subdomain_id() = 2;
      elem->set_node(0) = mesh->node_ptr(9);
      elem->set_node(1) = mesh->node_ptr(15);
      elem->set_node(2) = mesh->node_ptr(10);
      elem->set_node(3) = mesh->node_ptr(16);
      elem->set_node(4) = mesh->node_ptr(17);
      elem->set_node(5) = mesh->node_ptr(12);
            
      mesh->prepare_for_use();
      
      std::unique_ptr<libMesh::EquationSystems> _es;
      _es.reset( new libMesh::EquationSystems(*mesh) );

      GRINS::MultiphysicsSystem* _system;
      _system = &_es->add_system<GRINS::MultiphysicsSystem>( "GRINS-TEST" );
      _system->read_input_options( (*_input) );

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

      libMesh::UniquePtr<libMesh::PointLocatorBase> point_locator = mesh->sub_point_locator();

      GRINS::OverlappingFluidSolidMap mesh_overlap( (*_system), (*point_locator), solid_ids, fluid_ids, disp_vars );

      CPPUNIT_ASSERT_EQUAL(1,(int)mesh_overlap.solid_map().size());
      CPPUNIT_ASSERT_EQUAL(2,(int)(mesh_overlap.solid_map().find(0)->second).size());

      CPPUNIT_ASSERT_EQUAL(1,(int)(mesh_overlap.solid_map().find(0)->second).find(1)->first);
      CPPUNIT_ASSERT_EQUAL(2,(int)(mesh_overlap.solid_map().find(0)->second).find(2)->first);

      CPPUNIT_ASSERT_EQUAL(6,(int)((mesh_overlap.solid_map().find(0)->second).find(1)->second).size());
      CPPUNIT_ASSERT_EQUAL(3,(int)((mesh_overlap.solid_map().find(0)->second).find(2)->second).size());

      std::vector<unsigned int> elem1_qps(6);
      elem1_qps[0] = 0; elem1_qps[1] = 3; elem1_qps[2] = 4; elem1_qps[3] = 6; elem1_qps[4] = 7; elem1_qps[5] = 8;

      std::vector<unsigned int> elem2_qps(3);
      elem2_qps[0] = 1; elem2_qps[1] = 2; elem2_qps[2] = 5;

      for( unsigned int i = 0; i < 6; i++ )
        CPPUNIT_ASSERT_EQUAL(elem1_qps[i],((mesh_overlap.solid_map().find(0)->second).find(1)->second)[i]);

      for( unsigned int i = 0; i < 3; i++ )
        CPPUNIT_ASSERT_EQUAL(elem2_qps[i],((mesh_overlap.solid_map().find(0)->second).find(2)->second)[i]);


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
    
    //---------------------------------------------TRI-ON-TRI-TEST-------------------------------------------------
   
    void tri_on_tri_overlapping_test()
    {
      std::string input_file = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/overlap_variables.in";
      
      std::unique_ptr<GetPot> _input;
      _input.reset( new GetPot(input_file) );

      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      /*
        Overlapping element in the interior(TRI6): id = 0, node 0 = (0.0,0.5)
        Right interior element(TRI6): id = 1, node 0 = (1,1) 
        Left interior element(TRI6): id = 2, node 0 = (1,1)
            
      */

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.5),0 );
      mesh->add_point( libMesh::Point(-0.5,-0.5),1 );
      mesh->add_point( libMesh::Point(0.5,-0.5),2 );
      mesh->add_point( libMesh::Point(-0.25,0.0),3 );
      mesh->add_point( libMesh::Point(0.0,-0.5),4 );
      mesh->add_point( libMesh::Point(0.25,0.0),5 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Tri6 );
      elem->set_id(0);
      elem->subdomain_id() = 1;
      for (unsigned int n=0; n<6; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->add_point( libMesh::Point(1.0,1.0),6 );
      mesh->add_point( libMesh::Point(-1.0,-1.0),7 );
      mesh->add_point( libMesh::Point(1.0,-1.0),8 );
      mesh->add_point( libMesh::Point(0.0,0.0),9 );
      mesh->add_point( libMesh::Point(0.0,-1.0),10 );
      mesh->add_point( libMesh::Point(1.0,0.0),11 );

      elem = mesh->add_elem( new libMesh::Tri6 );
      elem->set_id(1);
      elem->subdomain_id() = 2;
      elem->set_node(0) = mesh->node_ptr(6);
      elem->set_node(1) = mesh->node_ptr(7);
      elem->set_node(2) = mesh->node_ptr(8);
      elem->set_node(3) = mesh->node_ptr(9);
      elem->set_node(4) = mesh->node_ptr(10);
      elem->set_node(5) = mesh->node_ptr(11);

      //mesh->add_point( libMesh::Point(1.0,1.0),6 );
      mesh->add_point( libMesh::Point(-1.0,1.0),12 );
      //mesh->add_point( libMesh::Point(-1.0,-1.0),7 );
      mesh->add_point( libMesh::Point(0.0,1.0),13 );
      mesh->add_point( libMesh::Point(-1.0,0.0),14 );
      //mesh->add_point( libMesh::Point(0.0,0.0),9 );
      
      elem = mesh->add_elem( new libMesh::Tri6 );
      elem->set_id(2);
      elem->subdomain_id() = 2;
      elem->set_node(0) = mesh->node_ptr(6);
      elem->set_node(1) = mesh->node_ptr(12);
      elem->set_node(2) = mesh->node_ptr(7);
      elem->set_node(3) = mesh->node_ptr(13);
      elem->set_node(4) = mesh->node_ptr(14);
      elem->set_node(5) = mesh->node_ptr(9);

      mesh->prepare_for_use();

      std::unique_ptr<libMesh::EquationSystems> _es;
      _es.reset( new libMesh::EquationSystems(*mesh) );

      GRINS::MultiphysicsSystem* _system;
      _system = &_es->add_system<GRINS::MultiphysicsSystem>( "GRINS-TEST" );
      _system->read_input_options( (*_input) );

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

      libMesh::UniquePtr<libMesh::PointLocatorBase> point_locator = mesh->sub_point_locator();

      GRINS::OverlappingFluidSolidMap mesh_overlap( (*_system), (*point_locator), solid_ids, fluid_ids, disp_vars );

      CPPUNIT_ASSERT_EQUAL(1,(int)mesh_overlap.solid_map().size());
      CPPUNIT_ASSERT_EQUAL(2,(int)(mesh_overlap.solid_map().find(0)->second).size());

      CPPUNIT_ASSERT_EQUAL(1,(int)(mesh_overlap.solid_map().find(0)->second).find(1)->first);
      CPPUNIT_ASSERT_EQUAL(2,(int)(mesh_overlap.solid_map().find(0)->second).find(2)->first);

      CPPUNIT_ASSERT_EQUAL(5,(int)((mesh_overlap.solid_map().find(0)->second).find(1)->second).size());
      CPPUNIT_ASSERT_EQUAL(2,(int)((mesh_overlap.solid_map().find(0)->second).find(2)->second).size());

      std::vector<unsigned int> elem1_qps(5);
      elem1_qps[0] = 0; elem1_qps[1] = 1; elem1_qps[2] = 3; elem1_qps[3] = 5; elem1_qps[4] = 6; 

      std::vector<unsigned int> elem2_qps(2);
      elem2_qps[0] = 2; elem2_qps[1] = 4; 

      for( unsigned int i = 0; i < 5; i++ )
        CPPUNIT_ASSERT_EQUAL(elem1_qps[i],((mesh_overlap.solid_map().find(0)->second).find(1)->second)[i]);

      for( unsigned int i = 0; i < 2; i++ )
        CPPUNIT_ASSERT_EQUAL(elem2_qps[i],((mesh_overlap.solid_map().find(0)->second).find(2)->second)[i]);


      CPPUNIT_ASSERT_EQUAL(2,(int)mesh_overlap.fluid_map().size());
      CPPUNIT_ASSERT_EQUAL(1,(int)(mesh_overlap.fluid_map().find(1)->second).size());
      CPPUNIT_ASSERT_EQUAL(1,(int)(mesh_overlap.fluid_map().find(2)->second).size());
      CPPUNIT_ASSERT_EQUAL(0,(int)(mesh_overlap.fluid_map().find(1)->second).find(0)->first);
      CPPUNIT_ASSERT_EQUAL(0,(int)(mesh_overlap.fluid_map().find(2)->second).find(0)->first);
      CPPUNIT_ASSERT_EQUAL(5,(int)((mesh_overlap.fluid_map().find(1)->second).find(0)->second).size());
      CPPUNIT_ASSERT_EQUAL(2,(int)((mesh_overlap.fluid_map().find(2)->second).find(0)->second).size());

      for( unsigned int i = 0; i < 5; i++ )
        CPPUNIT_ASSERT_EQUAL(elem1_qps[i],((mesh_overlap.fluid_map().find(1)->second).find(0)->second)[i]);

      for( unsigned int i = 0; i < 2; i++ )
        CPPUNIT_ASSERT_EQUAL(elem2_qps[i],((mesh_overlap.fluid_map().find(2)->second).find(0)->second)[i]);

      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

   private:


  };

  CPPUNIT_TEST_SUITE_REGISTRATION( OverlappingFluidSolidMeshTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
