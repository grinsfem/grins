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

#include "test_comm.h"
#include "grins_test_paths.h"

// GRINS
#include "grins/grins_enums.h"
#include "grins/mesh_builder.h"
#include "grins/rayfire_mesh.h"
#include "grins/math_constants.h"

#include "rayfire_test_base.h"

// libMesh
#include "libmesh/elem.h"
#include "libmesh/getpot.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_quad9.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/serial_mesh.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

#include <tuple>

namespace GRINSTesting
{
  class RayfireTestAMR3D : public CppUnit::TestCase,
                           public GRINSTesting::RayfireTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( RayfireTestAMR3D );

    CPPUNIT_TEST( single_hex );
    CPPUNIT_TEST( single_tet );
    CPPUNIT_TEST( through_vertex_postrefinment );
    CPPUNIT_TEST( near_vertex_postrefinement );
    CPPUNIT_TEST( refine_deformed_elem );
    CPPUNIT_TEST( refined_and_unrefined );
    CPPUNIT_TEST( two_refined );
    CPPUNIT_TEST( init_on_refined_mesh );

    CPPUNIT_TEST_SUITE_END();

  public:

    void single_hex()
    {
      libMesh::Point origin(0.0,0.1,0.0);
      libMesh::Point end_point(1.0,0.2,1.0);

      std::vector<unsigned int> children_in_rayfire;
      children_in_rayfire.push_back(0);
      children_in_rayfire.push_back(5);

      std::vector<unsigned int> children_not_in_rayfire;
      children_not_in_rayfire.push_back(1);
      children_not_in_rayfire.push_back(2);
      children_not_in_rayfire.push_back(3);
      children_not_in_rayfire.push_back(4);
      children_not_in_rayfire.push_back(6);
      children_not_in_rayfire.push_back(7);

      // HEX8
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_hex8_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);

      // HEX27
      mesh = this->build_square_hex27_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
    }

    void single_tet()
    {
      libMesh::Point origin(0.0,0.1,0.0);
      libMesh::Point end_point(0.25,0.25,0.5);

      std::vector<unsigned int> children_in_rayfire;
      children_in_rayfire.push_back(0);
      children_in_rayfire.push_back(4);
      children_in_rayfire.push_back(7);

      std::vector<unsigned int> children_not_in_rayfire;
      children_not_in_rayfire.push_back(1);
      children_not_in_rayfire.push_back(2);
      children_not_in_rayfire.push_back(3);
      children_not_in_rayfire.push_back(5);
      children_not_in_rayfire.push_back(6);

      // TET4
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_tet4_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);

      // TET10
      mesh = this->build_tet10_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
    }

    //! After refinement, the rayfire will travel through a vertex
    void through_vertex_postrefinment()
    {
      libMesh::Point origin(0.0,0.0,0.0);
      libMesh::Point end_point(1.0,1.0,1.0);

      std::vector<unsigned int> children_in_rayfire;
      children_in_rayfire.push_back(0);
      children_in_rayfire.push_back(7);

      std::vector<unsigned int> children_not_in_rayfire;
      children_not_in_rayfire.push_back(1);
      children_not_in_rayfire.push_back(2);
      children_not_in_rayfire.push_back(3);
      children_not_in_rayfire.push_back(4);
      children_not_in_rayfire.push_back(5);
      children_not_in_rayfire.push_back(6);
      // HEX8
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_hex8_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);

      // HEX27
      mesh = this->build_square_hex27_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
    }

    //! After refinement, the rayfire will travel very near (but not through) a vertex
    void near_vertex_postrefinement()
    {
      libMesh::Point origin(0.4999,0.0,0.0);
      libMesh::Point end_point(1.0,1.0,1.0);

      std::vector<unsigned int> children_in_rayfire;
      children_in_rayfire.push_back(0);
      children_in_rayfire.push_back(1);
      children_in_rayfire.push_back(7);

      std::vector<unsigned int> children_not_in_rayfire;
      children_not_in_rayfire.push_back(2);
      children_not_in_rayfire.push_back(3);
      children_not_in_rayfire.push_back(4);
      children_not_in_rayfire.push_back(5);
      children_not_in_rayfire.push_back(6);

      // HEX8
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_hex8_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);

      // HEX27
      mesh = this->build_square_hex27_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
    }

    //! HEX8 elem is deformed, rayfire goes within libMesh::TOLERANCE central node
    void refine_deformed_elem()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(3);

      mesh->add_point( libMesh::Point(0.0,0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.0,1.0,0.0),2 );
      mesh->add_point( libMesh::Point(0.0,1.0,0.0),3 );
      mesh->add_point( libMesh::Point(0.2,0.0,1.0),4 );
      mesh->add_point( libMesh::Point(1.1,0.0,1.1),5 );
      mesh->add_point( libMesh::Point(1.1,1.1,1.1),6 );
      mesh->add_point( libMesh::Point(0.2,1.1,1.0),7 );

      libMesh::Elem * e = mesh->add_elem( new libMesh::Hex8 );
      for (unsigned int n=0; n<8; n++)
        e->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      CPPUNIT_ASSERT_EQUAL( mesh->n_elem(), (libMesh::dof_id_type)1 );

      libMesh::Real z = 0.5250005;
      libMesh::Real x = z/5.0;
      libMesh::Real y = 0.25;

      libMesh::Point origin(x,y,z); // within libMesh::TOLERANCE of node 8
      libMesh::Point end_point(1.0477273181818181,y,z);

      std::vector<unsigned int> children_in_rayfire;
      children_in_rayfire.push_back(1);
      children_in_rayfire.push_back(4);
      children_in_rayfire.push_back(5);

      std::vector<unsigned int> children_not_in_rayfire;
      children_not_in_rayfire.push_back(0);
      children_not_in_rayfire.push_back(2);
      children_not_in_rayfire.push_back(3);
      children_not_in_rayfire.push_back(6);
      children_not_in_rayfire.push_back(7);

      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
    }

    //! 2 HEX8 elems, one is refined, one if not, rayfire travels boundary between them
    void refined_and_unrefined()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_2elem_hex8();

      libMesh::Elem * elem0 = mesh->elem_ptr(0);
      libMesh::Elem * elem1 = mesh->elem_ptr(1);

      libMesh::Point origin(0.0,0.75,1.0);
      libMesh::Real theta = 0.0;
      libMesh::Real phi = GRINS::Constants::pi/2.0;

      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta,phi) );
      rayfire->init(*mesh);

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(elem0->id()) );
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem1->id()) );

      elem0->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);

      libMesh::MeshRefinement mr(*mesh);
      mr.refine_elements();

      rayfire->reinit(*mesh);

      // make sure we didn't wander into the other elem
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem1->id()) );

      // we should not wander to other children in elem0
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem0->child_ptr(4)->id()) );
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem0->child_ptr(5)->id()) );

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(elem0->child_ptr(6)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(elem0->child_ptr(7)->id()) );
    }

    //! 2 HEX8 elems, both refined, rayfire travels boundary between them
    void two_refined()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_2elem_hex8();

      libMesh::Elem * elem0 = mesh->elem_ptr(0);
      libMesh::Elem * elem1 = mesh->elem_ptr(1);

      libMesh::Point origin(0.0,0.75,1.0);
      libMesh::Real theta = 0.0;
      libMesh::Real phi = GRINS::Constants::pi/2.0;

      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta,phi) );
      rayfire->init(*mesh);

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(elem0->id()) );
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem1->id()) );

      libMesh::MeshRefinement mr(*mesh);
      mr.uniformly_refine();

      rayfire->reinit(*mesh);

      // make sure we didn't wander into the other elem
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem1->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem1->child_ptr(1)->id()) );
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem1->child_ptr(2)->id()) );
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem1->child_ptr(3)->id()) );

      // we should not wander to other children in elem0
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem0->child_ptr(4)->id()) );
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem0->child_ptr(5)->id()) );

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(elem0->child_ptr(6)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(elem0->child_ptr(7)->id()) );
    }

    //! Call init() on an already refined mesh
    //!
    //! Mimics a run from restart
    void init_on_refined_mesh()
    {
      std::string input_string = this->mesh_3D("HEX27",3.0,3.0,3.0,3,3,3);
      std::stringstream ss;
      ss << input_string;

      GetPot input(ss);
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input);

      // refine elem 1 before rayfire init
      libMesh::MeshRefinement mr(*mesh);
      mesh->elem_ptr(1)->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);
      mr.refine_elements();

      CPPUNIT_ASSERT( !(mesh->elem_ptr(1)->active()) );
      CPPUNIT_ASSERT( mesh->elem_ptr(1)->has_children() );

      libMesh::Point origin(0.0,0.8,0.1);
      libMesh::Real theta = 0.0;
      libMesh::Real phi = GRINS::Constants::pi/2.0;

      // initialize the rayfire
      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta,phi) );
      rayfire->init(*mesh);

      // make sure we pick up the 2 unrefined elements
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(0) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(2) );

      // rayfire should not include elem 1 since it's INACTIVE
      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(1)) );

      // rayfire should contain children 0,1 of elem 1
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ref(1).child_ptr(2)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ref(1).child_ptr(3)->id()) );
    }

  private:

    std::shared_ptr<libMesh::UnstructuredMesh> build_2elem_hex8()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.0,1.0,0.0),2 );
      mesh->add_point( libMesh::Point(0.0,1.0,0.0),3 );
      
      mesh->add_point( libMesh::Point(0.0,0.0,1.0),4 );
      mesh->add_point( libMesh::Point(1.0,0.0,1.0),5 );
      mesh->add_point( libMesh::Point(1.0,1.0,1.0),6 );
      mesh->add_point( libMesh::Point(0.0,1.0,1.0),7 );
      
      mesh->add_point( libMesh::Point(0.0,0.0,2.0),8 );
      mesh->add_point( libMesh::Point(1.0,0.0,2.0),9 );
      mesh->add_point( libMesh::Point(1.0,1.0,2.0),10 );
      mesh->add_point( libMesh::Point(0.0,1.0,2.0),11 );

      libMesh::Elem * elem0 = mesh->add_elem( new libMesh::Hex8 );
      elem0->set_node(0) = mesh->node_ptr(0);
      elem0->set_node(1) = mesh->node_ptr(1);
      elem0->set_node(2) = mesh->node_ptr(2);
      elem0->set_node(3) = mesh->node_ptr(3);
      elem0->set_node(4) = mesh->node_ptr(4);
      elem0->set_node(5) = mesh->node_ptr(5);
      elem0->set_node(6) = mesh->node_ptr(6);
      elem0->set_node(7) = mesh->node_ptr(7);

      libMesh::Elem * elem1 = mesh->add_elem( new libMesh::Hex8 );
      elem1->set_node(0) = mesh->node_ptr(4);
      elem1->set_node(1) = mesh->node_ptr(5);
      elem1->set_node(2) = mesh->node_ptr(6);
      elem1->set_node(3) = mesh->node_ptr(7);
      elem1->set_node(4) = mesh->node_ptr(8);
      elem1->set_node(5) = mesh->node_ptr(9);
      elem1->set_node(6) = mesh->node_ptr(10);
      elem1->set_node(7) = mesh->node_ptr(11);

      mesh->prepare_for_use();

      return mesh;
    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( RayfireTestAMR3D );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
