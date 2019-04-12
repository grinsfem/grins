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
  class RayfireTestAMR2D : public CppUnit::TestCase,
                           public GRINSTesting::RayfireTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( RayfireTestAMR2D );

    CPPUNIT_TEST( single_elems );
    CPPUNIT_TEST( through_vertex_postrefinment );
    CPPUNIT_TEST( near_vertex_postrefinement );
    CPPUNIT_TEST( large_2D_mesh );
    CPPUNIT_TEST( refine_elem_not_on_rayfire );
    CPPUNIT_TEST( multiple_refinements );
    CPPUNIT_TEST( coarsen_elements );
    CPPUNIT_TEST( refine_and_coarsen );
    CPPUNIT_TEST( mixed_type_mesh );
    CPPUNIT_TEST( start_with_refined_mesh );
    CPPUNIT_TEST( refined_and_unrefined );
    CPPUNIT_TEST( refine_deformed_elem );
    CPPUNIT_TEST( refine_deformed_elem_near_tolerance );
    CPPUNIT_TEST( init_on_refined_mesh );

    CPPUNIT_TEST_SUITE_END();

  public:

    //! Refine a single element
    void single_elems()
    {
      libMesh::Point origin(0.0,0.1);
      libMesh::Point end_point(1.0,0.2);

      std::vector<unsigned int> children_in_rayfire;
      children_in_rayfire.push_back(0);
      children_in_rayfire.push_back(1);

      std::vector<unsigned int> children_not_in_rayfire;
      children_not_in_rayfire.push_back(2);
      children_not_in_rayfire.push_back(3);

      // QUAD4
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_quad4_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);

      // QUAD9
      mesh = this->build_square_quad9_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
    }

    //! After refinement, the rayfire will travel through a vertex
    void through_vertex_postrefinment()
    {
      libMesh::Point origin(0.0,0.0);
      libMesh::Point mid_point(0.5,0.5);
      libMesh::Point end_point(1.0,1.0);

      std::vector<unsigned int> children_in_rayfire;
      children_in_rayfire.push_back(0);
      children_in_rayfire.push_back(3);

      std::vector<unsigned int> children_not_in_rayfire;
      children_not_in_rayfire.push_back(1);
      children_not_in_rayfire.push_back(2);

      // QUAD4
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_quad4_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);

      // QUAD9
      mesh = this->build_square_quad9_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
    }

    //! After refinement, the rayfire will travel very near (but not through) a vertex
    void near_vertex_postrefinement()
    {
      libMesh::Point origin(0.4999,0.0);
      libMesh::Point end_point(1.0,1.0);

      std::vector<unsigned int> children_in_rayfire;
      children_in_rayfire.push_back(0);
      children_in_rayfire.push_back(1);
      children_in_rayfire.push_back(3);

      std::vector<unsigned int> children_not_in_rayfire;
      children_not_in_rayfire.push_back(2);

      // QUAD4
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_quad4_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);

      // QUAD9
      mesh = this->build_square_quad9_elem();
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
    }

    //! A 10x10 mesh with selectively refined elements
    void large_2D_mesh()
    {
      std::string quad4_string = this->mesh_2D("QUAD4",10.0,10.0,10,10);
      std::stringstream ss_quad4;
      ss_quad4 << quad4_string;

      GetPot input_quad4(ss_quad4);
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input_quad4);
      test_large_mesh(mesh);

      std::string quad9_string = this->mesh_2D("QUAD9",10.0,10.0,10,10);
      std::stringstream ss_quad9;
      ss_quad9 << quad9_string;

      GetPot input_quad9(ss_quad9);
      mesh = this->build_mesh(input_quad9);
      test_large_mesh(mesh);
    }

    //! Make sure refined elements not on the rayfire do not get added
    void refine_elem_not_on_rayfire()
    {
      libMesh::Point origin(0.0,6.5);
      libMesh::Point end_point(10.0,0.25);

      std::vector<unsigned int> children_in_rayfire;

      std::vector<unsigned int> children_not_in_rayfire;
      children_not_in_rayfire.push_back(0);
      children_not_in_rayfire.push_back(2);
      children_not_in_rayfire.push_back(1);
      children_not_in_rayfire.push_back(3);

      // QUAD 4
      std::string quad4_string = this->mesh_2D("QUAD4",10.0,10.0,10,10);
      std::stringstream ss_quad4;
      ss_quad4 << quad4_string;

      GetPot input_quad4(ss_quad4);
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input_quad4);
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
      
      // QUAD 9
      std::string quad9_string = this->mesh_2D("QUAD9",10.0,10.0,10,10);
      std::stringstream ss_quad9;
      ss_quad9 << quad9_string;

      GetPot input_quad9(ss_quad9);
      mesh = this->build_mesh(input_quad9);
      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
    }

    //! Test for multiple successive refinements
    void multiple_refinements()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_quad4_elem();
      test_multiple_refinements(mesh);

      mesh = this->build_square_quad9_elem();
      test_multiple_refinements(mesh);
    }

    //! Test of the coarsening functionality
    void coarsen_elements()
    {
      std::string quad4_string = this->mesh_2D("QUAD4",10.0,10.0,10,10);
      std::stringstream ss_quad4;
      ss_quad4 << quad4_string;

      GetPot input_quad4(ss_quad4);
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input_quad4);
      test_coarsen(mesh);

      std::string quad9_string = this->mesh_2D("QUAD9",10.0,10.0,10,10);
      std::stringstream ss_quad9;
      ss_quad9 << quad9_string;

      GetPot input_quad9(ss_quad9);
      mesh = this->build_mesh(input_quad9);
      test_coarsen(mesh);
    }

    //! Refine and coarsen before calling reinit()
    void refine_and_coarsen()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_quad4_elem();

      libMesh::Point origin(0.0,0.1);
      libMesh::Point end_point(1.0,0.1);
      libMesh::Real theta = this->calc_theta(origin,end_point);

      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta) );
      rayfire->init(*mesh);

      libMesh::MeshRefinement mr(*mesh);
      mr.uniformly_refine();
      rayfire->reinit(*mesh);
      mr.uniformly_refine();
      rayfire->reinit(*mesh);

      // refine elem(0)->child_ptr(1)->child_ptr(1)
      mesh->elem_ptr(0)->child_ptr(1)->child_ptr(1)->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);
      // coarsen elem(0)->child_ptr(0)
      for (unsigned int c=0; c<mesh->elem_ptr(0)->child_ptr(0)->n_children(); c++)
        mesh->elem_ptr(0)->child_ptr(0)->child_ptr(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);

      mr.refine_and_coarsen_elements();
      rayfire->reinit(*mesh);

      CPPUNIT_ASSERT( mesh->elem_ptr(0)->child_ptr(0)->active() );

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->child_ptr(1)->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->child_ptr(1)->child_ptr(1)->id()) );
    }

    //! Mesh contains 2 QUAD4 and 2 TRI3
    void mixed_type_mesh()
    {
      std::string input_string = this->mesh_mixed_quad_tri();
      std::stringstream ss;
      ss << input_string;

      GetPot input(ss);
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input);

      libMesh::Point origin(0.0,0.25);
      libMesh::Point end_point(1.0,0.25);
      libMesh::Real theta = this->calc_theta(origin,end_point);

      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta) );
      rayfire->init(*mesh);

      // check that the rayfire is through the right elems
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(0) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(6) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(9) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(8) );

      // and not through others
      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(2)) );
      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(3)) );
      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(4)) );
      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(5)) );
      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(1)) );
      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(7)) );

      // now do a uniform refinement
      libMesh::MeshRefinement mr(*mesh);
      mr.uniformly_refine();

      rayfire->reinit(*mesh);

      // check post-refinement
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(6)->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(9)->child_ptr(1)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(8)->child_ptr(2)->id()) );

      // coarsen two of the TRIs
      for (unsigned int c=0; c<mesh->elem_ptr(6)->n_children(); c++)
        {
          mesh->elem_ptr(6)->child_ptr(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
          mesh->elem_ptr(9)->child_ptr(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
        }

      mr.coarsen_elements();
      rayfire->reinit(*mesh);

      CPPUNIT_ASSERT( mesh->elem_ptr(6)->active() );
      CPPUNIT_ASSERT( mesh->elem_ptr(9)->active() );

      // check post-coarsen
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(8)->child_ptr(2)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(6) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(9) );

      for (unsigned int c=0; c<mesh->elem_ptr(6)->n_children(); c++)
        {
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem_ptr(6)->child_ptr(c)->id())) );
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem_ptr(9)->child_ptr(c)->id())) );
        }
    }

    //! Mesh given to init() is already refined
    void start_with_refined_mesh()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_quad4_elem();

      // uniform refinement
      libMesh::MeshRefinement mr(*mesh);
      mr.uniformly_refine();

      // now init the rayfire
      libMesh::Point origin(0.0,0.1);
      libMesh::Point end_point(1.0,0.1);
      libMesh::Real theta = this->calc_theta(origin,end_point);

      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta) );
      rayfire->init(*mesh);

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->id()) );

      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(2)->id())) );
      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(3)->id())) );

      // refine again
      mr.uniformly_refine();
      rayfire->reinit(*mesh);

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(0)->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(0)->child_ptr(1)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->child_ptr(1)->id()) );

      // uniformly coarsen twice to get back to a single elem
      mr.uniformly_coarsen();
      rayfire->reinit(*mesh);
      mr.uniformly_coarsen();
      rayfire->reinit(*mesh);

      // ensure the original elem is active
      CPPUNIT_ASSERT( mesh->elem_ptr(0)->active() );

      // the original elem should be in the rayfire
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(0) );

      // and no other elem
      for (unsigned int e=1; e<mesh->n_elem(); e++)
        CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(e)) );

    }

    //! 2 QUAD4 elems, one is refined, one if not, rayfire travels boundary between them
    void refined_and_unrefined()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.0,1.0),2 );
      mesh->add_point( libMesh::Point(0.0,1.0),3 );
      mesh->add_point( libMesh::Point(1.0,2.0),4 );
      mesh->add_point( libMesh::Point(0.0,2.0),5 );

      libMesh::Elem * elem0 = mesh->add_elem( new libMesh::Quad4 );
      elem0->set_node(0) = mesh->node_ptr(0);
      elem0->set_node(1) = mesh->node_ptr(1);
      elem0->set_node(2) = mesh->node_ptr(2);
      elem0->set_node(3) = mesh->node_ptr(3);

      libMesh::Elem * elem1 = mesh->add_elem( new libMesh::Quad4 );
      elem1->set_node(0) = mesh->node_ptr(3);
      elem1->set_node(1) = mesh->node_ptr(2);
      elem1->set_node(2) = mesh->node_ptr(4);
      elem1->set_node(3) = mesh->node_ptr(5);

      mesh->prepare_for_use();

      libMesh::Point origin(0.0,1.0);
      libMesh::Real theta = 0.0;

      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta) );
      rayfire->init(*mesh);

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(elem0->id()) );
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem1->id()) );

      elem0->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);

      libMesh::MeshRefinement mr(*mesh);
      mr.refine_elements();

      rayfire->reinit(*mesh);

      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem1->id()) );

      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem0->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem0->child_ptr(1)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(elem0->child_ptr(2)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(elem0->child_ptr(3)->id()) );
    }

    //! QUAD4 elem is deformed, rayfire goes within libMesh::TOLERANCE central node
    void refine_deformed_elem()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.1,1.1),2 );
      mesh->add_point( libMesh::Point(0.2,1.0),3 );

      libMesh::Elem * e = mesh->add_elem( new libMesh::Quad4 );
      for (unsigned int n=0; n<4; n++)
        e->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      CPPUNIT_ASSERT_EQUAL( mesh->n_elem(), (libMesh::dof_id_type)1 );

      libMesh::Real y = 0.5250005;
      libMesh::Real x = y/5.0;

      libMesh::Point origin(x,y); // within libMesh::TOLERANCE of node 8
      libMesh::Point end_point(1.0477273181818181,y);

      std::vector<unsigned int> children_in_rayfire;
      children_in_rayfire.push_back(1);
      children_in_rayfire.push_back(2);
      children_in_rayfire.push_back(3);

      std::vector<unsigned int> children_not_in_rayfire;
      children_not_in_rayfire.push_back(0);

      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
    }

    //! QUAD4 elem is slightly deformed, rayfire starts and ends very neay boundary of refined children.
    //!
    //! Based on an element from a larger mesh, hence the seemingly arbitrary coordinates
    void refine_deformed_elem_near_tolerance()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.26591876082146976,0.016285303166482301),2 );
      mesh->add_point( libMesh::Point(0.26539426622418877,0.016285726631637167),3 );
      mesh->add_point( libMesh::Point(0.26538091506882988,0.015714299294800886),0 );
      mesh->add_point( libMesh::Point(0.26590538328042834,0.015713833483130532),1 );

      libMesh::Elem * e = mesh->add_elem( new libMesh::Quad4 );
      for (unsigned int n=0; n<4; n++)
        e->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      CPPUNIT_ASSERT_EQUAL( mesh->n_elem(), (libMesh::dof_id_type)1 );

      libMesh::Real x = 0.26538759034362924;
      libMesh::Real y = 0.01600000000000000;

      libMesh::Point origin(x,y);
      libMesh::Point end_point(0.26591208215603906,y);

      std::vector<unsigned int> children_in_rayfire;
      children_in_rayfire.push_back(0);
      children_in_rayfire.push_back(2);
      children_in_rayfire.push_back(3);

      std::vector<unsigned int> children_not_in_rayfire;
      children_not_in_rayfire.push_back(1);

      this->amr_single_elem(mesh,origin,end_point,children_in_rayfire,children_not_in_rayfire);
    }

    //! Call init() on an already refined mesh
    //!
    //! Mimics a run from restart
    void init_on_refined_mesh()
    {
      std::string input_string = this->mesh_2D("QUAD4",3.0,3.0,3,3);
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

      libMesh::Point origin(0.0,0.1);
      libMesh::Real theta = 0.0;

      // initialize the rayfire
      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta) );
      rayfire->init(*mesh);

      // make sure we pick up the 2 unrefined elements
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(0) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(2) );

      // rayfire should not include elem 1 since it's INACTIVE
      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(1)) );

      // rayfire should contain children 0,1 of elem 1
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ref(1).child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ref(1).child_ptr(1)->id()) );
    }

  private:

    //! Refine specific elements on a large mesh
    void test_large_mesh(std::shared_ptr<libMesh::UnstructuredMesh> mesh)
    {
      libMesh::Point origin(0.0,6.5);
      libMesh::Point end_point(10.0,0.25);
      libMesh::Real theta = this->calc_theta(origin,end_point);

      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta) );
      rayfire->init(*mesh);

      // manually refine specific elements
      // map index 0: main_elem id to refine
      // map index 1: children that should be in the rayfire post-refinement
      std::map<unsigned int,std::vector<unsigned int> > refine_elems;
      refine_elems[60] = std::vector<unsigned int>();
      refine_elems[60].push_back(0);
      refine_elems[60].push_back(1);

      refine_elems[52] = std::vector<unsigned int>();
      refine_elems[52].push_back(0);

      refine_elems[34] = std::vector<unsigned int>();
      refine_elems[34].push_back(1);
      refine_elems[34].push_back(2);
      refine_elems[34].push_back(3);

      refine_elems[25] = std::vector<unsigned int>();
      refine_elems[25].push_back(3);

      refine_elems[26] = std::vector<unsigned int>();
      refine_elems[26].push_back(0);
      refine_elems[26].push_back(1);
      refine_elems[26].push_back(2);

      refine_elems[17] = std::vector<unsigned int>();
      refine_elems[17].push_back(2);
      refine_elems[17].push_back(3);

      refine_elems[9] = std::vector<unsigned int>();
      refine_elems[9].push_back(1);
      refine_elems[9].push_back(2);
      refine_elems[9].push_back(3);

      // check to make sure those elements are along the rayfire
      // and then set their refinement flags
      std::map<unsigned int,std::vector<unsigned int> >::iterator it = refine_elems.begin();
      for(; it != refine_elems.end(); it++)
        {
          libMesh::Elem * elem = mesh->elem_ptr( it->first );
          CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(elem->id()) );
          elem->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);
        }

      // refine the elements
      libMesh::MeshRefinement mr(*mesh);
      mr.refine_elements();

      // reinitialize the rayfire
      rayfire->reinit(*mesh);

      // check the refined elements and their children
      it = refine_elems.begin();
      for(; it != refine_elems.end(); it++)
        {
          libMesh::Elem * elem = mesh->elem_ptr( it->first );

          // make sure it was refined
          CPPUNIT_ASSERT_EQUAL( elem->refinement_flag(), libMesh::Elem::RefinementState::INACTIVE );

          // make sure it was deleted from the rayfire
          CPPUNIT_ASSERT( !( rayfire->map_to_rayfire_elem(elem->id()) ) );

          // check that the computed children are in the rayfire
          // and that the rest are not
          std::vector<unsigned int> children = it->second;
          unsigned int index = 0;
          for(unsigned int c=0; c<elem->n_children(); c++)
            {
              if (index<children.size())
                {
                  if (c==children[index])
                    {
                      index++;
                      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem( elem->child_ptr(c)->id() ) );
                    }
                }
              else
                CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem( elem->child_ptr(c)->id() )) );
            }
        }
    }

    //! Refine the mesh 3 times
    void test_multiple_refinements(std::shared_ptr<libMesh::UnstructuredMesh> mesh)
    {
      libMesh::Point origin(0.0,0.0);
      libMesh::Point end_point(1.0,0.25);
      libMesh::Real theta = this->calc_theta(origin,end_point);

      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta) );
      rayfire->init(*mesh);

      libMesh::MeshRefinement mr(*mesh);

      // do 3 successive uniform refinements
      for (unsigned int i=0; i<3; i++)
        {
          mr.uniformly_refine();
          rayfire->reinit(*mesh);
        }

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(0)->child_ptr(0)->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(0)->child_ptr(0)->child_ptr(1)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(0)->child_ptr(1)->child_ptr(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(0)->child_ptr(1)->child_ptr(1)->id()) );

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->child_ptr(0)->child_ptr(2)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->child_ptr(0)->child_ptr(3)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->child_ptr(1)->child_ptr(2)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem_ptr(0)->child_ptr(1)->child_ptr(1)->child_ptr(3)->id()) );
    }

    //! Test coarsening specific elements on a large mesh
    void test_coarsen(std::shared_ptr<libMesh::UnstructuredMesh> mesh)
    {
      libMesh::Point origin(0.0,6.5);
      libMesh::Point end_point(10.0,0.25);
      libMesh::Real theta = this->calc_theta(origin,end_point);

      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta) );
      rayfire->init(*mesh);

      // unifromly refine the elements
      libMesh::MeshRefinement mr(*mesh);
      mr.uniformly_refine();
      rayfire->reinit(*mesh);

      // coarsen specific elements along the rayfire
      for(unsigned int c=0; c<4; c++)
        {
          mesh->elem_ptr(60)->child_ptr(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
          mesh->elem_ptr(25)->child_ptr(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
          mesh->elem_ptr(26)->child_ptr(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
          mesh->elem_ptr(9)->child_ptr(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
        }

      mr.coarsen_elements();
      rayfire->reinit(*mesh);

      // check that the coarsened elems are in the rayfire
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(60) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(25) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(26) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(9) );

      // and that their children are not
      for(unsigned int c=0; c<4; c++)
        {
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem_ptr(60)->child_ptr(c)->id())) );
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem_ptr(25)->child_ptr(c)->id())) );
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem_ptr(26)->child_ptr(c)->id())) );
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem_ptr(9)->child_ptr(c)->id())) );
        }
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( RayfireTestAMR2D );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
