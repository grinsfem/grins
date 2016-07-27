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

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include "test_comm.h"
#include "grins_test_paths.h"

// GRINS
#include "grins/grins_enums.h"
#include "grins/mesh_builder.h"
#include "grins/rayfire_mesh.h"
#include "grins/math_constants.h"

// libMesh
#include "libmesh/elem.h"
#include "libmesh/getpot.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_quad9.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/serial_mesh.h"

namespace GRINSTesting
{
  class RayfireTestAMR : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( RayfireTestAMR );

    CPPUNIT_TEST( single_elems );
    CPPUNIT_TEST( through_vertex_postrefinment );
    CPPUNIT_TEST( large_2D_mesh );
    CPPUNIT_TEST( refine_elem_not_on_rayfire );
    CPPUNIT_TEST( multiple_refinements );
    CPPUNIT_TEST( coarsen_elements );
    CPPUNIT_TEST( refine_and_coarsen );
    CPPUNIT_TEST( mixed_type_mesh );
    CPPUNIT_TEST( start_with_refined_mesh );

    CPPUNIT_TEST_SUITE_END();

  public:

    //! Refine a single element
    void single_elems()
    {
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = build_quad4_elem();
      test_single_elem(mesh);

      mesh = build_quad9_elem();
      test_single_elem(mesh);
    }

    //! After refinement, the rayfire will travel through a vertex
    void through_vertex_postrefinment()
    {
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = build_quad4_elem();
      test_through_vertex(mesh);

      mesh = build_quad9_elem();
      test_through_vertex(mesh);
    }

    //! A 10x10 mesh with selectively refined elements
    void large_2D_mesh()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/mesh_quad4_100elem_2D.in";
      GetPot input_quad4(filename);
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input_quad4);
      test_large_mesh(mesh);

      filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/mesh_quad9_100elem_2D.in";
      GetPot input_quad9(filename);
      mesh = this->build_mesh(input_quad9);
      test_large_mesh(mesh);
    }

    //! Make sure refined elements not on the rayfire do not get added
    void refine_elem_not_on_rayfire()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/mesh_quad4_100elem_2D.in";
      GetPot input(filename);
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input);

      libMesh::Point origin(0.0,6.5);
      libMesh::Point end_point(10.0,0.25);
      libMesh::Real theta = calc_theta(origin,end_point);

      GRINS::SharedPtr<GRINS::RayfireMesh> rayfire = new GRINS::RayfireMesh(origin,theta);
      rayfire->init(*mesh);

      // elem 0 is not along the rayfire
      mesh->elem(0)->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);

      libMesh::MeshRefinement mr(*mesh);
      mr.refine_elements();

      rayfire->reinit(*mesh);

      // elem 0 should not be in the rayfire
      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem(0)->id()) ) );

      // no children of elem 0 should be in the rayfire
      for (unsigned int i=0; i<4; i++)
        CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem( mesh->elem(0)->child(i)->id() ) ) );
    }

    //! Test for multiple successive refinements
    void multiple_refinements()
    {
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = build_quad4_elem();
      test_multiple_refinements(mesh);

      mesh = build_quad9_elem();
      test_multiple_refinements(mesh);
    }

    //! Test of the coarsening functionality
    void coarsen_elements()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/mesh_quad4_100elem_2D.in";
      GetPot input_quad4(filename);
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input_quad4);
      test_coarsen(mesh);

      filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/mesh_quad9_100elem_2D.in";
      GetPot input_quad9(filename);
      mesh = this->build_mesh(input_quad9);
      test_coarsen(mesh);
    }

    //! Refine and coarsen before calling reinit()
    void refine_and_coarsen()
    {
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = build_quad4_elem();

      libMesh::Point origin(0.0,0.1);
      libMesh::Point end_point(1.0,0.1);
      libMesh::Real theta = calc_theta(origin,end_point);

      GRINS::SharedPtr<GRINS::RayfireMesh> rayfire = new GRINS::RayfireMesh(origin,theta);
      rayfire->init(*mesh);

      libMesh::MeshRefinement mr(*mesh);
      mr.uniformly_refine();
      rayfire->reinit(*mesh);
      mr.uniformly_refine();
      rayfire->reinit(*mesh);

      // refine elem(0)->child(1)->child(1)
      mesh->elem(0)->child(1)->child(1)->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);
      // coarsen elem(0)->child(0)
      for (unsigned int c=0; c<mesh->elem(0)->child(0)->n_children(); c++)
        mesh->elem(0)->child(0)->child(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);

      mr.refine_and_coarsen_elements();
      rayfire->reinit(*mesh);

      CPPUNIT_ASSERT( mesh->elem(0)->child(0)->active() );

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->child(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->child(1)->child(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->child(1)->child(1)->id()) );


    }

    //! Mesh contains 2 QUAD4 and 2 TRI3
    void mixed_type_mesh()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/mixed_quad_tri.in";
      GetPot input(filename);
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input);

      libMesh::Point origin(0.0,0.25);
      libMesh::Point end_point(1.0,0.25);
      libMesh::Real theta = calc_theta(origin,end_point);

      GRINS::SharedPtr<GRINS::RayfireMesh> rayfire = new GRINS::RayfireMesh(origin,theta);
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
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(6)->child(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(9)->child(1)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(8)->child(2)->id()) );

      // coarsen two of the TRIs
      for (unsigned int c=0; c<mesh->elem(6)->n_children(); c++)
        {
          mesh->elem(6)->child(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
          mesh->elem(9)->child(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
        }

      mr.coarsen_elements();
      rayfire->reinit(*mesh);

      CPPUNIT_ASSERT( mesh->elem(6)->active() );
      CPPUNIT_ASSERT( mesh->elem(9)->active() );

      // check post-coarsen
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(8)->child(2)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(6) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(9) );

      for (unsigned int c=0; c<mesh->elem(6)->n_children(); c++)
        {
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem(6)->child(c)->id())) );
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem(9)->child(c)->id())) );
        }
    }

    //! Mesh given to init() is already refined
    void start_with_refined_mesh()
    {
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = build_quad4_elem();

      // uniform refinement
      libMesh::MeshRefinement mr(*mesh);
      mr.uniformly_refine();

      // now init the rayfire
      libMesh::Point origin(0.0,0.1);
      libMesh::Point end_point(1.0,0.1);
      libMesh::Real theta = calc_theta(origin,end_point);

      GRINS::SharedPtr<GRINS::RayfireMesh> rayfire = new GRINS::RayfireMesh(origin,theta);
      rayfire->init(*mesh);

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->id()) );

      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem(0)->child(2)->id())) );
      CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem(0)->child(3)->id())) );

      // refine again
      mr.uniformly_refine();
      rayfire->reinit(*mesh);

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(0)->child(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(0)->child(1)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->child(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->child(1)->id()) );

      // uniformly coarsen twice to get back to a single elem
      mr.uniformly_coarsen();
      rayfire->reinit(*mesh);
      mr.uniformly_coarsen();
      rayfire->reinit(*mesh);

      // ensure the original elem is active
      CPPUNIT_ASSERT( mesh->elem(0)->active() );

      // the original elem should be in the rayfire
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(0) );

      // and no other elem
      for (unsigned int e=1; e<mesh->n_elem(); e++)
        CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(e)) );
    }

  private:

    libMesh::Real calc_theta(libMesh::Point& start, libMesh::Point end)
    {
      return std::atan2( (end(1)-start(1)), (end(0)-start(0)) );
    }

    GRINS::SharedPtr<libMesh::UnstructuredMesh> build_mesh( const GetPot& input )
    {
      GRINS::MeshBuilder mesh_builder;
      return mesh_builder.build( input, *TestCommWorld );
    }

    //! Build a single square QUAD4
    GRINS::SharedPtr<libMesh::UnstructuredMesh> build_quad4_elem()
    {
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = new libMesh::SerialMesh(*TestCommWorld);

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.0,1.0),2 );
      mesh->add_point( libMesh::Point(0.0,1.0),3 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Quad4 );
      for (unsigned int n=0; n<4; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      return mesh;
    }

    //! Build a single square QUAD9
    GRINS::SharedPtr<libMesh::UnstructuredMesh> build_quad9_elem()
    {
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = new libMesh::SerialMesh(*TestCommWorld);

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.0,1.0),2 );
      mesh->add_point( libMesh::Point(0.0,1.0),3 );
      mesh->add_point( libMesh::Point(0.5,0.0),4 );
      mesh->add_point( libMesh::Point(1.0,0.5),5 );
      mesh->add_point( libMesh::Point(0.5,1.0),6 );
      mesh->add_point( libMesh::Point(0.0,0.5),7 );
      mesh->add_point( libMesh::Point(0.5,0.5),8 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Quad9 );
      for (unsigned int n=0; n<9; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      return mesh;
    }

    //! Runs the test on single square QUADs of unit length and width
    void test_single_elem(GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh)
    {
      CPPUNIT_ASSERT_EQUAL( mesh->n_elem(), (libMesh::dof_id_type)1 );

      libMesh::Point origin(0.0,0.1);
      libMesh::Point end_point(1.0,0.2);
      libMesh::Real theta = calc_theta(origin,end_point);

      GRINS::SharedPtr<GRINS::RayfireMesh> rayfire = new GRINS::RayfireMesh(origin,theta);
      rayfire->init(*mesh);

      libMesh::Elem* elem = mesh->elem(0);
      CPPUNIT_ASSERT(elem);

      elem->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);

      libMesh::MeshRefinement mr(*mesh);
      mr.refine_elements();

      rayfire->reinit(*mesh);

      const libMesh::Elem* rayfire_elem0 = rayfire->map_to_rayfire_elem( elem->child(0)->id() );
      CPPUNIT_ASSERT(rayfire_elem0);

      const libMesh::Elem* rayfire_elem1 = rayfire->map_to_rayfire_elem( elem->child(1)->id() );
      CPPUNIT_ASSERT(rayfire_elem1);

      CPPUNIT_ASSERT( (rayfire_elem0->get_node(0))->absolute_fuzzy_equals(origin) );
      CPPUNIT_ASSERT( (rayfire_elem0->get_node(1))->absolute_fuzzy_equals(libMesh::Point(0.5,0.15)) );
      CPPUNIT_ASSERT( (rayfire_elem1->get_node(0))->absolute_fuzzy_equals(libMesh::Point(0.5,0.15)) );
      CPPUNIT_ASSERT( (rayfire_elem1->get_node(1))->absolute_fuzzy_equals(end_point) );
    }

    void test_through_vertex(GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh)
    {
      libMesh::Point origin(0.0,0.0);
      libMesh::Point end_point(1.0,1.0);
      libMesh::Real theta = calc_theta(origin,end_point);

      GRINS::SharedPtr<GRINS::RayfireMesh> rayfire = new GRINS::RayfireMesh(origin,theta);
      rayfire->init(*mesh);

      libMesh::Elem* elem = mesh->elem(0);
      CPPUNIT_ASSERT(elem);

      elem->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);

      libMesh::MeshRefinement mr(*mesh);
      mr.refine_elements();

      rayfire->reinit(*mesh);

      const libMesh::Elem* rayfire_elem0 = rayfire->map_to_rayfire_elem( elem->child(0)->id() );
      CPPUNIT_ASSERT(rayfire_elem0);

      const libMesh::Elem* rayfire_elem1 = rayfire->map_to_rayfire_elem( elem->child(3)->id() );
      CPPUNIT_ASSERT(rayfire_elem1);

      CPPUNIT_ASSERT( (rayfire_elem0->get_node(0))->absolute_fuzzy_equals(origin) );
      CPPUNIT_ASSERT( (rayfire_elem0->get_node(1))->absolute_fuzzy_equals(libMesh::Point(0.5,0.5)) );
      CPPUNIT_ASSERT( (rayfire_elem1->get_node(0))->absolute_fuzzy_equals(libMesh::Point(0.5,0.5)) );
      CPPUNIT_ASSERT( (rayfire_elem1->get_node(1))->absolute_fuzzy_equals(end_point) );
    }

    //! Refine specific elements on a large mesh
    void test_large_mesh(GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh)
    {
      libMesh::Point origin(0.0,6.5);
      libMesh::Point end_point(10.0,0.25);
      libMesh::Real theta = calc_theta(origin,end_point);

      GRINS::SharedPtr<GRINS::RayfireMesh> rayfire = new GRINS::RayfireMesh(origin,theta);
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
          libMesh::Elem* elem = mesh->elem( it->first );
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
          libMesh::Elem* elem = mesh->elem( it->first );

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
              if (c==children[index] && index<children.size())
                {
                  index++;
                  CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem( elem->child(c)->id() ) );
                }
              else
                CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem( elem->child(c)->id() )) );
            }
        }
    }

    //! Refine the mesh 3 times
    void test_multiple_refinements(GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh)
    {
      libMesh::Point origin(0.0,0.0);
      libMesh::Point end_point(1.0,0.25);
      libMesh::Real theta = calc_theta(origin,end_point);

      GRINS::SharedPtr<GRINS::RayfireMesh> rayfire = new GRINS::RayfireMesh(origin,theta);
      rayfire->init(*mesh);

      libMesh::MeshRefinement mr(*mesh);

      // do 3 successive uniform refinements
      for (unsigned int i=0; i<3; i++)
        {
          mr.uniformly_refine();
          rayfire->reinit(*mesh);
        }

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(0)->child(0)->child(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(0)->child(0)->child(1)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(0)->child(1)->child(0)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(0)->child(1)->child(1)->id()) );

      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->child(0)->child(2)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->child(0)->child(3)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->child(1)->child(2)->id()) );
      CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(mesh->elem(0)->child(1)->child(1)->child(3)->id()) );

    }

    //! Test coarsening specific elements on a large mesh
    void test_coarsen(GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh)
    {
      libMesh::Point origin(0.0,6.5);
      libMesh::Point end_point(10.0,0.25);
      libMesh::Real theta = calc_theta(origin,end_point);

      GRINS::SharedPtr<GRINS::RayfireMesh> rayfire = new GRINS::RayfireMesh(origin,theta);
      rayfire->init(*mesh);

      // unifromly refine the elements
      libMesh::MeshRefinement mr(*mesh);
      mr.uniformly_refine();
      rayfire->reinit(*mesh);

      // coarsen specific elements along the rayfire
      for(unsigned int c=0; c<4; c++)
        {
          mesh->elem(60)->child(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
          mesh->elem(25)->child(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
          mesh->elem(26)->child(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
          mesh->elem(9)->child(c)->set_refinement_flag(libMesh::Elem::RefinementState::COARSEN);
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
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem(60)->child(c)->id())) );
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem(25)->child(c)->id())) );
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem(26)->child(c)->id())) );
          CPPUNIT_ASSERT( !(rayfire->map_to_rayfire_elem(mesh->elem(9)->child(c)->id())) );
        }
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( RayfireTestAMR );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
