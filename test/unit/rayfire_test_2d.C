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
#include "libmesh/serial_mesh.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class RayfireTest2D : public CppUnit::TestCase,
                        public GRINSTesting::RayfireTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( RayfireTest2D );

    CPPUNIT_TEST( quad4_all_sides );
    CPPUNIT_TEST( quad9_all_sides );
    CPPUNIT_TEST( test_slanted_quad4 );
    CPPUNIT_TEST( test_vertical_fire );
    CPPUNIT_TEST( test_quad4_5elem );
    CPPUNIT_TEST( test_quad9_5elem );
    CPPUNIT_TEST( test_quad4_2D );
    CPPUNIT_TEST( test_quad9_2D );
    CPPUNIT_TEST( fire_through_vertex );
    CPPUNIT_TEST( origin_between_elems );
    CPPUNIT_TEST( quad4_off_origin );
    CPPUNIT_TEST( quadratic_top_quad9 );

    CPPUNIT_TEST_SUITE_END();

  public:

    void quad4_all_sides()
    {
      // vector of intersection points
      std::vector<libMesh::Point> pts(4);
      pts[0] = libMesh::Point(0.0,0.5);
      pts[1] = libMesh::Point(0.915243860856226,0.0);
      pts[2] = libMesh::Point(1.0,0.25);
      pts[3] = libMesh::Point(0.321046307967165,1.0);

      // create the mesh (single square QUAD4 element)
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.0,1.0),2 );
      mesh->add_point( libMesh::Point(0.0,1.0),3 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Quad4 );
      for (unsigned int n=0; n<4; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      this->run_test_on_all_point_combinations(pts,mesh);
    }

    void quad9_all_sides()
    {
      // vector of intersection points
      std::vector<libMesh::Point> pts(4);
      pts[0] = libMesh::Point(0.0,0.5);
      pts[1] = libMesh::Point(0.915243860856226,0.0);
      pts[2] = libMesh::Point(1.375,0.25);
      pts[3] = libMesh::Point(0.321046307967165,1.0);

      // create a non-rectangular QUAD9
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.0,1.0),2 );
      mesh->add_point( libMesh::Point(0.0,1.0),3 );
      mesh->add_point( libMesh::Point(0.5,0.0),4 );
      mesh->add_point( libMesh::Point(1.5,0.5),5 );
      mesh->add_point( libMesh::Point(0.5,1.0),6 );
      mesh->add_point( libMesh::Point(0.0,0.5),7 );
      mesh->add_point( libMesh::Point(0.5,0.5),8 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Quad9 );
      for (unsigned int n=0; n<9; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      this->run_test_on_all_point_combinations(pts,mesh);
    }

    void test_slanted_quad4()
    {
      // vector of intersection points
      std::vector<libMesh::Point> pts(4);
      pts[0] = libMesh::Point(0.5,0.0);
      pts[1] = libMesh::Point(1.75,0.75);
      pts[2] = libMesh::Point(1.5,0.9);
      pts[3] = libMesh::Point(-0.1,0.1);

      // create the mesh (single trapezoidal QUAD4 element)
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0),1 );
      mesh->add_point( libMesh::Point(2.0,1.0),2 );
      mesh->add_point( libMesh::Point(-0.5,0.5),3 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Quad4 );
      for (unsigned int n=0; n<4; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      this->run_test_on_all_point_combinations(pts,mesh);
    }

    void test_vertical_fire()
    {
      // create the mesh (single square QUAD4 element)
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.0,1.0),2 );
      mesh->add_point( libMesh::Point(0.0,1.0),3 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Quad4 );
      for (unsigned int n=0; n<4; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      libMesh::Point start_point(0.3,0.0);
      libMesh::Point end_point(0.3,1.0);

      libMesh::Real theta = GRINS::Constants::pi/2.0;

      this->run_test_with_mesh(mesh,start_point,theta,-1.0,end_point,0);
    }

    void test_quad4_5elem()
    {
      libMesh::Point origin = libMesh::Point(0,0.5);

      libMesh::Node calc_end_node_straight = libMesh::Node(5.0,0.5);
      this->run_test(origin,0.0,calc_end_node_straight,5,4,"quad4",1);

      libMesh::Node calc_end_node_angle = libMesh::Node(0.5/std::tan(0.15),1.0);
      this->run_test(origin,0.15,calc_end_node_angle,5,3,"quad4",1);

      libMesh::Node calc_end_node_neg_angle = libMesh::Node(0.5/std::tan(0.15),0.0);
      this->run_test(origin,-0.15,calc_end_node_neg_angle,5,3,"quad4",1);
    }

    void test_quad9_5elem()
    {
      libMesh::Point origin = libMesh::Point(0,0.5);

      libMesh::Node calc_end_node_straight = libMesh::Node(5.0,0.5);
      this->run_test(origin,0.0,calc_end_node_straight,5,4,"quad9",1);

      libMesh::Node calc_end_node_angle = libMesh::Node(0.5/std::tan(0.15),1.0);
      this->run_test(origin,0.15,calc_end_node_angle,5,3,"quad9",1);

      libMesh::Node calc_end_node_neg_angle = libMesh::Node(0.5/std::tan(0.15),0.0);
      this->run_test(origin,-0.15,calc_end_node_neg_angle,5,3,"quad9",1);
    }

    void test_quad4_2D()
    {
      libMesh::Point origin = libMesh::Point(0.0,1.5);

      libMesh::Node calc_end_node_straight = libMesh::Node(3.0,1.5);
      this->run_test(origin,0.0,calc_end_node_straight,9,5,"quad4",2);

      libMesh::Node calc_end_node_small_angle = libMesh::Node(3.0,1.5+3.0*std::tan(0.15));
      this->run_test(origin,0.15,calc_end_node_small_angle,9,5,"quad4",2);

      libMesh::Node calc_end_node_small_neg_angle = libMesh::Node(3.0,1.5+3.0*std::tan(-0.15));
      this->run_test(origin,-0.15,calc_end_node_small_neg_angle,9,5,"quad4",2);

      libMesh::Node calc_end_node_large_angle = libMesh::Node( (3.0-1.5)/std::tan(1.0), 3.0);
      this->run_test(origin,1.0,calc_end_node_large_angle,9,6,"quad4",2);

      libMesh::Node calc_end_node_large_neg_angle = libMesh::Node( (0.0-1.5)/std::tan(-1.0), 0.0);
      this->run_test(origin,-1.0,calc_end_node_large_neg_angle,9,0,"quad4",2);
    }

    void test_quad9_2D()
    {
      libMesh::Point origin = libMesh::Point(0.0,1.5);

      libMesh::Node calc_end_node_straight = libMesh::Node(3.0,1.5);
      this->run_test(origin,0.0,calc_end_node_straight,9,5,"quad9",2);

      libMesh::Node calc_end_node_small_angle = libMesh::Node(3.0,1.5+3.0*std::tan(0.15));
      this->run_test(origin,0.15,calc_end_node_small_angle,9,5,"quad9",2);

      libMesh::Node calc_end_node_small_neg_angle = libMesh::Node(3.0,1.5+3.0*std::tan(-0.15));
      this->run_test(origin,-0.15,calc_end_node_small_neg_angle,9,5,"quad9",2);

      libMesh::Node calc_end_node_large_angle = libMesh::Node( (3.0-1.5)/std::tan(1.0), 3.0);
      this->run_test(origin,1.0,calc_end_node_large_angle,9,6,"quad9",2);

      libMesh::Node calc_end_node_large_neg_angle = libMesh::Node( (0.0-1.5)/std::tan(-1.0), 0.0);
      this->run_test(origin,-1.0,calc_end_node_large_neg_angle,9,0,"quad9",2);
    }

    void fire_through_vertex()
    {
      libMesh::Point origin = libMesh::Point(0.0,0.0);

      libMesh::Node calc_end_node_straight = libMesh::Node(3.0,3.0);

      // 3x3 QUAD4 mesh
      this->run_test(origin,45.0*GRINS::Constants::pi/180.0,calc_end_node_straight,9,8,"quad4",2);

      // 3x3 QUAD9 mesh
      this->run_test(origin,45.0*GRINS::Constants::pi/180.0,calc_end_node_straight,9,8,"quad9",2);
    }

    void origin_between_elems()
    {
      libMesh::Point origin = libMesh::Point(0.0,1.0);
      libMesh::Node calc_end_node = libMesh::Node(3.0,3.0);
      libMesh::Real theta = calc_theta(origin,calc_end_node);

      this->run_test(origin,theta,calc_end_node,9,8,"quad4",2);
      this->run_test(origin,theta,calc_end_node,9,8,"quad9",2);
    }

    void quad4_off_origin()
    {
      // vector of intersection points
      std::vector<libMesh::Point> pts(4);
      pts[0] = libMesh::Point(1.4,0.8);
      pts[1] = libMesh::Point(2.2,1.3);
      pts[2] = libMesh::Point(2.0,2.25);
      pts[3] = libMesh::Point(1.3,1.6);

      // create the mesh (single trapezoidal QUAD4 element)
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(1.0,1.0),0 );
      mesh->add_point( libMesh::Point(2.0,0.5),1 );
      mesh->add_point( libMesh::Point(2.5,2.5),2 );
      mesh->add_point( libMesh::Point(1.5,2.0),3 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Quad4 );
      for (unsigned int n=0; n<4; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      this->run_test_on_all_point_combinations(pts,mesh);
    }

    void quadratic_top_quad9()
    {
      // vector of intersection points
      std::vector<libMesh::Point> pts(4);
      pts[0] = libMesh::Point(0.25,0.0);
      pts[1] = libMesh::Point(1.0,0.333);
      pts[2] = libMesh::Point(0.25,1.375);
      pts[3] = libMesh::Point(0.0,0.875);

      // create a non-rectangular QUAD9
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.0,1.0),2 );
      mesh->add_point( libMesh::Point(0.0,1.0),3 );
      mesh->add_point( libMesh::Point(0.5,0.0),4 );
      mesh->add_point( libMesh::Point(1.0,0.5),5 );
      mesh->add_point( libMesh::Point(0.5,1.5),6 );
      mesh->add_point( libMesh::Point(0.0,0.5),7 );
      mesh->add_point( libMesh::Point(0.5,0.5),8 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Quad9 );
      for (unsigned int n=0; n<9; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      this->run_test_on_all_point_combinations(pts,mesh);
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( RayfireTest2D );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
