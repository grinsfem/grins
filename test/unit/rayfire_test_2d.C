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
#include <string>

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
    CPPUNIT_TEST( test_5elem_inline );
    CPPUNIT_TEST( test_9elem_3x3 );
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

      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_quad4_elem();

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
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_quad9_elem();
      (mesh->node_ref(5))(0) = 1.5; // make right side non-linear
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
      // vector of intersection points
      std::vector<libMesh::Point> pts(2);
      pts[0] = libMesh::Point(0.3,0.0);
      pts[1] = libMesh::Point(0.3,1.0);
      
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_quad4_elem();

      this->run_test_on_all_point_combinations(pts,mesh);
    }

    void test_5elem_inline()
    {
      libMesh::Point origin = libMesh::Point(0,0.5);

      // end points
      std::vector<libMesh::Point> pts(3);
      pts[0] = libMesh::Point(5.0,0.5);
      pts[1] = libMesh::Point(0.5/std::tan(0.15),1.0);
      pts[2] = libMesh::Point(0.5/std::tan(0.15),0.0);

      // exit_elem IDs
      std::vector<unsigned int> exit_ids(3);
      exit_ids[0] = 4;
      exit_ids[1] = 3;
      exit_ids[2] = 3;

      std::string quad4_string = this->mesh_2D("QUAD4",5.0,1.0,5,1);
      std::string quad9_string = this->mesh_2D("QUAD9",5.0,1.0,5,1);

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,quad4_string);
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,quad9_string);
    }

    void test_9elem_3x3()
    {
      libMesh::Point origin = libMesh::Point(0.0,1.5);

       // end points
      std::vector<libMesh::Point> pts(5);
      pts[0] = libMesh::Point(3.0,1.5);
      pts[1] = libMesh::Point(3.0,1.5+3.0*std::tan(0.15));
      pts[2] = libMesh::Point(3.0,1.5+3.0*std::tan(-0.15));
      pts[3] = libMesh::Point((3.0-1.5)/std::tan(1.0),3.0);
      pts[4] = libMesh::Point((0.0-1.5)/std::tan(-1.0),0.0);

      // exit_elem IDs
      std::vector<unsigned int> exit_ids(5);
      exit_ids[0] = 5;
      exit_ids[1] = 5;
      exit_ids[2] = 5;
      exit_ids[3] = 6;
      exit_ids[4] = 0;

      std::string quad4_string = this->mesh_2D("QUAD4",3.0,3.0,3,3);
      std::string quad9_string = this->mesh_2D("QUAD9",3.0,3.0,3,3);

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,quad4_string);
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,quad9_string);
    }

    void fire_through_vertex()
    {
      libMesh::Point origin = libMesh::Point(0.0,0.0);

       // end points
      std::vector<libMesh::Point> pts(1);
      pts[0] = libMesh::Point(3.0,3.0);

      // exit_elem IDs
      std::vector<unsigned int> exit_ids(1);
      exit_ids[0] = 8;

      std::string quad4_string = this->mesh_2D("QUAD4",3.0,3.0,3,3);
      std::string quad9_string = this->mesh_2D("QUAD9",3.0,3.0,3,3);

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,quad4_string);
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,quad9_string);
    }

    void origin_between_elems()
    {
      
      libMesh::Point origin = libMesh::Point(0.0,1.0);

       // end points
      std::vector<libMesh::Point> pts(1);
      pts[0] = libMesh::Point(3.0,3.0);

      // exit_elem IDs
      std::vector<unsigned int> exit_ids(1);
      exit_ids[0] = 8;

      std::string quad4_string = this->mesh_2D("QUAD4",3.0,3.0,3,3);
      std::string quad9_string = this->mesh_2D("QUAD9",3.0,3.0,3,3);

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,quad4_string);
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,quad9_string);
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
