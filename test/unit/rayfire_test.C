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
  class RayfireTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( RayfireTest );

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

      run_test_on_all_point_combinations(pts,mesh);

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
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = new libMesh::SerialMesh(*TestCommWorld);

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

      run_test_on_all_point_combinations(pts,mesh);

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
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = new libMesh::SerialMesh(*TestCommWorld);

      mesh->set_mesh_dimension(2);

      mesh->add_point( libMesh::Point(0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0),1 );
      mesh->add_point( libMesh::Point(2.0,1.0),2 );
      mesh->add_point( libMesh::Point(-0.5,0.5),3 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Quad4 );
      for (unsigned int n=0; n<4; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      run_test_on_all_point_combinations(pts,mesh);

    }

    void test_vertical_fire()
    {
      // create the mesh (single square QUAD4 element)
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

      libMesh::Point start_point(0.3,0.0);
      libMesh::Point end_point(0.3,1.0);

      libMesh::Real theta = GRINS::Constants::pi/2.0;

      this->run_test_with_mesh(mesh,start_point,theta,end_point,0);
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


  private:

    void run_test(libMesh::Point& origin, libMesh::Real theta, libMesh::Node& calc_end_node, unsigned int n_elem, unsigned int exit_elem, std::string elem_type, unsigned int dim)
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/mesh_"+elem_type+"_"+std::to_string(n_elem)+"elem_"+std::to_string(dim)+"D.in";
      GetPot input(filename);
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input);

      // ensure the mesh has the desired number of elements
      CPPUNIT_ASSERT_EQUAL(n_elem,mesh->n_elem());

      run_test_with_mesh(mesh,origin,theta,calc_end_node,exit_elem);
    }

    void run_test_on_all_point_combinations(std::vector<libMesh::Point> pts, GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh)
    {
      // iterate over the starting points
      for(unsigned int i=0; i<pts.size(); i++)
      {
        libMesh::Point start_point = pts[i];

        // iterate over all the intersection points
        for(unsigned int j=0; j<pts.size(); j++)
        {
          if(j==i)
            continue;

          libMesh::Point end_point = pts[j];

          libMesh::Real theta = calc_theta(start_point,end_point);

          // run the test
          this->run_test_with_mesh(mesh,start_point,theta,end_point,0);
        }

      }

    }

    void run_test_with_mesh(GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh, libMesh::Point& origin, libMesh::Real theta, libMesh::Point& calc_end_point, unsigned int exit_elem)
    {
      GRINS::SharedPtr<GRINS::RayfireMesh> rayfire = new GRINS::RayfireMesh(origin,theta);
      rayfire->init(*mesh);

      const libMesh::Elem* original_elem = mesh->elem(exit_elem);

      const libMesh::Elem* rayfire_elem = rayfire->map_to_rayfire_elem(original_elem->id());

      if (!rayfire_elem)
        libmesh_error_msg("Attempted to map an element that is not in the Rayfire");

      CPPUNIT_ASSERT_DOUBLES_EQUAL(calc_end_point(0), (*(rayfire_elem->get_node(1)))(0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(calc_end_point(1), (*(rayfire_elem->get_node(1)))(1),libMesh::TOLERANCE);
    }

    GRINS::SharedPtr<libMesh::UnstructuredMesh> build_mesh( const GetPot& input )
    {
      GRINS::MeshBuilder mesh_builder;
      return mesh_builder.build( input, *TestCommWorld );
    }

    libMesh::Real calc_theta(libMesh::Point& start, libMesh::Point end)
    {
      return std::atan2( (end(1)-start(1)), (end(0)-start(0)) );
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( RayfireTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
