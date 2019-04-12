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
#include "libmesh/fe_interface.h"
#include "libmesh/serial_mesh.h"

#include "libmesh/cell_hex8.h"
#include "libmesh/cell_hex27.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_tet10.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class RayfireTest3D : public CppUnit::TestCase,
                        public GRINSTesting::RayfireTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( RayfireTest3D );

    CPPUNIT_TEST( hex8_all_sides );
    CPPUNIT_TEST( hex27_all_sides );
    CPPUNIT_TEST( tet4_all_sides );
    CPPUNIT_TEST( tet10_all_sides );
    CPPUNIT_TEST( hex_5elem_inline );
    CPPUNIT_TEST( hex_27elem_3x3x3 );
    CPPUNIT_TEST( fire_through_vertex );
    CPPUNIT_TEST( origin_between_elems );
    CPPUNIT_TEST( hex27_from_larger_mesh );
    CPPUNIT_TEST( another_hex27_from_larger_mesh );

    CPPUNIT_TEST_SUITE_END();

  public:

    void hex8_all_sides()
    {
      // vector of intersection points
      std::vector<libMesh::Point> pts(6);
      pts[0] = libMesh::Point(0.0,0.89,0.27); // left
      pts[1] = libMesh::Point(1.0,0.41,0.93); // right
      pts[2] = libMesh::Point(0.21,0.83,1.0); // top
      pts[3] = libMesh::Point(0.73,0.11,0.0); // bottom
      pts[4] = libMesh::Point(0.75,1.0,0.25); // back
      pts[5] = libMesh::Point(0.65,0.0,0.40); // front

      // create the mesh (single cube HEX8 element)
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_hex8_elem();
      this->run_test_on_all_point_combinations(pts,mesh);
    }

    void hex27_all_sides()
    {
      // vector of intersection points
      std::vector<libMesh::Point> pts(6);
      pts[0] = libMesh::Point(0.0,0.89,0.27);  // left
      pts[1] = libMesh::Point(1.375,0.25,0.5); // right
      pts[2] = libMesh::Point(0.21,0.83,1.0);  // top
      pts[3] = libMesh::Point(0.73,0.11,0.0);  // bottom
      pts[4] = libMesh::Point(0.75,1.0,0.25);  // back
      pts[5] = libMesh::Point(0.65,0.0,0.40);  // front

      // create the mesh (single HEX27 with non-linear side)
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_square_hex27_elem();
      (mesh->node_ref(22))(0) = 1.5; // make right face non-linear
      mesh->prepare_for_use();

      this->run_test_on_all_point_combinations(pts,mesh);
    }

    void tet4_all_sides()
    {
      // vector of intersection points
      std::vector<libMesh::Point> pts(4);
      pts[0] = libMesh::Point(0.11,0.17,0.0);
      pts[1] = libMesh::Point(0.55,0.0,0.21);
      pts[2] = libMesh::Point(0.0,0.83,0.04);
      pts[3] = libMesh::Point(0.25,0.25,0.5);

      // create the mesh (single TET4 element)
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_tet4_elem();

      this->run_test_on_all_point_combinations(pts,mesh);
    }

    void tet10_all_sides()
    {
      // vector of intersection points
      std::vector<libMesh::Point> pts(4);
      pts[0] = libMesh::Point(0.11,0.17,0.0);
      pts[1] = libMesh::Point(0.55,0.0,0.21);
      pts[2] = libMesh::Point(0.0,0.83,0.04);
      pts[3] = libMesh::Point(0.25,0.25,0.5);

      // create the mesh (single TET10 element)
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_tet10_elem();

      this->run_test_on_all_point_combinations(pts,mesh);
    }

    void hex_5elem_inline()
    {
      libMesh::Point origin = libMesh::Point(0,0.5,0.5);

      // end points
      std::vector<libMesh::Point> pts(5);
      pts[0] = libMesh::Point(5.0,0.5,0.5);    // right
      pts[1] = libMesh::Point(3.75,0.0,0.25);  // front
      pts[2] = libMesh::Point(1.25,0.77,1.0);  // top
      pts[3] = libMesh::Point(2.22,1.0,0.667); // back
      pts[4] = libMesh::Point(0.38,0.57,0.0);  // bottom

      // exit_elem IDs
      std::vector<unsigned int> exit_ids(5);
      exit_ids[0] = 4;
      exit_ids[1] = 3;
      exit_ids[2] = 1;
      exit_ids[3] = 2;
      exit_ids[4] = 0;

      std::string hex8_string = this->mesh_3D("HEX8",5.0,1.0,1.0,5,1,1);
      std::string hex27_string = this->mesh_3D("HEX27",5.0,1.0,1.0,5,1,1);

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,hex8_string);
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,hex27_string);
    }

    void hex_27elem_3x3x3()
    {
      libMesh::Point origin = libMesh::Point(0.0,1.5,1.5);

       // end points
      std::vector<libMesh::Point> pts(5);
      pts[0] = libMesh::Point(3.0,1.5,1.5);     // right
      pts[1] = libMesh::Point(1.875,0.0,2.584); // front
      pts[2] = libMesh::Point(1.1,2.67,3.0);    // top
      pts[3] = libMesh::Point(2.75,3.0,0.25);   // back
      pts[4] = libMesh::Point(2.23,2.59,0.0);   // bottom

      // exit_elem IDs
      std::vector<unsigned int> exit_ids(5);
      exit_ids[0] = 14;
      exit_ids[1] = 19;
      exit_ids[2] = 25;
      exit_ids[3] = 8;
      exit_ids[4] = 8;

      std::string hex8_string = this->mesh_3D("HEX8",3.0,3.0,3.0,3,3,3);
      std::string hex27_string = this->mesh_3D("HEX27",3.0,3.0,3.0,3,3,3);

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,hex8_string);
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,hex27_string);
    }

    void fire_through_vertex()
    {
      libMesh::Point origin = libMesh::Point(0.0,0.0,0.0);

       // end points
      std::vector<libMesh::Point> pts(1);
      pts[0] = libMesh::Point(3.0,3.0,3.0);

      // exit_elem IDs
      std::vector<unsigned int> exit_ids(1);
      exit_ids[0] = 26;

      std::string hex8_string = this->mesh_3D("HEX8",3.0,3.0,3.0,3,3,3);
      std::string hex27_string = this->mesh_3D("HEX27",3.0,3.0,3.0,3,3,3);

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,hex8_string);
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,hex27_string);
    }

    void origin_between_elems()
    {
      libMesh::Point origin = libMesh::Point(0.0,1.0,1.0);

       // end points
      std::vector<libMesh::Point> pts(1);
      pts[0] = libMesh::Point(3.0,3.0,3.0);

      // exit_elem IDs
      std::vector<unsigned int> exit_ids(1);
      exit_ids[0] = 26;

      std::string hex8_string = this->mesh_3D("HEX8",3.0,3.0,3.0,3,3,3);
      std::string hex27_string = this->mesh_3D("HEX27",3.0,3.0,3.0,3,3,3);

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,hex8_string);
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,hex27_string);
    }

    //! This test replicates a failure from a larger mesh, hence the seemingly arbitrary coordinates
    void hex27_from_larger_mesh()
    {
      libMesh::Point origin = libMesh::Point(0.0030875001992899573, 0.069253246753246747, -6.7698557104360308e-11);
      libMesh::Real theta = 1.57079632679;
      libMesh::Real phi = 1.57079632679;

      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(3);

      mesh->add_point( libMesh::Point(0.001122532015133644, 0.070792207792207781, 0.0011225320151336397),0 );
      mesh->add_point( libMesh::Point(0.0015874999999999999, 0.070792207792207795, -4.1403598545004364e-18),1 );
      mesh->add_point( libMesh::Point(0.0031562499999999993, 0.070792207792207781, -3.9482433878841953e-18),2 );
      mesh->add_point( libMesh::Point(0.0029363425824496629, 0.070792207792207795, 0.0013735659049402692),3 );

      mesh->add_point( libMesh::Point(0.001122532015133644, 0.069253246753246747, 0.0011225320151336397),4 );
      mesh->add_point( libMesh::Point(0.0015874999999999999, 0.069253246753246747, -4.0461256689816297e-18),5 );
      mesh->add_point( libMesh::Point(0.0031562499999999985, 0.069253246753246733, -3.8540092023653886e-18),6 );
      mesh->add_point( libMesh::Point(0.0029363425824496621, 0.069253246753246719, 0.0013735659049402692),7 );

      mesh->add_point( libMesh::Point(0.0014666587578616678, 0.070792207792207795, 0.00060750994887957582),8 );
      mesh->add_point( libMesh::Point(0.0023718749999999999, 0.070792207792207795, -4.0443016211923158e-18),9 );
      mesh->add_point( libMesh::Point(0.0030462962912248311, 0.070792207792207795, 0.00068678295247013264),10 );
      mesh->add_point( libMesh::Point(0.0020294372987916536, 0.070792207792207795, 0.0012480489600369543),11 );

      mesh->add_point( libMesh::Point(0.0011225320151336442, 0.070022727272727264, 0.0011225320151336399),12 );
      mesh->add_point( libMesh::Point(0.0015874999999999999, 0.070022727272727264, -4.0932427617410323e-18),13 );
      mesh->add_point( libMesh::Point(0.0031562499999999989, 0.070022727272727264, -3.9011262951247919e-18),14 );
      mesh->add_point( libMesh::Point(0.0029363425824496625, 0.070022727272727264, 0.0013735659049402692),15 );

      mesh->add_point( libMesh::Point(0.0014666587578616678, 0.069253246753246747, 0.00060750994887957582),16 );
      mesh->add_point( libMesh::Point(0.002371874999999999, 0.069253246753246733, -3.9500674356735084e-18),17 );
      mesh->add_point( libMesh::Point(0.0030462962912248303, 0.069253246753246733, 0.00068678295247013264),18 );
      mesh->add_point( libMesh::Point(0.0020294372987916531, 0.069253246753246733, 0.0012480489600369543),19 );

      mesh->add_point( libMesh::Point(0.0022564775245432493, 0.070792207792207795, 0.00064714645067485423),20 );
      mesh->add_point( libMesh::Point(0.0014666587578616675, 0.070022727272727264, 0.00060750994887957593),21 );
      mesh->add_point( libMesh::Point(0.0023718749999999994, 0.070022727272727264, -3.9971845284329117e-18),22 );
      mesh->add_point( libMesh::Point(0.0030462962912248307, 0.070022727272727264, 0.00068678295247013264),23 );
      mesh->add_point( libMesh::Point(0.0020294372987916531, 0.070022727272727264,  0.0012480489600369543),24 );
      mesh->add_point( libMesh::Point(0.0022564775245432489, 0.069253246753246733, 0.00064714645067485434),25 );

      mesh->add_point( libMesh::Point(0.0022564775245432489, 0.070022727272727264, 0.00064714645067485413),26 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Hex27 );
      for (unsigned int n=0; n<27; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta,phi) );

      rayfire->init(*mesh);

      // make sure we get a rayfire elem
      const libMesh::Elem * rayfire_elem = rayfire->map_to_rayfire_elem(0);

      CPPUNIT_ASSERT(rayfire_elem);
    }

    void another_hex27_from_larger_mesh()
    {
      libMesh::Point origin(0.00125,     0.65, 0.00237201);
      libMesh::Real theta = 1.57079632679;
      libMesh::Real phi = 1.57079632679;

      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(3);

      mesh->add_point( libMesh::Point(0.00164665,0.723,0.00393763),0);
      mesh->add_point( libMesh::Point(0.0030143,0.723,0.00300646),1);
      mesh->add_point( libMesh::Point(0.0028715,0.65,0.00285975),2);
      mesh->add_point( libMesh::Point(0.001577,0.65,0.0037506),3);
      mesh->add_point( libMesh::Point(0.00108207,0.723,0.00256845),4);
      mesh->add_point( libMesh::Point(0.00196242,0.723,0.0019608),5);
      mesh->add_point( libMesh::Point(0.00185348,0.65,0.00185105),6);
      mesh->add_point( libMesh::Point(0.00103306,0.65,0.0024282),7);
      mesh->add_point( libMesh::Point(0.00233048,0.723,0.00347205),8);
      mesh->add_point( libMesh::Point(0.0029429,0.6865,0.0029331),9);
      mesh->add_point( libMesh::Point(0.00222425,0.65,0.00330517),10);
      mesh->add_point( libMesh::Point(0.00161182,0.6865,0.00384412),11);
      mesh->add_point( libMesh::Point(0.00136436,0.723,0.00325304),12);
      mesh->add_point( libMesh::Point(0.00248836,0.723,0.00248363),13);
      mesh->add_point( libMesh::Point(0.00236249,0.65,0.0023554),14);
      mesh->add_point( libMesh::Point(0.00130503,0.65,0.0030894),15);
      mesh->add_point( libMesh::Point(0.00152225,0.723,0.00226463),16);
      mesh->add_point( libMesh::Point(0.00190795,0.6865,0.00190593),17);
      mesh->add_point( libMesh::Point(0.00144327,0.65,0.00213962),18);
      mesh->add_point( libMesh::Point(0.00105756,0.6865,0.00249832),19);
      mesh->add_point( libMesh::Point(0.00227736,0.6865,0.00338861),20);
      mesh->add_point( libMesh::Point(0.00192636,0.723,0.00286834),21);
      mesh->add_point( libMesh::Point(0.00242543,0.6865,0.00241952),22);
      mesh->add_point( libMesh::Point(0.00183376,0.65,0.0027224),23);
      mesh->add_point( libMesh::Point(0.00133469,0.6865,0.00317122),24);
      mesh->add_point( libMesh::Point(0.00148276,0.6865,0.00220212),25);
      mesh->add_point( libMesh::Point(0.00188006,0.6865,0.00279537),26);

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Hex27 );
      for (unsigned int n=0; n<27; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta,phi) );

      rayfire->init(*mesh);

      // make sure we get a rayfire elem
      const libMesh::Elem * rayfire_elem = rayfire->map_to_rayfire_elem(0);

      CPPUNIT_ASSERT(rayfire_elem);
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( RayfireTest3D );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
