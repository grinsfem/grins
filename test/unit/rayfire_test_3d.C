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
#include "libmesh/fe_interface.h"
#include "libmesh/serial_mesh.h"

#include "libmesh/cell_hex8.h"
#include "libmesh/cell_hex27.h"

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
    CPPUNIT_TEST( hex_5elem_inline );
    CPPUNIT_TEST( hex_27elem_3x3x3 );
    CPPUNIT_TEST( fire_through_vertex );
    CPPUNIT_TEST( origin_between_elems );

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
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );
      
      mesh->set_mesh_dimension(3);		

      mesh->add_point( libMesh::Point(0.0,0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.0,1.0,0.0),2 );
      mesh->add_point( libMesh::Point(0.0,1.0,0.0),3 );
      mesh->add_point( libMesh::Point(0.0,0.0,1.0),4 );
      mesh->add_point( libMesh::Point(1.0,0.0,1.0),5 );
      mesh->add_point( libMesh::Point(1.0,1.0,1.0),6 );
      mesh->add_point( libMesh::Point(0.0,1.0,1.0),7 );	

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Hex8 );		
      for (unsigned int n=0; n<8; n++)		
       elem->set_node(n) = mesh->node_ptr(n);		

      mesh->prepare_for_use();

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
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(3);

      mesh->add_point( libMesh::Point(0.0,0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0,0.0),1 );
      mesh->add_point( libMesh::Point(1.0,1.0,0.0),2 );
      mesh->add_point( libMesh::Point(0.0,1.0,0.0),3 );

      mesh->add_point( libMesh::Point(0.0,0.0,1.0),4 );
      mesh->add_point( libMesh::Point(1.0,0.0,1.0),5 );
      mesh->add_point( libMesh::Point(1.0,1.0,1.0),6 );
      mesh->add_point( libMesh::Point(0.0,1.0,1.0),7 );

      mesh->add_point( libMesh::Point(0.5,0.0,0.0),8 );
      mesh->add_point( libMesh::Point(1.0,0.5,0.0),9 );
      mesh->add_point( libMesh::Point(0.5,1.0,0.0),10 );
      mesh->add_point( libMesh::Point(0.0,0.5,0.0),11 );

      mesh->add_point( libMesh::Point(0.0,0.0,0.5),12 );
      mesh->add_point( libMesh::Point(1.0,0.0,0.5),13 );
      mesh->add_point( libMesh::Point(1.0,1.0,0.5),14 );
      mesh->add_point( libMesh::Point(0.0,1.0,0.5),15 );

      mesh->add_point( libMesh::Point(0.5,0.0,1.0),16 );
      mesh->add_point( libMesh::Point(1.0,0.5,1.0),17 );
      mesh->add_point( libMesh::Point(0.5,1.0,1.0),18 );
      mesh->add_point( libMesh::Point(0.0,0.5,1.0),19 );

      mesh->add_point( libMesh::Point(0.5,0.5,0.0),20 );
      mesh->add_point( libMesh::Point(0.5,0.0,0.5),21 );
      mesh->add_point( libMesh::Point(1.5,0.5,0.5),22 ); // non-linear
      mesh->add_point( libMesh::Point(0.5,1.0,0.5),23 );
      mesh->add_point( libMesh::Point(0.0,0.5,0.5),24 );
      mesh->add_point( libMesh::Point(0.5,0.5,1.0),25 );

      mesh->add_point( libMesh::Point(0.5,0.5,0.5),26 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Hex27 );
      for (unsigned int n=0; n<27; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

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

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,"mesh_hex8_5elem.in");
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,"mesh_hex27_5elem.in");
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

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,"mesh_hex8_27elem.in");
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,"mesh_hex27_27elem.in");
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

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,"mesh_hex8_27elem.in");
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,"mesh_hex27_27elem.in");
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

      this->run_test_with_mesh_from_file(origin,pts,exit_ids,"mesh_hex8_27elem.in");
      this->run_test_with_mesh_from_file(origin,pts,exit_ids,"mesh_hex27_27elem.in");
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( RayfireTest3D );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
