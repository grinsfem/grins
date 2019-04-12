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


#ifndef GRINS_RAYFIRE_TEST_BASE_H
#define GRINS_RAYFIRE_TEST_BASE_H

#ifdef GRINS_HAVE_CPPUNIT

#include "test_comm.h"
#include "grins_test_paths.h"

// GRINS
#include "grins/rayfire_mesh.h"
#include "grins/mesh_builder.h"

// libMesh
#include "libmesh/elem.h"
#include "libmesh/getpot.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/face_quad4.h"
#include "libmesh/face_quad9.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/cell_hex27.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_tet10.h"

#include <tuple>

namespace GRINSTesting
{
  class RayfireTestBase
  {
  protected:

    //! Given an origin and a vector of end points, we will iterate through them based on the mesh given in the supplied input file
    void run_test_with_mesh_from_file(libMesh::Point & origin, std::vector<libMesh::Point> & end_points, std::vector<unsigned int> & exit_elem_ids, const std::string & input_string)
    {
      std::stringstream ss;
      ss << input_string;

      GetPot input(ss);
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input);

      // iterate over the end points for the given origin
      for (unsigned int i=0; i<end_points.size(); ++i)
          this->run_test(mesh,origin,end_points[i],exit_elem_ids[i]);

    }

    //! Given a vector of points, we will test the rayfire on all possible combinations of them
    void run_test_on_all_point_combinations(std::vector<libMesh::Point> pts, std::shared_ptr<libMesh::UnstructuredMesh> mesh)
    {
      libmesh_assert(pts.size() >= 2);

      // iterate over the starting points
      for(unsigned int i=0; i<pts.size(); i++)
        {
          libMesh::Point origin = pts[i];

          // iterate over all the intersection points
          for(unsigned int j=0; j<pts.size(); j++)
            {
              if(j==i)
                continue;

              libMesh::Point end_point = pts[j];

              // run the test
              this->run_test(mesh,origin,end_point,0);
            }

        }

    }

    //! Calculate the angles from the given orign and end_point, then test the rayfire
    void run_test(std::shared_ptr<libMesh::UnstructuredMesh> mesh, libMesh::Point & origin, libMesh::Point & end_point, unsigned int exit_elem_id)
    {
      libMesh::Real theta = this->calc_theta(origin,end_point);

      std::shared_ptr<GRINS::RayfireMesh> rayfire;
      if (mesh->mesh_dimension() == 2)
        rayfire.reset( new GRINS::RayfireMesh(origin,theta) );
      else
        {
          libMesh::Real phi = this->calc_phi(origin,end_point);
          rayfire.reset( new GRINS::RayfireMesh(origin,theta,phi) );
        }

      rayfire->init(*mesh);

      // look at the elem where the rayfire should exit
      const libMesh::Elem * original_elem = mesh->elem_ptr(exit_elem_id);
      const libMesh::Elem * rayfire_elem = rayfire->map_to_rayfire_elem(original_elem->id());

      CPPUNIT_ASSERT(rayfire_elem);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(end_point(0), (*(rayfire_elem->node_ptr(1)))(0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(end_point(1), (*(rayfire_elem->node_ptr(1)))(1),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(end_point(2), (*(rayfire_elem->node_ptr(1)))(2),libMesh::TOLERANCE);
    }

    std::shared_ptr<libMesh::UnstructuredMesh> build_mesh( const GetPot & input )
    {
      GRINS::MeshBuilder mesh_builder;
      return mesh_builder.build( input, *TestCommWorld );
    }

    libMesh::Real calc_theta(libMesh::Point & start, libMesh::Point & end)
    {
      return std::atan2( (end(1)-start(1)), (end(0)-start(0)) );
    }

    libMesh::Real calc_phi(libMesh::Point & start, libMesh::Point & end)
    {
      libMesh::Real L = (end-start).norm();
      return std::acos( (end(2)-start(2))/L );
    }

    //! Build a single square QUAD4
    std::shared_ptr<libMesh::UnstructuredMesh> build_square_quad4_elem()
    {
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

      return mesh;
    }

    //! Build a single square QUAD9
    std::shared_ptr<libMesh::UnstructuredMesh> build_square_quad9_elem()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

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

    //! Build a single square HEX8
    std::shared_ptr<libMesh::UnstructuredMesh> build_square_hex8_elem()
    {
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

      return mesh;
    }

    //! Build a single square HEX27
    std::shared_ptr<libMesh::UnstructuredMesh> build_square_hex27_elem()
    {
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
      mesh->add_point( libMesh::Point(1.0,0.5,0.5),22 );
      mesh->add_point( libMesh::Point(0.5,1.0,0.5),23 );
      mesh->add_point( libMesh::Point(0.0,0.5,0.5),24 );
      mesh->add_point( libMesh::Point(0.5,0.5,1.0),25 );

      mesh->add_point( libMesh::Point(0.5,0.5,0.5),26 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Hex27 );
      for (unsigned int n=0; n<27; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      return mesh;
    }

    //! Build a single TET4
    std::shared_ptr<libMesh::UnstructuredMesh> build_tet4_elem()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(3);

      mesh->add_point( libMesh::Point(0.0,0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0,0.0),1 );
      mesh->add_point( libMesh::Point(0.0,1.0,0.0),2 );
      mesh->add_point( libMesh::Point(0.0,0.0,1.0),3 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Tet4 );
      for (unsigned int n=0; n<4; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      return mesh;
    }

    //! Build a single TET10
    std::shared_ptr<libMesh::UnstructuredMesh> build_tet10_elem()
    {
      std::shared_ptr<libMesh::UnstructuredMesh> mesh( new libMesh::SerialMesh(*TestCommWorld) );

      mesh->set_mesh_dimension(3);

      mesh->add_point( libMesh::Point(0.0,0.0,0.0),0 );
      mesh->add_point( libMesh::Point(1.0,0.0,0.0),1 );
      mesh->add_point( libMesh::Point(0.0,1.0,0.0),2 );
      mesh->add_point( libMesh::Point(0.0,0.0,1.0),3 );
      mesh->add_point( libMesh::Point(0.5,0.0,0.0),4 );
      mesh->add_point( libMesh::Point(0.5,0.5,0.0),5 );
      mesh->add_point( libMesh::Point(0.0,0.5,0.0),6 );
      mesh->add_point( libMesh::Point(0.0,0.0,0.5),7 );
      mesh->add_point( libMesh::Point(0.5,0.0,0.5),8 );
      mesh->add_point( libMesh::Point(0.0,0.5,0.5),9 );

      libMesh::Elem* elem = mesh->add_elem( new libMesh::Tet10 );
      for (unsigned int n=0; n<10; n++)
        elem->set_node(n) = mesh->node_ptr(n);

      mesh->prepare_for_use();

      return mesh;
    }

    void amr_single_elem( std::shared_ptr<libMesh::UnstructuredMesh> & mesh,
                          libMesh::Point & origin, libMesh::Point & end_point,
                          std::vector<unsigned int> & children_in_rayfire,
                          std::vector<unsigned int> & children_not_in_rayfire )
    {
      libMesh::Real theta = this->calc_theta(origin,end_point);

      std::shared_ptr<GRINS::RayfireMesh> rayfire;
      if (mesh->mesh_dimension() == 2)
        rayfire.reset( new GRINS::RayfireMesh(origin,theta) );
      else
        {
          libMesh::Real phi = this->calc_phi(origin,end_point);
          rayfire.reset( new GRINS::RayfireMesh(origin,theta,phi) );
        }

      rayfire->init(*mesh);

      libMesh::Elem * elem = mesh->elem_ptr(0);
      CPPUNIT_ASSERT(elem);

      elem->set_refinement_flag(libMesh::Elem::RefinementState::REFINE);

      libMesh::MeshRefinement mr(*mesh);
      mr.refine_elements();

      rayfire->reinit(*mesh);

      for (unsigned int c=0; c<children_in_rayfire.size(); c++)
        CPPUNIT_ASSERT( rayfire->map_to_rayfire_elem(elem->child_ptr(children_in_rayfire[c])->id()) );

      for (unsigned int c=0; c<children_not_in_rayfire.size(); c++)
        CPPUNIT_ASSERT( !rayfire->map_to_rayfire_elem(elem->child_ptr(children_not_in_rayfire[c])->id()) );

    }

    std::string mesh_2D(const std::string & elem_type, libMesh::Real x_max, libMesh::Real y_max, unsigned int nx, unsigned int ny)
    {
      std::string text = "[Mesh]\n";
                  text +=  "[./Generation]\n";
                  text +=    "dimension = '2'\n";
                  text +=    "element_type = '"+elem_type+"'\n";
                  text +=    "x_min = '0'\n";
                  text +=    "y_min = '0'\n";
                  text +=    "x_max = '"+std::to_string(x_max)+"'\n";
                  text +=    "y_max = '"+std::to_string(y_max)+"'\n";
                  text +=    "n_elems_x = '"+std::to_string(nx)+"'\n";
                  text +=    "n_elems_y = '"+std::to_string(ny)+"'";

      return text;
    }

    std::string mesh_mixed_quad_tri()
    {
      std::string text = "[Mesh]\n";
                  text += "[./Read]\n";
                  text +=   "filename = './grids/mixed_quad_tri_square_mesh.xda'";

      return text;
    }

    std::string mesh_3D(const std::string & elem_type, libMesh::Real x_max, libMesh::Real y_max, libMesh::Real z_max, unsigned int nx, unsigned int ny, unsigned int nz)
    {
      std::string text = "[Mesh]\n";
                  text +=  "[./Generation]\n";
                  text +=    "dimension = '3'\n";
                  text +=    "element_type = '"+elem_type+"'\n";
                  text +=    "x_min = '0'\n";
                  text +=    "y_min = '0'\n";
                  text +=    "z_min = '0'\n";
                  text +=    "x_max = '"+std::to_string(x_max)+"'\n";
                  text +=    "y_max = '"+std::to_string(y_max)+"'\n";
                  text +=    "z_max = '"+std::to_string(z_max)+"'\n";
                  text +=    "n_elems_x = '"+std::to_string(nx)+"'\n";
                  text +=    "n_elems_y = '"+std::to_string(ny)+"'\n";
                  text +=    "n_elems_z = '"+std::to_string(nz)+"'";

      return text;
    }
 
  };

} // namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT

#endif //GRINS_RAYFIRE_TEST_BASE_H
