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

namespace GRINSTesting
{
  class RayfireTestBase
  {
  protected:

    //! Given an origin and a vector of end points, we will iterate through them based on the mesh given in the supplied input file
    void run_test_with_mesh_from_file(libMesh::Point & origin, std::vector<libMesh::Point> & end_points, std::vector<unsigned int> & exit_elem_ids, const std::string & input_filename)
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/"+input_filename;
      GetPot input(filename);
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
    
  };

} // namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT

#endif //GRINS_RAYFIRE_TEST_BASE_H
