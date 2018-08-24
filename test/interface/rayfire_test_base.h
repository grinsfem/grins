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

    void run_test(libMesh::Point & origin, libMesh::Real theta, libMesh::Node & calc_end_node, unsigned int n_elem, unsigned int exit_elem, std::string elem_type, unsigned int dim)
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/mesh_"+elem_type+"_"+std::to_string(n_elem)+"elem_"+std::to_string(dim)+"D.in";
      GetPot input(filename);
      std::shared_ptr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input);

      // ensure the mesh has the desired number of elements
      CPPUNIT_ASSERT_EQUAL(n_elem,mesh->n_elem());

      this->run_test_with_mesh(mesh,origin,theta,calc_end_node,exit_elem);
    }

    void run_test_on_all_point_combinations(std::vector<libMesh::Point> pts, std::shared_ptr<libMesh::UnstructuredMesh> mesh)
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

    void run_test_with_mesh(std::shared_ptr<libMesh::UnstructuredMesh> mesh, libMesh::Point & origin, libMesh::Real theta, libMesh::Point & calc_end_point, unsigned int exit_elem)
    {
      std::shared_ptr<GRINS::RayfireMesh> rayfire( new GRINS::RayfireMesh(origin,theta) );
      rayfire->init(*mesh);

      const libMesh::Elem * original_elem = mesh->elem_ptr(exit_elem);
      const libMesh::Elem * rayfire_elem = rayfire->map_to_rayfire_elem(original_elem->id());

      CPPUNIT_ASSERT(rayfire_elem);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(calc_end_point(0), (*(rayfire_elem->node_ptr(1)))(0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(calc_end_point(1), (*(rayfire_elem->node_ptr(1)))(1),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(calc_end_point(2), (*(rayfire_elem->node_ptr(1)))(2),libMesh::TOLERANCE);
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
    
  };

} // namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT

#endif //GRINS_RAYFIRE_TEST_BASE_H
