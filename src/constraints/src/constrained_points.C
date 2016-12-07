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

// This class
#include "grins/constrained_points.h"

// GRINS
#include "grins/builder_helper.h"
#include "grins/variables_parsing.h"

// libMesh
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"
#include "libmesh/node.h"
#include "libmesh/system.h"

namespace GRINS
{
  ConstrainedPoints::ConstrainedPoints( const GetPot& input,
                                        libMesh::System& system )
    : _sys(system)
  {
    std::vector<std::string> var_sections;

    BuilderHelper::parse_var_sections_vector(input, var_sections);

    for (std::vector<std::string>::const_iterator
           it = var_sections.begin(); it != var_sections.end(); ++it)
      {
        const std::string & section_name = *it;
        const std::string constraint_section =
          VariablesParsing::variables_section() + '/' +
          section_name + '/' + "ConstrainedPoints";

        if (input.have_section (constraint_section))
          {
            std::vector<std::string> point_sections =
              input.get_subsection_names(constraint_section);

            for (std::vector<std::string>::const_iterator
                   pt_it = point_sections.begin();
                   pt_it != point_sections.end(); ++pt_it)
              {
                ConstrainedPoint pt;

                pt.name = *pt_it;

                const std::string point_section =
                  constraint_section + '/' + pt.name;

                // By default we constrain the first variable in our
                // variable names list
                std::string constrained_var_name =
                  input(VariablesParsing::variables_section() + '/' +
                        section_name + '/' + "names", std::string(),
                        0);

                // But we can override that
                constrained_var_name =
                  input(point_section + "/constraint_var",
                        constrained_var_name);

                pt.var = _sys.variable_number(constrained_var_name);

                const unsigned int constrained_loc_dim =
                  input.vector_variable_size(point_section+"/constraint_location");
                libmesh_assert_less_equal(constrained_loc_dim, LIBMESH_DIM);
                for (unsigned int d=0; d != constrained_loc_dim; ++d)
                  pt(d) = input(point_section+"/constraint_location", 0.0, d );

                const unsigned int n_constraining_points_x =
                  input.vector_variable_size(point_section+"/constraining_points_x");

#ifndef NDEBUG
                const unsigned int n_constraining_points_y =
                  input.vector_variable_size(point_section+"/constraining_points_y");
                if (n_constraining_points_y)
                  libmesh_assert_equal_to(n_constraining_points_x, n_constraining_points_y);

                const unsigned int n_constraining_points_z =
                  input.vector_variable_size(point_section+"/constraining_points_z");
                if (n_constraining_points_z)
                  libmesh_assert_equal_to(n_constraining_points_x, n_constraining_points_z);

                const unsigned int n_constraining_points_var =
                  input.vector_variable_size(point_section+"/constraining_points_var");
                libmesh_assert_equal_to(n_constraining_points_x, n_constraining_points_var);

                const unsigned int n_constraining_points_coeff =
                  input.vector_variable_size(point_section+"/constraining_points_coeff");
                if (n_constraining_points_coeff)
                  libmesh_assert_equal_to(n_constraining_points_x, n_constraining_points_coeff);
#endif // !NDEBUG

                for (unsigned int n=0; n != n_constraining_points_x; ++n)
                  {
                    ConstrainingPoint constraining_pt;
                    constraining_pt(0) = input(point_section+"/constraining_points_x", 0.0, n);
                    if (LIBMESH_DIM > 1)
                      constraining_pt(1) = input(point_section+"/constraining_points_y", 0.0, n);
                    if (LIBMESH_DIM > 2)
                      constraining_pt(2) = input(point_section+"/constraining_points_z", 0.0, n);

                    constraining_pt.coeff =
                      input(point_section+"/constraining_points_coeff", libMesh::Number(0), n);

                    const std::string var_name
                      (input(point_section+"/constraining_points_var", std::string(), n));

                    constraining_pt.var = _sys.variable_number(var_name);

                    pt.constrainers.push_back(constraining_pt);
                  }

                pt.rhs = input(point_section+"/constraint_rhs", 0.0);

                pt.forbid_overwrite = input(point_section+"/forbid_overwrite", true);

                _constrained_pts.push_back(pt);
              }
          }
      }
  }

  void ConstrainedPoints::constrain ()
  {
    libMesh::MeshBase & mesh = _sys.get_mesh();
    const unsigned int sys_num = _sys.number();

    libMesh::UniquePtr<libMesh::PointLocatorBase> locator =
      mesh.sub_point_locator();

    for (std::vector<ConstrainedPoint>::const_iterator
           it = _constrained_pts.begin(); it != _constrained_pts.end(); ++it)
      {
        const ConstrainedPoint & constrained_pt = *it;
        const libMesh::Point & pt = constrained_pt;

        // Find the node at the constrained point
        const libMesh::Node *constrained_node = locator->locate_node(pt);

        // If we're on a distributed mesh, not every processor might have found
        // the node, but someone needs to have found it.
        {
          bool found_node = constrained_node;
          _sys.comm().max(found_node);
          if (!found_node)
            libmesh_error_msg("Failed to find a Node at point " << pt);
        }

        // We'll only need to add constraints on processors
        // which have found the node.  But those processors may not be
        // able to find all the constraining nodes, so all processors
        // need to stay in play here.

        const libMesh::dof_id_type constrained_dof =
          constrained_node ? constrained_node->dof_number
            (sys_num, constrained_pt.var, 0) : 0;

        libmesh_assert
          (_sys.comm().semiverify (constrained_node ? &constrained_dof
                                   : libmesh_nullptr));

        libMesh::DofConstraintRow constraint_row;

        for (std::vector<ConstrainingPoint>::const_iterator
             jt = constrained_pt.constrainers.begin();
             jt != constrained_pt.constrainers.end(); ++jt)
          {
            const ConstrainingPoint & constraining_pt = *jt;
            const libMesh::Point & pt2 = constraining_pt;

            const libMesh::Node *constraining_node =
              locator->locate_node(pt2);

            // If we're on a distributed mesh, not every processor
            // might have found the node, but someone needs to have
            // found it.
            {
              bool found_node = constraining_node;
              _sys.comm().max(found_node);
              if (!found_node)
                libmesh_error_msg
                  ("Failed to find a Node at point " << pt);
            }

            libMesh::dof_id_type constraining_dof =
              constraining_node ? constraining_node->dof_number
                (sys_num, constraining_pt.var, 0) : 0;

            libmesh_assert
              (_sys.comm().semiverify (constraining_node ?
                                       &constraining_dof :
                                       libmesh_nullptr));

            _sys.comm().max(constraining_dof);

            if (constrained_node)
              constraint_row[constraining_dof] = constraining_pt.coeff;
          }

        if (constrained_node)
          _sys.get_dof_map().add_constraint_row
            (constrained_dof, constraint_row, constrained_pt.rhs,
             constrained_pt.forbid_overwrite);
      }
  }
} // end namespace GRINS
